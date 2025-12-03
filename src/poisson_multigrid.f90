module poisson_multigrid
  implicit none
  private
  public :: solve_poisson_multigrid

  integer, parameter :: dp = selected_real_kind(15)

contains

  subroutine solve_poisson_multigrid(phi, rhs, dx, dy, dz, nx, ny, nz, nlevels, npre, npost, tol)
    implicit none
    integer, intent(in) :: nx, ny, nz, nlevels, npre, npost
    real(kind=dp), intent(in) :: dx, dy, dz, tol
    real(kind=dp), intent(inout) :: phi(0:nx+1,0:ny+1,0:nz+1)
    real(kind=dp), intent(in) :: rhs(0:nx+1,0:ny+1,0:nz+1)

    call v_cycle(phi, rhs, dx, dy, dz, nx, ny, nz, nlevels, npre, npost, tol)
  end subroutine solve_poisson_multigrid

  !========================= MULTIGRID V-CYCLE ==========================
  recursive subroutine v_cycle(phi, rhs, dx, dy, dz, nx, ny, nz, level, npre, npost, tol)
    implicit none
    integer, intent(in) :: nx, ny, nz, level, npre, npost
    real(kind=dp), intent(in) :: dx, dy, dz, tol
    real(kind=dp), intent(in) :: rhs(0:nx+1,0:ny+1,0:nz+1)
    real(kind=dp), intent(inout) :: phi(0:nx+1,0:ny+1,0:nz+1)

    real(kind=dp), allocatable :: res(:,:,:), coarse_rhs(:,:,:), coarse_phi(:,:,:)
    integer :: i, j, k, ncx, ncy, ncz

    ! Pré-lissage (pré-smoothing)
    do i = 1, npre
       call gauss_seidel(phi, rhs, dx, dy, dz, nx, ny, nz)
    end do

    ! Calcul du résidu
    allocate(res(0:nx+1,0:ny+1,0:nz+1))
    call compute_residual(phi, rhs, res, dx, dy, dz, nx, ny, nz)

    ! Condition d'arrêt sur le résidu (norme max)
    if (level == 1) then
       call gauss_seidel(phi, rhs, dx, dy, dz, nx, ny, nz)
       deallocate(res)
       return
    end if

    ! Restriction
    ncx = nx/2
    ncy = ny/2
    ncz = nz/2
    allocate(coarse_rhs(0:ncx+1,0:ncy+1,0:ncz+1))
    allocate(coarse_phi(0:ncx+1,0:ncy+1,0:ncz+1))
    coarse_rhs = 0.0_dp
    coarse_phi = 0.0_dp

    call restrict_full_weighting(res, coarse_rhs, nx, ny, nz)

    ! Résolution récursive sur la grille grossière
    call v_cycle(coarse_phi, coarse_rhs, 2*dx, 2*dy, 2*dz, ncx, ncy, ncz, level-1, npre, npost, tol)

    ! Prolongation et correction
    call prolongation_add(phi, coarse_phi, nx, ny, nz)
    ! Impose les conditions de Dirichlet nulles explicitement
    phi(0,:,:) = 0.0_dp
    phi(nx+1,:,:) = 0.0_dp
    phi(:,0,:) = 0.0_dp
    phi(:,ny+1,:) = 0.0_dp
    phi(:,:,0) = 0.0_dp
    phi(:,:,nz+1) = 0.0_dp

    ! Post-lissage
    do i = 1, npost
       call gauss_seidel(phi, rhs, dx, dy, dz, nx, ny, nz)
    end do

    deallocate(res, coarse_rhs, coarse_phi)
  end subroutine v_cycle

  !========================= GAUSS-SEIDEL RED-BLACK ==========================
  subroutine gauss_seidel(phi, rhs, dx, dy, dz, nx, ny, nz)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=dp), intent(in) :: dx, dy, dz
    real(kind=dp), intent(inout) :: phi(0:nx+1,0:ny+1,0:nz+1)
    real(kind=dp), intent(in) :: rhs(0:nx+1,0:ny+1,0:nz+1)
    integer :: i, j, k
    real(kind=dp) :: dx2, dy2, dz2, denom

    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz
    denom = 2.0_dp*(1.0_dp/dx2 + 1.0_dp/dy2 + 1.0_dp/dz2)

    ! Red-black Gauss-Seidel
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (mod(i+j+k,2) == 0) then
                phi(i,j,k) = ( &
                     (phi(i+1,j,k) + phi(i-1,j,k))/dx2 + &
                     (phi(i,j+1,k) + phi(i,j-1,k))/dy2 + &
                     (phi(i,j,k+1) + phi(i,j,k-1))/dz2 - &
                     rhs(i,j,k)) / denom
             end if
          end do
       end do
    end do
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (mod(i+j+k,2) == 1) then
                phi(i,j,k) = ( &
                     (phi(i+1,j,k) + phi(i-1,j,k))/dx2 + &
                     (phi(i,j+1,k) + phi(i,j-1,k))/dy2 + &
                     (phi(i,j,k+1) + phi(i,j,k-1))/dz2 - &
                     rhs(i,j,k)) / denom
             end if
          end do
       end do
    end do
  end subroutine gauss_seidel

  !========================= CALCUL DU RÉSIDU ==========================
  subroutine compute_residual(phi, rhs, res, dx, dy, dz, nx, ny, nz)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=dp), intent(in) :: dx, dy, dz
    real(kind=dp), intent(in) :: phi(0:nx+1,0:ny+1,0:nz+1), rhs(0:nx+1,0:ny+1,0:nz+1)
    real(kind=dp), intent(out) :: res(0:nx+1,0:ny+1,0:nz+1)
    integer :: i, j, k
    real(kind=dp) :: dx2, dy2, dz2, lap

    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             lap = (phi(i+1,j,k) - 2*phi(i,j,k) + phi(i-1,j,k))/dx2 + &
                  (phi(i,j+1,k) - 2*phi(i,j,k) + phi(i,j-1,k))/dy2 + &
                  (phi(i,j,k+1) - 2*phi(i,j,k) + phi(i,j,k-1))/dz2
             res(i,j,k) = rhs(i,j,k) - lap
          end do
       end do
    end do
  end subroutine compute_residual

  !========================= RESTRICTION ==========================
  subroutine restrict_full_weighting(fine, coarse, nx, ny, nz)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=dp), intent(in) :: fine(0:nx+1,0:ny+1,0:nz+1)
    real(kind=dp), intent(out) :: coarse(0:nx/2+1,0:ny/2+1,0:nz/2+1)
    integer :: i, j, k, ic, jc, kc

    do kc = 1, nz/2
       k = 2*kc
       do jc = 1, ny/2
          j = 2*jc
          do ic = 1, nx/2
             i = 2*ic
             coarse(ic,jc,kc) = ( &
                  fine(i,j,k) * 0.5_dp + &
                  (fine(i+1,j,k) + fine(i-1,j,k) + fine(i,j+1,k) + fine(i,j-1,k) + fine(i,j,k+1) + fine(i,j,k-1)) * 0.0833_dp )
          end do
       end do
    end do
  end subroutine restrict_full_weighting

  !========================= PROLONGATION ==========================
  subroutine prolongation_add(fine, coarse, nx, ny, nz)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=dp), intent(inout) :: fine(0:nx+1,0:ny+1,0:nz+1)
    real(kind=dp), intent(in) :: coarse(0:nx/2+1,0:ny/2+1,0:nz/2+1)
    integer :: i, j, k, ic, jc, kc

    do kc = 1, nz/2
       do jc = 1, ny/2
          do ic = 1, nx/2
             i = 2*ic
             j = 2*jc
             k = 2*kc
             fine(i,j,k) = fine(i,j,k) + coarse(ic,jc,kc)
          end do
       end do
    end do
  end subroutine prolongation_add

end module poisson_multigrid

