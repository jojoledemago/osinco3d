module poisson

  implicit none

contains
  subroutine poisson_solver_0000(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax, idyn)
    !> Solve the Poisson equation using SOR method 
    !> with periodic bondary condition in the all directions.
    !
    !> INPUT:
    !> pp   : Initial guess and solution matrix for pressure
    !> rhs  : Right-hand side (source terms)
    !> dx   : Grid spacing in x-direction
    !> dy   : Grid spacing in y-direction
    !> dz   : Grid spacing in z-direction
    !> nx, ny, nz : Number of grid points in x, y, z direction
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence
    !> omega : Relaxation factor for SOR

    real(kind=8), intent(inout) :: pp(:,:,:), omega
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: eps
    integer, intent(in) :: nx, ny, nz, kmax, idyn

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
    real(kind=8), parameter :: factor = 1.05d0
    integer :: iter, i, j, k
    integer :: ip1, im1, jp1, jm1, kp1, km1

    print *, "* Poisson solver 0000 start"

    p_new = 0.d0
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz
    oneondx2 = 1.d0 / dx2
    oneondy2 = 1.d0 / dy2
    oneondz2 = 1.d0 / dz2
    twoondx2 = 2.d0 * oneondx2
    twoondy2 = 2.d0 * oneondy2
    twoondz2 = 2.d0 * oneondz2
    A = -(twoondx2 + twoondy2 + twoondz2)
    dmax_old = 1609.d0
    do iter = 1, kmax
       dmax = 0.d0
       do k = 1, nz
          ! Conditions aux limites périodiques en z
          if (k == 1) then
             km1 = nz
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = 1
          else
             km1 = k - 1
             kp1 = k + 1
          end if

          do j = 1, ny
             ! Conditions aux limites periodiques en y
             if (j == 1) then
                jm1 = ny
                jp1 = 2
             else if (j == ny) then
                jm1 = ny - 1
                jp1 = 1
             else
                jm1 = j - 1
                jp1 = j + 1
             end if

             do i = 1, nx
                ! Conditions aux limites périodiques en x
                if (i == 1) then
                   im1 = nx
                   ip1 = 2
                else if (i == nx) then
                   im1 = nx - 1
                   ip1 = 1
                else
                   im1 = i - 1
                   ip1 = i + 1
                end if

                ! Calculating the discretized Laplacian with SOR
                p_new(i,j,k) = (-oneondx2 * (pp(im1,j,k) + pp(ip1,j,k)) - &
                     oneondy2 * (pp(i,jm1,k) + pp(i,jp1,k)) - &
                     oneondz2 * (pp(i,j,km1) + pp(i,j,kp1)) + &
                     rhs(i,j,k)) / A
                ! Updating the maximum difference for convergence check
                d = abs(p_new(i, j, k) - pp(i, j, k))
                dmax = max(d, dmax)
                pp(i,j,k) = (1.d0 - omega) * pp(i,j,k) + omega * p_new(i,j,k)
             end do
          end do
       end do
       if (mod(iter, 100) == 0) then ! Print convergence every 100 iterations
          write(*, '(A14,I4,A1,E13.6)') ' Iteration : ', iter, ',', dmax
          if (idyn == 1) write(*,'(A23, F5.3)') ' With dynamic omega = ', omega
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       if (abs(dmax_old - dmax) < eps/1000.d0) then
          write(*, '(A29,E10.3)') " No converged at less than: ", eps
          exit
       end if
       if (iter > 1 .and. idyn == 1) then
          if (dmax > dmax_old) then
             omega = omega * (2.d0 - factor) ! Reduce omega if convergence deteriorates
          else if (dmax < 0.1d0 * dmax_old) then
             omega = min(omega * factor, 2.d0) ! Increase omega if convergence is good
          end if
       end if
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax
    if (idyn == 1) write(*,'(A23, F5.3)') ' With dynamic omega = ', omega

    return

  end subroutine poisson_solver_0000

  subroutine poisson_solver_0011(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax, idyn)
    !> Solve the Poisson equation using SOR method
    !> with free-slip boundary conditions in y directions and periodic
    !> in the other.
    !
    !> INPUT:
    !> pp   : Initial guess and solution matrix for pressure
    !> rhs  : Right-hand side (source terms)
    !> dx   : Grid spacing in x-direction
    !> dy   : Grid spacing in y-direction
    !> dz   : Grid spacing in z-direction
    !> nx, ny, nz : Number of grid points in x, y, z direction
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence
    !> omega : Relaxation factor for SOR

    real(kind=8), intent(inout) :: pp(:,:,:), omega
    real(kind=8), intent(in) :: rhs(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), intent(in) :: eps
    integer, intent(in) :: nx, ny, nz, kmax, idyn

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
    real(kind=8), parameter :: factor = 1.01d0

    integer :: iter, i, j, k
    integer :: ip1, im1, jp1, jm1, kp1, km1

    print *, "* Poisson solver 0011 start"

    dmax_old = 1609.d0
    p_new = 0.d0
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz
    oneondx2 = 1.d0 / dx2
    oneondy2 = 1.d0 / dy2
    oneondz2 = 1.d0 / dz2
    twoondx2 = 2.d0 * oneondx2
    twoondy2 = 2.d0 * oneondy2
    twoondz2 = 2.d0 * oneondz2
    A = -(twoondx2 + twoondy2 + twoondz2)
    do iter = 1, kmax
       dmax = 0.d0
       do k = 1, nz
          ! Conditions aux limites périodiques en z
          if (k == 1) then
             km1 = nz
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = 1
          else
             km1 = k - 1
             kp1 = k + 1
          end if
          do j = 1, ny
             ! Conditions aux limites libres en y
             if (j == 1) then
                jm1 = 2
                jp1 = 2
             else if (j == ny) then
                jm1 = ny - 1
                jp1 = ny - 1
             else
                jm1 = j - 1
                jp1 = j + 1
             end if
             do i = 1, nx
                ! Conditions aux limites périodiques en x
                if (i == 1) then
                   im1 = nx
                   ip1 = 2
                else if (i == nx) then
                   im1 = nx - 1
                   ip1 = 1
                else
                   im1 = i - 1
                   ip1 = i + 1
                end if
                ! Calculating the discretized Laplacian with SOR
                p_new(i,j,k) = (-oneondx2 * (pp(im1,j,k) + pp(ip1,j,k)) - &
                     oneondy2 * (pp(i,jm1,k) + pp(i,jp1,k)) - &
                     oneondz2 * (pp(i,j,km1) + pp(i,j,kp1)) + &
                     rhs(i,j,k)) / A
                ! Updating the maximum difference for convergence check
                d = abs(p_new(i, j, k) - pp(i, j, k))
                dmax = max(d, dmax)
                pp(i,j,k) = (1.d0 - omega) * pp(i,j,k) + omega * p_new(i,j,k)
             end do
          end do
       end do
       if (mod(iter, 100) == 0) then ! Print convergence every 100 iterations
          write(*, '(A14,I4,A1,E13.6)') ' Iteration : ', iter, ',', dmax
          if (idyn == 1) write(*,'(A23, F5.3)') ' With dynamic omega = ', omega
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       if (abs(dmax_old - dmax) < eps/1000.d0) then
          write(*, '(A29,E10.3)') " No converged at less than: ", eps
          exit
       end if
       if (iter > 1 .and. idyn == 1) then
          if (dmax > dmax_old) then
             omega = omega * (2.d0 - factor) ! Reduce omega if convergence deteriorates
          else if (dmax < 0.1d0 * dmax_old) then
             omega = min(omega * factor, 2.d0) ! Increase omega if convergence is good
          end if
       end if
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax
    if (idyn == 1) write(*,'(A23, F5.3)') ' With dynamic omega = ', omega

    return

  end subroutine poisson_solver_0011

  subroutine poisson_solver_111111(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax, idyn)
    !> Solve the Poisson equation using SOR method 
    !> with free-slip boundary conditions in x and y directions and periodic 
    !> in z.
    !
    !> INPUT:
    !> pp   : Initial guess and solution matrix for pressure
    !> rhs  : Right-hand side (source terms)
    !> dx   : Grid spacing in x-direction
    !> dy   : Grid spacing in y-direction
    !> dz   : Grid spacing in z-direction
    !> nx, ny, nz : Number of grid points in x, y, z direction
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !> idyn  : Flag for dynamic omega
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence
    !> omega : Relaxation factor for SOR

    real(kind=8), intent(inout) :: pp(:,:,:), omega  
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: eps
    integer, intent(in) :: nx, ny, nz, kmax, idyn

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
    real(kind=8), parameter :: factor = 1.05d0
    integer :: iter, i, j, k
    integer :: ip1, im1, jp1, jm1, kp1, km1

    print *, "* Poisson solver 111111 start"

    p_new = 0.d0
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz
    oneondx2 = 1.d0 / dx2
    oneondy2 = 1.d0 / dy2
    oneondz2 = 1.d0 / dz2
    twoondx2 = 2.d0 * oneondx2
    twoondy2 = 2.d0 * oneondy2
    twoondz2 = 2.d0 * oneondz2
    A = -(twoondx2 + twoondy2 + twoondz2)
    dmax_old = 1609.d0
    do iter = 1, kmax
       dmax = 0.d0
       do k = 1, nz
          ! Conditions aux limites free-slip en z
          if (k == 1) then
             km1 = 2
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = nz - 1
          else
             km1 = k - 1
             kp1 = k + 1
          end if

          do j = 1, ny
             ! Conditions aux limites libres en y
             if (j == 1) then
                jm1 = 2
                jp1 = 2
             else if (j == ny) then
                jm1 = ny - 1
                jp1 = ny - 1
             else
                jm1 = j - 1
                jp1 = j + 1
             end if
             do i = 1, nx
                ! Conditions aux limites glissement libre en x
                if (i == 1) then
                   im1 = 2 
                   ip1 = 2 
                else if (i == nx) then
                   im1 = nx - 1
                   ip1 = nx - 1
                else
                   im1 = i - 1
                   ip1 = i + 1
                end if
                ! Calculating the discretized Laplacian with SOR
                p_new(i,j,k) = (-oneondx2 * (pp(im1,j,k) + pp(ip1,j,k)) - &
                     oneondy2 * (pp(i,jm1,k) + pp(i,jp1,k)) - &
                     oneondz2 * (pp(i,j,km1) + pp(i,j,kp1)) + &
                     rhs(i,j,k)) / A
                ! Updating the maximum difference for convergence check
                d = abs(p_new(i, j, k) - pp(i, j, k))
                dmax = max(d, dmax)
                pp(i,j,k) = (1.d0 - omega) * pp(i,j,k) + omega * p_new(i,j,k)
             end do
          end do
       end do
       if (mod(iter, 100) == 0) then ! Print convergence every 100 iterations
          write(*, '(A14,I4,A1,E13.6)') ' Iteration : ', iter, ',', dmax
          if (idyn == 1) write(*,'(A23, F5.3)') ' With dynamic omega = ', omega
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       if (abs(dmax_old - dmax) < eps/1000.d0) then
          write(*, '(A29,E10.3)') " No converged at less than: ", eps
          exit
       end if
       if (iter > 1 .and. idyn == 1) then
          if (dmax > dmax_old) then
             omega = omega * (2.d0 - factor) ! Reduce omega if convergence deteriorates
          else if (dmax < 0.1d0 * dmax_old) then
             omega = min(omega * factor, 2.d0) ! Increase omega if convergence is good
          end if
       end if
       dmax_old = dmax
    end do

    write(*, '(A18,I4,A1,E13.6)') ' Exit iteration: ', iter, ',', dmax
    if (idyn == 1) write(*,'(A23, F5.3)') ' With dynamic omega = ', omega

    return

  end subroutine poisson_solver_111111

end module poisson
