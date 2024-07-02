module poisson

  use boundary_conditions

  implicit none

contains
  subroutine poisson_solver_0000(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax)
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
    !> omega : Relaxation factor for SOR
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence

    real(kind=8), intent(inout) :: pp(:,:,:)   
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: omega, eps
    integer, intent(in) :: nx, ny, nz, kmax

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
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
             km1 = nz - 1
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = 2
          else
             km1 = k - 1
             kp1 = k + 1
          end if

          do j = 1, ny
             ! Conditions aux limites periodiques en y
             if (j == 1) then
                jm1 = ny - 1
                jp1 = 2
             else if (j == ny) then
                jm1 = ny - 1
                jp1 = 2
             else
                jm1 = j - 1
                jp1 = j + 1
             end if

             do i = 1, nx
                ! Conditions aux limites périodiques en x
                if (i == 1) then
                   im1 = nx - 1
                   ip1 = 2
                else if (i == nx) then
                   im1 = nx - 1
                   ip1 = 2
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
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       if (abs(dmax_old - dmax) < eps/1000.d0) then
          print *, "* No converged at less than: ", eps
          exit
       end if
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax

    return

  end subroutine poisson_solver_0000

  subroutine poisson_solver_0011(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax)
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
    !> omega : Relaxation factor for SOR
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence

    real(kind=8), intent(inout) :: pp(:,:,:)   
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: omega, eps
    integer, intent(in) :: nx, ny, nz, kmax

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
    integer :: iter, i, j, k
    integer :: ip1, im1, jp1, jm1, kp1, km1

    print *, "* Poisson solver 0011 start"

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
             km1 = nz - 1
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = 2
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
                   im1 = nx - 1
                   ip1 = 2
                else if (i == nx) then
                   im1 = nx - 1
                   ip1 = 2
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
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       !if (abs(dmax_old - dmax) < 1.d-9) exit
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax

    return

  end subroutine poisson_solver_0011

  subroutine poisson_solver_111111(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax)
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
    !> omega : Relaxation factor for SOR
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence

    real(kind=8), intent(inout) :: pp(:,:,:)   
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: omega, eps
    integer, intent(in) :: nx, ny, nz, kmax

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
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
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       !if (abs(dmax_old - dmax) < 1.d-9) exit
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax

    return

  end subroutine poisson_solver_111111

  subroutine poisson_solver_2200(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax)
    !> Solve the Poisson equation using SOR method 
    !> with periodic boundary conditions in y and z directions and 
    !> inflow/outflow (Neuman dp/dx = 0) in x direction.
    !
    !> INPUT:
    !> pp   : Initial guess and solution matrix for pressure
    !> rhs  : Right-hand side (source terms)
    !> dx   : Grid spacing in x-direction
    !> dy   : Grid spacing in y-direction
    !> dz   : Grid spacing in z-direction
    !> nx, ny, nz : Number of grid points in x, y, z direction
    !> omega : Relaxation factor for SOR
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence

    real(kind=8), intent(inout) :: pp(:,:,:)   
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: omega, eps
    integer, intent(in) :: nx, ny, nz, kmax

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
    integer :: iter, i, j, k
    integer :: ip1, im1, jp1, jm1, kp1, km1

    print *, "* Poisson solver 2200 start"

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
             km1 = nz - 1
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = 2
          else
             km1 = k - 1
             kp1 = k + 1
          end if
          do j = 1, ny
             ! Conditions aux limites periodiques en y
             if (j == 1) then
                jm1 = ny - 1
                jp1 = 2
             else if (j == ny) then
                jm1 = ny - 1
                jp1 = 2
             else
                jm1 = j - 1
                jp1 = j + 1
             end if
             do i = 2, nx-1
                im1 = i - 1
                ip1 = i + 1
                p_new(i,j,k) = (-oneondx2 * (pp(im1,j,k) + pp(ip1,j,k)) - &
                     oneondy2 * (pp(i,jm1,k) + pp(i,jp1,k)) - &
                     oneondz2 * (pp(i,j,km1) + pp(i,j,kp1)) + &
                     rhs(i,j,k)) / A
                d = abs(p_new(i, j, k) - pp(i, j, k))
                dmax = max(d, dmax)
                pp(i,j,k) = (1.d0 - omega) * pp(i,j,k) + omega * p_new(i,j,k)
             end do
          end do
       end do
       ! Neuman bondaries condition in x-direction
       call pressure_neuman_bc_x(pp, nx)
       if (mod(iter, 100) == 0) then ! Print convergence every 100 iterations
          write(*, '(A14,I4,A1,E13.6)') ' Iteration : ', iter, ',', dmax
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       if (abs(dmax_old - dmax) < eps/1000.d0) then
          print *, "* No converged at less than: ", eps
          exit
       end if
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax

    return

  end subroutine poisson_solver_2200

  subroutine poisson_solver_2211(pp, rhs, dx, dy, dz, nx, ny, nz, &
       omega, eps, kmax)
    !> Solve the Poisson equation using SOR method 
    !> with free-slip boundary conditions in y directions, 
    !> periodic boundary conditions in z directions and 
    !> inflow/outflow in x-direction.
    !
    !> INPUT:
    !> pp   : Initial guess and solution matrix for pressure
    !> rhs  : Right-hand side (source terms)
    !> dx   : Grid spacing in x-direction
    !> dy   : Grid spacing in y-direction
    !> dz   : Grid spacing in z-direction
    !> nx, ny, nz : Number of grid points in x, y, z direction
    !> omega : Relaxation factor for SOR
    !> eps   : Convergence criterion
    !> kmax  : Maximum number of iterations
    !
    !> OUTPUT:
    !> pp   : Updated pressure field after convergence

    real(kind=8), intent(inout) :: pp(:,:,:)   
    real(kind=8), intent(in) :: rhs(:,:,:)     
    real(kind=8), intent(in) :: dx, dy, dz     
    real(kind=8), intent(in) :: omega, eps
    integer, intent(in) :: nx, ny, nz, kmax

    real(kind=8), dimension(nx,ny,nz) :: p_new
    real(kind=8) :: dx2, dy2, dz2, A, oneondx2, oneondy2, oneondz2
    real(kind=8) :: twoondx2, twoondy2, twoondz2
    real(kind=8) :: dmax, d, dmax_old
    integer :: iter, i, j, k
    integer :: ip1, im1, jp1, jm1, kp1, km1

    print *, "* Poisson solver 2211 start"

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
             km1 = nz - 1
             kp1 = 2
          else if (k == nz) then
             km1 = nz - 1
             kp1 = 2
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
             do i = 2, nx-1
                ip1 = i + 1
                im1 = i - 1
                ! Calculating the discretized Laplacian with SOR
                p_new(i,j,k) = (-oneondx2 * (pp(im1,j,k) + pp(ip1,j,k)) - &
                     oneondy2 * (pp(i,jm1,k) + pp(i,jp1,k)) - &
                     oneondz2 * (pp(i,j,km1) + pp(i,j,kp1)) + &
                     rhs(i,j,k)) / A
                ! Updating the maximum difference for convergence check
                d = abs(p_new(i, j, k) - pp(i, j, k))
                dmax = max(d, dmax)
                pp(i,j,k) = (1.d0 - omega) * pp(i,j,k)  + omega * p_new(i,j,k) 
             end do
          end do
       end do
       ! Neuman bondaries condition in x-direction
       call pressure_neuman_bc_x(pp, nx)
       if (mod(iter, 100) == 0) then ! Print convergence every 100 iterations
          write(*, '(A14,I4,A1,E13.6)') ' Iteration : ', iter, ',', dmax
       end if
       if (dmax < eps) exit ! Condition d'arrêt
       !if (abs(dmax_old - dmax) < 1.d-9) exit
       dmax_old = dmax
    end do

    write(*, '(A19,I4,A1,E13.6)') ' Exit iteration : ', iter, ',', dmax

    return

  end subroutine poisson_solver_2211

end module poisson
