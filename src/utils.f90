module utils
  use derivation
  use functions
  use IOfunctions
  implicit none

contains

  subroutine calcul_u_base(u_base, u, dy)

    real(kind=8), dimension(:), intent(out) :: u_base
    real(kind=8), dimension(:), intent(in) :: u
    real(kind=8), intent(in) :: dy

    call dery1D(u_base, u, dy)
    call normalize1D(u_base, ny)

    return
  end subroutine calcul_u_base

  subroutine gene_pink_noise(pink_noise, n)

    integer, intent(in) :: n
    real(kind=8), dimension(n), intent(out) :: pink_noise

    integer, parameter :: m = 100

    real(kind=8), parameter :: pi = acos(-1.d0)
    real(kind=8) :: a, b, sum

    integer :: i, k
    real(kind=8) :: theta
    real(kind=8), dimension(m) :: phase

    ! Parameters of the Weierstrass fractal
    a = 0.5
    b = 2.0

    ! Generate random phases
    call random_seed()
    call random_number(phase)
    phase = 2.0 * pi * phase  ! Random phases between 0 and 2pi

    ! Initialize the pink noise array
    pink_noise = 0.0

    ! Generate pink noise using the Weierstrass series
    do i = 1, n
       sum = 0.0
       do k = 1, m
          theta = b**k * real(i) / real(n)
          sum = sum + a**k * cos(theta + phase(k))
       end do
       pink_noise(i) = sum
    end do

    ! Normalize the pink noise to be between -1 and 1
    pink_noise = 2.0 * (pink_noise - minval(pink_noise)) / &
         (maxval(pink_noise) - minval(pink_noise)) - 1.0

    return
  end subroutine gene_pink_noise

  subroutine normalize1D(f, n1)
    !> Normalize the values of a given array to the range [0, 1] based on the 
    !> absolute values.
    !> INPUT:
    !> n1     : Number of grid points
    !> OUTPUT: 
    !> f      : array containing the values to be normalized
    real(kind=8), intent(inout) :: f(:)
    integer, intent(in) :: n1
    real(kind=8) :: f_min, f_max, tiny_value
    real(kind=8), dimension(:), allocatable :: f_abs
    integer :: i

    ! Allocate a temporary array for absolute values
    allocate(f_abs(n1))
    tiny_value = tiny(1.d0)

    ! Compute the absolute values
    do i = 1, n1
       f_abs(i) = abs(f(i))
    end do

    ! Find the minimum and maximum values in the absolute array
    f_min = minval(f_abs)
    f_max = maxval(f_abs)

    ! Normalize the original array values based on absolute values
    if ((f_max - f_min) > tiny_value) then
       do i = 1, n1
          f(i) = (abs(f(i)) - f_min) / (f_max - f_min)
       end do
    else
       ! If all values are the same, set them all to 0
       f = 0.0d0
    end if

    ! Deallocate the temporary array
    deallocate(f_abs)

  end subroutine normalize1D

  subroutine normalize2D(f, n1, n2)
    !> Normalize the values of a given array to the range [0, 1] based on the 
    !> absolute values.
    !> INPUT:
    !> n1     : Number of grid points
    !> n2     : Number of grid points
    !> OUTPUT: 
    !> f      : array containing the values to be normalized
    real(kind=8), intent(inout) :: f(:,:)
    integer, intent(in) :: n1, n2
    real(kind=8) :: f_min, f_max, tiny_value
    real(kind=8), dimension(:,:), allocatable :: f_abs
    integer :: i, j

    ! Allocate a temporary array for absolute values
    allocate(f_abs(n1,n2))
    tiny_value = tiny(1.d0)

    ! Compute the absolute values
    do j = 1, n2
       do i = 1, n1
          f_abs(i,j) = abs(f(i,j))
       end do
    end do

    ! Find the minimum and maximum values in the absolute array
    f_min = minval(f_abs)
    f_max = maxval(f_abs)

    ! Normalize the original array values based on absolute values
    if ((f_max - f_min) > tiny_value) then
       do j = 1, n2
          do i = 1, n1
             f(i, j) = (abs(f(i, j)) - f_min) / (f_max - f_min)
          end do
       end do
    else
       ! If all values are the same, set them all to 0
       f = 0.0d0
    end if

    ! Deallocate the temporary array
    deallocate(f_abs)

  end subroutine normalize2D

  subroutine calculate_residuals(u, v, w, old_u, old_v, old_w, dt, t_ref, u_ref, nx, ny, nz, it)
    integer, intent(in) :: nx, ny, nz, it
    real(kind=8), dimension(nx,ny,nz), intent(in) :: u, v, w, old_u, old_v, old_w
    real(kind=8), intent(in) :: dt, t_ref, u_ref
    real(kind=8) :: res_u, res_v, res_w, aa, bb, cc, linf_u, linf_v, linf_w
    real(kind=8), dimension(nx,ny,nz) :: a, b, c
    real(kind=8), parameter :: smallest_value = tiny(1.d0)
    integer :: i, j, k, ia, ja, ka, ib, jb, kb, ic, jc, kc

    print *, "* Residual on velocity"

    res_u = 0.d0
    res_v = 0.d0
    res_w = 0.d0
    linf_u = 0.d0
    linf_v = 0.d0
    linf_w = 0.d0

    do k = 2, nz-1
       do j = 2, ny-1
          do i = 2, nx-1
             res_u = res_u + (abs(old_u(i,j,k) - u(i,j,k)) / (2.D0 * dt))**2.d0
             res_v = res_v + (abs(old_v(i,j,k) - v(i,j,k)) / (2.D0 * dt))**2.d0
             res_w = res_w + (abs(old_w(i,j,k) - w(i,j,k)) / (2.D0 * dt))**2.d0
             a(i,j,k) = abs(old_u(i,j,k) - u(i,j,k)) / (2.D0 * dt)
             b(i,j,k) = abs(old_v(i,j,k) - v(i,j,k)) / (2.D0 * dt)
             c(i,j,k) = abs(old_w(i,j,k) - w(i,j,k)) / (2.D0 * dt)
             linf_u = max(a(i,j,k), linf_u)
             linf_v = max(b(i,j,k), linf_v)
             linf_w = max(c(i,j,k), linf_w)
          end do
       end do
    end do
    do k = 2, nz-1
       do j = 2, ny-1
          do i = 2, nx-1
             if (abs(linf_u - a(i,j,k)) <= smallest_value) then
                ia = i
                ja = j
                ka = k
             end if
             if (abs(linf_v - b(i,j,k)) <= smallest_value) then
                ib = i
                jb = j
                kb = k
             end if
             if (abs(linf_w - c(i,j,k)) <= smallest_value) then
                ic = i
                jc = j
                kc = k
             end if
          end do
       end do
    end do

    aa = (t_ref / u_ref) * linf_u
    bb = (t_ref / u_ref) * linf_v
    cc = (t_ref / u_ref) * linf_w
    res_u = (t_ref / u_ref) * sqrt((1.d0 / real(nx*ny*nz)) * res_u)
    res_v = (t_ref / u_ref) * sqrt((1.d0 / real(nx*ny*nz)) * res_v)
    res_w = (t_ref / u_ref) * sqrt((1.d0 / real(nx*ny*nz)) * res_w)

    call print_residuals(res_u, res_v, res_w, aa, ia, ja, ka, bb, ib, jb, kb, cc, ic, jc, kc)
    call save_residu(it, res_u, res_v, res_w)

    if (res_u > 1.d6 .or. res_u > 1.d6 .or. res_u > 1.d6) then
       print *, "residue too high"
       stop
    end if

  end subroutine calculate_residuals

  subroutine old_values(u, v, w, old_u, old_v, old_w, nx, ny, nz)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx,ny,nz), intent(in) :: u, v, w
    real(kind=8), dimension(nx,ny,nz), intent(inout) :: old_u, old_v, old_w

    old_u = u
    old_v = v
    old_w = w

    return
  end subroutine old_values

  subroutine compute_cfl(cflx, cfly, cflz, ux, uy, uz, dx, dy, dz, dt)

    !> Subroutine to compute the CFL numbers in the x, y, and z directions.
    !> The CFL numbers are computed using the velocity components
    !> and the cell sizes in each direction.
    !>
    !> INPUT
    !>   ux(:,:,:):  Velocity field in the x direction
    !>   uy(:,:,:):  Velocity field in the y direction
    !>   uz(:,:,:):  Velocity field in the z direction
    !>   dx       :  Cell size in the x direction
    !>   dy       :  Cell size in the y direction
    !>   dz       :  Cell size in the z direction
    !>   dt       :  Time step size
    !>
    !> OUTPUT
    !>   cflx     :  CFL number in the x direction
    !>   cfly     :  CFL number in the y direction
    !>   cflz     :  CFL number in the z direction

    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz, dt
    real(kind=8), intent(out) :: cflx, cfly, cflz

    cflx = maxval(abs(ux)) * dt /dx
    cfly = maxval(abs(uy)) * dt /dy
    cflz = maxval(abs(uz)) * dt /dz

    return

  end subroutine compute_cfl

  subroutine calcul_cpu_time(go_start, start_time, end_time, it, nt, &
       sum_elapsed_time)
    integer, intent(in) :: it, nt
    real(kind=8), intent(in) :: go_start, start_time, end_time
    real(kind=8), intent(inout) :: sum_elapsed_time
    real(kind=8) :: elapsed_time, remaining_cpu_time, time_since_start, mean_elapsed_time
    integer :: remaining_it
    character(len=7) :: remaining_cpu_hm, time_since_start_hm

    elapsed_time = end_time - start_time
    sum_elapsed_time = sum_elapsed_time + elapsed_time
    time_since_start = end_time - go_start
    remaining_it = nt - it
    mean_elapsed_time = sum_elapsed_time / real(it, kind=8)
    remaining_cpu_time = mean_elapsed_time * &
         real(remaining_it, kind=8)

    remaining_cpu_hm = seconds_to_hm(int(remaining_cpu_time))
    time_since_start_hm = seconds_to_hm(int(time_since_start))

    call print_cpu_time(elapsed_time, remaining_cpu_hm, time_since_start_hm, mean_elapsed_time)

    return

  end subroutine calcul_cpu_time

end module utils

