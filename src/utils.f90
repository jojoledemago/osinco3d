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

    return
  end subroutine calcul_u_base

  subroutine normalize1D(f, base)

    !> Normalize the values of a given array to the range [base, 1]
    !> INPUT:
    !>   f      : array containing the values to be normalized
    !>   base   : lower bound of the normalized range (must be < 1)
    !> OUTPUT:
    !>   f      : array with normalized values in [base, 1]

    implicit none
    real(kind=8), intent(in) :: base
    real(kind=8), intent(inout) :: f(:)
    real(kind=8) :: max_val, min_val, range_val
    real(kind=8), parameter :: epsilon = 1.0d-12

    max_val = maxval(f)
    min_val = minval(f)
    range_val = max_val - min_val

    if (abs(range_val) < epsilon) then
       f = base
    else
       f = base + (1.d0 - base) * (f - min_val) / range_val
    end if

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
    character(len=13) :: remaining_cpu_hm, time_since_start_hm

    elapsed_time = end_time - start_time
    sum_elapsed_time = sum_elapsed_time + elapsed_time
    time_since_start = end_time - go_start

    remaining_it = max(0, nt - it)

    if (it > 0) then
       mean_elapsed_time = sum_elapsed_time / real(it, kind=8)
    else
       mean_elapsed_time = 0.d0
    end if

    remaining_cpu_time = mean_elapsed_time * real(remaining_it, kind=8)

    remaining_cpu_hm     = seconds_to_hm(int(remaining_cpu_time))
    time_since_start_hm  = seconds_to_hm(int(time_since_start))

    call print_cpu_time(elapsed_time, remaining_cpu_hm, time_since_start_hm, mean_elapsed_time)

    return

  end subroutine calcul_cpu_time

  subroutine statistics_calc(ux, uy, uz, nx, ny, nz, &
       dx, dy, dz, re, t)
    !> Calculate various statistical quantities of the flow field, including kinetic energy,
    !> dissipation rate, and enstrophy.
    !>
    !> INPUT:
    !>   ux, uy, uz  - Velocity components (3D arrays)
    !>   nx, ny, nz  - Number of grid points in the x, y, and z directions
    !>   dx, dy, dz  - Mesh size steps in the x, y, and z directions
    !>   dt          - Time step
    !>   re          - Reynolds number
    !>   t           - time 

    implicit none

    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: ux(nx,ny,nz), uy(nx,ny,nz), uz(nx,ny,nz)
    real(kind=8), intent(in) :: dx, dy, dz, re, t

    integer :: i, j, k
    real(kind=8) :: ux_mean, uy_mean, uz_mean
    real(kind=8) :: duxdx_mean, duxdy_mean, duxdz_mean
    real(kind=8) :: duydx_mean, duydy_mean, duydz_mean
    real(kind=8) :: duzdx_mean, duzdy_mean, duzdz_mean
    real(kind=8) :: e_k, eps, eps2, dzeta, xnu, t1
    real(kind=8), allocatable :: duxdx(:,:,:), duxdy(:,:,:), duxdz(:,:,:)
    real(kind=8), allocatable :: duydx(:,:,:), duydy(:,:,:), duydz(:,:,:)
    real(kind=8), allocatable :: duzdx(:,:,:), duzdy(:,:,:), duzdz(:,:,:)

    ! Allocate intermediate arrays
    allocate(duxdx(nx,ny,nz), duxdy(nx,ny,nz), duxdz(nx,ny,nz))
    allocate(duydx(nx,ny,nz), duydy(nx,ny,nz), duydz(nx,ny,nz))
    allocate(duzdx(nx,ny,nz), duzdy(nx,ny,nz), duzdz(nx,ny,nz))

    ! Calculate mean velocity components
    ux_mean = average_3d_array(ux**2)
    uy_mean = average_3d_array(uy**2)
    uz_mean = average_3d_array(uz**2)

    ! Calculate derivatives of velocity components
    call derxi(duxdx, ux, dx)
    call deryp(duxdy, ux, dy)
    call derzp(duxdz, ux, dz)
    call derxp(duydx, uy, dx)
    call deryi(duydy, uy, dy)
    call derzp(duydz, uy, dz)
    call derxp(duzdx, uz, dx)
    call deryp(duzdy, uz, dy)
    call derzi(duzdz, uz, dz)

    ! Calculate mean values of the derivatives
    duxdx_mean = average_3d_array(duxdx**2)
    duxdy_mean = average_3d_array(duxdy**2)
    duxdz_mean = average_3d_array(duxdz**2)
    duydx_mean = average_3d_array(duydx**2)
    duydy_mean = average_3d_array(duydy**2)
    duydz_mean = average_3d_array(duydz**2)
    duzdx_mean = average_3d_array(duzdx**2)
    duzdy_mean = average_3d_array(duzdy**2)
    duzdz_mean = average_3d_array(duzdz**2)

    ! Calculate kinematic viscosity
    xnu = 1.0d0 / re

    dzeta = 0.d0
    eps = 0.d0
    e_k = 0.d0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ! Calculate kinetic energy
             e_k = e_k + 0.5d0 * ( &
                  ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)
             ! Calculate dissipation rate
             eps = eps + 0.5d0 * xnu * ( &
                  (2.d0 * duxdx(i,j,k))**2 + &
                  (2.d0 * duydy(i,j,k))**2 + &
                  (2.d0 * duzdz(i,j,k))**2 + &
                  2.d0 * (duxdy(i,j,k) + duydx(i,j,k))**2 + &
                  2.d0 * (duxdz(i,j,k) + duzdx(i,j,k))**2 + &
                  2.d0 * (duydz(i,j,k) + duzdy(i,j,k))**2)

             ! Calculate enstrophy (related to dissipation)
             dzeta = dzeta + 0.5d0 * ( &
                  (duzdy(i,j,k) - duydz(i,j,k))**2 + &
                  (duxdz(i,j,k) - duzdx(i,j,k))**2 + &
                  (duydx(i,j,k) - duxdy(i,j,k))**2)

          end do
       end do
    end do
    e_k = e_k / real(nx*ny*nz, kind=8)
    eps = eps / real(nx*ny*nz, kind=8)
    dzeta = dzeta / real(nx*ny*nz, kind=8)

    ! Calculate 2nd derivatives of velocity components
    call derxxi(duxdx, ux, dx)
    call deryyp(duxdy, ux, dy)
    call derzzp(duxdz, ux, dz)
    call derxxp(duydx, uy, dx)
    call deryyi(duydy, uy, dy)
    call derzzp(duydz, uy, dz)
    call derxxp(duzdx, uz, dx)
    call deryyp(duzdy, uz, dy)
    call derzzi(duzdz, uz, dz)
    eps2 = 0.d0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ! Calculate turbulent dissipation rate
             t1 = (-xnu) * ( & 
                  ux(i,j,k) * (duxdx(i,j,k)+duxdy(i,j,k)+duxdz(i,j,k)) + &
                  uy(i,j,k) * (duydx(i,j,k)+duydy(i,j,k)+duydz(i,j,k)) + &
                  uz(i,j,k) * (duzdx(i,j,k)+duzdy(i,j,k)+duzdz(i,j,k)))
             eps2 = eps2 + t1
          end do
       end do
    end do
    eps2 = eps2 / real(nx*ny*nz, kind=8)

    call write_statistics(t, e_k, eps, eps2, dzeta, &
         ux_mean, uy_mean, uz_mean, &
         duxdx_mean, duxdy_mean, duxdz_mean, &
         duydx_mean, duydy_mean, duydz_mean, &
         duzdx_mean, duzdy_mean, duzdz_mean)

    ! Deallocate intermediate arrays
    deallocate(duxdx, duxdy, duxdz)
    deallocate(duydx, duydy, duydz)
    deallocate(duzdx, duzdy, duzdz)

    return
  end subroutine statistics_calc

  subroutine calc_visu_data_size(datasize, n1, n2, n3, itstop, itstart, nfre, num_var)
    real(kind=8), intent(out) :: datasize
    integer, intent(in) :: itstop, itstart, nfre, num_var, n1, n2, n3

    integer :: num_enr, num_elements
    real(kind=8) :: element_size

    num_elements = n1 * n2 * n3
    element_size = 7.45d-9
    num_enr = (itstop - itstart + 1) /  nfre + 1
    datasize = real(num_elements, kind=8) * &
         element_size * real(num_var, kind=8) * real(num_enr, kind=8)

    return
  end subroutine calc_visu_data_size

  subroutine smooth_field(f, nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(inout) :: f(:,:,:)
    real(kind=8) :: tmp(nx, ny, nz)
    integer :: i,j,k
    integer :: ip, im, jp, jm, kp, km
    tmp = f
    do k = 1, nz
       kp = modulo(k, nz) + 1        ! k+1 périodique
       km = modulo(k - 2, nz) + 1    ! k-1 périodique
       do j = 1, ny
          jp = modulo(j, ny) + 1
          jm = modulo(j - 2, ny) + 1
          do i = 1, nx
             ip = modulo(i, nx) + 1
             im = modulo(i - 2, nx) + 1

             f(i,j,k) = ( &
                  6.0d0 * tmp(i,j,k) + &
                  tmp(ip,j,k) + tmp(im,j,k) + &
                  tmp(i,jp,k) + tmp(i,jm,k) + &
                  tmp(i,j,kp) + tmp(i,j,km) ) / 12.0d0
          end do
       end do
    end do

    return
  end subroutine smooth_field

  subroutine apply_gaussian_filter(field, kernel, nx, ny, nz, radius)
    integer, intent(in) :: nx, ny, nz, radius
    real(kind=8), dimension(nx,ny,nz), intent(inout) :: field
    real(kind=8),  dimension(-radius:radius), intent(in) :: kernel

    integer :: i,j,k,l,ii
    real(kind=8), dimension(nx,ny,nz) :: tmp1, tmp2

    tmp1 = 0.0d0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             do l = -radius, radius
                ii = mod(i - 1 + l + nx, nx) + 1
                tmp1(i,j,k) = tmp1(i,j,k) + kernel(l) * field(ii,j,k)
             end do
          end do
       end do
    end do

    tmp2 = 0.0d0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             do l = -radius, radius
                ii = mod(j - 1 + l + ny, ny) + 1
                tmp2(i,j,k) = tmp2(i,j,k) + kernel(l) * tmp1(i,ii,k)
             end do
          end do
       end do
    end do

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             do l = -radius, radius
                ii = mod(k - 1 + l + nz, nz) + 1
                tmp1(i,j,k) = tmp1(i,j,k) + kernel(l) * tmp2(i,j,ii)
             end do
          end do
       end do
    end do

    field = tmp1

    return

  end subroutine apply_gaussian_filter

  subroutine init_kernel(radius, kernel, sigma, kernel_type)
    integer, intent(in) :: radius
    real(kind=8), dimension(-radius:radius), intent(out) :: kernel
    real(kind=8), intent(in) :: sigma
    character(len=*), intent(in) :: kernel_type

    integer :: r
    real(kind=8) :: xk, sum_k
    real(kind=8), parameter :: tiny_value = 1.0d-12

    select case (trim(kernel_type))
    case ("gaussian")
       do r = -radius, radius
          xk = real(r, kind=8)
          kernel(r) = exp(-0.5d0 * (xk / sigma) ** 2)
       end do

    case ("dog")  ! Difference of Gaussians (band-pass)
       do r = -radius, radius
          xk = real(r, kind=8)
          kernel(r) = exp(-0.5d0 * (xk / sigma)**2) - 0.5d0 * exp(-0.5d0 * (xk / (2.0d0 * sigma))**2)
       end do

    case ("mexican")
       do r = -radius, radius
          xk = real(r, kind=8)
          kernel(r) = (1.d0 - (xk/sigma)**2) * exp(-0.5d0 * (xk / sigma)**2)
       end do

    case default
       print*, "Unknown kernel type: ", kernel_type
       stop
    end select

    sum_k = sum(kernel)
    if (abs(sum_k) > tiny_value) then
       kernel = kernel / sum_k
    else
       print *, "Warning: sum_k too small for normalization!"
    end if

    return
  end subroutine init_kernel

end module utils

