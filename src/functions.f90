module functions
  implicit none
contains
  function seconds_to_hm(seconds) result(hm)
    implicit none
    integer, intent(in) :: seconds
    integer :: h, m, d
    character(len=13) :: hm  ! Long enough for "d d hh:mm"

    ! Convert seconds to hours and minutes
    h = seconds / 3600
    m = mod(seconds, 3600) / 60

    ! Check if the number of hours exceeds 72
    if (h <= 72) then
       ! Format as hh:mm
       write(hm, '(I2.2, A, I2.2)') h, ':', m
    else
       ! Convert to days and hours: d d hh:mm
       d = h / 24
       h = mod(h, 24)
       write(hm, '(I1, A, I2.2, A, I2.2)') d, ' d ', h, ':', m
    end if
  end function seconds_to_hm

  function function_stats(f, nx, ny, nz) result(stats)
    implicit none
    real(kind=8), dimension(:,:,:), intent(in) :: f
    integer, intent(in) :: nx, ny, nz
    real(kind=8) :: stats(6)
    real(kind=8) :: sum, f_min, f_max
    integer :: i, j, k
    integer :: imax, jmax, kmax

    sum = 0.0d0
    f_min = huge(f_min)
    f_max = -huge(f_max)
    imax = 1
    jmax = 1
    kmax = 1
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             sum = sum + f(i, j, k)
             if (f(i, j, k) < f_min) f_min = f(i, j, k)
             if (f(i, j, k) > f_max) then
                f_max = f(i, j, k)
                imax = i
                jmax = j
                kmax = k
             endif
          end do
       end do
    end do

    stats(1) = f_min
    stats(2) = f_max
    stats(3) = sum / real(nx * ny * nz, kind=8)
    stats(4) = imax
    stats(5) = jmax
    stats(6) = kmax
  end function function_stats

  function compute_velocity_magnitude(u, v, w, nx, ny, nz) result(magnitude)
    !> Compute the magnitude of the velocity vector.
    !>
    !> INPUT:
    !> u(:,:,:)       : X-component of velocity
    !> v(:,:,:)       : Y-component of velocity
    !> w(:,:,:)       : Z-component of velocity
    !> nx, ny, nz     : Number of grid points in each direction
    !> OUTPUT:
    !> magnitude(:,:,:) : Magnitude of the velocity vector
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)
    real(kind=8) :: magnitude(nx, ny, nz)

    integer :: i, j, k

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             magnitude(i, j, k) = sqrt(u(i,j,k)**2 + &
                  v(i,j,k)**2 + w(i,j,k)**2)
          end do
       end do
    end do

  end function compute_velocity_magnitude

  function average_2d_array(tab) result(avg)
    !> Calculate the average value of a 2D array
    !>
    !> INPUT:
    !> tab(ny, nz) : 2D array of real numbers
    !>
    !> OUTPUT:
    !> avg         : Average value of the 2D array

    real(kind=8), intent(in) :: tab(:,:)
    real(kind=8) :: avg
    integer :: ny, nz, i, j
    real(kind=8) :: sum

    ny = size(tab, 1)
    nz = size(tab, 2)

    sum = 0.0d0

    ! Calculate the sum of all elements in the 2D array
    do i = 1, ny
       do j = 1, nz
          sum = sum + tab(i, j)
       end do
    end do

    ! Calculate the average
    avg = sum / (ny * nz)

  end function average_2d_array

  function average_3d_array(tab) result(avg)
    !> Calculate the average value of a 3D array
    !>
    !> INPUT:
    !> tab(nx, ny, nz) : 3D array of real numbers
    !>
    !> OUTPUT:
    !> avg         : Average value of the 3D array

    real(kind=8), intent(in) :: tab(:,:,:)
    real(kind=8) :: avg
    integer :: nx, ny, nz, i, j, k
    real(kind=8) :: sum

    nx = size(tab, 1)
    ny = size(tab, 2)
    nz = size(tab, 3)

    sum = 0.0d0

    ! Calculate the sum of all elements in the 3D array
    do i = 1, nx
       do j = 1, ny
          do k = 1, nz
             sum = sum + tab(i, j, k)
          end do
       end do
    end do

    ! Calculate the average
    avg = sum / (nx * ny * nz)

  end function average_3d_array

  function contains_nan(tab) result(has_nan)
    !> Check if the 3D array contains any NaN values
    !>
    !> INPUT:
    !> tab(n1, n2, n3) : 3D array of real numbers
    !>
    !> OUTPUT:
    !> has_nan        : Logical value indicating presence of NaN
    !
    implicit none
    real(kind=8), dimension(:,:,:), intent(in) :: tab
    logical :: has_nan
    integer :: i, j, k
    integer :: n1, n2, n3

    ! Initialize result to false
    has_nan = .false.

    ! Get the dimensions of the array
    n1 = size(tab, 1)
    n2 = size(tab, 2)
    n3 = size(tab, 3)

    ! Loop through each element to check for NaN
    do i = 1, n1
       do j = 1, n2
          do k = 1, n3
             if (isnan(tab(i, j, k))) then
                has_nan = .true.
                return
             end if
          end do
       end do
    end do

  end function contains_nan

end module functions
