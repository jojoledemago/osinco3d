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
    if (h <= 24) then
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

  ! Function that generates a random integer between min and max
  function random_between(min, max) result(r)
    implicit none
    integer, intent(in) :: min, max
    integer :: r
    real(kind=8) :: random_value

    ! Initialize the random number generator seed
    call random_seed()

    ! Generate a random real number between 0 and 1
    call random_number(random_value)

    ! Convert the real number to an integer between min and max
    r = min + int(random_value * (max - min + 1))
  end function random_between

  function cumsum(arr) result(res)
    !> Cumulative sum of the elements in an array.
    !>
    !> INPUT:
    !> arr         : Input array to compute the cumulative sum
    !>
    !> OUTPUT:
    !> res         : Array of cumulative sums

    real(kind=8), intent(in) :: arr(:)  ! Input array
    real(kind=8) :: res(size(arr))      ! Cumulative sum result
    integer :: i

    res(1) = arr(1)
    do i = 2, size(arr)
       res(i) = res(i-1) + arr(i)
    end do
  end function cumsum

  pure function compute_wave_number(i, n, L) result(k)
    !> Compute physical wave number for FFT-based spectral methods.
    !> Handles FFT frequency mapping: i in 0:n-1 → [-n/2, n/2)
    integer, intent(in) :: i      ! Index (0 to n-1)
    integer, intent(in) :: n      ! Grid size
    real(kind=8), intent(in) :: L ! Domain length
    real(kind=8) :: k             ! Resulting wave number

    k = 2.d0 * 3.141592653589793d0 * (i - n * merge(1, 0, i >= n/2)) / L
  end function compute_wave_number

end module functions
