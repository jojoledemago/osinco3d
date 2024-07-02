module functions
  implicit none
contains
  function seconds_to_hm(seconds) result(hm)
    implicit none
    integer, intent(in) :: seconds
    integer :: h, m
    character(len=5) :: hm

    ! Convert seconds to hours, minutes, and remaining seconds
    h = seconds / 3600
    m = mod(seconds, 3600) / 60

    ! Convert integers to character strings and format as h:mm:ss
    write(hm, '(I02.2, A, I2.2)') h, ':', m
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

end module functions
