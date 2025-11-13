module derivation

  implicit none
contains

  subroutine derx_00(df, f, dx)

    !> X first derivation of function f with O(dx6) scheme
    !> Case for Periodic boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dx            : mesh size step
    !> OUTPUT: df(nx, ny, nz): X derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydx
    real(kind=8) :: a, b, c
    integer :: nx, i

    nx = size(f, 1)

    ! Coefficients for derivation
    sixtydx  = 60.d0 * dx

    a  =  1.d0 / sixtydx
    b  =  9.d0 / sixtydx
    c  = 45.d0 / sixtydx

    ! Implementation of periodic boundary conditions for X-axis

    df(1, :, :) = a * (f(4, :, :) - f(nx-2, :, :)) - &
         b * (f(3, :, :) - f(nx-1, :, :)) + &
         c * (f(2, :, :) - f(nx, :, :))
    df(2, :, :) = a * (f(5, :, :) - f(nx-1, :, :)) - &
         b * (f(4, :, :) - f(nx, :, :)) + &
         c * (f(3, :, :) - f(1   , :, :))
    df(3, :, :) = a * (f(6, :, :) - f(nx, :, :)) - &
         b * (f(5, :, :) - f(1   , :, :)) + &
         c * (f(4, :, :) - f(2   , :, :))
    do i = 4, nx-3
       df(i, :, :) = a * (f(i+3, :, :) - f(i-3, :, :)) - &
            b * (f(i+2, :, :) - f(i-2, :, :)) + &
            c * (f(i+1, :, :) - f(i-1, :, :))
    end do
    df(nx-2, : ,:) = a * (f(1   , :, :) - f(nx-5, :, :)) - &
         b * (f(nx  , :, :) - f(nx-4, :, :)) + &
         c * (f(nx-1, :, :) - f(nx-3, :, :))
    df(nx-1, : ,:) = a * (f(2   , :, :) - f(nx-4, :, :)) - &
         b * (f(1   , :, :) - f(nx-3, :, :)) + &
         c * (f(nx  , :, :) - f(nx-2, :, :))
    df(nx  , :, :) = a * (f(3   , :, :) - f(nx-3, :, :)) - &
         b * (f(2   , :, :) - f(nx-2, :, :)) + &
         c * (f(1   , :, :) - f(nx-1, :, :))

    return

  end subroutine derx_00

  subroutine derxp_11(df, f, dx)
    !> X first derivation of even function f with O(dx6) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dx            : mesh size step
    !> OUTPUT: df(nx, ny, nz): X derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydx
    real(kind=8) :: a, b, c
    integer :: nx, i

    nx = size(f, 1)

    ! Coefficients for derivation
    sixtydx = 60.d0 * dx
    a = 1.d0 / sixtydx
    b = 9.d0 / sixtydx
    c = 45.d0 / sixtydx

    ! Boundary handling for periodic conditions
    df(1, :, :) = 0.d0
    df(2, :, :) = a * (f(5, :, :) - f(3, :, :)) - &
         b * (f(4, :, :) - f(2, :, :)) + &
         c * (f(3, :, :) - f(1, :, :))
    df(3, :, :) = a * (f(6, :, :) - f(2, :, :)) - &
         b * (f(5, :, :) - f(1, :, :)) + &
         c * (f(4, :, :) - f(2, :, :))
    do i = 4, nx-3
       df(i, :, :) = a * (f(i+3, :, :) - f(i-3, :, :)) - &
            b * (f(i+2, :, :) - f(i-2, :, :)) + &
            c * (f(i+1, :, :) - f(i-1, :, :))
    end do
    df(nx-2, :, :) = a * (f(nx-1, :, :) - f(nx-5, :, :)) - &
         b * (f(nx, :, :) - f(nx-4, :, :)) + &
         c * (f(nx-1, :, :) - f(nx-3, :, :))
    df(nx-1, :, :) = a * (f(nx-2, :, :) - f(nx-4, :, :)) - &
         b * (f(nx-1, :, :) - f(nx-3, :, :)) + &
         c * (f(nx, :, :) - f(nx-2, :, :))
    df(nx, :, :) = 0.d0

    return

  end subroutine derxp_11

  subroutine derxi_11(df, f, dx)

    !> X first derivative of even function f with O(dx6) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dx            : mesh size step
    !> OUTPUT: df(nx, ny, nz): X derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydx
    real(kind=8) :: a, b, c
    integer :: nx, i

    nx = size(f, 1)

    ! Coefficients for derivation
    sixtydx = 60.d0 * dx
    a = 1.d0 / sixtydx
    b = 9.d0 / sixtydx
    c = 45.d0 / sixtydx

    ! Boundary handling for free-slip conditions
    df(1, :, :) = a * (f(4, :, :) + f(4, :, :)) - &
         b * (f(3, :, :) + f(3, :, :)) + &
         c * (f(2, :, :) + f(2, :, :))
    df(2, :, :) = a * (f(5, :, :) + f(3, :, :)) - &
         b * (f(4, :, :) + f(2, :, :)) + &
         c * (f(3, :, :) - f(1, :, :))
    df(3, :, :) = a * (f(6, :, :) + f(2, :, :)) - &
         b * (f(5, :, :) - f(1, :, :)) + &
         c * (f(4, :, :) - f(2, :, :))
    do i = 4, nx-3
       df(i, :, :) = a * (f(i+3, :, :) - f(i-3, :, :)) - &
            b * (f(i+2, :, :) - f(i-2, :, :)) + &
            c * (f(i+1, :, :) - f(i-1, :, :))
    end do
    df(nx-2, :, :) = a * (-f(nx-1, :, :) - f(nx-5, :, :)) - &
         b * (f(nx, :, :) - f(nx-4, :, :)) + &
         c * (f(nx-1, :, :) - f(nx-3, :, :))
    df(nx-1, :, :) = a * (-f(nx-2, :, :) - f(nx-4, :, :)) - &
         b * (-f(nx-1, :, :) - f(nx-3, :, :)) + &
         c * (f(nx, :, :) - f(nx-2, :, :))
    df(nx, :, :) = a * (-f(nx-3, :, :) - f(nx-3, :, :)) - &
         b * (-f(nx-2, :, :) - f(nx-2, :, :)) + &
         c * (-f(nx-1, :, :) - f(nx-1, :, :))

    return

  end subroutine derxi_11

  subroutine dery_00(df, f, dy)

    !> Y first derivation of function f with O(dy6) scheme
    !> Case for Periodic boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dy            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Y derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydy
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f, 2)

    ! Coefficients for derivation
    sixtydy = 60.d0 * dy
    a = 1.d0 / sixtydy
    b = 9.d0 / sixtydy
    c = 45.d0 / sixtydy

    ! Boundary handling for periodic conditions
    df(:, 1, :) = a * (f(:, 4, :) - f(:, ny-2, :)) - &
         b * (f(:, 3, :) - f(:, ny-1, :)) + &
         c * (f(:, 2, :) - f(:, ny, :))
    df(:, 2, :) = a * (f(:, 5, :) - f(:, ny-1, :)) - &
         b * (f(:, 4, :) - f(:, ny, :)) + &
         c * (f(:, 3, :) - f(:, 1   , :))
    df(:, 3, :) = a * (f(:, 6, :) - f(:, ny, :)) - &
         b * (f(:, 5, :) - f(:, 1   , :)) + &
         c * (f(:, 4, :) - f(:, 2   , :))
    do j = 4, ny-3
       df(:, j, :) = a * (f(:, j+3, :) - f(:, j-3, :)) - &
            b * (f(:, j+2, :) - f(:, j-2, :)) + &
            c * (f(:, j+1, :) - f(:, j-1, :))
    end do
    df(:, ny-2, :) = a * (f(:, 1, :) - f(:, ny-5, :)) - &
         b * (f(:, ny, :) - f(:, ny-4, :)) + &
         c * (f(:, ny-1, :) - f(:, ny-3, :))
    df(:, ny-1, :) = a * (f(:, 2, :) - f(:, ny-4, :)) - &
         b * (f(:, 1, :) - f(:, ny-3, :)) + &
         c * (f(:, ny, :) - f(:, ny-2, :))
    df(:, ny, :) = a * (f(:, 3, :) - f(:, ny-3, :)) - &
         b * (f(:, 2, :) - f(:, ny-2, :)) + &
         c * (f(:, 1, :) - f(:, ny-1, :))

    return

  end subroutine dery_00

  subroutine deryp_11(df, f, dy)

    !> Y first derivation of odd function f with O(dy6) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dy            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Y derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydy
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f, 2)

    ! Coefficients for derivation
    sixtydy = 60.d0 * dy
    a = 1.d0 / sixtydy
    b = 9.d0 / sixtydy
    c = 45.d0 / sixtydy

    ! Boundary handling for periodic conditions
    df(:, 1, :) = 0.d0
    df(:, 2, :) = a * (f(:, 5, :) - f(:, 3, :)) - &
         b * (f(:, 4, :) - f(:, 2, :)) + &
         c * (f(:, 3, :) - f(:, 1, :))
    df(:, 3, :) = a * (f(:, 6, :) - f(:, 2, :)) - &
         b * (f(:, 5, :) - f(:, 1, :)) + &
         c * (f(:, 4, :) - f(:, 2, :))
    do j = 4, ny-3
       df(:, j, :) = a * (f(:, j+3, :) - f(:, j-3, :)) - &
            b * (f(:, j+2, :) - f(:, j-2, :)) + &
            c * (f(:, j+1, :) - f(:, j-1, :))
    end do
    df(:, ny-2, :) = a * (f(:, ny-1, :) - f(:, ny-5, :)) - &
         b * (f(:, ny, :) - f(:, ny-4, :)) + &
         c * (f(:, ny-1, :) - f(:, ny-3, :))
    df(:, ny-1, :) = a * (f(:, ny-2, :) - f(:, ny-4, :)) - &
         b * (f(:, ny-1, :) - f(:, ny-3, :)) + &
         c * (f(:, ny, :) - f(:, ny-2, :))
    df(:, ny, :) = 0.d0

    return

  end subroutine deryp_11

  subroutine deryi_11(df, f, dy)

    !> Y first derivation of even function f with O(dy6) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dy            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Y derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydy
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f, 2)

    ! Coefficients for derivation
    sixtydy = 60.d0 * dy
    a = 1.d0 / sixtydy
    b = 9.d0 / sixtydy
    c = 45.d0 / sixtydy

    ! Boundary handling for periodic conditions
    df(:, 1, :) = a * (f(:, 4, :) + f(:, 4, :)) - &
         b * (f(:, 3, :) + f(:, 3, :)) + &
         c * (f(:, 2, :) + f(:, 2, :))
    df(:, 2, :) = a * (f(:, 5, :) + f(:, 3, :)) - &
         b * (f(:, 4, :) + f(:, 2, :)) + &
         c * (f(:, 3, :) - f(:, 1, :))
    df(:, 3, :) = a * (f(:, 6, :) + f(:, 2, :)) - &
         b * (f(:, 5, :) - f(:, 1, :)) + &
         c * (f(:, 4, :) - f(:, 2, :))
    do j = 4, ny-3
       df(:, j, :) = a * (f(:, j+3, :) - f(:, j-3, :)) - &
            b * (f(:, j+2, :) - f(:, j-2, :)) + &
            c * (f(:, j+1, :) - f(:, j-1, :))
    end do
    df(:, ny-2, :) = a * (-f(:, ny-1, :) - f(:, ny-5, :)) - &
         b * (f(:, ny, :) - f(:, ny-4, :)) + &
         c * (f(:, ny-1, :) - f(:, ny-3, :))
    df(:, ny-1, :) = a * (-f(:, ny-2, :) - f(:, ny-4, :)) - &
         b * (-f(:, ny-1, :) - f(:, ny-3, :)) + &
         c * (f(:, ny, :) - f(:, ny-2, :))
    df(:, ny, :) = a * (-f(:, ny-3, :) - f(:, ny-3, :)) - &
         b * (-f(:, ny-2, :) - f(:, ny-2, :)) + &
         c * (-f(:, ny-1, :) - f(:, ny-1, :))

    return

  end subroutine deryi_11

  subroutine derz_00(df, f, dz)

    !> Z first derivation of function f with O(dz6) scheme
    !> Case for Periodic boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dz            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Z derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydz
    real(kind=8) :: a, b, c
    integer :: nz, k

    nz = size(f, 3)

    ! Coefficients for derivation
    sixtydz = 60.d0 * dz
    a = 1.d0 / sixtydz
    b = 9.d0 / sixtydz
    c = 45.d0 / sixtydz

    ! Boundary handling for periodic conditions
    df(:, :, 1) = a * (f(:, :, 4) - f(:, :, nz-2)) - &
         b * (f(:, :, 3) - f(:, :, nz-1)) + &
         c * (f(:, :, 2) - f(:, :, nz))
    df(:, :, 2) = a * (f(:, :, 5) - f(:, :, nz-1)) - &
         b * (f(:, :, 4) - f(:, :, nz)) + &
         c * (f(:, :, 3) - f(:, :, 1   ))
    df(:, :, 3) = a * (f(:, :, 6) - f(:, :, nz)) - &
         b * (f(:, :, 5) - f(:, :, 1   )) + &
         c * (f(:, :, 4) - f(:, :, 2   ))
    do k = 4, nz-3
       df(:, :, k) = a * (f(:, :, k+3) - f(:, :, k-3)) - &
            b * (f(:, :, k+2) - f(:, :, k-2)) + &
            c * (f(:, :, k+1) - f(:, :, k-1))
    end do
    df(:, :, nz-2) = a * (f(:, :, 1) - f(:, :, nz-5)) - &
         b * (f(:, :, nz) - f(:, :, nz-4)) + &
         c * (f(:, :, nz-1) - f(:, :, nz-3))
    df(:, :, nz-1) = a * (f(:, :, 2) - f(:, :, nz-4)) - &
         b * (f(:, :, 1) - f(:, :, nz-3)) + &
         c * (f(:, :, nz) - f(:, :, nz-2))
    df(:, :, nz) = a * (f(:, :, 3) - f(:, :, nz-3)) - &
         b * (f(:, :, 2) - f(:, :, nz-2)) + &
         c * (f(:, :, 1) - f(:, :, nz-1))

    return

  end subroutine derz_00

  subroutine derzp_11(df, f, dz)

    !> Z first derivation of even function f with O(dz6) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dz            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Z derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydz
    real(kind=8) :: a, b, c
    integer :: nz, k

    nz = size(f, 3)

    ! Coefficients for derivation
    sixtydz = 60.d0 * dz
    a = 1.d0 / sixtydz
    b = 9.d0 / sixtydz
    c = 45.d0 / sixtydz

    ! Boundary handling for free-slip conditions
    df(:, :, 1) = 0.d0
    df(:, :, 2) = a * (f(:, :, 5) - f(:, :, 3)) - &
         b * (f(:, :, 4) - f(:, :, 2)) + &
         c * (f(:, :, 3) - f(:, :, 1))
    df(:, :, 3) = a * (f(:, :, 6) - f(:, :, 2)) - &
         b * (f(:, :, 5) - f(:, :, 1)) + &
         c * (f(:, :, 4) - f(:, :, 2))
    do k = 4, nz-3
       df(:, :, k) = a * (f(:, :, k+3) - f(:, :, k-3)) - &
            b * (f(:, :, k+2) - f(:, :, k-2)) + &
            c * (f(:, :, k+1) - f(:, :, k-1))
    end do
    df(:, :, nz-2) = a * (f(:, :, nz-1) - f(:, :, nz-5)) - &
         b * (f(:, :, nz) - f(:, :, nz-4)) + &
         c * (f(:, :, nz-1) - f(:, :, nz-3))
    df(:, :, nz-1) = a * (f(:, :, nz-2) - f(:, :, nz-4)) - &
         b * (f(:, :, nz-1) - f(:, :, nz-3)) + &
         c * (f(:, :, nz) - f(:, :, nz-2))
    df(:, :, nz) = 0.d0

    return

  end subroutine derzp_11

  subroutine derzi_11(df, f, dz)

    !> Z first derivative of even function f with O(dz6) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dz            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Z derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: sixtydz
    real(kind=8) :: a, b, c
    integer :: nz, k

    nz = size(f, 3)

    ! Coefficients for derivation
    sixtydz = 60.d0 * dz
    a = 1.d0 / sixtydz
    b = 9.d0 / sixtydz
    c = 45.d0 / sixtydz

    ! Boundary handling for free-slip conditions
    df(:, :, 1) = a * (f(:, :, 4) + f(:, :, 4)) - &
         b * (f(:, :, 3) + f(:, :, 3)) + &
         c * (f(:, :, 2) + f(:, :, 2))
    df(:, :, 2) = a * (f(:, :, 5) + f(:, :, 3)) - &
         b * (f(:, :, 4) + f(:, :, 2)) + &
         c * (f(:, :, 3) - f(:, :, 1))
    df(:, :, 3) = a * (f(:, :, 6) + f(:, :, 2)) - &
         b * (f(:, :, 5) - f(:, :, 1)) + &
         c * (f(:, :, 4) - f(:, :, 2))
    do k = 4, nz-3
       df(:, :, k) = a * (f(:, :, k+3) - f(:, :, k-3)) - &
            b * (f(:, :, k+2) - f(:, :, k-2)) + &
            c * (f(:, :, k+1) - f(:, :, k-1))
    end do
    df(:, :, nz-2) = a * (-f(:, :, nz-1) - f(:, :, nz-5)) - &
         b * (f(:, :, nz) - f(:, :, nz-4)) + &
         c * (f(:, :, nz-1) - f(:, :, nz-3))
    df(:, :, nz-1) = a * (-f(:, :, nz-2) - f(:, :, nz-4)) - &
         b * (-f(:, :, nz-1) - f(:, :, nz-3)) + &
         c * (f(:, :, nz) - f(:, :, nz-2))
    df(:, :, nz) = a * (-f(:, :, nz-3) - f(:, :, nz-3)) - &
         b * (-f(:, :, nz-2) - f(:, :, nz-2)) + &
         c * (-f(:, :, nz-1) - f(:, :, nz-1))

    return

  end subroutine derzi_11

  subroutine derz_2dsim(df, f, dz)
    ! Z first derivation of function f for a 2D simulation
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dz            : mesh size step
    ! OUTPUT: df(nx, ny, nz): Z derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in):: dz
    real(kind=8), intent(out) :: df(:,:,:)

    df = 0.d0 * f * dz

    return
  end subroutine derz_2dsim

  subroutine derxx_00(df, f, dx)

    ! X second derivation of function f with O(dx4) scheme
    ! Case for Periodic boundary conditions
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dx            : mesh size step
    ! OUTPUT: df(nx, ny, nz): X derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedxsquare
    real(kind=8) :: a, b, c
    integer :: nx, i

    nx = size(f, 1)

    ! Coefficients for derivation
    twelvedxsquare = 12.d0 * dx * dx

    a =  1.d0 / twelvedxsquare
    b = 16.d0 / twelvedxsquare
    c = 30.d0 / twelvedxsquare

    df(1, :, :) = -a * (f(nx-1, :, :) + f(3, :, :)) + &
         b * (f(nx, :, :) + f(2, :, :)) - &
         c * f(1, :, :)
    df(2, :, :) = -a * (f(nx, :, :) + f(4, :, :)) + &
         b * (f(1, :, :) + f(3, :, :)) - &
         c * f(2, :, :)
    do i = 3, nx - 2
       df(i, :, :) = -a * (f(i-2, :, :) + f(i+2, :, :)) + &
            b * (f(i-1, :, :) + f(i+1, :, :)) - &
            c * f(i, :, :)
    end do
    df(nx-1, :, :) = -a * (f(nx-3, :, :) + f(1, :, :)) + &
         b * (f(nx-2, :, :) + f(nx, :, :)) - &
         c * f(nx-1, :, :)
    df(nx, :, :) = -a * (f(nx-2, :, :) + f(2, :, :)) + &
         b * (f(nx-1, :, :) + f(1, :, :)) - &
         c * f(nx, :, :)

    return
  end subroutine derxx_00

  subroutine derxxi_11(df, f, dx)

    !> X second derivative of even function f with O(dx4) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dx            : mesh size step
    !> OUTPUT: df(nx, ny, nz): X second derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedxsquare
    real(kind=8) :: a, b, c
    integer :: nx, i

    nx = size(f, 1)

    ! Coefficients for derivation
    twelvedxsquare = 12.d0 * dx * dx

    a =  1.d0 / twelvedxsquare
    b = 16.d0 / twelvedxsquare
    c = 30.d0 / twelvedxsquare

    ! Boundary handling for free-slip conditions
    df(1, :, :) = -a * (-f(3, :, :) + f(3, :, :)) + &
         b * (-f(2, :, :) + f(2, :, :)) - &
         c * f(1, :, :)
    df(2, :, :) = -a * (-f(2, :, :) + f(4, :, :)) + &
         b * (f(1, :, :) + f(3, :, :)) - &
         c * f(2, :, :)
    do i = 3, nx - 2
       df(i, :, :) = -a * (f(i-2, :, :) + f(i+2, :, :)) + &
            b * (f(i-1, :, :) + f(i+1, :, :)) - &
            c * f(i, :, :)
    end do
    df(nx-1, :, :) = -a * (f(nx-3, :, :) - f(nx-1, :, :)) + &
         b * (f(nx-2, :, :) + f(nx, :, :)) - &
         c * f(nx-1, :, :)
    df(nx, :, :) = -a * (f(nx-2, :, :) - f(nx-2, :, :)) + &
         b * (f(nx-1, :, :) - f(nx-1, :, :)) - &
         c * f(nx, :, :)

    return

  end subroutine derxxi_11

  subroutine derxxp_11(df, f, dx)

    !> X second derivative of odd function f with O(dx4) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dx            : mesh size step
    !> OUTPUT: df(nx, ny, nz): X second derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedxsquare
    real(kind=8) :: a, b, c
    integer :: nx, i

    nx = size(f, 1)

    ! Coefficients for derivation
    twelvedxsquare = 12.d0 * dx * dx

    a =  1.d0 / twelvedxsquare
    b = 16.d0 / twelvedxsquare
    c = 30.d0 / twelvedxsquare

    ! Boundary handling for free-slip conditions
    df(1, :, :) = -a * (f(3, :, :) + f(3, :, :)) + &
         b * (f(2, :, :) + f(2, :, :)) - &
         c * f(1, :, :)
    df(2, :, :) = -a * (f(2, :, :) + f(4, :, :)) + &
         b * (f(1, :, :) + f(3, :, :)) - &
         c * f(2, :, :)
    do i = 3, nx - 2
       df(i, :, :) = -a * (f(i-2, :, :) + f(i+2, :, :)) + &
            b * (f(i-1, :, :) + f(i+1, :, :)) - &
            c * f(i, :, :)
    end do
    df(nx-1, :, :) = -a * (f(nx-3, :, :) + f(nx-1, :, :)) + &
         b * (f(nx-2, :, :) + f(nx, :, :)) - &
         c * f(nx-1, :, :)
    df(nx, :, :) = -a * (f(nx-2, :, :) + f(nx-2, :, :)) + &
         b * (f(nx-1, :, :) + f(nx-1, :, :)) - &
         c * f(nx, :, :)

    return

  end subroutine derxxp_11

  subroutine deryy_00(df, f, dy)

    ! Y second derivation of function f with O(dy4) scheme
    ! Case for Periodic boundary conditions
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dy            : mesh size step
    !         odd           : true if f is odd
    ! OUTPUT: df(nx, ny, nz): Y derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedysquare
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f, 2)

    ! Coefficients for derivation
    twelvedysquare = 12.d0 * dy * dy

    a =  1.d0 / twelvedysquare
    b = 16.d0 / twelvedysquare
    c = 30.d0 / twelvedysquare

    df(:, 1, :) = -a * (f(:, ny-1, :) + f(:, 3, :)) + &
         b * (f(:, ny, :) + f(:, 2, :)) - &
         c * f(:, 1, :)
    df(:, 2, :) = -a * (f(:, ny, :) + f(:, 4, :)) + &
         b * (f(:, 1, :) + f(:, 3, :)) - &
         c * f(:, 2, :)
    do j = 3, ny - 2
       df(:, j, :) = -a * (f(:, j-2, :) + f(:, j+2, :)) + &
            b * (f(:, j-1, :) + f(:, j+1, :)) - &
            c * f(:, j, :)
    end do
    df(:, ny-1, :) = -a * (f(:, ny-3, :) + f(:, 1, :)) + &
         b * (f(:, ny-2, :) + f(:, ny, :)) - &
         c * f(:, ny-1, :)
    df(:, ny, :) = -a * (f(:, ny-2, :) + f(:, 2, :)) + &
         b * (f(:, ny-1, :) + f(:, 1, :)) - &
         c * f(:, ny, :)

    return

  end subroutine deryy_00

  subroutine deryyp_11(df, f, dy)

    ! Y second derivation of odd function f with O(dy4) scheme
    ! Case for Free-slip boundary conditions
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dy            : mesh size step
    !         odd           : true if f is odd
    ! OUTPUT: df(nx, ny, nz): Y derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedysquare
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f, 2)

    ! Coefficients for derivation
    twelvedysquare = 12.d0 * dy * dy

    a =  1.d0 / twelvedysquare
    b = 16.d0 / twelvedysquare
    c = 30.d0 / twelvedysquare

    df(:, 1, :) = -a * (f(:, 3, :) + f(:, 3, :)) + &
         b * (f(:, 2, :) + f(:, 2, :)) - &
         c * f(:, 1, :)
    df(:, 2, :) = -a * (f(:, 2, :) + f(:, 4, :)) + &
         b * (f(:, 1, :) + f(:, 3, :)) - &
         c * f(:, 2, :)
    do j = 3, ny - 2
       df(:, j, :) = -a * (f(:, j-2, :) + f(:, j+2, :)) + &
            b * (f(:, j-1, :) + f(:, j+1, :)) - &
            c * f(:, j, :)
    end do
    df(:, ny-1, :) = -a * (f(:, ny-3, :) + f(:, ny-1, :)) + &
         b * (f(:, ny-2, :) + f(:, ny, :)) - &
         c * f(:, ny-1, :)
    df(:, ny, :) = -a * (f(:, ny-2, :) + f(:, ny-2, :)) + &
         b * (f(:, ny-1, :) + f(:, ny-1, :)) - &
         c * f(:, ny, :)

    return

  end subroutine deryyp_11

  subroutine deryyi_11(df, f, dy)

    ! Y second derivation of even function f with O(dy4) scheme
    ! Case for Free-slip boundary conditions
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dy            : mesh size step
    !         odd           : true if f is odd
    ! OUTPUT: df(nx, ny, nz): Y derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedysquare
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f, 2)

    ! Coefficients for derivation
    twelvedysquare = 12.d0 * dy * dy

    a =  1.d0 / twelvedysquare
    b = 16.d0 / twelvedysquare
    c = 30.d0 / twelvedysquare

    df(:, 1, :) = -a * (-f(:, 3, :) + f(:, 3, :)) + &
         b * (-f(:, 2, :) + f(:, 2, :)) - &
         c * f(:, 1, :)
    df(:, 2, :) = -a * (-f(:, 2, :) + f(:, 4, :)) + &
         b * (f(:, 1, :) + f(:, 3, :)) - &
         c * f(:, 2, :)
    do j = 3, ny - 2
       df(:, j, :) = -a * (f(:, j-2, :) + f(:, j+2, :)) + &
            b * (f(:, j-1, :) + f(:, j+1, :)) - &
            c * f(:, j, :)
    end do
    df(:, ny-1, :) = -a * (f(:, ny-3, :) - f(:, ny-1, :)) + &
         b * (f(:, ny-2, :) + f(:, ny, :)) - &
         c * f(:, ny-1, :)
    df(:, ny, :) = -a * (f(:, ny-2, :) - f(:, ny-2, :)) + &
         b * (f(:, ny-1, :) - f(:, ny-1, :)) - &
         c * f(:, ny, :)

    return

  end subroutine deryyi_11

  subroutine derzz_00(df, f, dz)

    ! Z second derivation of function f with O(dz4) scheme
    ! Case for Periodic boundary conditions
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dz            : mesh size step
    ! OUTPUT: df(nx, ny, nz): Z derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedzsquare
    real(kind=8) :: a, b, c
    integer :: nz, k

    nz = size(f, 3)

    ! Coefficients for derivation
    twelvedzsquare = 12.d0 * dz * dz

    a =  1.d0 / twelvedzsquare
    b = 16.d0 / twelvedzsquare
    c = 30.d0 / twelvedzsquare

    df(:, :, 1) = -a * (f(:, :, nz-1) + f(:, :, 3)) + &
         b * (f(:, :, nz) + f(:, :, 2)) - &
         c * f(:, :, 1)
    df(:, :, 2) = -a * (f(:, :, nz) + f(:, :, 4)) + &
         b * (f(:, :, 1) + f(:, :, 3)) - &
         c * f(:, :, 2)
    do k = 3, nz - 2
       df(:, :, k) = -a * (f(:, :, k-2) + f(:, :, k+2)) + &
            b * (f(:, :, k-1) + f(:, :, k+1)) - &
            c * f(:, :, k)
    end do
    df(:, :, nz-1) = -a * (f(:, :, nz-3) + f(:, :, 1)) + &
         b * (f(:, :, nz-2) + f(:, :, nz)) - &
         c * f(:, :, nz-1)
    df(:, :, nz) = -a * (f(:, :, nz-2) + f(:, :, 2)) + &
         b * (f(:, :, nz-1) + f(:, :, 1)) - &
         c * f(:, :, nz)
    return

  end subroutine derzz_00

  subroutine derzzi_11(df, f, dz)

    !> Z second derivative of even function f with O(dz4) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dz            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Z second derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedzsquare
    real(kind=8) :: a, b, c
    integer :: nz, k

    nz = size(f, 3)

    ! Coefficients for derivation
    twelvedzsquare = 12.d0 * dz * dz

    a =  1.d0 / twelvedzsquare
    b = 16.d0 / twelvedzsquare
    c = 30.d0 / twelvedzsquare

    ! Boundary handling for free-slip conditions
    df(:, :, 1) = -a * (-f(:, :, 3) + f(:, :, 3)) + &
         b * (-f(:, :, 2) + f(:, :, 2)) - &
         c * f(:, :, 1)
    df(:, :, 2) = -a * (-f(:, :, 2) + f(:, :, 4)) + &
         b * (f(:, :, 1) + f(:, :, 3)) - &
         c * f(:, :, 2)
    do k = 3, nz - 2
       df(:, :, k) = -a * (f(:, :, k-2) + f(:, :, k+2)) + &
            b * (f(:, :, k-1) + f(:, :, k+1)) - &
            c * f(:, :, k)
    end do
    df(:, :, nz-1) = -a * (f(:, :, nz-3) - f(:, :, nz-1)) + &
         b * (f(:, :, nz-2) + f(:, :, nz)) - &
         c * f(:, :, nz-1)
    df(:, :, nz) = -a * (f(:, :, nz-2) - f(:, :, nz-2)) + &
         b * (f(:, :, nz-1) - f(:, :, nz-1)) - &
         c * f(:, :, nz)

    return

  end subroutine derzzi_11

  subroutine derzzp_11(df, f, dz)

    !> Z second derivative of odd function f with O(dz4) scheme
    !> Case for Free-slip boundary conditions
    !>
    !> INPUT:  f(nx, ny, nz) : flow variable
    !>         dz            : mesh size step
    !> OUTPUT: df(nx, ny, nz): Z second derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    real(kind=8) :: twelvedzsquare
    real(kind=8) :: a, b, c
    integer :: nz, k

    nz = size(f, 3)

    ! Coefficients for derivation
    twelvedzsquare = 12.d0 * dz * dz

    a =  1.d0 / twelvedzsquare
    b = 16.d0 / twelvedzsquare
    c = 30.d0 / twelvedzsquare

    ! Boundary handling for free-slip conditions
    df(:, :, 1) = -a * (f(:, :, 3) + f(:, :, 3)) + &
         b * (f(:, :, 2) + f(:, :, 2)) - &
         c * f(:, :, 1)
    df(:, :, 2) = -a * (f(:, :, 2) + f(:, :, 4)) + &
         b * (f(:, :, 1) + f(:, :, 3)) - &
         c * f(:, :, 2)
    do k = 3, nz - 2
       df(:, :, k) = -a * (f(:, :, k-2) + f(:, :, k+2)) + &
            b * (f(:, :, k-1) + f(:, :, k+1)) - &
            c * f(:, :, k)
    end do
    df(:, :, nz-1) = -a * (f(:, :, nz-3) + f(:, :, nz-1)) + &
         b * (f(:, :, nz-2) + f(:, :, nz)) - &
         c * f(:, :, nz-1)
    df(:, :, nz) = -a * (f(:, :, nz-2) + f(:, :, nz-2)) + &
         b * (f(:, :, nz-1) + f(:, :, nz-1)) - &
         c * f(:, :, nz)

    return

  end subroutine derzzp_11

  subroutine derzz_2dsim(df, f, dz)
    ! Z first derivation of function f for a 2D simulation
    !
    ! INPUT:  f(nx, ny, nz) : flow variable
    !         dz            : mesh size step
    ! OUTPUT: df(nx, ny, nz): Z derivative of the flow variable

    real(kind=8), intent(in) :: f(:,:,:)
    real(kind=8), intent(in) :: dz
    real(kind=8), intent(out) :: df(:,:,:)

    df = 0.d0 * f * dz

    return
  end subroutine derzz_2dsim

  subroutine dery1D(df, f, dy)
    !> Y first derivation of function f with O(dy6) scheme
    !> Case for modified boundary conditions
    !>
    !> INPUT:  f(ny)  : flow variable
    !>         dy     : mesh size step
    !> OUTPUT: df(ny) : Y derivative of the f function

    real(kind=8), intent(in) :: f(:)
    real(kind=8), intent(in) :: dy
    real(kind=8), intent(out) :: df(:)

    real(kind=8) :: sixtydy
    real(kind=8) :: a, b, c
    integer :: ny, j

    ny = size(f)

    ! Coefficients for derivation
    sixtydy = 60.d0 * dy
    a = 1.d0 / sixtydy
    b = 9.d0 / sixtydy
    c = 45.d0 / sixtydy

    ! Boundary handling for modified conditions
    df(1) = (-f(3) + 4.d0 * f(2) - 3.d0 * f(1)) / (2.d0 * dy)
    df(2) = (f(3) - f(1)) / (2.d0 * dy)
    df(3) = (-f(5) + 8.d0 * f(4) - 8.d0 * f(2) + f(1)) / (12.d0 * dy)

    do j = 4, ny-3
       df(j) = a * (f(j+3) - f(j-3)) - &
            b * (f(j+2) - f(j-2)) + &
            c * (f(j+1) - f(j-1))
    end do

    ! Boundary handling for modified conditions
    df(ny-2) = (-f(ny) + 8.d0 * f(ny-1) - 8.d0 * f(ny-3) + f(ny-4)) &
         / (12.d0 * dy)
    df(ny-1) = (f(ny) - f(ny-2)) / (2.d0 * dy)
    df(ny  ) = (3.d0 * f(ny) - 4.d0 * f(ny-1) + f(ny-2)) / (2.d0 * dy)

    return
  end subroutine dery1D

end module derivation
