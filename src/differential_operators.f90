module diffoper
  use initialization
  use derivation
  implicit none

contains
  subroutine divergence(divf, fx, fy, fz, dx, dy, dz, nx, ny, nz, odd)

    !> INPUT:
    !> fx, fy, fz: vector field composantes
    !> dx, dy, dz: mesh step size
    !> nx, ny, nz: number of grid points
    !> odd       : odd = 1 -> fy odd, odd = 0 -> fy even
    !>
    !> OUTPUT:
    !> divf      : divergence of the vector field (fx, fy, fz)

    integer, intent(in) :: nx, ny, nz, odd
    real(kind=8), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), intent(out) :: divf(:,:,:)

    real(kind=8), dimension(nx,ny,nz) :: dfx, dfy, dfz

    if (odd == 0) then
       call derxp(dfx, fx, dx)
       call deryp(dfy, fy, dy)
       call derzp(dfz, fz, dz)
    else
       call derxi(dfx, fx, dx)
       call deryi(dfy, fy, dy)
       call derzi(dfz, fz, dz)
    end if

    divf = dfx + dfy + dfz

    return
  end subroutine divergence

  subroutine rotational(rotx, roty, rotz, ux, uy, uz, dx, dy, dz, nx, ny, nz)
    !> Calculate the rotational of a vector field U = (ux, uy, uz)
    !> on a 3D grid with resolutions dx, dy, dz.
    !
    !> INPUT:
    !> ux(nx, ny, nz) : x-component of velocity
    !> uy(nx, ny, nz) : y-component of velocity
    !> uz(nx, ny, nz) : z-component of velocity
    !> dx, dy, dz     : grid spacing in each direction
    !> nx, ny, nz     : number of grid points in each direction
    !
    !> OUTPUT:
    !> rotx(nx, ny, nz) : x-component of the curl
    !> roty(nx, ny, nz) : y-component of the curl
    !> rotz(nx, ny, nz) : z-component of the curl

    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), intent(out) :: rotx(:,:,:), roty(:,:,:), rotz(:,:,:)
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz) :: duzdy, duydz, duxdz
    real(kind=8), dimension(nx, ny, nz) :: duzdx, duydx, duxdy

    ! Calculate partial derivatives necessary for the curl
    call deryp(duzdy, uz, dy)  ! Derivative of uz with respect to y
    call derzp(duydz, uy, dz)  ! Derivative of uy with respect to z
    rotx(:,:,:) = duzdy - duydz

    call derzp(duxdz, ux, dz)  ! Derivative of ux with respect to z
    call derxp(duzdx, uz, dx)  ! Derivative of uz with respect to x
    roty(:,:,:) = duxdz - duzdx

    call derxp(duydx, uy, dx)  ! Derivative of uy with respect to x
    call deryp(duxdy, ux, dy)  ! Derivative of ux with respect to y
    rotz(:,:,:) = duydx - duxdy

    return
  end subroutine rotational

  subroutine calculate_Q_criterion(Q, ux, uy, uz, dx, dy, dz, nx, ny, nz)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), intent(out) :: Q(:,:,:)

    real(kind=8) :: duxdx(nx,ny,nz), duydy(nx,ny,nz), duzdz(nx,ny,nz)
    real(kind=8) :: duxdy(nx,ny,nz), duxdz(nx,ny,nz), duydx(nx,ny,nz)
    real(kind=8) :: duydz(nx,ny,nz), duzdx(nx,ny,nz), duzdy(nx,ny,nz)

    ! Compute derivatives
    call derxi(duxdx, ux, dx)
    call deryi(duydy, uy, dy)
    call derzi(duzdz, uz, dz)

    call deryp(duxdy, ux, dy)
    call derzp(duxdz, ux, dz)
    call derxp(duydx, uy, dx)

    call derzp(duydz, uy, dz)
    call derxp(duzdx, uz, dx)
    call deryp(duzdy, uz, dy)

    ! Calculate Q-criterion at each point
    Q = - 0.5d0 * (duxdx ** 2 + duydy ** 2 + duzdz ** 2) - &
         duxdy * duydx - duxdz * duzdx - duydz * duzdy

    return
  end subroutine calculate_Q_criterion

end module diffoper
