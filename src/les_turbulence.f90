module les_turbulence
  use IOfunctions
  use derivation
  use initialization
  implicit none

contains

  subroutine calculate_nu_t(nu_t, ux, uy, uz, dx, dy, dz, cs, delta)
    !> Calculates the sub-grid scale turbulent viscosity (nu_t) using the 
    !> Smagorinsky model.
    !>
    !> INPUT:
    !> ux(nx, ny, nz)    : Current x-component of velocity
    !> uy(nx, ny, nz)    : Current y-component of velocity
    !> uz(nx, ny, nz)    : Current z-component of velocity
    !> dx                : Grid spacing in the x-direction
    !> dy                : Grid spacing in the y-direction
    !> dz                : Grid spacing in the z-direction
    !> cs                : Smagorinsky constant
    !> delta             : filter size
    !>
    !> OUTPUT:
    !> nu_t(nx, ny, nz)  : Calculated sub-grid scale turbulent viscosity

    ! Define input and output variables
    real(kind=8), dimension(:,:,:), intent(out) :: nu_t
    real(kind=8), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(kind=8), intent(in) :: dx, dy, dz, cs, delta

    integer :: nx, ny, nz
    real(kind=8), allocatable, dimension(:,:,:) :: duxdx, duxdy, duxdz
    real(kind=8), allocatable, dimension(:,:,:) :: duydx, duydy, duydz
    real(kind=8), allocatable, dimension(:,:,:) :: duzdx, duzdy, duzdz
    real(kind=8), allocatable, dimension(:,:,:) :: sij
    real(kind=8) :: stats(6)

    ! Grid dimensions
    nx = size(ux, 1)
    ny = size(ux, 2)
    nz = size(ux, 3)

    ! Allocate arrays for partial derivatives
    allocate(duxdx(nx,ny,nz), duxdy(nx,ny,nz), duxdz(nx,ny,nz))
    allocate(duydx(nx,ny,nz), duydy(nx,ny,nz), duydz(nx,ny,nz))
    allocate(duzdx(nx,ny,nz), duzdy(nx,ny,nz), duzdz(nx,ny,nz))
    allocate(sij(nx,ny,nz))

    ! Compute the partial derivatives of ux
    call derxi(duxdx, ux, dx)
    call deryp(duxdy, ux, dy)
    call derzp(duxdz, ux, dz)

    ! Compute the partial derivatives of uy
    call derxp(duydx, uy, dx)
    call deryi(duydy, uy, dy)
    call derzp(duydz, uy, dz)

    ! Compute the partial derivatives of uz
    call derxp(duzdx, uz, dx)
    call deryp(duzdy, uz, dy)
    call derzi(duzdz, uz, dz)

    ! Calculate the deformation rate tensor Sij
    sij = sqrt(2.0 * (duxdx**2 + duydy**2 + duzdz**2 + &
         0.5 * ((duxdy + duydx)**2 + (duxdz + duzdx)**2 + (duydz + duzdy)**2)))

    ! Compute the sub-grid scale turbulent viscosity nu_t
    nu_t = (cs * delta)**2.d0 * sij
    stats = function_stats(nu_t, nx, ny, nz)
    call print_nu_t_statistics(stats)

    ! Deallocate arrays
    deallocate(duxdx, duxdy, duxdz)
    deallocate(duydx, duydy, duydz)
    deallocate(duzdx, duzdy, duzdz)
    deallocate(sij)

    return
  end subroutine calculate_nu_t

  subroutine calculate_tau_ij(tau_ij, ux, uy, uz, dx, dy, dz, nu_t)
    !> Calculates the sub-grid scale stress tensor (tau_ij) for a 
    !> Large Eddy Simulation (LES).
    !>
    !> This subroutine computes the components of the sub-grid stress tensor 
    !> using the filtered velocity field and the sub-grid scale turbulent 
    !> viscosity (nu_t). 
    !> The stress tensor represents the effects of unresolved small-scale 
    !> motions on the resolved velocity field.
    !>
    !> INPUT: 
    !> ux(nx, ny, nz)    : x-component of the velocity field
    !> uy(nx, ny, nz)    : y-component of the velocity field
    !> uz(nx, ny, nz)    : z-component of the velocity field
    !> nu_t(nx, ny, nz)  : Sub-grid scale turbulent viscosity
    !> dx                : Grid spacing in the x-direction
    !> dy                : Grid spacing in the y-direction
    !> dz                : Grid spacing in the z-direction
    !>
    !> OUTPUT: 
    !> tau_ij(nx, ny, nz, 6) : Sub-grid scale stress tensor, where:
    !>     tau_ij(:,:,:,1) = tau_xx (normal stress in x-direction)
    !>     tau_ij(:,:,:,2) = tau_yy (normal stress in y-direction)
    !>     tau_ij(:,:,:,3) = tau_zz (normal stress in z-direction)
    !>     tau_ij(:,:,:,4) = tau_xy (shear stress between x and y)
    !>     tau_ij(:,:,:,5) = tau_xz (shear stress between x and z)
    !>     tau_ij(:,:,:,6) = tau_yz (shear stress between y and z)
    !>
    !> The stress tensor is symmetric: 
    !> tau_xy = tau_yx, tau_xz = tau_zx, tau_yz = tau_zy.

    real(kind=8), dimension(:,:,:,:), intent(out) :: tau_ij
    real(kind=8), dimension(:,:,:), intent(in) :: ux, uy, uz, nu_t
    real(kind=8), intent(in) :: dx, dy, dz

    integer :: nx, ny, nz
    real(kind=8), allocatable, dimension(:,:,:) :: duxdx, duxdy, duxdz
    real(kind=8), allocatable, dimension(:,:,:) :: duydx, duydy, duydz
    real(kind=8), allocatable, dimension(:,:,:) :: duzdx, duzdy, duzdz

    ! Grid dimensions
    nx = size(ux, 1)
    ny = size(ux, 2)
    nz = size(ux, 3)

    ! Allocate arrays for partial derivatives
    allocate(duxdx(nx,ny,nz), duxdy(nx,ny,nz), duxdz(nx,ny,nz))
    allocate(duydx(nx,ny,nz), duydy(nx,ny,nz), duydz(nx,ny,nz))
    allocate(duzdx(nx,ny,nz), duzdy(nx,ny,nz), duzdz(nx,ny,nz))

    ! Compute the partial derivatives of ux
    call derxi(duxdx, ux, dx)
    call deryp(duxdy, ux, dy)
    call derzp(duxdz, ux, dz)

    ! Compute the partial derivatives of uy
    call derxp(duydx, uy, dx)
    call deryi(duydy, uy, dy)
    call derzp(duydz, uy, dz)

    ! Compute the partial derivatives of uz
    call derxp(duzdx, uz, dx)
    call deryp(duzdy, uz, dy)
    call derzi(duzdz, uz, dz)

    tau_ij(:,:,:,1) = -2.d0 * nu_t(:,:,:) * duxdx(:,:,:)    !> tau_xx
    tau_ij(:,:,:,2) = -2.d0 * nu_t(:,:,:) * duydy(:,:,:)    !> tau_yy
    tau_ij(:,:,:,3) = -2.d0 * nu_t(:,:,:) * duzdz(:,:,:)    !> tau_zz
    tau_ij(:,:,:,4) = -nu_t(:,:,:) * (duxdy(:,:,:) + duydx) !> tau_xy
    tau_ij(:,:,:,5) = -nu_t(:,:,:) * (duxdz(:,:,:) + duzdx) !> tau_xz
    tau_ij(:,:,:,6) = -nu_t(:,:,:) * (duydz(:,:,:) + duzdy) !> tau_yz

    ! Deallocate arrays
    deallocate(duxdx, duxdy, duxdz)
    deallocate(duydx, duydy, duydz)
    deallocate(duzdx, duzdy, duzdz)

    return
  end subroutine calculate_tau_ij

  subroutine calculate_dtau_ij_dxj(dtau_dx, dtau_dy, dtau_dz, &
       tau_ij, dx, dy, dz)
    !> Calculates the divergence of the sub-grid scale stress tensor for a 
    !> Large Eddy Simulation (LES).
    !>
    !> This subroutine computes the partial derivatives of the sub-grid scale 
    !> stress tensor (tau_ij)
    !> with respect to the spatial coordinates x, y, and z. These derivatives 
    !> represent the contribution
    !> of the sub-grid scale stresses to the momentum equation.
    !>
    !> INPUT:
    !> tau_ij(nx, ny, nz, 6) : Sub-grid scale stress tensor, where:
    !>     tau_ij(:,:,:,1) = tau_xx (normal stress in x-direction)
    !>     tau_ij(:,:,:,2) = tau_yy (normal stress in y-direction)
    !>     tau_ij(:,:,:,3) = tau_zz (normal stress in z-direction)
    !>     tau_ij(:,:,:,4) = tau_xy (shear stress between x and y)
    !>     tau_ij(:,:,:,5) = tau_xz (shear stress between x and z)
    !>     tau_ij(:,:,:,6) = tau_yz (shear stress between y and z)
    !> dx                 : Grid spacing in the x-direction
    !> dy                 : Grid spacing in the y-direction
    !> dz                 : Grid spacing in the z-direction
    !>
    !> OUTPUT:
    !> dtau_dx(nx, ny, nz) : 
    !>          Divergence of tau_ij in the x-direction (d(tau_ij)/dx_j)
    !> dtau_dy(nx, ny, nz) : 
    !>          Divergence of tau_ij in the y-direction (d(tau_ij)/dy_j)
    !> dtau_dz(nx, ny, nz) : 
    !>          Divergence of tau_ij in the z-direction (d(tau_ij)/dz_j)
    !>
    !> The divergence components are computed for each of the six terms in the 
    !> stress tensor using central finite differences.

    implicit none
    real(kind=8), dimension(:,:,:,:), intent(in) :: tau_ij
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(:,:,:), intent(out) :: dtau_dx, dtau_dy, dtau_dz

    integer :: nx, ny, nz
    real(kind=8), allocatable, dimension(:,:,:) :: &
         dtauxx_dx, dtauxy_dx, dtauxz_dx
    real(kind=8), allocatable, dimension(:,:,:) :: &
         dtauxx_dy, dtauxy_dy, dtauxz_dy
    real(kind=8), allocatable, dimension(:,:,:) :: &
         dtauxx_dz, dtauxy_dz, dtauxz_dz

    ! Get grid dimensions
    nx = size(tau_ij, 1)
    ny = size(tau_ij, 2)
    nz = size(tau_ij, 3)

    ! Allocate arrays for partial derivatives
    allocate(dtauxx_dx(nx,ny,nz), dtauxy_dx(nx,ny,nz), dtauxz_dx(nx,ny,nz))
    allocate(dtauxx_dy(nx,ny,nz), dtauxy_dy(nx,ny,nz), dtauxz_dy(nx,ny,nz))
    allocate(dtauxx_dz(nx,ny,nz), dtauxy_dz(nx,ny,nz), dtauxz_dz(nx,ny,nz))

    ! Compute derivatives in the x-direction
    call derxp(dtauxx_dx, tau_ij(:,:,:,1), dx)
    call derxp(dtauxy_dx, tau_ij(:,:,:,4), dx)
    call derxp(dtauxz_dx, tau_ij(:,:,:,5), dx)

    ! Compute derivatives in the y-direction
    call deryp(dtauxx_dy, tau_ij(:,:,:,4), dy)
    call deryp(dtauxy_dy, tau_ij(:,:,:,2), dy)
    call deryp(dtauxz_dy, tau_ij(:,:,:,6), dy)

    ! Compute derivatives in the z-direction
    call derzp(dtauxx_dz, tau_ij(:,:,:,5), dz)
    call derzp(dtauxy_dz, tau_ij(:,:,:,6), dz)
    call derzp(dtauxz_dz, tau_ij(:,:,:,3), dz)

    ! Sum up contributions for each component
    dtau_dx = dtauxx_dx + dtauxy_dx + dtauxz_dx  ! d(tau_ij)/dxj for x-component
    dtau_dy = dtauxx_dy + dtauxy_dy + dtauxz_dy  ! d(tau_ij)/dyj for y-component
    dtau_dz = dtauxx_dz + dtauxy_dz + dtauxz_dz  ! d(tau_ij)/dzj for z-component

    ! Deallocate arrays
    deallocate(dtauxx_dx, dtauxy_dx, dtauxz_dx)
    deallocate(dtauxx_dy, dtauxy_dy, dtauxz_dy)
    deallocate(dtauxx_dz, dtauxy_dz, dtauxz_dz)

    return
  end subroutine calculate_dtau_ij_dxj

end module les_turbulence
