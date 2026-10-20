module les_turbulence
  use IOfunctions
  use derivation
  use initialization
  use functions
  implicit none

  !> LES statistics and printouts (nu_t, tau_ij, dtau_ij/dxj) at every call.
  !> Set to .false. (or drive it from the main program with
  !> mod(itime, ioutput) == 0) to skip the global reductions and the I/O.
  logical :: les_verbose = .true.

  !> Work array (allocated once) for calculate_dtau_ij_dxj
  real(kind=8), allocatable, private :: les_scratch(:,:,:)

contains

  subroutine compute_velocity_gradients(g, ux, uy, uz, dx, dy, dz)
    !> Computes the 9 velocity gradients once, so that they can be shared
    !> between the LES model and the convective/viscous terms.
    !>
    !> OUTPUT: g(nx,ny,nz,9), with
    !>   g(:,:,:,1) = dux/dx   g(:,:,:,2) = dux/dy   g(:,:,:,3) = dux/dz
    !>   g(:,:,:,4) = duy/dx   g(:,:,:,5) = duy/dy   g(:,:,:,6) = duy/dz
    !>   g(:,:,:,7) = duz/dx   g(:,:,:,8) = duz/dy   g(:,:,:,9) = duz/dz

    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), intent(out) :: g(:,:,:,:)

    call derxi(g(:,:,:,1), ux, dx)
    call deryp(g(:,:,:,2), ux, dy)
    call derzp(g(:,:,:,3), ux, dz)

    call derxp(g(:,:,:,4), uy, dx)
    call deryi(g(:,:,:,5), uy, dy)
    call derzp(g(:,:,:,6), uy, dz)

    call derxp(g(:,:,:,7), uz, dx)
    call deryp(g(:,:,:,8), uz, dy)
    call derzi(g(:,:,:,9), uz, dz)

  end subroutine compute_velocity_gradients

  subroutine calculate_nu_t_from_grad(nu_t, g, cs, delta, nu)
    !> Smagorinsky model: nu_t = (cs*delta)**2 * |S|, computed from the
    !> already available gradients (see compute_velocity_gradients).
    !>
    !> INPUT:  g(nx,ny,nz,9) : velocity gradients
    !>         cs, delta     : Smagorinsky constant, filter size
    !>         nu            : molecular viscosity (used for statistics output)
    !> OUTPUT: nu_t(nx,ny,nz)

    real(kind=8), intent(in) :: g(:,:,:,:)
    real(kind=8), intent(in) :: cs, delta, nu
    real(kind=8), intent(out) :: nu_t(:,:,:)

    integer :: i, j, k, nx, ny, nz
    real(kind=8) :: s11, s22, s33, s12, s13, s23, s_mag, cs2
    real(kind=8) :: stats(6)

    nx = size(g, 1)
    ny = size(g, 2)
    nz = size(g, 3)

    cs2 = (cs * delta)**2

    !$omp parallel do collapse(2) default(shared) private(i, s11, s22, s33, s12, s13, s23, s_mag)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ! Symmetric strain rate tensor components
             s11 = g(i,j,k,1)
             s22 = g(i,j,k,5)
             s33 = g(i,j,k,9)
             s12 = 0.5d0 * (g(i,j,k,2) + g(i,j,k,4))
             s13 = 0.5d0 * (g(i,j,k,3) + g(i,j,k,7))
             s23 = 0.5d0 * (g(i,j,k,6) + g(i,j,k,8))

             ! |S| = sqrt(2 * S_ij * S_ij)
             s_mag = sqrt(2.0d0 * (s11**2 + s22**2 + s33**2 + &
                  2.0d0 * (s12**2 + s13**2 + s23**2)))

             nu_t(i,j,k) = cs2 * s_mag
          end do
       end do
    end do

    if (les_verbose) then
       stats = function_stats(nu_t, nx, ny, nz)
       call print_nu_t_statistics(stats, nu)
    end if

  end subroutine calculate_nu_t_from_grad

  subroutine calculate_tau_ij_from_grad(tau_ij, g, nu_t)
    !> Sub-grid scale stress tensor tau_ij from the already computed gradients.
    !>
    !> OUTPUT: tau_ij(nx,ny,nz,6) :
    !>     1 = tau_xx, 2 = tau_yy, 3 = tau_zz,
    !>     4 = tau_xy, 5 = tau_xz, 6 = tau_yz

    real(kind=8), intent(out) :: tau_ij(:,:,:,:)
    real(kind=8), intent(in) :: g(:,:,:,:)
    real(kind=8), intent(in) :: nu_t(:,:,:)

    integer :: i, j, k, nx, ny, nz
    real(kind=8) :: nut

    nx = size(g, 1)
    ny = size(g, 2)
    nz = size(g, 3)

    !$omp parallel do collapse(2) default(shared) private(i, nut)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             nut = nu_t(i,j,k)
             tau_ij(i,j,k,1) = -2.d0 * nut * g(i,j,k,1)                 ! tau_xx
             tau_ij(i,j,k,2) = -2.d0 * nut * g(i,j,k,5)                 ! tau_yy
             tau_ij(i,j,k,3) = -2.d0 * nut * g(i,j,k,9)                 ! tau_zz
             tau_ij(i,j,k,4) = -nut * (g(i,j,k,2) + g(i,j,k,4))         ! tau_xy
             tau_ij(i,j,k,5) = -nut * (g(i,j,k,3) + g(i,j,k,7))         ! tau_xz
             tau_ij(i,j,k,6) = -nut * (g(i,j,k,6) + g(i,j,k,8))         ! tau_yz
          end do
       end do
    end do

    if (les_verbose) call print_tau_ij_max_values(tau_ij)

  end subroutine calculate_tau_ij_from_grad

  subroutine calculate_nu_t(nu_t, ux, uy, uz, dx, dy, dz, cs, delta, nu)
    !> Calculates the sub-grid scale turbulent viscosity (nu_t) using the
    !> Smagorinsky model.
    !>
    !> Interface kept for backward compatibility. If the gradients are
    !> already available, prefer compute_velocity_gradients followed by
    !> calculate_nu_t_from_grad (avoids recomputing 9 derivatives).
    !>
    !> INPUT:
    !> ux, uy, uz (nx, ny, nz) : velocity components
    !> dx, dy, dz              : grid spacings
    !> cs                      : Smagorinsky constant
    !> delta                   : filter size
    !> nu                      : molecular viscosity
    !>
    !> OUTPUT:
    !> nu_t(nx, ny, nz)        : sub-grid scale turbulent viscosity

    real(kind=8), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(kind=8), intent(in) :: dx, dy, dz, cs, delta, nu
    real(kind=8), dimension(:,:,:), intent(out) :: nu_t

    real(kind=8), allocatable :: g(:,:,:,:)

    allocate(g(size(ux,1), size(ux,2), size(ux,3), 9))
    call compute_velocity_gradients(g, ux, uy, uz, dx, dy, dz)
    call calculate_nu_t_from_grad(nu_t, g, cs, delta, nu)
    deallocate(g)

  end subroutine calculate_nu_t

  subroutine calculate_tau_ij(tau_ij, ux, uy, uz, dx, dy, dz, nu_t)
    !> Calculates the sub-grid scale stress tensor (tau_ij) for a
    !> Large Eddy Simulation (LES).
    !>
    !> Interface kept for backward compatibility (see calculate_tau_ij_from_grad).
    !>
    !> INPUT:
    !> ux, uy, uz (nx, ny, nz) : velocity components
    !> nu_t(nx, ny, nz)        : sub-grid scale turbulent viscosity
    !> dx, dy, dz              : grid spacings
    !>
    !> OUTPUT:
    !> tau_ij(nx, ny, nz, 6) :
    !>     tau_ij(:,:,:,1) = tau_xx
    !>     tau_ij(:,:,:,2) = tau_yy
    !>     tau_ij(:,:,:,3) = tau_zz
    !>     tau_ij(:,:,:,4) = tau_xy
    !>     tau_ij(:,:,:,5) = tau_xz
    !>     tau_ij(:,:,:,6) = tau_yz
    !> The tensor is symmetric (tau_yx = tau_xy, etc.).

    real(kind=8), dimension(:,:,:,:), intent(out) :: tau_ij
    real(kind=8), dimension(:,:,:), intent(in) :: ux, uy, uz, nu_t
    real(kind=8), intent(in) :: dx, dy, dz

    real(kind=8), allocatable :: g(:,:,:,:)

    allocate(g(size(ux,1), size(ux,2), size(ux,3), 9))
    call compute_velocity_gradients(g, ux, uy, uz, dx, dy, dz)
    call calculate_tau_ij_from_grad(tau_ij, g, nu_t)
    deallocate(g)

  end subroutine calculate_tau_ij

  subroutine calculate_dtau_ij_dxj(dtau_dx, dtau_dy, dtau_dz, &
       tau_ij, dx, dy, dz)
    !> Calculates the divergence of the sub-grid scale stress tensor for a
    !> Large Eddy Simulation (LES).
    !>
    !> INPUT:
    !> tau_ij(nx, ny, nz, 6) : sub-grid scale tensor (order: xx, yy, zz,
    !>                         xy, xz, yz)
    !> dx, dy, dz            : grid spacings
    !>
    !> OUTPUT:
    !> dtau_dx(nx, ny, nz) : d(tau_xj)/dxj
    !> dtau_dy(nx, ny, nz) : d(tau_yj)/dxj
    !> dtau_dz(nx, ny, nz) : d(tau_zj)/dxj
    !>
    !> A single work array is used (instead of 9); the summation order is
    !> identical to the previous version.

    real(kind=8), dimension(:,:,:,:), intent(in) :: tau_ij
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(:,:,:), intent(out) :: dtau_dx, dtau_dy, dtau_dz

    integer :: nx, ny, nz

    nx = size(tau_ij, 1)
    ny = size(tau_ij, 2)
    nz = size(tau_ij, 3)

    if (allocated(les_scratch)) then
       if (any(shape(les_scratch) /= [nx, ny, nz])) deallocate(les_scratch)
    end if
    if (.not. allocated(les_scratch)) allocate(les_scratch(nx, ny, nz))

    ! tau_ij(nx, ny, nz, 6) :      ! 
    !     tau_ij(:,:,:,1) = tau_xx !
    !     tau_ij(:,:,:,2) = tau_yy !
    !     tau_ij(:,:,:,3) = tau_zz !
    !     tau_ij(:,:,:,4) = tau_xy !
    !     tau_ij(:,:,:,5) = tau_xz !
    !     tau_ij(:,:,:,6) = tau_yz !

    ! x-component : d(tau_xx)/dx + d(tau_xy)/dy + d(tau_xz)/dz
    call derxp(dtau_dx, tau_ij(:,:,:,1), dx)
    call deryi(les_scratch, tau_ij(:,:,:,4), dy)
    dtau_dx = dtau_dx + les_scratch
    call derzi(les_scratch, tau_ij(:,:,:,5), dz)
    dtau_dx = dtau_dx + les_scratch

    ! y-component : d(tau_xy)/dx + d(tau_yy)/dy + d(tau_yz)/dz
    call derxi(dtau_dy, tau_ij(:,:,:,4), dx)
    call deryp(les_scratch, tau_ij(:,:,:,2), dy)
    dtau_dy = dtau_dy + les_scratch
    call derzi(les_scratch, tau_ij(:,:,:,6), dz)
    dtau_dy = dtau_dy + les_scratch

    ! z-component : d(tau_xz)/dx + d(tau_yz)/dy + d(tau_zz)/dz
    call derxi(dtau_dz, tau_ij(:,:,:,5), dx)
    call deryi(les_scratch, tau_ij(:,:,:,6), dy)
    dtau_dz = dtau_dz + les_scratch
    call derzp(les_scratch, tau_ij(:,:,:,3), dz)
    dtau_dz = dtau_dz + les_scratch

    if (les_verbose) call print_dtau_ij_dxj_max_values(dtau_dx, dtau_dy, dtau_dz)

    return
  end subroutine calculate_dtau_ij_dxj

end module les_turbulence
