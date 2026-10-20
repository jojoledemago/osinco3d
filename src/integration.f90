module integration
  use functions
  use poisson
  use multigrid_poisson
  use fft_poisson
  use derivation
  use initialization
  use diffoper
  use les_turbulence
  use IOfunctions, only: print_in_outflow_rate, write_velocity_diverged
  implicit none

  !> Persistent work arrays (allocated on first call, reused afterwards):
  !> they avoid dozens of 3D allocations at every time step.
  real(kind=8), allocatable :: ws_g(:,:,:,:)    ! (nx,ny,nz,9) velocity gradients
  real(kind=8), allocatable :: ws_d2(:,:,:,:)   ! (nx,ny,nz,3) second derivatives / scratch
  real(kind=8), allocatable :: ws_tau(:,:,:,:)  ! (nx,ny,nz,6) sub-grid scale tensor (LES)
  real(kind=8), allocatable :: ws_dtau(:,:,:,:) ! (nx,ny,nz,3) div(tau) (LES)

  private :: ws_g, ws_d2, ws_tau, ws_dtau
  private :: ensure_workspace, select_time_scheme
  private :: update_momentum_component, check_diverged, bc_to_int

contains

  subroutine ensure_workspace(nx, ny, nz, need_tau)
    !> Allocates the work arrays (and reallocates them if the grid size changes).
    integer, intent(in) :: nx, ny, nz
    logical, intent(in) :: need_tau

    if (allocated(ws_g)) then
       if (size(ws_g,1) /= nx .or. size(ws_g,2) /= ny .or. &
            size(ws_g,3) /= nz) then
          deallocate(ws_g, ws_d2, ws_dtau)
          if (allocated(ws_tau)) deallocate(ws_tau)
       end if
    end if

    if (.not. allocated(ws_g)) then
       allocate(ws_g(nx,ny,nz,9), ws_d2(nx,ny,nz,3), ws_dtau(nx,ny,nz,3))
    end if
    if (need_tau .and. .not. allocated(ws_tau)) then
       allocate(ws_tau(nx,ny,nz,6))
    end if

  end subroutine ensure_workspace

  subroutine select_time_scheme(itscheme, itime, adt, bdt, cdt, &
       adu, bdu, cdu, ischeme)
    !> Selects the time integration coefficients.
    !> ischeme : 1 = Euler, 2 = Adams-Bashforth 2, 3 = Adams-Bashforth 3
    integer, intent(in) :: itscheme, itime
    real(kind=8), dimension(3), intent(in) :: adt, bdt, cdt
    real(kind=8), intent(out) :: adu, bdu, cdu
    integer, intent(out) :: ischeme

    if (itscheme == 1 .or. itime == 1) then
       ischeme = 1
    elseif (itscheme == 2 .or. itime == 2) then
       ischeme = 2
    elseif (itscheme == 3) then
       ischeme = 3
    else
       print *, 'itscheme:', itscheme, ' unrecognized'
       stop
    end if

    adu = adt(ischeme)
    bdu = bdt(ischeme)
    cdu = cdt(ischeme)

  end subroutine select_time_scheme

  subroutine update_momentum_component(fu, u_pred, u, d2, ux, uy, uz, gc, &
       dtau, onere, with_sgs, adu, bdu, cdu, itscheme)
    !> For one velocity component u: computes in a single pass
    !>   fu(:,:,:,1) = onere * lap(u) - (u . grad u) [- dtau]
    !>   u_pred      = u + adu*fu1 + bdu*fu2 + cdu*fu3
    !> then shifts the flux history (AB2/AB3).
    !>
    !> d2(:,:,:,1:3) : second derivatives d2u/dx2, d2u/dy2, d2u/dz2
    !> gc(:,:,:,1:3) : gradient of u (du/dx, du/dy, du/dz)

    real(kind=8), intent(inout) :: fu(:,:,:,:)
    real(kind=8), intent(out) :: u_pred(:,:,:)
    real(kind=8), intent(in) :: u(:,:,:), ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: d2(:,:,:,:), gc(:,:,:,:), dtau(:,:,:)
    real(kind=8), intent(in) :: onere, adu, bdu, cdu
    logical, intent(in) :: with_sgs
    integer, intent(in) :: itscheme

    integer :: i, j, k, nx, ny, nz
    real(kind=8) :: f1, f2, f3

    nx = size(u, 1)
    ny = size(u, 2)
    nz = size(u, 3)

    !$omp parallel do collapse(2) default(shared) private(i, f1, f2, f3)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             f1 = onere * ((d2(i,j,k,1) + d2(i,j,k,2)) + d2(i,j,k,3)) - &
                  (ux(i,j,k) * gc(i,j,k,1) + uy(i,j,k) * gc(i,j,k,2) + &
                  uz(i,j,k) * gc(i,j,k,3))
             if (with_sgs) f1 = f1 - dtau(i,j,k)

             f2 = fu(i,j,k,2)
             f3 = fu(i,j,k,3)

             u_pred(i,j,k) = u(i,j,k) + adu * f1 + bdu * f2 + cdu * f3

             fu(i,j,k,1) = f1
             if (itscheme == 3) then
                fu(i,j,k,3) = f2
                fu(i,j,k,2) = f1
             else if (itscheme == 2) then
                fu(i,j,k,2) = f1
             end if
          end do
       end do
    end do

  end subroutine update_momentum_component

  subroutine predict_velocity(ux_pred, uy_pred, uz_pred, ux, uy, uz, &
       fux, fuy, fuz, re, adt, bdt, cdt, itime, itscheme, &
       dx, dy, dz, nx, ny, nz, iles, cs, delta, nu_t)

    !> Predicts the velocity field for a fluid dynamics simulation using
    !> various time integration schemes.
    !>
    !> INPUT:
    !> ux(nx, ny, nz)    : Current x-component of velocity
    !> uy(nx, ny, nz)    : Current y-component of velocity
    !> uz(nx, ny, nz)    : Current z-component of velocity
    !> fux(nx, ny, nz, 3): Storage for velocity ux fluxes or calculations
    !> across time steps
    !> fuy(nx, ny, nz, 3): Storage for velocity uy fluxes or calculations
    !> across time steps
    !> fuz(nx, ny, nz, 3): Storage for velocity uz fluxes or calculations
    !> across time steps
    !> re                : Reynolds number
    !> adt(3)            : Coefficients for current time step integration
    !> bdt(3)            : Coefficients for previous time step integration
    !> cdt(3)            : Coefficients for two time steps back integration
    !> itime             : Current time index
    !> itscheme          : Integration scheme indicator (1=Euler, 2=AB2, 3=AB3)
    !> dx                : Grid spacing in x-direction
    !> dy                : Grid spacing in y-direction
    !> dz                : Grid spacing in z-direction
    !> nx                : Number of grid points in x-direction
    !> ny                : Number of grid points in y-direction
    !> nz                : Number of grid points in z-direction
    !> iles              : Integer for apply LES
    !> cs                : Smagoranski constant
    !> delta             : filter size
    !>
    !> OUTPUT:
    !> ux_pred(nx, ny, nz): Predicted x-component of velocity for next time step
    !> uy_pred(nx, ny, nz): Predicted y-component of velocity for next time step
    !> uz_pred(nx, ny, nz): Predicted z-component of velocity for next time step
    !> nu_t(nx, ny, nz)   : Turbulent viscosity
    !>
    !> NOTES:
    !> The 9 velocity gradients are computed only once and shared between
    !> the LES model and the convective terms.

    integer, intent(in) :: itscheme, itime
    integer, intent(in) :: nx, ny, nz, iles
    real(kind=8), intent(in) :: re, cs, delta
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(3), intent(in) :: adt, bdt, cdt
    real(kind=8), dimension(:,:,:), intent(inout) :: ux, uy, uz
    real(kind=8), dimension(:,:,:), intent(inout) :: nu_t
    real(kind=8), dimension(:,:,:,:), intent(inout) :: fux, fuy, fuz
    real(kind=8), dimension(:,:,:), intent(out) :: ux_pred, uy_pred, uz_pred

    real(kind=8) :: onere, adu, bdu, cdu
    integer :: ischeme
    logical :: with_sgs

    print *, "* Predict velocity"

    with_sgs = (iles == 1)
    call ensure_workspace(nx, ny, nz, with_sgs)

    call select_time_scheme(itscheme, itime, adt, bdt, cdt, &
         adu, bdu, cdu, ischeme)
    print*, cscheme(ischeme)

    onere = 1.d0 / re

    ! The 9 velocity gradients (computed once)
    call compute_velocity_gradients(ws_g, ux, uy, uz, dx, dy, dz)

    if (with_sgs) then
       print *, " LES model calculation"
       call calculate_nu_t_from_grad(nu_t, ws_g, cs, delta, onere)
       call calculate_tau_ij_from_grad(ws_tau, ws_g, nu_t)
       call calculate_dtau_ij_dxj(ws_dtau(:,:,:,1), ws_dtau(:,:,:,2), &
            ws_dtau(:,:,:,3), ws_tau, dx, dy, dz)
    else
       nu_t = 0.0d0
    end if

    ! --- ux ---
    call derxxi(ws_d2(:,:,:,1), ux, dx)
    call deryyp(ws_d2(:,:,:,2), ux, dy)
    call derzzp(ws_d2(:,:,:,3), ux, dz)
    call update_momentum_component(fux, ux_pred, ux, ws_d2, ux, uy, uz, &
         ws_g(:,:,:,1:3), ws_dtau(:,:,:,1), onere, with_sgs, &
         adu, bdu, cdu, itscheme)

    ! --- uy ---
    call derxxp(ws_d2(:,:,:,1), uy, dx)
    call deryyi(ws_d2(:,:,:,2), uy, dy)
    call derzzp(ws_d2(:,:,:,3), uy, dz)
    call update_momentum_component(fuy, uy_pred, uy, ws_d2, ux, uy, uz, &
         ws_g(:,:,:,4:6), ws_dtau(:,:,:,2), onere, with_sgs, &
         adu, bdu, cdu, itscheme)

    ! --- uz ---
    call derxxp(ws_d2(:,:,:,1), uz, dx)
    call deryyp(ws_d2(:,:,:,2), uz, dy)
    call derzzi(ws_d2(:,:,:,3), uz, dz)
    call update_momentum_component(fuz, uz_pred, uz, ws_d2, ux, uy, uz, &
         ws_g(:,:,:,7:9), ws_dtau(:,:,:,3), onere, with_sgs, &
         adu, bdu, cdu, itscheme)

    print *, ""

    return
  end subroutine predict_velocity

  integer function bc_to_int(bc_str, axis)
    !> Converts the boundary condition string ('00' periodic, '11' free-slip)
    !> into an integer code for the multigrid solver.
    character(len=*), intent(in) :: bc_str
    character(len=1), intent(in) :: axis

    if (bc_str == '00') then
       bc_to_int = 0
    else if (bc_str == '11') then
       bc_to_int = 1
    else
       write(*,*) 'multigrid: unsupported boundary condition "', &
            trim(bc_str), '" in direction ', axis
       stop
    end if
  end function bc_to_int

  subroutine correct_pression(pp, ux_pred, uy_pred, uz_pred, dx, dy, dz, &
       nx, ny, nz, dt, omega, eps, kmax, idyn, psolver)

    !> This subroutine corrects the pressure field using the divergence of the predicted velocity field
    !> and solves the Poisson equation to enforce incompressibility (divergence-free condition).
    !>
    !> INPUT:
    !> ux_pred(nx, ny, nz) : Current x-component of predicted velocity
    !> uy_pred(nx, ny, nz) : Current y-component of predicted velocity
    !> uz_pred(nx, ny, nz) : Current z-component of predicted velocity
    !> dx                   : Grid spacing in the x-direction
    !> dy                   : Grid spacing in the y-direction
    !> dz                   : Grid spacing in the z-direction
    !> nx, ny, nz           : Number of grid points in the x, y, z directions
    !> dt                   : Time step size
    !> eps                  : Convergence criterion for poisson solver
    !> kmax                 : Maximum number of iterations for poisson solver
    !> idyn                 : Flg for omega dynamic
    !> multigrid            : Flag for process multigrid solver for poisson
    !>
    !> OUTPUT:
    !> pp(nx, ny, nz)       : Updated pressure field after convergence
    !> omega                : Relaxation factor for SOR solver

    real(kind=8), intent(inout) :: pp(:,:,:), omega
    real(kind=8), intent(in) :: ux_pred(:,:,:), uy_pred(:,:,:), uz_pred(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz, dt, eps
    integer, intent(in) :: nx, ny, nz, kmax, idyn, psolver

    real(kind=8), allocatable :: rhs(:,:,:)
    integer :: bc(3)

    print *, "* Correction pression"
    allocate(rhs(nx, ny, nz))

    ! Divergence of the predicted velocity, written directly into rhs
    call divergence(rhs, ux_pred, uy_pred, uz_pred, &
         dx, dy, dz, nx, ny, nz, 1)  ! Assume 1 for odd conditions

    ! Right-hand side of the Poisson equation (in place)
    rhs = rhs / dt

    if (psolver == 0) then
       ! >>> Use standard SOR solver
       call poisson_solver(pp, rhs, dx, dy, dz, nx, ny, nz, omega, eps, kmax, idyn)
    elseif (psolver == 1) then
       ! >> Use multigrid solver
       bc(1) = bc_to_int(bc_x, 'x')
       bc(2) = bc_to_int(bc_y, 'y')
       bc(3) = bc_to_int(bc_z, 'z')
       call poisson_multigrid_solver(pp, rhs, dx, dy, dz, nx, ny, nz, &
            bc, eps, kmax, verbose=.true.)!, iter_out, res_out, verbose, order)
    elseif(psolver == 2) then
       ! >> Use sepectral method
       bc(1) = bc_to_int(bc_x, 'x')
       bc(2) = bc_to_int(bc_y, 'y')
       bc(3) = bc_to_int(bc_z, 'z')
       call poisson_fft_solver(pp, rhs, dx, dy, dz, nx, ny, nz, &
            bc, verbose=.true.)
    else
       write(*,*) 'multigrid must be 0 or 1 in parametars.o3d'
    end if

    deallocate(rhs)
    print *,""

    return
  end subroutine correct_pression

  subroutine check_diverged(u)
    !> Stops the simulation if the field contains a NaN or if |u| > 1000.
    real(kind=8), intent(in) :: u(:,:,:)

    if (contains_nan(u) .or. maxval(abs(u)) > 1000.d0) then
       call write_velocity_diverged()
       stop
    end if
  end subroutine check_diverged

  subroutine correct_velocity(ux, uy, uz, ux_pred, uy_pred, uz_pred, &
       pp, dt, dx, dy, dz, nx, ny, nz)
    !> Correct the velocity components based on the pressure gradient.
    !> This step is part of the fractional step method to solve
    !> incompressible Navier-Stokes equations.
    !> The function also includes a stability check to ensure that
    !> the corrected velocity components do not contain NaN values
    !> and stay bounded (|u| <= 1000).
    !
    !> INPUT:
    !> ux_pred(nx, ny, nz): Predicted x-component of velocity
    !> uy_pred(nx, ny, nz): Predicted y-component of velocity
    !> uz_pred(nx, ny, nz): Predicted z-component of velocity
    !> pp(nx, ny, nz)     : Current pressure field
    !> dt                 : Time step size
    !> dx, dy, dz         : Grid spacing in x, y, z directions
    !> nx, ny, nz         : Number of grid points in x, y, z directions
    !
    !> OUTPUT:
    !> ux(nx, ny, nz)     : Corrected x-component of velocity
    !> uy(nx, ny, nz)     : Corrected y-component of velocity
    !> uz(nx, ny, nz)     : Corrected z-component of velocity
    !>
    !> NOTES:
    !> If a NaN (or |u| > 1000) is found in any component, the simulation
    !> is halted after writing the diverged fields.

    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: ux_pred(:,:,:), uy_pred(:,:,:), uz_pred(:,:,:)
    real(kind=8), intent(in) :: pp(:,:,:)
    real(kind=8), intent(in) :: dt, dx, dy, dz
    integer, intent(in) :: nx, ny, nz

    print *, "* Correct velocity"

    call ensure_workspace(nx, ny, nz, .false.)

    ! A single work array for the pressure gradient.
    ! Subtracting because gradient points from high to low pressure
    call derxp(ws_d2(:,:,:,1), pp, dx)
    ux = ux_pred - dt * ws_d2(:,:,:,1)
    call deryp(ws_d2(:,:,:,1), pp, dy)
    uy = uy_pred - dt * ws_d2(:,:,:,1)
    call derzp(ws_d2(:,:,:,1), pp, dz)
    uz = uz_pred - dt * ws_d2(:,:,:,1)

    ! Stability check
    call check_diverged(ux)
    call check_diverged(uy)
    call check_diverged(uz)

    print*, ""

    return
  end subroutine correct_velocity

  subroutine transeq(phi, ux, uy, uz, src, fphi, re, sc, adt, bdt, cdt, &
       itime, itscheme, dx, dy, dz, nx, ny, nz, iles, nu_t)
    !> Solves the 3D transport equation for a passive scalar.
    !>
    !> The equation solved is:
    !>     d(phi)/dt + u . grad(phi) = div(alpha_eff * grad(phi)) + src
    !>
    !> The subroutine supports Euler, Adams-Bashforth 2 and 3 schemes.
    !> A conservative post-processing step is applied to clip phi
    !> within [0,1] and redistribute any excess/deficit to preserve
    !> the global mean.
    !>
    !> INPUT:
    !>   phi(nx,ny,nz)   : Scalar field to be transported
    !>   ux, uy, uz       : Velocity components (x, y, z)
    !>   src(nx,ny,nz)    : Source term
    !>   re               : Reynolds number
    !>   sc               : Schmidt number
    !>   adt(3), bdt(3), cdt(3) : Time integration coefficients
    !>   itime            : Current time step index
    !>   itscheme         : Integration scheme (1=Euler, 2=AB2, 3=AB3)
    !>   dx, dy, dz       : Grid spacings in x, y, z
    !>   nx, ny, nz       : Grid dimensions
    !>   iles             : Flag for LES model (1=on, 0=off)
    !>   nu_t(nx,ny,nz)   : Turbulent viscosity (LES)
    !>
    !> OUTPUT:
    !>   phi(nx,ny,nz)    : Updated scalar field after one time step
    !>   fphi(nx,ny,nz,3) : Time-step fluxes for AB2/AB3 integration

    implicit none

    integer, intent(in) :: itime, itscheme, nx, ny, nz, iles
    real(kind=8), intent(in) :: re, sc, dx, dy, dz
    real(kind=8), dimension(3), intent(in) :: adt, bdt, cdt
    real(kind=8), dimension(nx,ny,nz), intent(in) :: ux, uy, uz, src, nu_t
    real(kind=8), dimension(nx,ny,nz), intent(inout) :: phi
    real(kind=8), dimension(nx,ny,nz,3), intent(inout) :: fphi

    character(len=20), parameter :: scheme_name(3) = &
         [character(len=20) :: " Euler", " Adams-Bashforth 2", &
         " Adams-Bashforth 3"]

    real(kind=8) :: adu, bdu, cdu, alpha0, alpha, f1, f2, f3
    real(kind=8) :: ncells, phi_old_avg, phi_new_avg, excess, wsum
    integer :: i, j, k, ischeme
    logical :: with_les

    print *, "* Transport equation"

    with_les = (iles == 1)
    call ensure_workspace(nx, ny, nz, .false.)

    !--- Select time integration coefficients ---
    call select_time_scheme(itscheme, itime, adt, bdt, cdt, &
         adu, bdu, cdu, ischeme)
    print *, trim(scheme_name(ischeme))

    !--- Molecular diffusion coefficient (turbulent part added per point) ---
    alpha0 = 1.d0 / (re * sc)

    !--- First derivatives (order 6 centered) ---
    call derxp(ws_g(:,:,:,1), phi, dx)
    call deryp(ws_g(:,:,:,2), phi, dy)
    call derzp(ws_g(:,:,:,3), phi, dz)

    !--- Second derivatives (order 4 centered) ---
    call derxxp(ws_d2(:,:,:,1), phi, dx)
    call deryyp(ws_d2(:,:,:,2), phi, dy)
    call derzzp(ws_d2(:,:,:,3), phi, dz)

    !--- Fluxes (diffusion + convection + source), time advance and
    !    history shift, all in one pass ---
    !$omp parallel do collapse(2) default(shared) private(i, alpha, f1, f2, f3)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             alpha = alpha0
             if (with_les) alpha = alpha0 + nu_t(i,j,k) / sc

             f1 = alpha * ((ws_d2(i,j,k,1) + ws_d2(i,j,k,2)) + ws_d2(i,j,k,3)) - &
                  (ux(i,j,k) * ws_g(i,j,k,1) + uy(i,j,k) * ws_g(i,j,k,2) + &
                  uz(i,j,k) * ws_g(i,j,k,3)) + src(i,j,k)

             f2 = fphi(i,j,k,2)
             f3 = fphi(i,j,k,3)

             phi(i,j,k) = phi(i,j,k) + adu * f1 + bdu * f2 + cdu * f3

             fphi(i,j,k,1) = f1
             if (itscheme == 3) then
                fphi(i,j,k,3) = f2
                fphi(i,j,k,2) = f1
             else if (itscheme == 2) then
                fphi(i,j,k,2) = f1
             end if
          end do
       end do
    end do

    !===============================
    ! Post-processing: conservative clipping
    !===============================
    ncells = real(nx, kind=8) * real(ny, kind=8) * real(nz, kind=8)

    ! 1) Global mean before clipping
    phi_old_avg = sum(phi) / ncells

    ! 2) Clip phi to [0,1]
    phi = max(0.d0, min(1.d0, phi))

    ! 3) Excess/deficit after clipping
    phi_new_avg = sum(phi) / ncells
    excess = phi_old_avg - phi_new_avg

    ! 4-6) Redistribute only if clipping actually changed the mean.
    !      (If nothing was clipped, excess is exactly 0 and this step
    !      would have no effect.)
    if (abs(excess) > 0.d0) then
       ! Weights proportional to distance from the bounds
       wsum = sum(min(phi, 1.d0 - phi))
       if (wsum > tiny(1.d0)) then
          phi = phi + excess * (min(phi, 1.d0 - phi) / wsum)
          ! Reclip after redistribution for safety
          phi = max(0.d0, min(1.d0, phi))
       else
          print *, " transeq: cannot redistribute clipping excess (phi is 0/1 everywhere)"
       end if
    end if

    print *, ""
    return
  end subroutine transeq

end module integration
