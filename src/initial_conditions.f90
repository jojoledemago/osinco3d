module initial_conditions
  use functions
  use utils
  use noise_module
  use initialization
  use derivation
  use visualization
  use IOfunctions
  use poisson
  use diffoper

  implicit none

  !> Abstract interface for initializing flow conditions.
  interface
     subroutine initialize_type(ux, uy, uz, pp, phi, x, y, z, nx, ny, nz, &
          l0, ratio, nscr)
       real(kind=8), intent(in) :: x(:), y(:), z(:)
       real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
       real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
       real(kind=8), intent(in) :: l0, ratio 
       integer, intent(in) :: nx, ny, nz, nscr
     end subroutine initialize_type
  end interface

  !> Pointer to a function that initializes the flow conditions.
  procedure(initialize_type), pointer :: init_condition => null()

contains

  subroutine initialize_vortex(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr)    
    !> Initialize the flow field with a vortex centered at the middle of the domain.
    !> This subroutine sets the velocity field (ux, uy, uz) and the pressure (pp)
    !> based on a Gaussian profile for a vortex.
    !>
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Characteristic length scale (radius of the vortex)
    !> ratio     : Ratio of maximum velocity to the characteristic velocity
    !> nscr     : Flag for scalar transport equation (0 or 1)
    !>
    !> OUTPUT:
    !> ux(nx, ny, nz)  : X-component of velocity initialized with the vortex
    !> uy(nx, ny, nz)  : Y-component of velocity initialized with the vortex
    !> uz(nx, ny, nz)  : Z-component of velocity initialized (zero for this vortex)
    !> pp(nx, ny, nz)  : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized to a constant value
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in) :: l0, ratio 
    integer, intent(in) :: nx, ny, nz, nscr

    real(kind=8) :: dx, dy
    real(kind=8) :: cv, rv, xc, yc, zc, r, A
    real(kind=8), dimension(nx, ny, nz) :: ksi, dksidx, dksidy
    integer :: i, j, k

    print *, "* Condition Initiale for a convected vortex"
    ! Vortex parameters (could be inputs)
    cv = u0   ! Max velocity (example value)
    rv = l0   ! Vortex radius (example value)

    A = init_noise_x * init_noise_y * init_noise_z ! Unused for no warning
    ! during compilation
    xc = 0.5d0 * (x(nx) + x(1))
    yc = 0.5d0 * (y(ny) + y(1))
    zc = 0.5d0 * (z(ny) + z(1))

    ! Initialize the vortex pattern
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             r = sqrt((x(i) - xc)**2 + (y(j) - yc)**2)
             ksi(i,j,k) = cv * exp(-0.5 * (r/rv)**2)
          end do
       end do
    end do

    ! Compute derivatives to find velocity components
    dx = (x(nx) - x(1)) / real(nx - 1, kind=8)
    dy = (y(ny) - y(1)) / real(ny - 1, kind=8)
    call derxi(dksidx, ksi, dx)
    call deryi(dksidy, ksi, dy)

    ux(:,:,:) = ratio * u0 + dksidy(:,:,:)  ! Create swirl
    uy(:,:,:) = - dksidx(:,:,:)
    uz(:,:,:) = 0.0d0
    pp = 0.0d0 ! -10.0566540d0 

    if (nscr == 1) then
       phi = 0.d0
    end if

  end subroutine initialize_vortex

  subroutine initialize_taylor_green_vortex(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr)
    !> Initialize the flow field with Taylor-Green vortex.
    !> This subroutine sets the velocity field (ux, uy, uz) based on the exact
    !> solution for the Taylor-Green vortex in a periodic domain.
    !>
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Characteristic length scale for the problem
    !> ratio    : Ratio of maximum velocity to the characteristic velocity
    !> nscr     : Flag for scalar transport equation (0 or 1)
    !>
    !> OUTPUT:
    !> ux(nx, ny, nz)  : X-component of velocity initialized with Taylor-Green vortex
    !> uy(nx, ny, nz)  : Y-component of velocity initialized with Taylor-Green vortex
    !> uz(nx, ny, nz)  : Z-component of velocity initialized (zero for this vortex)
    !> pp(nx, ny, nz)  : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized to a constant value
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr

    integer :: i, j, k 
    real(kind=8) :: twopi, twox, twoy, twoz
    real(kind=8), parameter :: smooth_width = 0.05d0
    real(kind=8) :: R, dist, cx, cy, cz, delt

    print *, "* Condition Initiale for a the Taylor-Green vortex"
    if (ici == 1) print *, "* No excitation of the initial condition in the case of TGV"

    ! during compilation
    twopi = 2.d0 * acos(-1.d0)
    do k = 1, nz
       twoz = 2.d0 * z(k)
       do j = 1, ny
          twoy = 2.d0 * y(j)
          do i = 1, nx
             twox = 2.d0 * x(i)
             ux(i,j,k) =  ratio * u0/l0 * sin(x(i)) * cos(y(j))* cos(z(k))
             uy(i,j,k) = -ratio * u0/l0 * cos(x(i)) * sin(y(j))* cos(z(k))
             uz(i,j,k) = 0.d0
             pp(i,j,k) = 0.d0 !0.0625d0 * (cos(twox) + cos(twoy)) * (cos(twoz) + 2.d0) !https://cfd.ku.edu/hiocfd/case_c3.5.pdf
          end do
       end do
    end do
    if (nscr == 1) then
       print*, "* Initialize Scalar"
       R = 0.25d0 * (x(nx) - x(1))
       cx = 0.5d0 * (x(1) + x(nx))
       cy = 0.5d0 * (y(1) + y(ny))
       cz = 0.5d0 * (z(1) + z(nz))
       delt = smooth_width * R 
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dist = sqrt( (x(i) - cx)**2 + (y(j) - cy)**2 + (z(k) - cz)**2 )
                phi(i,j,k) = 0.5d0 * (1.d0 - tanh((dist - R) / delt))
             end do
          end do
       end do
    end if

    return

  end subroutine initialize_taylor_green_vortex

  subroutine initialize_planar_jet(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr)
    !> Initialize the flow field for a planar jet.
    !> This subroutine sets the velocity field (ux, uy, uz) for multiple jets
    !> that are coplanar, allowing for the simulation of interactions between
    !> the jets.
    !>
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Characteristic length scale for the jets
    !> ratio    : Ratio of maximum velocity to the characteristic velocity
    !> nscr     : Flag for scalar transport equation (0 or 1)
    !>
    !> OUTPUT:
    !> ux(nx, ny, nz)  : X-component of velocity initialized for the coplanar jet
    !> uy(nx, ny, nz)  : Y-component of velocity initialized for the coplanar jet
    !> uz(nx, ny, nz)  : Z-component of velocity initialized (zero for this jet)
    !> pp(nx, ny, nz)  : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized to a constant value
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in) :: l0, ratio 
    integer, intent(in) :: nx, ny, nz, nscr

    real(kind=8) :: u1, u2, t1, t2, t3, t4, theta_o
    real(kind=8) :: phi1, phi2
    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: pi
    integer :: i, j, k

    print *, "* Condition Initiale for a Planar Jet simulation"

    pi = acos(-1.d0)
    u1 = u0 / (1-ratio)
    u2 = u1 * ratio
    phi1 = 1.d0
    phi2 = 0.d0
    theta_o = 1./30. * l0
    t1 = 0.5d0 * (u1 + u2) 
    t2 = 0.5d0 * (u1 - u2) 
    t3 = l0 / (4.d0 * theta_o)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             t4 = (2.d0 * abs(y(j))) / l0
             ux(i,j,k) = x(i) * 0.d0 + t1 + t2 * tanh(t3 * (1.d0 - t4))
             uy(i,j,k) = y(j) * 0.d0
             uz(i,j,k) = z(k) * 0.d0
             if (nscr == 1) then
                phi(i,j,k) = 0.5d0 * (phi1 + phi2) + 0.5d0 * (phi1 - phi2) * &
                     tanh(t3 * (1.d0 - t4))
             end if
          end do
       end do
    end do
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 0.d0 !+ magnitude * magnitude / 2.d0

    return

  end subroutine initialize_planar_jet

  subroutine initialize_coplanar_jet(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr)
    !> Initialize the flow field for a coplanar jet.
    !> This subroutine sets the velocity field (ux, uy, uz) for multiple jets
    !> that are coplanar, allowing for the simulation of interactions between
    !> the jets.
    !>
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Characteristic length scale for the jets
    !> ratio    : Ratio of maximum velocity to the characteristic velocity
    !> nscr     : Flag for scalar transport equation (0 or 1)
    !>
    !> OUTPUT:
    !> ux(nx, ny, nz)  : X-component of velocity initialized for 
    !>                   the coplanar jet
    !> uy(nx, ny, nz)  : Y-component of velocity initialized for 
    !>                   the coplanar jet
    !> uz(nx, ny, nz)  : Z-component of velocity initialized (zero for this jet)
    !> pp(nx, ny, nz)  : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized to a constant value
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in) :: l0, ratio 
    integer, intent(in) :: nx, ny, nz, nscr

    real(kind=8) :: u1, u2, u3, theta_1, theta_2, d1, d2, h1, h2, hm
    real(kind=8) :: phi1, phi2
    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: pi
    integer :: i, j, k

    print *, "* Condition Initiale for a Co-planar Jet simulation"

    pi = acos(-1.d0)
    u2 = u0
    u1 = u2 / ratio
    u3 = 0.0d0 * u2
    d1 = l0
    d2 = 2.d0 * d1
    h1 = 0.5d0 * d1
    h2 = 0.5d0 * d2
    theta_1 = h1 / 10.d0
    theta_2 = h2 / 25.d0
    hm = 0.5d0 * (h1 + h2) 
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (abs(y(j)) < hm) then
                ux(i,j,k) = x(i) * 0.d0 + 0.5d0 * (u1 + u2) + 0.5d0 * (u2 - u1) * &
                     tanh((abs(y(j)) - h1) / (2.d0 * theta_1))
                if (nscr == 1) then
                   phi1 = 0.d0
                   phi2 = 1.d0
                   phi(i,j,k) = 0.5d0 * (phi1 + phi2) + 0.5d0 * (phi2 - phi1) * &
                        tanh((abs(y(j)) - h1) / (2.d0 * theta_1))
                end if
             else
                ux(i,j,k) = 0.5d0 * (u2 + u3) + 0.5d0 * (u3 - u2) * &
                     tanh((abs(y(j)) - h2) / (2.d0 * theta_2))
                if (nscr == 1) then
                   phi1 = 1.d0
                   phi2 = 0.d0
                   phi(i,j,k) = 0.5d0 * (phi1 + phi2) + 0.5d0 * (phi2 - phi1) * &
                        tanh((abs(y(j)) - h2) / (2.d0 * theta_2))
                end if
             end if
             uy(i,j,k) = y(j) * 0.d0
             uz(i,j,k) = z(k) * 0.d0
          end do
       end do
    end do
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 0.d0 !+ magnitude * magnitude / 2.d0

    return

  end subroutine initialize_coplanar_jet

  subroutine initialize_mixing_layer(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr)
    !> Initialize the flow field for a mixing layer.
    !> This subroutine sets the velocity field (ux, uy, uz) to simulate the
    !> interaction of two streams of fluid with different velocities, creating
    !> a mixing layer.
    !>
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Characteristic length scale for the mixing layer
    !> ratio    : Ratio of maximum velocity to the characteristic velocity
    !> nscr     : Flag for scalar transport equation (0 or 1)
    !>
    !> OUTPUT:
    !> ux(nx, ny, nz)  : X-component of velocity initialized for the mixing layer
    !> uy(nx, ny, nz)  : Y-component of velocity initialized for the mixing layer
    !> uz(nx, ny, nz)  : Z-component of velocity initialized (zero for this case)
    !> pp(nx, ny, nz)  : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized to a constant value
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in) :: l0, ratio 
    integer, intent(in) :: nx, ny, nz, nscr

    real(kind=8), parameter :: tiny_value = 1.d-12
    real(kind=8) :: u1, u2, t1, t2, t3, theta_o
    real(kind=8) :: pi
    integer :: i, j, k

    print *, "* Condition Initiale for a Mixing Layer simulation"

    pi = acos(-1.d0)
    if (ratio < tiny_value .and. ratio > -tiny_value) then
       u2 = u0
       u1 = 0.d0
    else
       u2 = u0 / (1.d0 - ratio)
       u1 = u2 * ratio
    end if
    ! theta_o = 1.d0 / (13.d0 * l0) 
    theta_o = log(2.d0) / (4.d0 * l0)
    t1 = 0.5d0 * (u2 + u1) 
    t2 = 0.5d0 * (u1 - u2)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             t3 = y(j) * log(2.d0) / (2.0d0 * theta_o)
             ux(i,j,k) = x(i) * 0.d0 + t1 - t2 * tanh(t3)
             uy(i,j,k) = y(j) * 0.d0
             uz(i,j,k) = z(k) * 0.d0
             if (nscr == 1) then
                phi(i,j,k) = 0.5d0 - 0.5d0 * tanh(t3)
             end if
          end do
       end do
    end do

    pp = 0.d0

    return
  end subroutine initialize_mixing_layer

  subroutine initialize_hit(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr)
    !> Initialize the flow field with Homogeneous Isotropic Turbulence (HIT).
    !>
    !> Builds a random, solenoidal (div u = 0 in the spectral sense), zero-mean
    !> velocity field with a prescribed Passot-Pouquet energy spectrum in a
    !> periodic domain, using the noise generator of noise_module.
    !> The turbulence intensity is set by an rms normalization, and the
    !> amplitude is then capped so that |u_i| <= u0 everywhere.
    !>
    !> Passot-Pouquet spectrum:
    !>   E(k) = (32/3) sqrt(2/pi) (u'^2 / kp) (k/kp)^4 exp(-2 (k/kp)^2)
    !> which gives the Taylor microscale lambda = 2/kp, hence
    !>   Re_lambda = u' * lambda / nu = 2 * Re / (ratio * kp * l0)
    !> before the simulation is even started.
    !>
    !> Grid convention: dx = xlx/(nx-1), x(i) = x0 + (i-1)*dx, so x(nx) = x0 + xlx.
    !>
    !> INPUT:
    !> x(nx), y(ny), z(nz) : Coordinates of the grid points
    !> nx, ny, nz          : Number of grid points in each direction
    !> l0                  : Reference length scale
    !> re_in               : Reynolds number, Re = u0 * l0 / nu
    !>                       (hence nu = u0 * l0 / re)
    !> u0_in               : Reference velocity, upper bound of |u_i|
    !> ratio               : u0 / u'rms (typically 4 to 5; larger = weaker turbulence)
    !> nscr                : Flag for scalar transport equation (0 or 1)
    !>
    !> OUTPUT:
    !> ux, uy, uz (nx,ny,nz) : Random solenoidal velocity field
    !> pp  (nx,ny,nz)        : Pressure field, initialized to zero
    !> phi (nx,ny,nz)        : Scalar field (smoothed sphere if nscr == 1, else zero)
    !>
    !> Host-associated variables (from modules): xlx, yly, zlz, ici.
    !>
    !> NOTE: the discretization is assumed to be 6th-order centered finite
    !> differences (explicit stencil 3/4, -3/20, 1/60). The cutoff wavenumber
    !> and the resolution warnings are calibrated accordingly.

    real(kind=8), intent(in)  :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in)  :: l0, ratio
    integer, intent(in)       :: nx, ny, nz, nscr

    integer, parameter :: iseed = 78145389           ! seed (reproducibility)
    real(kind=8), parameter :: kpeak_modes = 4.d0    ! spectrum peak, in box modes
    real(kind=8), parameter :: kcut_div    = 6.d0    ! kcut = N/kcut_div * kfund
    ! (6th-order FD: theta_cut ~ 1.05)
    real(kind=8), parameter :: smooth_width = 0.05d0 ! scalar interface width / R

    real(kind=8) :: twopi, kfund, kpeak, kcut
    real(kind=8) :: nu, urms_target, urms, urms_eff, umax, fac
    real(kind=8) :: mx, my, mz, dxm
    real(kind=8) :: lambda, relam, eps, eta, teddy
    real(kind=8) :: theta, kmod_err
    real(kind=8) :: R, dist, cx, cy, cz, delt
    real(kind=8), allocatable :: tx(:,:,:), ty(:,:,:), tz(:,:,:)
    integer :: i, j, k, nxu, nyu, nzu, nmin

    print *, "* Initial condition for Homogeneous Isotropic Turbulence (HIT)"
    if (ici == 1) print *, "* No excitation of the initial condition in the case of HIT"

    twopi = 2.d0 * acos(-1.d0)

    ! Unique points (grid is built with dx = xlx/(nx-1): last plane = periodic image)
    nxu = nx - 1; nyu = ny - 1; nzu = nz - 1
    nmin = min(nxu, nyu, nzu)
    nmin = min(nxu, nyu, nzu)

    ! ------------------------------------------------------------------
    ! Physical parameters
    ! ------------------------------------------------------------------
    nu = u0 * l0 / re
    urms_target = u0 * l0 / ratio

    ! ------------------------------------------------------------------
    ! Wavenumbers
    ! ------------------------------------------------------------------
    kfund = min(twopi / xlx, twopi / yly, twopi / zlz)   ! fundamental wavenumber
    kpeak = kpeak_modes * kfund
    kcut  = real(nmin, 8) / kcut_div * kfund             ! no energy beyond kcut:
    ! scales poorly resolved by
    ! the FD scheme are excluded
    if (kcut < 2.d0 * kpeak) then
       print *, "WARNING: kcut < 2*kpeak -> truncated spectrum; increase the ", &
            "resolution or decrease kpeak_modes"
    end if

    ! ------------------------------------------------------------------
    ! Raw random solenoidal field (amplitude is renormalized below)
    ! ------------------------------------------------------------------
    call init_noise_seed(iseed)
    allocate(tx(nx,ny,nz), ty(nx,ny,nz), tz(nx,ny,nz))
    call generate_solenoidal_noise(tx, ty, tz, nx, ny, nz, &
         xlx, yly, zlz, kpeak, kcut, 1.d0)

    ! Remove residual mean (over unique points)
    mx = sum(tx(1:nxu,1:nyu,1:nzu)) / real(nxu*nyu*nzu, 8)
    my = sum(ty(1:nxu,1:nyu,1:nzu)) / real(nxu*nyu*nzu, 8)
    mz = sum(tz(1:nxu,1:nyu,1:nzu)) / real(nxu*nyu*nzu, 8)
    tx = tx - mx
    ty = ty - my
    tz = tz - mz

    ! rms normalization: u'rms = sqrt( <ux^2 + uy^2 + uz^2> / 3 )
    urms = sqrt( ( sum(tx(1:nxu,1:nyu,1:nzu)**2) &
         + sum(ty(1:nxu,1:nyu,1:nzu)**2) &
         + sum(tz(1:nxu,1:nyu,1:nzu)**2) ) &
         / (3.d0 * real(nxu*nyu*nzu, 8)) )
    fac = urms_target / urms

    ! Safeguard: enforce |u_i| <= u0
    umax = fac * max(maxval(abs(tx)), maxval(abs(ty)), maxval(abs(tz)))
    if (umax > u0) then
       print *, "WARNING: umax/u0 =", umax / u0, &
            " > 1 -> amplitude reduced so that umax = u0"
       fac  = fac * u0 / umax
       umax = u0
    end if
    urms_eff = fac * urms

    ux = fac * tx
    uy = fac * ty
    uz = fac * tz

    pp  = 0.d0
    phi = 0.d0

    ! ------------------------------------------------------------------
    ! Theoretical turbulence diagnostics (Passot-Pouquet: lambda = 2/kpeak)
    ! ------------------------------------------------------------------
    lambda = 2.d0 / kpeak
    relam  = urms_eff * lambda / nu
    eps    = 15.d0 * nu * urms_eff**2 / lambda**2        ! isotropic dissipation
    eta    = (nu**3 / eps)**0.25d0                       ! Kolmogorov scale
    teddy  = (twopi / kpeak) / urms_eff                  ! approx. eddy turnover time
    dxm    = min(xlx/real(nxu,8), yly/real(nyu,8), zlz/real(nzu,8))

    ! Modified wavenumber error of the 6th-order centered scheme at kcut:
    !   k'dx = 3/2 sin(t) - 3/10 sin(2t) + 1/30 sin(3t),  t = k*dx
    theta    = kcut * dxm
    kmod_err = abs( (1.5d0*sin(theta) - 0.3d0*sin(2.d0*theta) &
         + sin(3.d0*theta)/30.d0) / theta - 1.d0 )

    call print_hit_values(nu, urms_eff, umax, u0, relam, lambda, eta, eps, &
         teddy, dxm, kcut, theta, kmod_err)

    call hit_check_les_dns(eta, dxm, kpeak, iles)


    ! ------------------------------------------------------------------
    ! Scalar field: smoothed sphere (tanh profile)
    ! ------------------------------------------------------------------
    if (nscr == 1) then
       print *, "* Initialize Scalar"
       R    = 0.25d0 * min(xlx, yly, zlz)
       cx   = 0.5d0 * (x(1) + x(nx))
       cy   = 0.5d0 * (y(1) + y(ny))
       cz   = 0.5d0 * (z(1) + z(nz))
       delt = smooth_width * R
       if (delt < 3.d0 * dxm) &
            print *, "WARNING: scalar interface thinner than 3*dx -> ", &
            "under-resolved, expect oscillations"
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dist = sqrt((x(i)-cx)**2 + (y(j)-cy)**2 + (z(k)-cz)**2)
                phi(i,j,k) = 0.5d0 * (1.d0 - tanh((dist - R) / delt))
             end do
          end do
       end do
    end if

    deallocate(tx, ty, tz)

  end subroutine initialize_hit

  subroutine initialize_vortex_merge(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr)
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr


    real(kind=8) :: dx, dy
    real(kind=8) :: cv, rv, xc1, yc1, xc2, yc2, r1, r2, A
    real(kind=8), dimension(nx, ny, nz) :: ksi1, dksi1dx, dksi1dy
    real(kind=8), dimension(nx, ny, nz) :: ksi2, dksi2dx, dksi2dy
    integer :: i, j, k

    print *, "* Condition Initiale for a vortex merge"
    ! Vortex parameters (could be inputs)
    cv = u0   ! Max velocity (example value)
    rv = 0.5d0 * l0   ! Vortex radius (example value)

    A = z(1) ! Unused for no warning
    ! during compilation
    xc1 = 0.5d0 * (x(nx) + x(1)) - 1.5d0 * rv
    xc2 = 0.5d0 * (x(nx) + x(1)) + 1.5d0 * rv
    yc1 = 0.5d0 * (y(ny) + y(1)) 
    yc2 = yc1

    ! Initialize the vortex pattern

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             r1 = sqrt((x(i) - xc1)**2 + (y(j) - yc1)**2)
             ksi1(i,j,k) = cv * exp(-0.5 * (r1/rv)**2)
             r2 = sqrt((x(i) - xc2)**2 + (y(j) - yc2)**2)
             ksi2(i,j,k) = cv * exp(-0.5 * (r2/rv)**2)
          end do
       end do
    end do

    ! Compute derivatives to find velocity components
    dx = (x(nx) - x(1)) / real(nx - 1, kind=8)
    dy = (y(ny) - y(1)) / real(ny - 1, kind=8)
    call derxi(dksi1dx, ksi1, dx)
    call deryi(dksi1dy, ksi1, dy)
    call derxi(dksi2dx, ksi2, dx)
    call deryi(dksi2dy, ksi2, dy)

    ux = u0 * ( dksi1dy + ratio * dksi2dy )
    uy = -u0 * ( dksi1dx + ratio * dksi2dx )
    uz(:,:,:) = 0.0d0
    pp = 0.0d0 

    if (nscr == 1) then
       do k=1,nz; do j=1,ny; do i=1,nx
          phi(i,j,k) = exp( - ((x(i)-(xc1+xc2)/2.d0)**2 + (y(j)-yc1)**2) / (2.d0*rv**2) )
       end do; end do; end do
    end if
  end subroutine initialize_vortex_merge

  subroutine set_initialization_type(typesim)
    !> Sets the function pointer to the appropriate initialization function
    !> based on the simulation type identifier.
    !> INPUT:
    !> typesim : Integer representing the type of initialization to use
    integer, intent(in) :: typesim

    select case(typesim)
    case(1)
       init_condition => initialize_vortex 
    case(2)
       init_condition => initialize_taylor_green_vortex
    case(3)
       init_condition => initialize_planar_jet
    case(4)
       init_condition => initialize_coplanar_jet
    case(5)
       init_condition => initialize_mixing_layer
    case(6)
       init_condition => initialize_hit
    case default
       print *, "Warning: Invalid initialization method selected. Using default."
       init_condition => initialize_vortex  ! Fallback to default method
    end select

    return
  end subroutine set_initialization_type

  subroutine add_turbulent_init(ux, uy, uz, nx, ny, nz, dy, &
       u0, init_noise_x, init_noise_y, init_noise_z)
    !> Add a divergence-free perturbation to the velocity field:
    !>   (1) a random 3D noise, isotropic in physical wavenumber (x,z), and
    !>   (2) a coherent 2D (spanwise-uniform) Kelvin-Helmholtz forcing made of
    !>       the fundamental mode m_f and its sub-harmonic m_f/2 (pairing).
    !>
    !> Both parts are u' = curl(A), A_c = g(y) * a_c(x,z), with g(y) taken from
    !> the base flow (|dU/dy| normalized). The x and z derivatives of the
    !> potentials use the FD6 MODIFIED wavenumbers, so u' is divergence-free for
    !> the FD6 operators of the solver (centered stencil 3/4, -3/20, 1/60),
    !> up to round-off.
    !>
    !> Noise: von Karman-like spectrum peaked at kp = 2*pi/(7*l0), cut at
    !>        kc = min(kc_fac*kp, 0.8/max(dx,dz)); max |u'_i| = u0*max(init_noise_c).
    !> Coherent: A = (0, 0, g(y) a(x)), u' = (dA_z/dy, -dA_z/dx, 0);
    !>        fundamental m_f = nint(kp*nx*dx/(2*pi)) (= 4 for lx = 4*lambda_KH),
    !>        sub-harmonic m_f/2 with relative velocity amplitude coh_sub;
    !>        the sum is scaled so that max |u'_i| = coh_amp*u0.
    !>
    !> INPUT:
    !> ux, uy, uz       : base velocity field (before any perturbation)
    !> dy               : grid spacing in y
    !> u0               : reference velocity
    !> init_noise_c     : intensity of the random potential A_c (c = x,y,z);
    !>                    all <= 0 disables the noise
    !> Host-associated (module) variables used: x(:), z(:), dx, dz, l0
    !> Tunable parameters (below): kc_fac, coh_amp (0 disables the coherent part),
    !>                    coh_sub, coh_phase_f, coh_phase_s
    !> OUTPUT:
    !> ux, uy, uz       : velocity with the perturbations added
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: u0, dy
    real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
    real(kind=8), intent(inout) :: ux(nx,ny,nz), uy(nx,ny,nz), uz(nx,ny,nz)
 
    ! ---- Tunable parameters ----
    real(kind=8), parameter :: kc_fac      = 4.d0    ! noise cutoff = kc_fac * kp
    real(kind=8), parameter :: coh_amp     = 0.0d0  ! max |u'|/u0 of the coherent part
    real(kind=8), parameter :: coh_sub     = 0.5d0   ! sub-harmonic / fundamental amplitude
    real(kind=8), parameter :: coh_phase_f = 0.d0    ! phase of the fundamental
    real(kind=8), parameter :: coh_phase_s = 0.d0    ! phase of the sub-harmonic
 
    real(kind=8), parameter :: s = 11.d0 / 3.d0
    real(kind=8), parameter :: fk = 1.d0
    integer, parameter :: nsmooth = 0
 
    real(kind=8) :: pi, kp, kc, kx, kz, kxm, kzm, kapp
    real(kind=8) :: amp, phase, arg, e_k, w(3), r, wq, phq
    real(kind=8) :: g(ny), dg(ny), umax, scale, umaxc, scalec
    real(kind=8) :: maxpux, maxpuy, maxpuz
    real(kind=8), allocatable :: a(:,:,:), dax(:,:,:), daz(:,:,:)
    real(kind=8), allocatable :: pux(:,:,:), puy(:,:,:), puz(:,:,:)
    real(kind=8), allocatable :: ac(:), dac(:), pcx(:,:), pcy(:,:)
    integer :: i, j, k, m, n, c, mmax, nmax, mf, mq, q
    integer :: nseed
    integer, allocatable :: seed(:)
    logical :: noise_on, coh_on
 
    pi = acos(-1.d0)
    w = (/ init_noise_x, init_noise_y, init_noise_z /)
    noise_on = maxval(w) > 0.d0
    coh_on   = coh_amp > 0.d0
    if (.not. noise_on .and. .not. coh_on) return
 
    ! Physical peak wavenumber (most unstable KH mode, lambda = 7*l0)
    kp = max(fk * 2.d0 * pi / (7.d0 * l0), 2.d0 * pi / (nx * dx))
 
    ! ---- Envelope in y, from the BASE flow (before any perturbation) ----
    call compute_envelope(g, dg, ux(1,:,1), dy, ny, nsmooth)
 
    ! =====================================================================
    ! (1) Random 3D noise
    ! =====================================================================
    if (noise_on) then
       print *, "* Add a noise perturbution"
 
       ! Reproducible random numbers
       call random_seed(size=nseed)
       allocate(seed(nseed))
       seed = 19850223
       call random_seed(put=seed)
       deallocate(seed)
 
       kc   = min(kc_fac * kp, 0.8d0 / max(dx, dz))
       mmax = min(int(kc * nx * dx / (2.d0 * pi)), nx/2 - 1)
       nmax = 0
       if (nz > 1) nmax = min(int(kc * nz * dz / (2.d0 * pi)), nz/2 - 1)  ! 2D: no spanwise modes
       print '(A,2F8.3)',  "  perturbation kp, kc      = ", kp, kc
       print '(A,2I5)',    "  perturbation mmax, nmax  = ", mmax, nmax
 
       allocate(a(nx,nz,3), dax(nx,nz,3), daz(nx,nz,3))
       a = 0.d0; dax = 0.d0; daz = 0.d0
 
       ! Random potential a_c(x,z) and its derivatives
       do c = 1, 3
          if (w(c) <= 0.d0) cycle
          do m = 0, mmax
             do n = -nmax, nmax
                if (m == 0 .and. n <= 0) cycle   ! half-plane (real field)
                kx = 2.d0 * pi * real(m, kind=8) / (nx * dx)
                kz = 0.d0
                if (nz > 1) kz = 2.d0 * pi * real(n, kind=8) / (nz * dz)
                kapp = sqrt(kx*kx + kz*kz)
                if (kapp > kc) cycle
 
                ! FD6 modified wavenumbers (centered stencil 3/4, -3/20, 1/60)
                kxm = (1.5d0*sin(kx*dx) - 0.3d0*sin(2.d0*kx*dx) &
                     + sin(3.d0*kx*dx)/30.d0) / dx
                kzm = 0.d0
                if (nz > 1) kzm = (1.5d0*sin(kz*dz) - 0.3d0*sin(2.d0*kz*dz) &
                     + sin(3.d0*kz*dz)/30.d0) / dz
 
                ! von Karman-like spectrum for the velocity, A ~ u/|k|
                e_k = (kapp/kp)**4 / (1.d0 + (kapp/kp)**2)**s
                amp = sqrt(e_k) / kapp
                call random_number(r)
                phase = 2.d0 * pi * r
                do k = 1, nz
                   do i = 1, nx
                      arg = kx * x(i) + kz * z(k) + phase
                      a(i,k,c)   = a(i,k,c)   + w(c) * amp * cos(arg)
                      dax(i,k,c) = dax(i,k,c) - w(c) * amp * kxm * sin(arg)
                      daz(i,k,c) = daz(i,k,c) - w(c) * amp * kzm * sin(arg)
                   end do
                end do
             end do
          end do
       end do
 
       ! u' = curl(A)
       allocate(pux(nx,ny,nz), puy(nx,ny,nz), puz(nx,ny,nz))
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                pux(i,j,k) = dg(j) * a(i,k,3) - g(j) * daz(i,k,2)
                puy(i,j,k) = g(j) * (daz(i,k,1) - dax(i,k,3))
                puz(i,j,k) = g(j) * dax(i,k,2) - dg(j) * a(i,k,1)
             end do
          end do
       end do
 
       ! Global normalization (one scalar: keeps div = 0)
       maxpux = maxval(abs(pux))
       maxpuy = maxval(abs(puy))
       maxpuz = maxval(abs(puz))
       umax = max(maxpux, maxpuy, maxpuz)
       scale = 0.d0
       if (umax > 0.d0) then
          scale = u0 * maxval(w) / umax
          ux = ux + scale * pux
          uy = uy + scale * puy
          uz = uz + scale * puz
       end if
 
       ! Max of the perturbation actually ADDED (after scaling)
       call print_turbulent_perturb(scale*maxpux, scale*maxpuy, scale*maxpuz)
       deallocate(a, dax, daz, pux, puy, puz)
    end if
 
    ! =====================================================================
    ! (2) Coherent 2D forcing: fundamental + sub-harmonic
    !     A = (0, 0, g(y) ac(x)),  u' = (dA_z/dy, -dA_z/dx, 0)
    !     discrete divergence: D6x(dg*ac) + D6y(-g*dac) = dg*dac - dac*dg = 0
    ! =====================================================================
    if (coh_on) then
       mf = max(nint(kp * nx * dx / (2.d0 * pi)), 1)
       allocate(ac(nx), dac(nx), pcx(nx,ny), pcy(nx,ny))
       ac = 0.d0; dac = 0.d0
 
       do q = 1, 2
          if (q == 1) then
             mq = mf;    wq = 1.d0;     phq = coh_phase_f
          else
             if (mf < 2 .or. mod(mf, 2) /= 0 .or. coh_sub <= 0.d0) exit
             mq = mf/2;  wq = coh_sub;  phq = coh_phase_s
          end if
          kx  = 2.d0 * pi * real(mq, kind=8) / (nx * dx)
          kxm = (1.5d0*sin(kx*dx) - 0.3d0*sin(2.d0*kx*dx) &
               + sin(3.d0*kx*dx)/30.d0) / dx
          do i = 1, nx
             arg    = kx * x(i) + phq
             ac(i)  = ac(i)  + wq * cos(arg) / kx
             dac(i) = dac(i) - wq * (kxm / kx) * sin(arg)
          end do
       end do
 
       do j = 1, ny
          do i = 1, nx
             pcx(i,j) =  dg(j) * ac(i)
             pcy(i,j) = -g(j)  * dac(i)
          end do
       end do
 
       umaxc = max(maxval(abs(pcx)), maxval(abs(pcy)))
       if (umaxc > 0.d0) then
          scalec = coh_amp * u0 / umaxc
          do k = 1, nz
             ux(:,:,k) = ux(:,:,k) + scalec * pcx
             uy(:,:,k) = uy(:,:,k) + scalec * pcy
          end do
       end if
 
       print *, "* Add a coherent 2D Kelvin-Helmholtz forcing"
       print '(A,I0)',    "  fundamental mode m        = ", mf
       if (mf >= 2 .and. mod(mf, 2) == 0 .and. coh_sub > 0.d0) then
          print '(A,I0)',    "  sub-harmonic mode m       = ", mf/2
          print '(A,F6.3)',  "  sub / fundamental         = ", coh_sub
       else
          print *, "  no sub-harmonic (m_f odd, < 2, or coh_sub <= 0)"
       end if
       print '(A,ES10.3)', "  coherent max |u'|/u0      = ", coh_amp
       deallocate(ac, dac, pcx, pcy)
    end if
 
  end subroutine add_turbulent_init

  subroutine add_oscillations_init(ux, uy, uz, &
       nx, ny, nz, dy, u0, init_noise_x, init_noise_y, init_noise_z, &
       typesim)
    !> Add oscillatory initialization to the velocity components.
    !>
    !> INPUT:
    !> nx           : Number of grid points in the x-direction
    !> ny           : Number of grid points in the y-direction
    !> nz           : Number of grid points in the z-direction
    !> dy           : Grid spacing in the y-direction
    !> u0           : Reference Velocity of the flow
    !> init_noise_x : Amplitude of the x-direction oscillations
    !> init_noise_y : Amplitude of the y-direction oscillations
    !> init_noise_z : Amplitude of the z-direction oscillations
    !> typesim      : Type of simulation for normalization of u_base
    !> OUTPUT:
    !> ux           : X-component of velocity field with oscillations
    !> uy           : Y-component of velocity field with oscillations
    !> uz           : Z-component of velocity field with oscillations
    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dy, u0, init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, typesim
    real(kind=8), dimension(ny) :: u_base, du_basedy
    real(kind=8) :: phase_x, alpha, phi, alpha_3d, beta_3d
    integer :: i, j, k
    integer, parameter :: kx = 9, n = 5
    real(kind=8), parameter :: pi = acos(-1.d0)

    ! Compute base velocity profile
    call calcul_u_base(u_base, ux(1,:,1), dy)
    call dery1D(du_basedy, u_base, dy)
    if (typesim == 5) then
       call normalize1D(u_base, 0.0d0)
       call normalize1D(du_basedy, 0.0d0)
    else
       call normalize1D(u_base, -1.d0)
    end if

    ! Add oscillations to velocity components
    if (typesim == 5 .or. typesim == 3) then

       alpha = real(n, kind=8) * 2.d0 * pi / xlx 
       phi = 0.5d0 * pi
       alpha_3d = 4.d0 * 2.d0 * pi / xlx
       beta_3d  = 2.d0 * pi / zlz

       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                !ux(i,j,k) = ux(i,j,k) + &
                !     u0 * init_noise_x * u_base(j) * &
                !     (sin(8.d0 * pi * x(i) / xlx) &
                !     + sin(4.d0 * pi * x(i) / xlx) / 8.d0 &
                !     + sin(2.d0 * pi * x(i) / xlx) / 16.d0)
                !uy(i,j,k) = uy(i,j,k) + &
                !     u0 * init_noise_y * u_base(j) * &
                !     (cos(8.d0 * pi * x(i) / xlx) &
                !     + cos(4.d0 * pi * x(i) / xlx) / 8.d0 &
                !     + cos(2.d0 * pi * x(i) / xlx) / 16.d0)
                ux(i,j,k) = ux(i,j,k) + &
                     u0 * u_base(j) * &
                     (init_noise_y * cos(alpha * x(i)) / alpha + &
                     0.5d0 * init_noise_y * &
                     cos(2.d0 * alpha * x(i) + phi) / (2.d0 * alpha))
                uy(i,j,k) = uy(i,j,k) + &
                     u0 * u_base(j) * &
                     (init_noise_y * sin(alpha * x(i)) + &
                     0.5d0 * init_noise_y * sin(2.d0 * alpha * x(i) + phi))
                uz(i,j,k) = uz(i,j,k) 
                ! 3D Perturbation
                ux(i,j,k) = ux(i,j,k) + &
                     u0 * init_noise_z * &
                     (beta_3d / alpha_3d) * &
                     u_base(j) * &
                     cos(alpha_3d * x(i)) * &
                     sin(beta_3d * z(k))
                uz(i,j,k) = uz(i,j,k) - &
                     u0 * init_noise_z * &
                     u_base(j) * &
                     sin(alpha_3d * x(i)) * &
                     cos(beta_3d * z(k))
             end do
          end do
       end do
    else
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                phase_x = 2.d0 * pi * kx * x(i) / xlx

                ux(i,j,k) = ux(i,j,k) + &
                     u0 * init_noise_x * u_base(j) * sin(phase_x)
                uy(i,j,k) = uy(i,j,k) + &
                     u0 * init_noise_y * u_base(j) * sin(phase_x)
                uz(i,j,k) = uz(i,j,k) + &
                     u0 * init_noise_z * u_base(j) * sin(phase_x)
             end do
          end do
       end do
    end if

    return

  end subroutine add_oscillations_init

  subroutine apply_2dsim(uz)
    !> Apply 0. to uz velocity compoment in a 2d case
    !>
    !> INPUT:
    !> uz :  Z-component of velocity field

    real(kind=8), intent(inout) :: uz(:,:,:)

    uz = 0.d0

    return
  end subroutine apply_2dsim

end module initial_conditions
