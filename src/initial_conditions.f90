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
    real(kind=8), parameter :: smooth_width = 0.2d0
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
             pp(i,j,k) = 0.0625d0 * (cos(twox) + cos(twoy)) * (cos(twoz) + 2.d0) !https://cfd.ku.edu/hiocfd/case_c3.5.pdf
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
                phi(i,j,k) = 0.5d0 * (1.d0 - tanh((dist - R) / delta))
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
    u1 = u0 / (1. - ratio)
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
    pp = 1.d0 !+ magnitude * magnitude / 2.d0

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
    theta_o = 1.d0 / (13.d0 * l0)
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
    pp = 1.d0

    return
  end subroutine initialize_mixing_layer

  subroutine initialize_homogeneous_isotropic_turbulence(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr)
    implicit none
    integer, intent(in) :: nx, ny, nz, nscr
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(in) :: l0, ratio
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)

    integer :: i,j,k
    integer :: radius ! kernel size = 2*radius+1
    real(kind=8) :: sigma
    real(kind=8), dimension(nx,ny,nz) :: ax, ay, az
    real(kind=8), dimension(nx,ny,nz) :: dazdx, daydx
    real(kind=8), dimension(nx,ny,nz) :: dazdy, daxdy
    real(kind=8), dimension(nx,ny,nz) :: daydz, daxdz
    real(kind=8), dimension(:), allocatable :: kernel
    real(kind=8) :: maxval_ux, maxval_uy, maxval_uz, max_u
    real(kind=8), parameter :: tiny_value = 1.0d-12
    logical, parameter :: use_filter = .true.
    ! Declaration for scalar initialization
    real(kind=8), parameter :: smooth_width = 0.2d0
    real(kind=8) :: R, dist, cx, cy, cz, delt

    integer, dimension(8) :: seed, clock_vals

    radius = int(l0)
    sigma = u0
    allocate(kernel(-radius:radius))

    call date_and_time(values = clock_vals)
    do i = 1, 8
       seed(i) = clock_vals(mod(i, 8) + 1) + i * 17
    end do

    call random_seed(put = seed)
    call random_number(ax)
    call random_number(ay)
    call random_number(az)

    ax = 2.0d0 * (ax-0.5d0)
    ay = 2.0d0 * (ay-0.5d0)
    az = 2.0d0 * (az-0.5d0)

    if (use_filter) then
       call init_kernel(radius, kernel, sigma, 'mexican')
       call apply_gaussian_filter(ax, kernel, nx, ny, nz, radius)
       call apply_gaussian_filter(ay, kernel, nx, ny, nz, radius)
       call apply_gaussian_filter(az, kernel, nx, ny, nz, radius)
    end if

    call derxi(daydx, ay, dx)
    call derxi(dazdx, az, dx)
    call deryi(dazdy, az, dy)
    call deryi(daxdy, ax, dy)
    call derzi(daydz, ay, dz)
    call derzi(daxdz, ax, dz)

    ux = dazdy - daydz
    uy = daxdz - dazdx
    uz = daydx - daxdy

    maxval_ux = maxval(abs(ux))
    maxval_uy = maxval(abs(uy))
    maxval_uz = maxval(abs(uz))

    max_u = max(maxval_ux, maxval_uy, maxval_uz)

    if (max_u > tiny_value) then
       ux = ux / max_u
       uy = uy / max_u
       uz = uz / max_u
    end if

    pp = 1.d0

    if (nscr == 1) then
       print*, "* Initialize Scalar"
       R = ratio ! 0.25d0 * (x(nx) - x(1))
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
  end subroutine initialize_homogeneous_isotropic_turbulence

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
       init_condition => initialize_homogeneous_isotropic_turbulence
    case default
       print *, "Warning: Invalid initialization method selected. Using default."
       init_condition => initialize_vortex  ! Fallback to default method
    end select

    return
  end subroutine set_initialization_type

  subroutine add_turbulent_init(ux, uy, uz, &
       nx, ny, nz, dy, u0, init_noise_x, init_noise_y, init_noise_z, &
       typesim)
    !> Add noise initialization to the velocity components.
    !> 
    !> INPUT:
    !> nx           : Number of grid points in the x-direction
    !> ny           : Number of grid points in the y-direction
    !> nz           : Number of grid points in the z-direction
    !> dy           : Grid spacing in the y-direction
    !> u0           : Reference Velocity of the flow
    !> init_noise_x : Intensity of the initialization noise
    !> init_noise_y : Intensity of the initialization noise
    !> init_noise_z : Intensity of the initialization noise
    !> typesim      : Type of simulation for normalistaion of u_base
    !> OUTPUT:
    !> ux         : X-component of velocity field with noise
    !> uy         : Y-component of velocity field with noise 
    !> uz         : Z-component of velocity field with noise
    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dy, u0, init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, typesim
    real(kind=8), dimension(nx) :: ux_noise 
    real(kind=8), dimension(nx) :: uy_noise
    real(kind=8), dimension(nx) :: uz_noise
    real(kind=8), dimension(ny) :: u_base
    real(kind=8) :: k0, s
    integer :: i, j, k

    call calcul_u_base(u_base, ux(1,:,1), dy)
    if (typesim == 5) then
       call normalize1D(u_base, 0.0d0)
    else 
       call normalize1D(u_base, -1.d0)
    end if
    if (init_noise_x > 0.d0) then
       call print_noise_gene("Ux")
    end if
    if (init_noise_y > 0.d0) then
       call print_noise_gene("Uy")
    end if
    if (init_noise_z > 0.d0) then
       call print_noise_gene("Uz")
    end if

    do k = 1, nz
       k0 = 10.d0
       s = 11.d0 / 3.d0
       call noise_generator_1D(nx, k0, s, ux_noise, 12345)
       call noise_generator_1D(nx, k0, s, uy_noise, 54321)
       call noise_generator_1D(nx, k0, s, uz_noise, 98765)
       call normalize1D(ux_noise, -1.0d0)
       call normalize1D(uy_noise, -1.0d0)
       call normalize1D(uz_noise, -1.0d0)
       do j = 1, ny
          do i = 1, nx 
             ux(i,j,k) = ux(i,j,k) + &
                  u0 * init_noise_x * u_base(j) * ux_noise(i)
             uy(i,j,k) = uy(i,j,k) + &
                  u0 * init_noise_y * u_base(j) * uy_noise(i)
             uz(i,j,k) = uz(i,j,k) + &
                  u0 * init_noise_z * u_base(j) * uz_noise(i)
          end do
       end do
    end do

    return

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
    real(kind=8), dimension(ny) :: u_base
    real(kind=8) :: phase_x
    integer :: i, j, k
    integer, parameter :: kx = 3 ! Number of oscillations in x-direction
    real(kind=8), parameter :: pi = acos(-1.d0)

    ! Compute base velocity profile
    call calcul_u_base(u_base, ux(1,:,1), dy)
    if (typesim == 5) then
       call normalize1D(u_base, 0.0d0)
    else
       call normalize1D(u_base, -1.d0)
    end if

    ! Add oscillations to velocity components
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
