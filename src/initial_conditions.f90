module initial_conditions
  use functions
  use utils
  use noise_mod
  use initialization
  use derivation
  use visualization
  use IOfunctions
  implicit none

  !> Abstract interface for initializing flow conditions.
  interface
     subroutine initialize_type(ux, uy, uz, pp, phi, x, y, z, nx, ny, nz, &
          l0, ratio, nscr, ici, init_noise_x, init_noise_y, init_noise_z)
       real(kind=8), intent(in) :: x(:), y(:), z(:)
       real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
       real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
       real(kind=8), intent(in) :: l0, ratio 
       real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
       integer, intent(in) :: nx, ny, nz, nscr, ici
     end subroutine initialize_type
  end interface

  !> Pointer to a function that initializes the flow conditions.
  procedure(initialize_type), pointer :: init_condition => null()

contains

  subroutine initialize_vortex(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr, ici, &
       init_noise_x, init_noise_y, init_noise_z)
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
    !> ici      : Initial condition index (specifies the type of initialization)
    !> init_noise_x, init_noise_y, init_noise_z : Noise parameters (not used here)
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
    real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: dx, dy
    real(kind=8) :: cv, rv, xc, yc, zc, r, A
    real(kind=8), dimension(nx, ny, nz) :: ksi, dksidx, dksidy
    integer :: i, j, k

    print *, "* Condition Initiale for a convected vortex"
    if (ici == 1) print *, "* No excitation of the initial condition in the case of TGV"
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
       x, y, z, nx, ny, nz, l0, ratio, nscr, ici, &
       init_noise_x, init_noise_y, init_noise_z)
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
    !> ici      : Initial condition index (specifies the type of initialization)
    !> init_noise_x, init_noise_y, init_noise_z : Noise parameters (not used here)
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
    real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, nscr, ici

    integer :: i, j, k 
    integer, parameter :: n_jets_y=11, n_jets_z=11
    real(kind=8) :: twopi, twox, twoy, twoz, y_length, z_length, A

    print *, "* Condition Initiale for a the Taylor-Green vortex"
    if (ici == 1) print *, "* No excitation of the initial condition in the case of TGV"

    A = init_noise_x * init_noise_y * init_noise_z ! Unused for no warning 
    ! during compilation
    twopi = 2.d0 * acos(-1.d0)
    y_length = y(ny) - y(1)
    z_length = z(nz) - z(1)
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
             if (nscr == 1) then
                phi(i,j,k) = 0.5d0 * &
                     (cos(twopi * n_jets_y * (y(j) - y(1)) / y_length) * &                  
                     cos(twopi * n_jets_z * (z(k) - z(1)) / z_length) + 1.d0)
             end if
          end do
       end do
    end do

    return

  end subroutine initialize_taylor_green_vortex

  subroutine initialize_planar_jet(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr, ici, &
       init_noise_x, init_noise_y, init_noise_z)
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
    !> ici      : Initial condition index (specifies the type of initialization)
    !> init_noise_x, init_noise_y, init_noise_z : Noise parameters
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
    real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: u1, u2, t1, t2, t3, t4, theta_o
    real(kind=8) :: phi1, phi2
    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: A, x_disturb
    real(kind=8) :: pi
    integer :: i, j, k, kx

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
    if (ici == 1 .or. ici == 3) then
       A = init_noise_y * u0
       kx = 3
       dy = y(2) - y(1)
       call calcul_u_base(u_base, ux(1,:,1), dy)
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                x_disturb = u_base(j) * sin(2.d0 * pi * kx * x(i) / xlx)
                uy(i,j,k) = uy(i,j,k) + A * x_disturb 
             end do
          end do
       end do
    end if
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 1.d0 !+ magnitude * magnitude / 2.d0

    return

  end subroutine initialize_planar_jet

  subroutine initialize_coplanar_jet(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr, ici, &
       init_noise_x, init_noise_y, init_noise_z)
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
    !> ici      : Initial condition index (specifies the type of initialization)
    !> init_noise_x, init_noise_y, init_noise_z : Noise parameters
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
    real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: u1, u2, u3, theta_1, theta_2, d1, d2, h1, h2, hm
    real(kind=8) :: phi1, phi2
    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: A, y_disturb
    real(kind=8) :: pi
    real(kind=8) :: dy
    integer :: i, j, k, ky

    print *, "* Condition Initiale for a Co-planar Jet simulation"

    pi = acos(-1.d0)
    u2 = u0
    u1 = u2 / ratio
    u3 = 0.091d0 * u2
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
                ux(i,j,k) = 0.5d0 * (u1 + u2) + 0.5d0 * (u2 - u1) * &
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
    if (ici == 1 .or. ici == 3) then
       A = 0.03 * u0
       ky = 3
       dy = y(2) - y(1)
       call calcul_u_base(u_base, ux(1,:,1), dy)
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                y_disturb = cos(2.0 * pi * ky * x(i) / xlx) * u_base(j)
                uy(i,j,k) = uy(i,j,k) + A * y_disturb
             end do
          end do
       end do
    end if
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 0.d0 !+ magnitude * magnitude / 2.d0

    return

  end subroutine initialize_coplanar_jet

  subroutine initialize_mixing_layer(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr, ici, &
       init_noise_x, init_noise_y, init_noise_z)
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
    !> ici      : Initial condition index (specifies the type of initialization)
    !> init_noise_x, init_noise_y, init_noise_z : Noise parameters (not used here)
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
    real(kind=8), intent(in) :: init_noise_x, init_noise_y, init_noise_z
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: u1, u2, t1, t2, t3, t4, theta_o
    real(kind=8) :: A, x_disturb
    real(kind=8) :: pi
    integer :: i, j, k, kx

    print *, "* Condition Initiale for a Mixing Layer simulation"

    pi = acos(-1.d0)
    u2 = u0 / (1.d0 - ratio)
    u1 = u2 * ratio
    theta_o = 1.d0 / (30.d0 * l0)
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
    if (ici == 1 .or. ici == 3) then
       A = init_noise_y * u0
       kx = 1
       dy = y(2) - y(1)
       call calcul_u_base(u_base, ux(1,:,1), dy)
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                x_disturb = u_base(j) * sin(2.d0 * pi * kx * x(i) / xlx)
                uy(i,j,k) = uy(i,j,k) + A * x_disturb 
             end do
          end do
       end do
    end if
    pp = 1.d0

    return
  end subroutine initialize_mixing_layer

  subroutine get_inflow(inflow, ux, uy, uz, pp)
    !> Extract the inflow conditions from the flow field variables.
    !>
    !> This subroutine extracts the inflow conditions from the flow field 
    !> variables, including velocity components (ux, uy, uz) and pressure (pp), 
    !> at the inflow boundary (x = 0).
    !>
    !> INPUT:
    !> ux(nx,ny,nz) : X-component of velocity field
    !> uy(nx,ny,nz) : Y-component of velocity field
    !> uz(nx,ny,nz) : Z-component of velocity field
    !> pp(nx,ny,nz) : Pressure field
    !>
    !> OUTPUT:
    !> inflow(nx,ny,4): Array containing the inflow conditions extracted 
    !>                  from the velocity and pressure fields. The last 
    !>                  dimension contains the following components:
    !>                  - inflow(:,:,1) : X-component of velocity at 
    !>                    the inflow boundary
    !>                  - inflow(:,:,2) : Y-component of velocity at 
    !>                    the inflow boundary
    !>                  - inflow(:,:,3) : Z-component of velocity at 
    !>                    the inflow boundary
    !>                  - inflow(:,:,4) : Pressure at the inflow boundary
    real(kind=8), intent(out) :: inflow(:,:,:)
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:), pp(:,:,:)

    inflow(:,:,1) = ux(1,:,:)
    inflow(:,:,2) = uy(1,:,:)
    inflow(:,:,3) = uz(1,:,:)
    inflow(:,:,4) = pp(1,:,:)

    return
  end subroutine get_inflow

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
       call noise_generator_1D(nx, k0, s, ux_noise)
       call noise_generator_1D(nx, k0, s, uy_noise)
       call noise_generator_1D(nx, k0, s, uz_noise)
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

