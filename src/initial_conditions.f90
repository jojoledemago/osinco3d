module initial_conditions
  use functions
  use utils
  use initialization
  use derivation
  use visualization
  implicit none

  !> Abstract interface for initializing flow conditions.
  interface
     subroutine initialize_type(ux, uy, uz, pp, phi, x, y, z, nx, ny, nz, l0, ratio, nscr, ici)
       real(kind=8), intent(in) :: x(:), y(:), z(:)
       real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
       real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
       real(kind=8), optional, intent(in) :: l0, ratio
       integer, intent(in) :: nx, ny, nz, nscr, ici
     end subroutine initialize_type
  end interface

  !> Pointer to a function that initializes the flow conditions.
  procedure(initialize_type), pointer :: init_condition => null()

contains

  subroutine initialize_vortex(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr, ici)
    !> Initialize the flow field with a vortex centered at the middle of the domain.
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> nscr     : Flag for the scalar transport equation
    !> OUTPUT:
    !> ux(nx, ny, nz)  : X-component of velocity initialized with the vortex
    !> uy(nx, ny, nz)  : Y-component of velocity initialized with the vortex
    !> uz(nx, ny, nz)  : Z-component of velocity initialized (zero for this vortex)
    !> pp(nx, ny, nz)  : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized to a constant value

    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), optional, intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: dx, dy
    real(kind=8) :: cv, rv, xc, yc, zc, r
    real(kind=8), dimension(nx, ny, nz) :: ksi, dksidx, dksidy
    integer :: i, j, k

    print *, "* Condition Initiale for a convected vortex"
    if (ici == 1) print *, "* No excitation of the initial condition in the case of TGV"
    ! Vortex parameters (could be inputs)
    cv = u0   ! Max velocity (example value)
    rv = l0   ! Vortex radius (example value)

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
       x, y, z, nx, ny, nz, l0, ratio, nscr, ici)
    !> Initialize the flow field with the Taylor-Green vortex
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> nscr      : Flag for the scalar transport equation
    !> OUTPUT:
    !> ux(nx, ny, nz) : X-component of velocity initialized with the vortex
    !> uy(nx, ny, nz) : Y-component of velocity initialized with the vortex
    !> uz(nx, ny, nz) : Z-component of velocity initialized (zero for this vortex)
    !> pp(nx, ny, nz) : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field

    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), optional, intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr, ici

    integer :: i, j, k 
    integer, parameter :: n_jets_y=11, n_jets_z=11
    real(kind=8) :: twopi, y_length, z_length

    print *, "* Condition Initiale for a the Taylor-Green vortex"
    if (ici == 1) print *, "* No excitation of the initial condition in the case of TGV"

    twopi = 2.d0 * acos(-1.d0)
    y_length = y(ny) - y(1)
    z_length = z(nz) - z(1)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ux(i,j,k) =  ratio * u0/l0 * sin(x(i)) * cos(y(j))* cos(z(k))
             uy(i,j,k) = -ratio * u0/l0 * cos(x(i)) * sin(y(j))* cos(z(k))
             uz(i,j,k) = 0.d0
             pp(i,j,k) = 0.d0
             if (nscr == 1) then
                phi(i,j,k) = 0.5d0 * &
                     (sin(twopi * n_jets_y * (y(j) - y(1)) / y_length) * &                  
                     sin(twopi * n_jets_z * (z(k) - z(1)) / z_length) + 1.d0)
             end if
          end do
       end do
    end do

    return
  end subroutine initialize_taylor_green_vortex

  subroutine initialize_planar_jet(ux, uy, uz, pp, phi, &
       x, y, z, nx, ny, nz, l0, ratio, nscr, ici)
    !> Initialize the flow field with a planar jet simulation
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Height of the planar jet
    !> ratio    : Ratio between u1 and u2
    !> nscr      : Flag for the scalar transport equation
    !> OUTPUT:
    !> ux(nx, ny, nz) : X-component of velocity initialized with the vortex
    !> uy(nx, ny, nz) : Y-component of velocity initialized with the vortex
    !> uz(nx, ny, nz) : Z-component of velocity initialized (zero for this      vortex)
    !> pp(nx, ny, nz) : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Scalar field initialized 

    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), optional, intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: u1, u2, t1, t2, t3, t4, theta_o
    real(kind=8) :: phi1, phi2
    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: A, x_disturb, y_disturb, z_disturb
    real(kind=8) :: pi, scale_factor
    integer :: i, j, k, h, kx, ky, kz, nharmonics

    print *, "* Condition Initiale for a Planar Jet simulation"

    pi = acos(-1.d0)
    u1 = u0 / (1 - ratio)
    u2 = u1 * ratio
    phi1 = 1.d0
    phi2 = 0.d0
    theta_o = 1./20. * l0
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
       A = 0.01 * u0
       kx = 4
       ky = 2
       kz = 1
       nharmonics = 8
       dy = y(2) - y(1)
       ! Calcul du facteur de normalisation pour garantir que la somme des harmoniques est entre -1 et 1
       scale_factor = 0.d0
       do h = 1, nharmonics
          scale_factor = scale_factor + 1.d0 / h
       end do
       call calcul_u_base(u_base, ux(1,:,1), dy)
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                x_disturb = 0.d0
                y_disturb = 0.d0
                z_disturb = 0.d0
                do h = 1, nharmonics
                   x_disturb = x_disturb + &
                        sin(2.0 * pi * kx * h * x(i) / xlx) * u_base(j)
                   y_disturb = y_disturb + &
                        cos(2.0 * pi * ky * h * x(i) / xlx) * u_base(j)
                   z_disturb = z_disturb + &
                        sin(2.0 * pi * kz * h * x(i) / xlx) * u_base(j)
                end do
                ux(i,j,k) = ux(i,j,k) + A * x_disturb / scale_factor
                uy(i,j,k) = uy(i,j,k) + A * y_disturb / scale_factor
                uz(i,j,k) = uz(i,j,k) + A * z_disturb / scale_factor
             end do
          end do
       end do
    end if
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 1.d0 !+ magnitude * magnitude / 2.d0

    return

  end subroutine initialize_planar_jet

  subroutine initialize_coplanar_jet(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr, ici)
    !> Initialize the flow field with a co-planar jet simulation
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Height of the planar jet
    !> ratio    : Ratio between u1 and u2
    !> nscr      : Flag for the scalar transport equation
    !> OUTPUT:
    !> ux(nx, ny, nz) : X-component of velocity initialized with the vortex
    !> uy(nx, ny, nz) : Y-component of velocity initialized with the vortex
    !> uz(nx, ny, nz) : Z-component of velocity initialized (zero for this      vortex)
    !> pp(nx, ny, nz) : Pressure field initialized to a constant value
    !> phi(nx, ny, nz) : Pressure field

    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), optional, intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8) :: u1, u2, u3, theta_1, theta_2, d1, d2, h1, h2, hm
    real(kind=8) :: phi1, phi2
    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: A, x_disturb, y_disturb
    real(kind=8) :: pi, scale_factor
    real(kind=8) :: dy
    integer :: i, j, k, h, kx, ky, nharmonics

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
                ux(i,j,k) = 0.5d0 * (u2 + u3) + 0.5d0 * (u3- u2) * &
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
       A = 0.01 * u0
       kx = 5
       ky = 3
       nharmonics = 6
       dy = y(2) - y(1)
       ! Calcul du facteur de normalisation pour garantir que la somme des harmoniques est entre -1 et 1
       scale_factor = 0.d0
       do h = 1, nharmonics
          scale_factor = scale_factor + 1.d0 / h
       end do
       call calcul_u_base(u_base, ux(1,:,1), dy)
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                x_disturb = 0.d0
                y_disturb = 0.d0
                do h = 1, nharmonics
                   x_disturb = x_disturb + &
                        sin(2.0 * pi * kx * h * x(i) / xlx) * u_base(j)
                   y_disturb = y_disturb + &
                        cos(2.0 * pi * ky * h * x(i) / xlx) * u_base(j)
                end do
                ux(i,j,k) = ux(i,j,k) + A * x_disturb / scale_factor
                uy(i,j,k) = uy(i,j,k) + A * y_disturb / scale_factor
             end do
          end do
       end do
    end if
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 0.d0 !+ magnitude * magnitude / 2.d0

    return

  end subroutine initialize_coplanar_jet

  subroutine initialize_mixing_layer(ux, uy, uz, pp, phi, x, y, z, &
       nx, ny, nz, l0, ratio, nscr, ici)
    !> Initialize the flow field with a mixing simulation
    !> INPUT:
    !> x(nx)    : X-coordinates of the domain
    !> y(ny)    : Y-coordinates of the domain
    !> z(nz)    : Z-coordinates of the domain (not used in 2D cases)
    !> nx       : Number of grid points in the x-direction
    !> ny       : Number of grid points in the y-direction
    !> nz       : Number of grid points in the z-direction
    !> l0       : Mixing layer thickness
    !> ratio    : Ratio between u1 and u2
    !> nscr     : Flag for the scalar transport equation
    !> ici      : Flag for sub-harmonic disturb
    !> OUTPUT:
    !> ux(nx, ny, nz) : X-component of velocity initialized with the vortex
    !> uy(nx, ny, nz) : Y-component of velocity initialized with the vortex
    !> uz(nx, ny, nz) : Z-component of velocity initialized 
    !> pp(nx, ny, nz) : Dynamic pressure
    !> phi(nx, ny, nz) : Scalar fields

    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(out) :: pp(:,:,:), phi(:,:,:)
    real(kind=8), optional, intent(in) :: l0, ratio
    integer, intent(in) :: nx, ny, nz, nscr, ici

    real(kind=8), dimension(nx,ny,nz) :: magnitude
    real(kind=8) :: u1, u2, t1, t2, t3, theta_o
    real(kind=8) :: A, x_disturb, y_disturb, z_disturb
    real(kind=8) :: sigma, pi, scale_factor
    integer :: i, j, k, h, kx, ky, kz, nharmonics

    print *, "* Condition Initiale for a Mixing Layer simulation"

    pi = acos(-1.d0)
    u1 = (2.d0 * u0 * ratio) / (1.d0 + abs(ratio))
    u2 = u1 / ratio
    theta_o = 1.d0 / 1.d0 * l0
    sigma = theta_o
    t1 = 0.5d0 * (u2 + u1) 
    t2 = 0.5d0 * (u1 - u2)
    A = 0.01 * u0
    kx = 8
    ky = 4
    kz = 2
    nharmonics = 8
    dy = y(2) - y(1)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             t3 = 2.d0 * y(j) / theta_o
             ux(i,j,k) = x(i) * 0.d0 + t1 - t2 * tanh(t3)
             uy(i,j,k) = y(j) * 0.d0
             uz(i,j,k) = z(k) * 0.d0
             if (nscr == 1) then
                phi(i,j,k) = 1.d0 - exp(-2.d0 * sigma * ((y(j) / theta_o)**2)) ! 0.5d0 - 0.5d0 * tanh(t3)
             end if
          end do
       end do
    end do
    if (ici == 1 .or. ici == 3) then
       ! Calcul du facteur de normalisation pour garantir que la somme des harmoniques est entre -1 et 1
       scale_factor = 0.d0
       do h = 1, nharmonics
          scale_factor = scale_factor + 1.d0 / h
       end do
       call calcul_u_base(u_base, ux(1,:,1), dy)
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                x_disturb = 0.d0
                y_disturb = 0.d0
                z_disturb = 0.d0
                do h = 1, nharmonics
                   x_disturb = x_disturb + sin(2.0 * pi * kx * h * x(i) / xlx) * u_base(j)
                   y_disturb = y_disturb + cos(2.0 * pi * ky * h * x(i) / xlx) * u_base(j)
                   z_disturb = z_disturb + sin(2.0 * pi * kz * h * x(k) / xlx) * u_base(j)
                end do
                ux(i,j,k) = ux(i,j,k) + A * x_disturb / scale_factor
                uy(i,j,k) = uy(i,j,k) + A * y_disturb / scale_factor
                uz(i,j,k) = uz(i,j,k) + A * z_disturb / scale_factor
             end do
          end do
       end do
    end if
    magnitude = compute_velocity_magnitude(ux, uy, uz, nx, ny, nz) 
    pp = 1.d0 !+ magnitude * magnitude / 2.d0

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
       nx, ny, nz, dy, u0, init_noise)
    !> Add noise initialization to the velocity components.
    !> 
    !> INPUT:
    !> nx         : Number of grid points in the x-direction
    !> ny         : Number of grid points in the y-direction
    !> nz         : Number of grid points in the z-direction
    !> dy         : Grid spacing in the y-direction
    !> u0         : Reference Velocity of the flow
    !> init_noise : Intensity of the initialization noise
    !> OUTPUT:
    !> ux         : X-component of velocity field with noise
    !> uy         : Y-component of velocity field with noise 
    !> uz         : Z-component of velocity field with noise
    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dy, u0, init_noise
    integer, intent(in) :: nx, ny, nz

    real(kind=8), dimension(ny, nz) :: ux_noise, uy_noise, uz_noise, white_noise
    real(kind=8), dimension(ny) :: u_base
    integer :: i, j, k
    ! Parameters 

    call calcul_u_base(u_base, ux(1,:,1), dy)

    do i = 1, nx
       call generate_white_noise(white_noise, ny, nz)
       call apply_energy_spectrum(ux_noise, white_noise, ny, nz)
       call generate_white_noise(white_noise, ny, nz)
       call apply_energy_spectrum(uy_noise, white_noise, ny, nz)
       call generate_white_noise(white_noise, ny, nz)
       call apply_energy_spectrum(uz_noise, white_noise, ny, nz)
       do j = 1, ny
          do k = 1, nz
             ux(i,j,k) = ux(i,j,k) + u0 * init_noise * u_base(j) * ux_noise(j,k)
             uy(i,j,k) = uy(i,j,k) + u0 * init_noise * u_base(j) * uy_noise(j,k)
             uz(i,j,k) = uz(i,j,k) + u0 * init_noise * u_base(j) * uz_noise(j,k)
          end do
       end do
    end do

    return

  end subroutine add_turbulent_init

end module initial_conditions

