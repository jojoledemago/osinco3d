module initialization
  use poisson
  use derivation
  implicit none

  !> This module initializes and configures the simulation parameters, variables, and methods.
  !> It sets up the grid, time stepping, and derivative calculation methods based on boundary conditions and simulation type.

  real(kind=8) :: divu_stats(6)
  character(len=10) :: response
  ! Domain definition variables
  integer :: nx, ny, nz           !> Grid points in the x, y, and z directions
  real(kind=8) :: xlx, yly, zlz   !> Physical dimensions of the domain in x, y, and z directions
  real(kind=8) :: x0, y0, z0      !> Origin of the x, y and z direction
  real(kind=8) :: dx, dy, dz      !> Grid spacing in x, y, and z directions
  real(kind=8) ,allocatable, dimension(:) :: x, y, z !> coordonne of the domain

  ! Time integration variables
  integer :: itstart, itstop, itime, irestart !> Start, stop, current time step indices, flag for the restart case
  integer :: num, numx
  real(kind=8) :: time, time0               !> Time of simulation
  integer :: itscheme                !> Numerical scheme for time integration: 1-Euler, 2-AB2, 3-AB3
  real(kind=8) :: dt, cfl, cflx, cfly, cflz !> Time step size and CFL number
  real(kind=8), dimension(3) :: adt, bdt, cdt  !> Coefficients for the current and previous time steps in AB schemes

  ! Poisson resolution 
  real(kind=8) :: omega  !> relaxation parameter (SOR)
  real(kind=8) :: eps    !> Convergence criteria
  integer :: kmax        !> Maximum iterations

  ! Flow parameters
  real(kind=8) :: re, sc, u0, l0, u_ref, t_ref
  integer :: typesim ! 0=Read fields.bin file, 1=Convected vorte
  real(kind=8) :: ratio !> ratio bteween u1 and u2 for the planar jet simulation
  ! Boundary conditions for x and y directions
  integer :: nbcx1, nbcxn, nbcy1, nbcyn, nbcz1, nbczn  !> Boundary condition types for the start and end in x and y directions
  !> Inflow fileds
  real(kind=8), allocatable, dimension(:,:,:) :: inflow
  real(kind=8), allocatable, dimension(:) :: u_base 
  !real(kind=8), allocatable, dimension(:,:) :: u_noise

  ! Flag for simulation dimensionality
  integer :: sim2d  !> Flag to indicate if the simulation is 2D (1) or 3D (0)
  ! Test Table with the domain size
  real(kind=8), allocatable, dimension(:,:,:) :: tt(:,:,:)
  ! Flow variables
  real(kind=8), allocatable, dimension(:,:,:) :: ux, uy, uz  !> Velocity components
  real(kind=8), allocatable, dimension(:,:,:) :: phi, src  !> passive scalar field and source terme
  real(kind=8), allocatable, dimension(:,:,:) :: old_ux, old_uy, old_uz  !> Velocity components at itime-1
  real(kind=8), allocatable, dimension(:,:,:) :: pp          !> Pressure field
  real(kind=8), allocatable, dimension(:,:,:,:) :: fux, fuy, fuz  !> Flux arrays for velocity components
  real(kind=8), allocatable, dimension(:,:,:,:) :: fphi  !> Flux arrays for passive scalar
  real(kind=8), allocatable, dimension(:,:,:) :: ux_pred, uy_pred, uz_pred  !> Predicted velocities for the next time step
  real(kind=8), allocatable, dimension(:,:,:) :: rotx, roty, rotz !> Composantes of the rotationnal vector
  real(kind=8), allocatable, dimension(:,:,:) :: divu, divu_pred !> Composantes of the divergence velocity
  real(kind=8), allocatable, dimension(:,:,:) :: q_criterion !> Q criterion

  ! Flag for transport equation of scalar 
  integer :: nscr
  ! Save data parameters
  integer :: nfre, nsve 
  real(kind=8) :: xpro, ypro, zpro !> x, y, z absisses of the profil for visu
  integer :: ipro, jpro, kpro !> indices of the profil for visu

  ! Turbulence Inflow
  integer :: iin !> Inflow conditions (0: classic, 1: turbinit)
  real(kind=8) :: inflow_noise !> Turbulence intensity (1=100%) !! Inflow conditiona
  integer :: ici !> Initial conditions ((0: classic, 1: turbulent)
  real(kind=8) :: init_noise !> Turbulence intensity (1=100%) !! Inflow condition

  ! Parameters for x and y BC
  integer, parameter :: PERIODIC = 0
  integer, parameter :: FREE_SLIP = 1
  integer, parameter :: INFLOW_OUTFLOW = 2

  ! Variables for CPU time evaluation
  real(kind=8) :: go_time, start_time, end_time, sum_elapsed_time

  interface
     subroutine der_type(df, f, d)
       real(kind=8), intent(in) :: f(:,:,:), d
       real(kind=8), intent(out) :: df(:,:,:)
     end subroutine der_type
  end interface

  interface
     subroutine poi_type(pp, rhs, dx, dy, dz, nx, ny, nz, &
          omega, eps, kmax)
       real(kind=8), intent(inout) :: pp(:,:,:)
       real(kind=8), intent(in) :: rhs(:,:,:)
       real(kind=8), intent(in) :: dx, dy, dz
       real(kind=8), intent(in) :: omega, eps
       integer, intent(in) :: nx, ny, nz, kmax
     end subroutine poi_type
  end interface

  procedure(der_type), pointer :: derxp => null(), derxxp => null()
  procedure(der_type), pointer :: derxi => null(), derxxi => null()
  procedure(der_type), pointer :: deryp => null(), deryyp => null()
  procedure(der_type), pointer :: deryi => null(), deryyi => null()
  procedure(der_type), pointer :: derzp => null(), derzzp => null()
  procedure(der_type), pointer :: derzi => null(), derzzi => null()
  procedure(poi_type), pointer :: poisson_solver => null()

contains

  subroutine parameters()
    integer :: io_status
    ! Namelists
    namelist /Domain/ nx, ny, nz, xlx, yly, zlz, x0, y0, z0
    namelist /BoundaryConditions/ nbcx1, nbcxn, nbcy1, nbcyn, nbcz1, nbczn, sim2d
    namelist /AdvanceTime/ itscheme, cfl, itstart, itstop, irestart
    namelist /FlowParam/ u0, l0, re, typesim, ratio
    namelist /PoissonEq/ omega, eps, kmax
    namelist /Scalar/ nscr, sc
    namelist /VisuParameters/ nfre, nsve, xpro, ypro, zpro
    namelist /InitInflow/ iin, inflow_noise, ici, init_noise

    open(unit=10, file='paramaters.o3d', status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
       print *, 'Error opening configuration file.'
       stop
    end if

    ! Lire les namelists en fonction des sections du fichier de configuration
    read(10, FlowParam, iostat=io_status)
    if (io_status /= 0) then
       print *, 'Error reading parameters.'
       stop
    end if

    read(10, Domain, iostat=io_status)
    read(10, BoundaryConditions, iostat=io_status)
    read(10, AdvanceTime, iostat=io_status)
    read(10, PoissonEq, iostat=io_status)
    read(10, Scalar, iostat=io_status)
    read(10, VisuParameters, iostat=io_status)
    read(10, InitInflow, iostat=io_status)

    close(10)

  end subroutine parameters

  subroutine variables()
    integer :: i
    real(kind=8) :: dmin, cnu

    allocate(x(nx), y(ny), z(nz))
    allocate(ux(nx,ny,nz), uy(nx,ny,nz), uz(nx,ny,nz), pp(nx,ny,nz))
    allocate(phi(nx,ny,nz), src(nx,ny,nz))
    allocate(ux_pred(nx,ny,nz), uy_pred(nx,ny,nz), uz_pred(nx,ny,nz))
    allocate(old_ux(nx,ny,nz), old_uy(nx,ny,nz), old_uz(nx,ny,nz))
    allocate(fux(nx,ny,nz,3), fuy(nx,ny,nz,3), fuz(nx,ny,nz,3))
    allocate(fphi(nx,ny,nz,3))
    allocate(rotx(nx,ny,nz), roty(nx,ny,nz), rotz(nx,ny,nz))
    allocate(q_criterion(nx,ny,nz))
    allocate(divu(nx,ny,nz), divu_pred(nx,ny,nz))
    allocate(tt(nx,ny,nz))
    allocate(inflow(ny,nz,4), u_base(ny))

    ux = 0.d0
    uy = 0.d0
    uz = 0.d0
    pp = 0.d0
    phi = 0.d0
    fux = 0.d0
    fuy = 0.d0
    fuz = 0.d0
    fphi = 0.d0

    num = 0
    numx = 1
    dx = xlx / real(nx - 1, kind=8)
    dy = yly / real(ny - 1, kind=8)
    dz = zlz / real(nz - 1, kind=8)
    dmin = huge(1.d0)
    dmin = min(dmin, dx)
    dmin = min(dmin, dy)
    dmin = min(dmin, dz)
    dt = cfl * dmin / u0
    u_ref = u0
    t_ref = xlx / u0
    adt(1) = dt
    bdt(1) = 0.d0
    cdt(1) = 0.d0
    adt(2) =  3.d0 * dt / 2.d0
    bdt(2) = -1.d0 * dt / 2.d0
    cdt(2) =  0.d0
    adt(3) =  23.d0 * dt / 12.d0
    bdt(3) = -16.d0 * dt / 12.d0
    cdt(3) =   5.d0 * dt / 12.d0
    time = 0.d0

    cnu = 1.d0 / re * (dt / (dmin*dmin))
    if (cnu > cfl) then
       dt = 0.1d0 * dmin * dmin * re 
       print *, "Predominant viscous terms, cnu =", cnu
    end if

    do i = 1, nx
       x(i) = x0 + real(i-1, kind=8) * dx
    end do
    do i = 1, ny
       y(i) = y0 + real(i-1, kind=8) * dy
    end do
    do i = 1, nz
       z(i) = z0 + real(i-1, kind=8) * dz
    end do

    call set_profil_position(ipro, jpro, kpro, xpro, ypro, zpro, &
         x, y, z, nx, ny, nz)

  end subroutine variables

  subroutine schemes()

    if (nbcx1 == PERIODIC .and. nbcxn == PERIODIC) then
       derxp => derx_00
       derxxp => derxx_00
       derxi => derx_00
       derxxi => derxx_00
    elseif (nbcx1 == FREE_SLIP .and. nbcxn == FREE_SLIP) then
       derxp => derxp_11
       derxxp => derxxp_11
       derxi => derxi_11
       derxxi => derxxi_11
    else
       print *, "Unrecognized boundarx laxer txpes for Y-axis."
       print *, "nbcx1: ", nbcx1, "nbcxn: ", nbcxn
       stop
    end if

    if (nbcy1 == PERIODIC .and. nbcyn == PERIODIC) then
       deryp => dery_00
       deryyp => deryy_00
       deryi => dery_00
       deryyi => deryy_00
    elseif (nbcy1 == FREE_SLIP .and. nbcyn == FREE_SLIP) then
       deryp => deryp_11
       deryyp => deryyp_11
       deryi => deryi_11
       deryyi => deryyi_11
    else
       print *, "Unrecognized boundary layer types for Y-axis."
       print *, "nbcy1: ", nbcy1, "nbcyn: ", nbcyn
       stop
    end if

    if (sim2d == 0) then
       if (nbcz1 == PERIODIC .and. nbczn == PERIODIC) then
          derzp => derz_00
          derzzp => derzz_00
          derzi => derz_00
          derzzi => derzz_00
       elseif (nbcz1 == FREE_SLIP .and. nbczn == FREE_SLIP) then
          derzp => derzp_11
          derzzp => derzzp_11
          derzi => derzi_11
          derzzi => derzzi_11
       else
          print *, "Unrecognized boundarz lazer tzpes for Y-axis."
          print *, "nbcz1: ", nbcz1, "nbczn: ", nbczn
          stop
       end if
    elseif (sim2d == 1) then
       derzi => derz_2dsim
       derzzi => derzz_2dsim
       derzp => derz_2dsim
       derzzp => derzz_2dsim
    end if

    if (nbcx1 == PERIODIC .and. nbcxn == PERIODIC .and. &
         nbcy1 == PERIODIC .and. nbcyn == PERIODIC) then
       poisson_solver => poisson_solver_0000
    else if (nbcx1 == PERIODIC .and. nbcxn == PERIODIC .and. &
         nbcy1 == FREE_SLIP .and. nbcyn == FREE_SLIP) then
       poisson_solver => poisson_solver_0011
    else if (nbcx1 == FREE_SLIP .and. nbcxn == FREE_SLIP .and. &
         nbcy1 == FREE_SLIP .and. nbcyn == FREE_SLIP) then
       poisson_solver => poisson_solver_111111
    else if (nbcx1 == INFLOW_OUTFLOW .and. nbcxn == INFLOW_OUTFLOW .and. &
         nbcy1 == PERIODIC .and. nbcyn == PERIODIC) then
       poisson_solver => poisson_solver_2200
    else if (nbcx1 == INFLOW_OUTFLOW .and. nbcxn == INFLOW_OUTFLOW .and. &
         nbcy1 == FREE_SLIP .and. nbcyn == FREE_SLIP) then
       poisson_solver => poisson_solver_2211
    end if

    return 
  end subroutine schemes

  subroutine set_profil_position(ipro, jpro, kpro, xpro, ypro, zpro, &
       x, y, z, nx, ny, nz)
    integer, intent(out) :: ipro, jpro, kpro
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: x(:), y(:), z(:)
    real(kind=8), intent(in) :: xpro, ypro, zpro

    integer :: i, j, k

    i = 1
    j = 1
    k = 1
    do while (x(i) < xpro .and. i <= nx)
       i = i+1
    end do
    ipro = i
    do while (y(j) < ypro .and. j <= ny)
       j = j+1
    end do
    jpro = j
    do while (z(k) < zpro .and. k <= nz)
       k = k+1
    end do
    kpro = k

    return
  end subroutine set_profil_position

end module initialization
