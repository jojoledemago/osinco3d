module boundary_conditions
  use derivation
  use functions
  implicit none

contains

  subroutine pressure_neuman_bc_x(pp, nx)
    !> Apply Neumann boundary conditions in the x-direction: dp(1)/dx = 0 and
    !> dp(nx)/dx=0
    !>
    !> This subroutine sets the Neumann boundary conditions for the pressure
    !> field pp at the x-boundaries. It ensures that the gradient of the
    !> pressure in the x-direction at the boundaries is zero, which corresponds
    !> to a no-flux condition.
    !>
    !> INPUT:
    !> nx          : Number of grid points in the x-direction
    !> OUTPUT:
    !> pp(nx,ny,nz): Pressure field with updated boundary values
    !>               in the x-direction

    integer, intent(in) :: nx
    real(kind=8), intent(inout) :: pp(:,:,:)

    real(kind=8) :: ust, qst

    ! Coefficients for the Neumann boundary conditions (O(dx^2) discretization)
    ust = -1.d0 / 3.d0
    qst =  4.d0 / 3.d0

    ! Apply Neumann BC at the first x-boundary
    pp(1 , :, :) = ust * pp(3   , :, :) + qst * pp(2   , :, :)

    ! Apply Neumann BC at the last x-boundary
    pp(nx, :, :) = ust * pp(nx-2, :, :) + qst * pp(nx-1, :, :)

    return
  end subroutine pressure_neuman_bc_x

  subroutine velocity_dirichlet_bc_x0(ux, uy, uz, inflow)
    !> Apply Dirichlet boundary conditions for velocity components at x=0.
    !>
    !> This subroutine sets the Dirichlet boundary conditions for the velocity 
    !> components (ux, uy, uz) at the x=0 boundary. The velocity components are
    !> specified based on the inflow conditions provided.
    !>
    !> INPUT:
    !> ux(nx,ny,nz)    : X-component of velocity field
    !> uy(nx,ny,nz)    : Y-component of velocity field
    !> uz(nx,ny,nz)    : Z-component of velocity field
    !> inflow(nx,ny,4) : Array containing the inflow conditions extracted from 
    !>                   the velocity and pressure fields. The last dimension 
    !>                   contains the following components:
    !>                   - inflow(:,:,1) : X-component of velocity at 
    !>                   the inflow boundary
    !>                   - inflow(:,:,2) : Y-component of velocity at 
    !>                   the inflow boundary
    !>                   - inflow(:,:,3) : Z-component of velocity at 
    !>                   the inflow boundary
    !> iin             : Flag for the noise at the inlet boundary
    !>
    !> OUTPUT:
    !> ux(:,:,:)     : X-component of velocity field with updated boundary 
    !>                 values 
    !> uy(:,:,:)     : Y-component of velocity field with updated boundary 
    !>                 values 
    !> uz(:,:,:)     : Z-component of velocity field with updated boundary 
    !>                 values 

    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: inflow(:,:,:)

    ux(1,:,:) = inflow(:,:,1) 
    uy(1,:,:) = inflow(:,:,2)
    uz(1,:,:) = inflow(:,:,3)

    return
  end subroutine velocity_dirichlet_bc_x0

  subroutine velocity_neuman_bc_x0(ux, uy, uz)
    !> Apply Neumann boundary conditions for velocity components at the
    !> inflow boundary (x = 0).
    !>
    !> This subroutine sets the Neumann boundary conditions for the velocity
    !> components (ux, uy, uz) at the x = nx boundary. The Neumann condition
    !> specifies that the derivative of the velocity components with respect
    !> to x is zero (du/dx = 0), using a second-order backward scheme.
    !>
    !> INPUT:
    !> ux(nx,ny,nz) : X-component of the velocity field
    !> uy(nx,ny,nz) : Y-component of the velocity field
    !> uz(nx,ny,nz) : Z-component of the velocity field
    !>
    !> OUTPUT:
    !> ux(:,:,:)    : X-component of the velocity field with updated boundary
    !>                values
    !> uy(:,:,:)    : Y-component of the velocity field with updated boundary
    !>                values
    !> uz(:,:,:)    : Z-component of the velocity field with updated boundary
    !>                values

    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8) :: ust, qst
    integer :: nx

    nx = size(ux, 1)
    ! Coefficients for the Neumann boundary conditions using a second-order 
    ! backward scheme
    ust = 1.d0 / 3.d0 
    qst = 4.d0 / 3.d0
    ! Apply Neumann boundary conditions at the outflow boundary
    ux(1,:,:) = - ust * ux(3,:,:) + qst * ux(2,:,:)
    uy(1,:,:) = - ust * uy(3,:,:) + qst * uy(2,:,:)
    uz(1,:,:) = - ust * uz(3,:,:) + qst * uz(2,:,:)

  end subroutine velocity_neuman_bc_x0

  subroutine add_u_noise(ux, uy, uz, inflow_noise, u_base)
    ! Add noise to the velocity components at the inlet boundary
    ! This subroutine modifies the inlet boundary values of the velocity
    ! components by adding turbulent noise. The noise is generated and scaled
    ! based on the given inflow noise intensity and the velocity profile.
    ! 
    ! INPUT:
    ! - inflow_noise   : Intensity of the noise at the inlet boundary (input)
    ! - u_base         : Velocity profile between 0 and 1 to apply noise at the 
    !                    maximum velocity gradient (input)
    ! OUTPUT:
    ! - ux, uy, uz     : Velocity components (input/output)


    ! Declaration of arguments
    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: inflow_noise, u_base(:)

    ! Local variables
    integer :: ny, nz
    real(kind=8), dimension(:, :), allocatable :: u_noise

    ! Determine the size of the second and third dimensions
    ny = size(ux, 2)
    nz = size(ux, 3)

    ! Allocate array for the noise
    allocate(u_noise(ny, nz))

    ! Add noise to the x-component of the velocity
    call inflow_turbulent(u_noise, u_base, inflow_noise, ny, nz)
    ux(1,:,:) = ux(1,:,:) + u_noise

    ! Add noise to the y-component of the velocity
    call inflow_turbulent(u_noise, u_base, inflow_noise, ny, nz)
    uy(1,:,:) = uy(1,:,:) + u_noise

    ! Add noise to the z-component of the velocity
    call inflow_turbulent(u_noise, u_base, inflow_noise, ny, nz)
    uz(1,:,:) = uz(1,:,:) + u_noise

    ! Deallocate the noise array
    deallocate(u_noise)

    return
  end subroutine add_u_noise


  subroutine velocity_neuman_bc_xn(ux, uy, uz)
    !> Apply Neumann boundary conditions for velocity components at the
    !> outflow boundary (x = nx).
    !>
    !> This subroutine sets the Neumann boundary conditions for the velocity
    !> components (ux, uy, uz) at the x = nx boundary. The Neumann condition
    !> specifies that the derivative of the velocity components with respect
    !> to x is zero (du/dx = 0), using a second-order backward scheme.
    !>
    !> INPUT:
    !> ux(nx,ny,nz) : X-component of the velocity field
    !> uy(nx,ny,nz) : Y-component of the velocity field
    !> uz(nx,ny,nz) : Z-component of the velocity field
    !>
    !> OUTPUT:
    !> ux(:,:,:)    : X-component of the velocity field with updated boundary
    !>                values
    !> uy(:,:,:)    : Y-component of the velocity field with updated boundary
    !>                values
    !> uz(:,:,:)    : Z-component of the velocity field with updated boundary
    !>                values

    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8) :: ust, qst
    integer :: nx

    nx = size(ux, 1)
    ! Coefficients for the Neumann boundary conditions using a second-order 
    ! backward scheme
    ust = 1.d0 / 3.d0 
    qst = 4.d0 / 3.d0
    ! Apply Neumann boundary conditions at the outflow boundary
    ux(nx,:,:) = - ust * ux(nx-2,:,:) + qst * ux(nx-1,:,:)
    uy(nx,:,:) = - ust * uy(nx-2,:,:) + qst * uy(nx-1,:,:)
    uz(nx,:,:) = - ust * uz(nx-2,:,:) + qst * uz(nx-1,:,:)

  end subroutine velocity_neuman_bc_xn

  subroutine velocity_dirichlet_bc_xn(bux, buy, buz, ux, uy, uz, &
       u0, dx, dy, dt, met)
    !> Apply Dirichlet boundary conditions for velocity components at the
    !> outflow boundary (x = nx).
    !> 
    !> Outflow non-reflective computation by resolving 1D convective equation.
    !> 
    !> This subroutine updates the boundary values for the velocity components
    !> (ux, uy, uz) at the outflow boundary (x = nx) using a Dirichlet boundary
    !> condition. The outflow condition is implemented by resolving the 1D
    !> convective equation to ensure non-reflective behavior at the boundary.
    !> 
    !> INPUT:
    !> ux(nx,ny,nz) : X-component of the velocity field.
    !> uy(nx,ny,nz) : Y-component of the velocity field.
    !> uz(nx,ny,nz) : Z-component of the velocity field.    
    !> u0           : Reference veolcity 
    !> dx           : Mesh size step in the x-direction.
    !> dy           : Mesh size step in the y-direction.
    !> dt           : Time step size.
    !> met          : Methode to evaluate convection coefficient cx
    !>
    !> OUTPUT:
    !> bux(:,:)   : X-component of the velocity field with updated boundary 
    !> values.
    !> buy(:,:)   : Y-component of the velocity field with updated boundary 
    !> values.
    !> buz(:,:)   : Z-component of the velocity field with updated boundary 
    !> values.

    real(kind=8), intent(inout) :: bux(:,:), buy(:,:), buz(:,:)
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:), dx, dy, dt, u0
    integer, intent(in) :: met

    real(kind=8) :: udx, cx, ux_interp, u_mean, distance_dx, distance_dxy
    integer :: nx, ny, nz, j, k

    nx = size(ux, 1)
    ny = size(uy, 2)
    nz = size(uz, 3)
    udx = 1.d0 / dx

    if (met == 1) then
       ! ! Use the reference value of velocity
       cx = u0 * dt * udx
       bux = ux(nx, :, :) - cx * (ux(nx, :, :) - ux(nx-1, :, :))
       buy = uy(nx, :, :) - cx * (uy(nx, :, :) - uy(nx-1, :, :))
       buz = uz(nx, :, :) - cx * (uz(nx, :, :) - uz(nx-1, :, :))
    else if (met == 2) then
       ! Use the mean value of velocity at the boundary
       u_mean = average_2d_array(ux(nx-1, :, :))
       cx = u_mean * dt * udx
       bux = ux(nx, :, :) - cx * (ux(nx, :, :) - ux(nx-1, :, :))
       buy = uy(nx, :, :) - cx * (uy(nx, :, :) - uy(nx-1, :, :))
       buz = uz(nx, :, :) - cx * (uz(nx, :, :) - uz(nx-1, :, :))
    else if (met == 3) then
       ! Use the discret value of velocity at the boundary (nx-1)
       do k = 1, nz
          do j = 1, ny
             cx = ux(nx, j, k) * dt * udx
             bux(j, k) = ux(nx, j, k) - cx * (ux(nx, j, k) - ux(nx-1, j, k))
             buy(j, k) = uy(nx, j, k) - cx * (uy(nx, j, k) - uy(nx-1, j, k))
             buz(j, k) = uz(nx, j, k) - cx * (uz(nx, j, k) - uz(nx-1, j, k))
          end do
       end do
    else if (met == 4) then
       ! Interpolate the velocity field at the boundary
       distance_dx = dx
       distance_dxy = sqrt(dx**2 + dy**2)
       do k = 1, nz
          j = 1
          ux_interp = ux(nx-1, j, k)
          cx = ux_interp * dt * udx
          bux(j, k) = ux(nx, j, k) - cx * (ux(nx, j, k) - ux(nx-1, j, k))
          buy(j, k) = uy(nx, j, k) - cx * (uy(nx, j, k) - uy(nx-1, j, k))
          buz(j, k) = uz(nx, j, k) - cx * (uz(nx, j, k) - uz(nx-1, j, k))
          do j = 2, ny-1
             ux_interp = (ux(nx-1, j, k) / distance_dx + &
                  ux(nx-1, j-1, k) / distance_dxy + &
                  ux(nx-1, j+1, k) / distance_dxy) / &
                  (1.0 / distance_dx + 2.0 / distance_dxy)
             cx = ux_interp * dt * udx
             bux(j, k) = ux(nx, j, k) - cx * (ux(nx, j, k) - ux(nx-1, j, k))
             buy(j, k) = uy(nx, j, k) - cx * (uy(nx, j, k) - uy(nx-1, j, k))
             buz(j, k) = uz(nx, j, k) - cx * (uz(nx, j, k) - uz(nx-1, j, k))
          end do
          j = ny
          ux_interp = ux(nx-1, j, k)
          cx = ux_interp * dt * udx
          bux(j, k) = ux(nx, j, k) - cx * (ux(nx, j, k) - ux(nx-1, j, k))
          buy(j, k) = uy(nx, j, k) - cx * (uy(nx, j, k) - uy(nx-1, j, k))
          buz(j, k) = uz(nx, j, k) - cx * (uz(nx, j, k) - uz(nx-1, j, k))
       end do
    end if
    return

  end subroutine velocity_dirichlet_bc_xn

  subroutine inflow_turbulent(u_noise, u_base, inflow_noise, ny, nz)
    integer, intent(in) :: ny, nz
    real(kind=8), dimension(ny, nz), intent(out) :: u_noise
    real(kind=8), dimension(ny), intent(in) :: u_base
    real(kind=8), intent(in) :: inflow_noise
    real(kind=8), dimension(ny, nz) :: modulated_noise, white_noise
    integer :: j, k


    call generate_white_noise(white_noise, ny, nz)
    call apply_energy_spectrum(modulated_noise, white_noise, ny, nz)

    do k = 1, nz
       do j = 1, ny
          u_noise(j, k) = inflow_noise * modulated_noise(j, k) * u_base(j)
       end do
    end do
  end subroutine inflow_turbulent

  subroutine apply_energy_spectrum(modulated_noise, white_noise, n1, n2)
    integer, intent(in) :: n1, n2
    real(kind=8), dimension(n1, n2), intent(in) :: white_noise
    real(kind=8), dimension(n1, n2), intent(out) :: modulated_noise
    integer :: i, j
    real(kind=8) :: kx, ky, energy_spectrum, tiny_value, max_value

    tiny_value = tiny(1.d0)
    modulated_noise = 0.d0

    do j = 1, n2
       ky = real(j-1, kind=8)
       if (j > n2/2) ky = ky - n2  ! Gérer les fréquences négatives
       do i = 1, n1
          kx = real(i-1, kind=8)
          if (i > n1/2) kx = kx - n1  ! Gérer les fréquences négatives
          ! Éviter la division par zéro à l'origine des fréquences
          if (sqrt(kx**2 + ky**2) < tiny_value) then
             energy_spectrum = 0.d0
          else
             ! E(k) ~ k^(-5/3) (Spectre de Kolmogorov)
             energy_spectrum = sqrt(kx**2 + ky**2)**(-5.0d0/6.0d0)
          end if

          modulated_noise(i, j) = white_noise(i, j) * energy_spectrum
       end do
    end do

    max_value = maxval(abs(modulated_noise))
    if (max_value > 0.d0) then
       modulated_noise = 28.d0 * modulated_noise / max_value
    end if
    return
  end subroutine apply_energy_spectrum

  subroutine generate_white_noise(noise, n1, n2)
    integer, intent(in) :: n1, n2
    real(kind=8), dimension(n1, n2), intent(out) :: noise
    integer :: seed, i
    integer, dimension(:), allocatable :: seed_array
    integer :: n_seed

    ! Use the system clock to get a unique seed
    call random_seed(size = n_seed)
    allocate(seed_array(n_seed))
    call system_clock(count=seed)

    ! Initialize seed array with different values
    do i = 1, n_seed
       seed_array(i) = seed + i
    end do
    ! Set the random seed
    call random_seed(put=seed_array)
    call random_number(noise)

    noise = 2.0d0 * noise - 1.0d0  ! Scale to range [-1, 1]
    return
  end subroutine generate_white_noise

end module boundary_conditions
