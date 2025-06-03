module boundary_conditions
  use derivation
  use functions
  implicit none

contains

  subroutine pressure_neuman_bc_x0(pp)
    !> Apply Neumann boundary conditions in the x-direction: dp(1)/dx = 0 
    !>
    !> This subroutine sets the Neumann boundary conditions for the pressure
    !> field pp at the x-boundaries. It ensures that the gradient of the
    !> pressure in the x-direction at the boundaries is zero, which corresponds
    !> to a no-flux condition.
    !>
    !> INPUT:
    !> OUTPUT:
    !> pp(nx,ny,nz): Pressure field with updated boundary values
    !>               in the x-direction

    real(kind=8), intent(inout) :: pp(:,:,:)

    real(kind=8) :: ust, qst

    ! Coefficients for the Neumann boundary conditions (O(dx^2) discretization)
    ust = -1.d0 / 3.d0
    qst =  4.d0 / 3.d0

    ! Apply Neumann BC at the first x-boundary
    pp(1 , :, :) = ust * pp(3   , :, :) + qst * pp(2   , :, :)

    return
  end subroutine pressure_neuman_bc_x0

  subroutine pressure_neuman_bc_xn(pp)
    !> Apply Neumann boundary conditions in the x-direction: dp(nx)/dx=0
    !>
    !> This subroutine sets the Neumann boundary conditions for the pressure
    !> field pp at the x-boundaries. It ensures that the gradient of the
    !> pressure in the x-direction at the boundaries is zero, which corresponds
    !> to a no-flux condition.
    !>
    !> INPUT:
    !> OUTPUT:
    !> pp(nx,ny,nz): Pressure field with updated boundary values
    !>               in the x-direction

    real(kind=8), intent(inout) :: pp(:,:,:)

    real(kind=8) :: ust, qst
    integer :: nx

    nx = size(pp, 1)

    ! Coefficients for the Neumann boundary conditions (O(dx^2) discretization)
    ust = -1.d0 / 3.d0
    qst =  4.d0 / 3.d0

    ! Apply Neumann BC at the first x-boundary
    pp(nx , :, :) = ust * pp(nx-2   , :, :) + qst * pp(nx-1   , :, :)

    return
  end subroutine pressure_neuman_bc_xn

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
       u0, dx, dt, met)
    !> Apply non-reflective Dirichlet boundary conditions for velocity
    !> components at the outflow boundary (x = nx).
    !>
    !> This subroutine applies the 1D convective Orlanski-type outflow condition
    !> to reduce reflections at the boundary.
    !>
    !> INPUT:
    !> ux(nx,ny,nz) : X-component of the velocity field.
    !> uy(nx,ny,nz) : Y-component of the velocity field.
    !> uz(nx,ny,nz) : Z-component of the velocity field.
    !> u0           : Reference velocity (used in fixed or target damping).
    !> dx           : Grid spacing in the x-direction.
    !> dt           : Time step.
    !> met          : Method for estimating the convective velocity:
    !>                1 = constant reference velocity (u0),
    !>                2 = mean velocity at boundary,
    !>                3 = local velocity at nx-1.
    !>
    !> OUTPUT:
    !> bux(:,:)     : X-component of the velocity field with updated boundary values.
    !> buy(:,:)     : Y-component of the velocity field with updated boundary values.
    !> buz(:,:)     : Z-component of the velocity field with updated boundary values.

    real(kind=8), intent(inout) :: bux(:,:), buy(:,:), buz(:,:)
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: dx, dt, u0
    integer, intent(in) :: met

    integer :: nx, ny, nz, j, k
    real(kind=8) :: cx, u_mean, udx

    nx = size(ux,1)
    ny = size(uy,2)
    nz = size(uz,3)
    udx = 1.d0 / dx

    do k = 1, nz
       do j = 1, ny
             ! Estimate convective velocity coefficient cx depending on method
             if (met == 1) then
                cx = u0 * dt * udx
             else if (met == 2) then
                u_mean = 0.5d0 * (ux(nx-2,j,k) + ux(nx-1,j,k))
                cx = u_mean * dt * udx
             else if (met == 3) then
                cx = ux(nx-1,j,k) * dt * udx
             else
                cx = u0 * dt * udx
             end if

             ! Orlanski convective update (1st-order)
             bux(j,k) = ux(nx,j,k) - cx * (ux(nx,j,k) - ux(nx-1,j,k))
             buy(j,k) = uy(nx,j,k) - cx * (uy(nx,j,k) - uy(nx-1,j,k))
             buz(j,k) = uz(nx,j,k) - cx * (uz(nx,j,k) - uz(nx-1,j,k))
       end do
    end do

    return
  end subroutine velocity_dirichlet_bc_xn

  subroutine inflow_turbulent(u_noise, u_base, inflow_noise, ny, nz)
    integer, intent(in) :: ny, nz
    real(kind=8), dimension(ny, nz), intent(out) :: u_noise
    real(kind=8), dimension(ny), intent(in) :: u_base
    real(kind=8), intent(in) :: inflow_noise
    real(kind=8), dimension(1, ny, nz) :: pink_noise
    integer :: j, k

    call generate_pink_noise(1, ny, nz, pink_noise, 256)
    do k = 1, nz
       do j = 1, ny
          u_noise(j, k) = inflow_noise * pink_noise(1, j, k) * u_base(j)
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
    ! call random_seed(size = n_seed)
    n_seed = random_between(256, 512)
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

  subroutine generate_pink_noise(nx, ny, nz, noise, n_rand)
    integer, intent(in) :: nx, ny, nz, n_rand
    real(kind=8), intent(out) :: noise(nx, ny, nz)
    integer :: i, j, k, l
    real(kind=8), dimension(:), allocatable :: white_noise
    real(kind=8), dimension(nx, ny, nz) :: filtered_noise
    integer, parameter :: num_filters = 128
    real(kind=8), dimension(num_filters) :: filters
    real(kind=8) :: sum_filters, noise_val
    integer :: seed(8)
    integer :: ind
    real(kind=8) :: min_val, max_val, sca
    real(kind=8), parameter :: tol = 1.0e-10

    ! Check input dimensions
    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
       print *, "Error: Dimensions must be positive."
       return
    end if

    ! Initialize white noise
    allocate(white_noise(nx*ny*nz))
    ! Initialize random seed based on the system clock
    call system_clock(count=seed(1), count_rate=seed(2), count_max=seed(3))
    seed = [(mod(seed(1) + i + n_rand, seed(3)), i = 1, 8)]
    call random_seed(put=seed)
    call random_number(white_noise)

    ! Initialize filters
    filters = [(1.0 / (i + 1)**2, i = 1, num_filters)]

    ! Apply Voss-McCartney filter to the white noise
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ind = i + (j - 1) * nx + (k - 1) * nx * ny
             noise_val = 0.0
             sum_filters = 0.0
             do l = 1, num_filters
                sum_filters = sum_filters + filters(l)
                noise_val = noise_val + white_noise(ind) * filters(l)
             end do
             filtered_noise(i, j, k) = noise_val / sum_filters
          end do
       end do
    end do

    ! Normalize the filtered noise to be between -1 and 1
    min_val = minval(filtered_noise)
    max_val = maxval(filtered_noise)
    if (abs(max_val - min_val) > tol) then
       sca = 2.0 / (max_val - min_val)
       noise = (filtered_noise - min_val) * sca - 1.0
    else
       noise = filtered_noise  ! No normalization needed if max_val == min_val
    end if

    ! Clean up
    deallocate(white_noise)

  end subroutine generate_pink_noise

end module boundary_conditions
