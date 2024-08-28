module integration
  use functions
  use poisson
  use derivation
  use initialization, only: inflow
  use diffoper
  use boundary_conditions
  use IOfunctions, only: print_in_outflow_rate, write_velocity_diverged 
  implicit none

contains

  subroutine predict_velocity(ux_pred, uy_pred, uz_pred, ux, uy, uz, &
       fux, fuy, fuz, re, adt, bdt, cdt, itime, itscheme, inflow, &
       dx, dy, dz, nx, ny, nz)

    !> Predicts the velocity field for a fluid dynamics simulation using various time integration schemes.
    !>
    !> INPUT: 
    !> ux(nx, ny, nz)    : Current x-component of velocity
    !> uy(nx, ny, nz)    : Current y-component of velocity
    !> uz(nx, ny, nz)    : Current z-component of velocity
    !> fux(nx, ny, nz, 3): Storage for velocity ux fluxes or calculations across time steps
    !> fuy(nx, ny, nz, 3): Storage for velocity uy fluxes or calculations across time steps
    !> fuz(nx, ny, nz, 3): Storage for velocity uz fluxes or calculations across time steps
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
    !>
    !> OUTPUT: 
    !> ux_pred(nx, ny, nz): Predicted x-component of velocity for next time step
    !> uy_pred(nx, ny, nz): Predicted y-component of velocity for next time step
    !> uz_pred(nx, ny, nz): Predicted z-component of velocity for next time step

    integer, intent(in) :: itscheme, itime
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: re
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(3), intent(in) :: adt, bdt, cdt
    real(kind=8), dimension(:,:,:), intent(inout) :: ux, uy, uz
    real(kind=8), dimension(:,:,:), intent(in) :: inflow
    real(kind=8), dimension(:,:,:,:), intent(inout) :: fux, fuy, fuz
    real(kind=8), dimension(:,:,:), intent(out) :: ux_pred, uy_pred, uz_pred


    real(kind=8) :: onere, adu, bdu, cdu
    real(kind=8), allocatable, dimension(:,:,:) :: duxdx, duxdy, duxdz
    real(kind=8), allocatable, dimension(:,:,:) :: duydx, duydy, duydz
    real(kind=8), allocatable, dimension(:,:,:) :: duzdx, duzdy, duzdz
    real(kind=8), allocatable, dimension(:,:,:) :: d2uxdx2, d2uxdy2, d2uxdz2
    real(kind=8), allocatable, dimension(:,:,:) :: d2uydx2, d2uydy2, d2uydz2
    real(kind=8), allocatable, dimension(:,:,:) :: d2uzdx2, d2uzdy2, d2uzdz2

    print *, "* Predict velocity"

    ! Allocation de mémoire pour les dérivées
    allocate(duxdx(nx,ny,nz), duxdy(nx,ny,nz), duxdz(nx,ny,nz))
    allocate(duydx(nx,ny,nz), duydy(nx,ny,nz), duydz(nx,ny,nz))
    allocate(duzdx(nx,ny,nz), duzdy(nx,ny,nz), duzdz(nx,ny,nz))
    allocate(d2uxdx2(nx,ny,nz), d2uxdy2(nx,ny,nz), d2uxdz2(nx,ny,nz))
    allocate(d2uydx2(nx,ny,nz), d2uydy2(nx,ny,nz), d2uydz2(nx,ny,nz))
    allocate(d2uzdx2(nx,ny,nz), d2uzdy2(nx,ny,nz), d2uzdz2(nx,ny,nz))

    if (itscheme == 1 .or. itime == 1) then
       print *, " Euler"
       adu = adt(1)
       bdu = bdt(1)
       cdu = cdt(1)
    elseif (itscheme == 2 .or. itime == 2) then
       print *, " Adams-Bashforth 2"
       adu = adt(2)
       bdu = bdt(2)
       cdu = cdt(2)
    elseif (itscheme == 3) then
       print *, " Adams-Bashforth 3"
       adu = adt(3)
       bdu = bdt(3)
       cdu = cdt(3)
    else
       adu = 0.d0
       bdu = 0.d0
       cdu = 0.d0
       print *, 'itscheme:', itscheme, ' unrecognized'
       stop 
    end if
    onere = 1.d0 / re

    ! Compute partial first derivatives of ux
    ! with respect to x, y, z axis
    call derxi(duxdx, ux, dx)
    call deryp(duxdy, ux, dy)
    call derzp(duxdz, ux, dz)
    ! Compute of the partial first derivatives of ux
    ! with respect to x, y, z axis
    call derxxi(d2uxdx2, ux, dx)
    call deryyp(d2uxdy2, ux, dy)
    call derzzp(d2uxdz2, ux, dz)

    ! Update fux for the new time step
    ! using viscosity and convective termes
    fux(:, :, :, 1) = onere * (d2uxdx2 + d2uxdy2 + d2uxdz2) - &
         (ux * duxdx + uy * duxdy + uz * duxdz)
    ! Compute predicted velocities ux using the selected time integration scheme
    ux_pred = ux + adu * fux(:, :, :, 1) + &
         bdu * fux(:, :, :, 2) + &
         cdu * fux(:, :, :, 3)

    ! Compute partial first derivatives of uy
    ! with respect to x, y, z axis
    call derxp(duydx, uy, dx)
    call deryi(duydy, uy, dy)
    call derzp(duydz, uy, dz)
    ! Compute of the partial first derivatives of uy
    ! with respect to x, y, z axis
    call derxxp(d2uydx2, uy, dx)
    call deryyi(d2uydy2, uy, dy)
    call derzzp(d2uydz2, uy, dz)

    ! Update fuy for the new time step
    ! using viscosity and convective termes
    fuy(:, :, :, 1) = onere * (d2uydx2 + d2uydy2 + d2uydz2) - &
         (ux * duydx + uy * duydy + uz * duydz)
    ! Compute predicted velocities uy using the selected time integration scheme
    uy_pred = uy + adu * fuy(:, :, :, 1) + &
         bdu * fuy(:, :, :, 2) + &
         cdu * fuy(:, :, :, 3)

    ! Compute partial first derivatives of uz
    ! with respect to x, y, z axis
    call derxp(duzdx, uz, dx)
    call deryp(duzdy, uz, dy)
    call derzi(duzdz, uz, dz)
    ! Compute of the partial first derivatives of uz
    ! with respect to x, y, z axis
    call derxxp(d2uzdx2, uz, dx)
    call deryyp(d2uzdy2, uz, dy)
    call derzzi(d2uzdz2, uz, dz)

    ! Update fuz for the new time step
    ! using viscosity and convective termes
    fuz(:, :, :, 1) = onere * (d2uzdx2 + d2uzdy2 + d2uzdz2) - &
         (ux * duzdx + uy * duzdy + uz * duzdz)
    ! Compute predicted velocities uz using the selected time integration scheme
    uz_pred = uz + adu * fuz(:, :, :, 1) + &
         bdu * fuz(:, :, :, 2) + &
         cdu * fuz(:, :, :, 3)

    ! Case inflow/outflow x-boundary-condition for u_pred: apply neuman inflow 
    ! for x0 and neuman condition for xmax
    if (nbcx1 == INFLOW_OUTFLOW) then
       !call velocity_neuman_bc_x0(ux_pred, uy_pred, uz_pred)
       call velocity_dirichlet_bc_x0(ux, uy, uz, inflow)
    end if
    if (nbcxn == INFLOW_OUTFLOW) then
       call velocity_neuman_bc_xn(ux_pred, uy_pred, uz_pred)
    end if

    select case(itscheme)
    case(2)
       fux(:,:,:,2) = fux(:,:,:,1)
       fuy(:,:,:,2) = fuy(:,:,:,1)
       fuz(:,:,:,2) = fuz(:,:,:,1)
    case(3)
       fux(:,:,:,3) = fux(:,:,:,2)
       fuy(:,:,:,3) = fuy(:,:,:,2)
       fuz(:,:,:,3) = fuz(:,:,:,2)
       fux(:,:,:,2) = fux(:,:,:,1)
       fuy(:,:,:,2) = fuy(:,:,:,1)
       fuz(:,:,:,2) = fuz(:,:,:,1)
    end select

    deallocate(duxdx, duxdy, duxdz, d2uxdx2, d2uxdy2, d2uxdz2)
    deallocate(duydx, duydy, duydz, d2uydx2, d2uydy2, d2uydz2)
    deallocate(duzdx, duzdy, duzdz, d2uzdx2, d2uzdy2, d2uzdz2)
    print *, ""

    return
  end subroutine predict_velocity

  subroutine correct_pression(pp, ux_pred, uy_pred, uz_pred, dx, dy, dz, &
       nx, ny, nz, dt, omega, eps, kmax)

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
    !> omega                : Relaxation factor for SOR solver
    !> eps                  : Convergence criterion for SOR solver
    !> kmax                 : Maximum number of iterations for SOR solver
    !>
    !> OUTPUT:
    !> pp(nx, ny, nz)       : Updated pressure field after convergence

    real(kind=8), intent(inout) :: pp(:,:,:)
    real(kind=8), intent(in) :: ux_pred(:,:,:), uy_pred(:,:,:), uz_pred(:,:,:)
    real(kind=8), intent(in) :: dx, dy, dz, dt, omega, eps
    integer, intent(in) :: nx, ny, nz, kmax

    real(kind=8), allocatable :: divu_pred(:,:,:), rhs(:,:,:)

    print *, "* Correction pression"
    ! Allocate memory for the divergence and RHS arrays
    allocate(divu_pred(nx, ny, nz))
    allocate(rhs(nx, ny, nz))

    ! Compute the divergence of the predicted velocity field
    call divergence(divu_pred, ux_pred, uy_pred, uz_pred, &
         dx, dy, dz, nx, ny, nz, 1)  ! Assume 1 for odd conditions

    ! Prepare the right-hand side term for the Poisson equation
    rhs = divu_pred / dt
    if (nbcx1 == INFLOW_OUTFLOW) then
       call pressure_neuman_bc_x(rhs, nx)
    end if

    ! Solve the Poisson equation to update the pressure field
    call poisson_solver(pp, rhs, dx, dy, dz, nx, ny, nz, &
         omega, eps, kmax)

    ! Free allocated memory
    deallocate(divu_pred, rhs)
    print *,""

    return
  end subroutine correct_pression

  subroutine correct_velocity(ux, uy, uz, ux_pred, uy_pred, uz_pred, &
       pp, dt, dx, dy, dz, nx, ny, nz)
    !> Correct the velocity components based on the pressure gradient.
    !> This step is part of the fractional step method to solve
    !> incompressible Navier-Stokes equations.
    !> Periodic boundary conditions are applied in the x-direction
    !> to handle boundary cells.
    !> The function also includes a stability check to ensure that
    !> the corrected velocity components do not contain NaN values.
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
    !> This subroutine performs a stability check on the corrected
    !> velocity fields to detect any NaN values. If a NaN is found
    !> in any component of the velocity field, the simulation is
    !> halted and an error message is output.

    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: ux_pred(:,:,:), uy_pred(:,:,:), uz_pred(:,:,:)
    real(kind=8), intent(in) :: pp(:,:,:)
    real(kind=8), intent(in) :: dt, dx, dy, dz
    integer, intent(in) :: nx, ny, nz

    real(kind=8), dimension(nx, ny, nz) :: dpdx, dpdy, dpdz
    logical :: has_nan

    print *, "* Correct velocity"
    ! Calculate the pressure gradients
    call derxp(dpdx, pp, dx)
    call deryp(dpdy, pp, dy)
    call derzp(dpdz, pp, dz)

    ! Update the velocities based on the pressure gradient
    ! Subtracting because gradient points from high to low pressure
    ux = ux_pred - dt * dpdx
    uy = uy_pred - dt * dpdy
    uz = uz_pred - dt * dpdz

    ! Stability check: Verify no NaN values in the velocity fields
    has_nan = contains_nan(ux)
    if (has_nan .or. maxval(ux) > 1000.) then
       call write_velocity_diverged()
       stop
    end if

    has_nan = contains_nan(uy)
    if (has_nan .or. maxval(uy) > 1000.) then
       call write_velocity_diverged()
       stop
    end if

    has_nan = contains_nan(uz)
    if (has_nan .or. maxval(uz) > 1000.) then
       call write_velocity_diverged()
       stop
    end if

    print*, ""

    return
  end subroutine correct_velocity

  subroutine transeq(phi, ux, uy, uz, src, fphi, re, sc, adt, bdt, cdt, &
       itime, itscheme, dx, dy, dz, nx, ny, nz)
    !> Solves the transport equation for a scalar variable.
    !> INPUT:
    !> - phi(nx,ny,nz)   : Scalar variable to be transported
    !> - ux(nx,ny,nz)    : x-component of velocity
    !> - uy(nx,ny,nz)    : y-component of velocity
    !> - uz(nx,ny,nz)    : z-component of velocity
    !> - src(nx,ny,nz)   : source terme of the scalar
    !> - re              : Reynolds number
    !> - sc              : Schmidt number
    !> - adt(3)          : Coefficients for the current time step integration
    !> - bdt(3)          : Coefficients for the previous time step integration
    !> - cdt(3)          : Coefficients for two time steps back integration
    !> - itime           : Current time index
    !> - itscheme        : Integration scheme indicator (1=Euler, 2=AB2, 3=AB3)
    !> - dx, dy, dz      : Mesh spacing in x, y, z directions
    !> - nx, ny, nz      : Number of grid points in x, y, z directions
    !> OUTPUT:
    !> - phi(nx,ny,nz)   : Scalar variable updated after integration
    !> - fphi(nx,ny,nz,3): Flux computed for 3 time step

    integer, intent(in) :: itime, itscheme, nx, ny, nz
    real(kind=8), dimension(nx,ny,nz), intent(in) :: ux, uy, uz, src
    real(kind=8), dimension(3), intent(in) :: adt, bdt, cdt
    real(kind=8), intent(in) :: re, sc, dx, dy, dz

    real(kind=8), dimension(nx,ny,nz), intent(inout) :: phi
    real(kind=8), dimension(nx,ny,nz,3), intent(inout) :: fphi

    real(kind=8) :: adu, bdu, cdu, oneresc
    real(kind=8), dimension(:,:,:), allocatable :: dphidx, dphidy, dphidz
    real(kind=8), dimension(:,:,:), allocatable :: dphidx2, dphidy2, dphidz2


    print *, "* Transport equation"

    allocate(dphidx(nx,ny,nz), dphidy(nx,ny,nz), dphidz(nx,ny,nz))
    allocate(dphidx2(nx,ny,nz), dphidy2(nx,ny,nz), dphidz2(nx,ny,nz))

    if (itscheme == 1 .or. itime == 1) then
       print *, " Euler"
       adu = adt(1)
       bdu = bdt(1)
       cdu = cdt(1)
    elseif (itscheme == 2 .or. itime == 2) then
       print *, " Adams-Bashforth 2"
       adu = adt(2)
       bdu = bdt(2)
       cdu = cdt(2)
    elseif (itscheme == 3) then
       print *, " Adams-Bashforth 3"
       adu = adt(3)
       bdu = bdt(3)
       cdu = cdt(3)
    else
       adu = 0.d0
       bdu = 0.d0
       cdu = 0.d0
       print *, 'itscheme:', itscheme, ' unrecognized'
       stop
    end if

    oneresc = 1.d0 / (re * sc)

    call derxp(dphidx, phi, dx)
    call deryp(dphidy, phi, dy)
    call derzp(dphidz, phi, dz)

    call derxxp(dphidx2, phi, dx)
    call deryyp(dphidy2, phi, dy)
    call derzzp(dphidz2, phi, dz)

    fphi(:,:,:,1) = oneresc * (dphidx2 + dphidy2 + dphidz2) - &
         (ux * dphidx + uy * dphidy + uz * dphidz) + src

    phi = phi + adu * fphi(:, :, :, 1) + &
         bdu * fphi(:, :, :, 2) + cdu * fphi(:, :, :, 3)

    select case(itscheme)

    case(2)
       fphi(:,:,:,2) = fphi(:,:,:,1)
    case(3)
       fphi(:,:,:,3) = fphi(:,:,:,2)
       fphi(:,:,:,2) = fphi(:,:,:,1)
    end select


    deallocate(dphidx, dphidy, dphidz)
    deallocate(dphidx2, dphidy2, dphidz2)

    print *, ""
    return

  end subroutine transeq

  subroutine apply_outflow_condition(ux, uy, uz, pp, old_ux, old_uy, old_uz, &
       u0, dx, dy, dz, dt)
    !> Apply outflow conditions for velocity components and pressure at the 
    !> outflow boundary (x = nx).
    !>
    !> This subroutine applies outflow boundary conditions to the velocity
    !> components (ux, uy, uz) and pressure (pp) at the outflow boundary
    !> (x = nx). It first updates the velocity components using a Dirichlet
    !> condition, then adjusts the velocity components based on the pressure
    !> gradients in the y and z directions to ensure a smooth outflow.
    !>
    !> INPUT:
    !> ux(nx,ny,nz)     : X-component of the velocity field.
    !> uy(nx,ny,nz)     : Y-component of the velocity field.
    !> uz(nx,ny,nz)     : Z-component of the velocity field.
    !> pp(nx,ny,nz)     : Pressure field.
    !> old_ux(nx,ny,nz) : X-component of the velocity field at itime-1.
    !> old_uy(nx,ny,nz) : Y-component of the velocity field at itime-1.
    !> old_uz(nx,ny,nz) : Z-component of the velocity field at itime-1.
    !> u0               : Reference veolcity 
    !> dx               : Mesh size step in the x-direction.
    !> dy               : Mesh size step in the y-direction.
    !> dz               : Mesh size step in the z-direction.
    !> dt               : Time step size.
    !>
    !> OUTPUT:
    !> ux(:,:,:)    : Updated X-component of the velocity field with outflow
    !> conditions applied.
    !> uy(:,:,:)    : Updated Y-component of the velocity field with outflow
    !> conditions applied.
    !> uz(:,:,:)    : Updated Z-component of the velocity field with outflow
    !> conditions applied.

    real(kind=8), intent(inout) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: old_ux(:,:,:), old_uy(:,:,:), old_uz(:,:,:)
    real(kind=8), intent(in) :: pp(:,:,:)
    real(kind=8), intent(in) :: u0, dx, dy, dz, dt

    real(kind=8), dimension(size(ux,2),size(ux,3)) :: bux, buy, buz
    real(kind=8), dimension(size(pp,1), size(pp,2),size(pp,3)) :: dpdy
    real(kind=8), dimension(size(pp,1), size(pp,2),size(pp,3)) :: dpdz
    real(kind=8) :: um_in, um_out
    integer :: nx, ny, nz
    integer :: j, k

    nx = size(ux,1)
    ny = size(ux,2)
    nz = size(ux,3)
    call velocity_dirichlet_bc_xn(bux, buy, buz, old_ux, old_uy, old_uz, &
         u0, dx, dy, dt, 1) !> 1 -> use u0 
    !> 2 -> use velocity mean value of ux(nx-1, :, :)
    !> 3 -> use discret velocity value of ux(nx-1, :, :)
    !> 4 -> use ux_interp for cx evaluation
    call deryp(dpdy, pp, dy)
    call derzp(dpdz, pp, dz)

    !Computation of the flow rate Inflow/Outflow
    um_in = 0.d0
    do k = 1, nz
       do j = 1, ny
          um_in = um_in + ux(1,j,k)
       end do
    end do
    um_in = um_in / real(ny * nz, kind=8)
    um_out = 0.d0
    do k = 1, nz
       do j = 1, ny
          um_out = um_out + bux(j,k)
       end do
    end do
    um_out = um_out / real(ny * nz, kind=8)

    call print_in_outflow_rate(um_in,um_out)
    do k = 1, nz
       do j = 1, ny
          bux(j,k) = bux(j,k) - um_out + um_in
       end do
    end do

    ux(nx,:,:) = bux(:,:)
    uy(nx,:,:) = buy(:,:) + dpdy(nx,:,:) * dt
    uz(nx,:,:) = buz(:,:) + dpdz(nx,:,:) * dt

    return
  end subroutine apply_outflow_condition

end module integration

