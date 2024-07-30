program osinco3d
  use IOfunctions
  use initialization
  use initial_conditions
  use integration
  use visualization
  use boundary_conditions
  implicit none

  call print_osinco3d_title()
  call check_directories()
  call header_stats()
  call parameters()
  call print_parameters()
  call variables()
  call print_variables()
  call schemes()
  ! Calculation visu module size
  call calc_visu_data_size(datasize, nx, ny, nz, itstop, &
       itstart, nfre, 7)
  call write_visu_data_size(datasize)
  ! Check out dir existences
  ! Set the simulation type
  if (typesim > 0) then
     call set_initialization_type(typesim)
     call init_condition(ux, uy, uz, pp, phi, x, y, z, &
          nx, ny, nz, l0, ratio, nscr, ici)
     time0 = 0.d0
  else
     call read_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time0)
     phi = 0.d0
  end if
  if (nbcx1 == INFLOW_OUTFLOW) then
     call get_inflow(inflow, ux, uy, uz, pp)
  end if
  if (ici == 2 .or. ici == 3) then
     call add_turbulent_init(ux, uy, uz, &
          nx, ny, nz, dy, u0, init_noise)
  end if
  if (iin == 1) then
     call calcul_u_base(u_base, ux(1,:,kpro), dy)
  end if
  call rotational(rotx, roty, rotz, ux, uy, uz, dx, dy, dz, nx, ny, nz)
  call calculate_Q_criterion(q_criterion, ux, uy, uz, dx, dy, dz, nx, ny, nz)
  call divergence(divu, ux, uy, uz, &
       dx, dy, dz, nx, ny, nz, 0)
  divu_stats = function_stats(divu, nx, ny, nz)
  call print_divu_statistics(divu_stats, .false.)
  if (nscr > 0) then
     call visu(phi(:,:,kpro), x, y, nx, ny, num)
  else 
     call visu(rotz(:,:,kpro), x, y, nx, ny, num)
  end if
  call visualize_2d(x, y, nx, ny, &
       ux(:,:,kpro), uy(:,:,kpro), uz(:,:,kpro), pp(:,:,kpro), &
       rotx(:,:,kpro), roty(:,:,kpro), rotz(:,:,kpro), q_criterion(:,:,kpro), &
       divu(:,:,kpro), phi(:,:,kpro), time, "outputs/solution_xy.dat")
  call visualize_2d(y, z, ny, nz, &
       ux(ipro,:,:), uy(ipro,:,:), uz(ipro,:,:), pp(ipro,:,:), &
       rotx(ipro,:,:), roty(ipro,:,:), rotz(ipro,:,:), q_criterion(ipro,:,:), &
       divu(ipro,:,:), phi(ipro,:,:), time, "outputs/solution_yz.dat")
  call visualize_2d(x, z, nx, nz, &
       ux(:,jpro,:), uy(:,jpro,:), uz(:,jpro,:), pp(:,jpro,:), &
       rotx(:,jpro,:), roty(:,jpro,:), rotz(:,jpro,:), q_criterion(:,jpro,:), &
       divu(:,jpro,:), phi(:,jpro,:), time, "outputs/solution_xz.dat")
  call write_profile(x, nx, ux(:,jpro,kpro), uy(:,jpro,kpro), &
       uz(:,jpro,kpro), pp(:,jpro,kpro), "outputs/profil_x.dat")
  call write_profile(y, ny, ux(ipro,:,kpro), uy(ipro,:,kpro), &
       uz(ipro,:,kpro), pp(ipro,:,kpro), "outputs/profil_y.dat")
  call write_profile(z, nz, ux(ipro,jpro,:), uy(ipro,jpro,:), &
       uz(ipro,jpro,:), pp(ipro,jpro,:), "outputs/profil_z.dat")
  call statistics_calc(ux, uy, uz, nx, ny, nz, &
       dx, dy, dz, re, 0.d0)
  print *, ""
  print *, "Do you want to start the loop? (yes/no)"
  read(*, '(A3)') response
  if (response /= 'yes') then
     stop
  end if
  call CPU_TIME(go_time)
  do itime = itstart, itstop
     call CPU_TIME(start_time)
     time = time0 + itime * dt
     write(*,*) "========================"
     write(*,'(A13, I6, A1, I6)') "Itération: ", itime, "/", itstop
     write(*,'(A8, F10.3, A1, F6.0)') "TIME = ", time, "/", itstop * dt
     write(*,*) "========================"
     call old_values(ux, uy, uz, old_ux, old_uy, old_uz, nx, ny, nz)
     if (iin == 1) then
        call add_u_noise(ux, uy, uz, inflow_noise, u_base)
     end if
     call predict_velocity(ux_pred, uy_pred, uz_pred, ux, uy, uz, &
          fux, fuy, fuz, re, adt, bdt, cdt, itime, itscheme, inflow, &
          dx, dy, dz, nx, ny, nz)
     call correct_pression(pp, ux_pred, uy_pred, uz_pred, dx, dy, dz, &
          nx, ny, nz, dt, omega, eps, kmax)
     call correct_velocity(ux, uy, uz, ux_pred, uy_pred, uz_pred, &
          pp, dt, dx, dy, dz, nx, ny, nz)
     if (nscr == 1) then
        call transeq(phi, ux, uy, uz, src, fphi, re, sc, adt, bdt, cdt, &
             itime, itscheme, dx, dy, dz, nx, ny, nz)
     end if
     call divergence(divu_pred, ux_pred, uy_pred, uz_pred, &
          dx, dy, dz, nx, ny, nz, 0)
     divu_stats = function_stats(divu_pred, nx, ny, nz)
     call print_divu_statistics(divu_stats, .true.)
     call divergence(divu, ux, uy, uz, &
          dx, dy, dz, nx, ny, nz, 0)
     divu_stats = function_stats(divu, nx, ny, nz)
     call print_divu_statistics(divu_stats, .false.)

     if (nbcx1 == INFLOW_OUTFLOW) then
        call velocity_dirichlet_bc_x0(ux, uy, uz, inflow)
     end if
     if (nbcxn == INFLOW_OUTFLOW) then
        call apply_outflow_condition(ux, uy, uz, pp, &
             old_ux, old_uy, old_uz, u0, dx, dy, dz, dt)
     end if
     call print_velocity_values(ux, uy, uz)
     if (nscr == 1) call print_scalar_values(phi)
     call compute_cfl(cflx, cfly, cflz, ux, uy, uz, dx, dy, dz, dt)
     call print_cfl(cflx, cfly, cflz)
     if (modulo(itime, nfre) == 0) then
        call rotational(rotx, roty, rotz, ux, uy, uz, dx, dy, dz, nx, ny, nz)
        call calculate_Q_criterion(q_criterion, &
             ux, uy, uz, dx, dy, dz, nx, ny, nz)
        call visualize_2d(x, y, nx, ny, &
             ux(:,:,kpro), uy(:,:,kpro), uz(:,:,kpro), pp(:,:,kpro), &
             rotx(:,:,kpro), roty(:,:,kpro), rotz(:,:,kpro), q_criterion(:,:,kpro), &
             divu(:,:,kpro), phi(:,:,kpro), time, "outputs/solution_xy.dat")
        call visualize_2d(y, z, ny, nz, &
             ux(ipro,:,:), uy(ipro,:,:), uz(ipro,:,:), pp(ipro,:,:), &
             rotx(ipro,:,:), roty(ipro,:,:), rotz(ipro,:,:), q_criterion(ipro,:,:), &
             divu(ipro,:,:), phi(ipro,:,:), time, "outputs/solution_yz.dat")
        call visualize_2d(x, z, nx, nz, &
             ux(:,jpro,:), uy(:,jpro,:), uz(:,jpro,:), pp(:,jpro,:), &
             rotx(:,jpro,:), roty(:,jpro,:), rotz(:,jpro,:), q_criterion(:,jpro,:), &
             divu(:,jpro,:), phi(:,jpro,:), time, "outputs/solution_xz.dat")
        call write_all_data(ux, uy, uz, rotx, roty, rotz, q_criterion, pp, phi, numx, nscr)
        call write_xdmf(nx, ny, nz, dx, dy, dz, x0, y0, z0, numx, nscr)

        if (nscr > 1) then
           call visu(phi(:,:,kpro), x, y, nx, ny, num)
        else 
           call visu(rotz(:,:,kpro), x, y, nx, ny, num)
        end if
     end if
     if (modulo(itime, 25) == 0) then
        call calculate_residuals(ux, uy, uz, &
             old_ux, old_uy, old_uz, dt, &
             t_ref, u_ref, nx, ny, nz, itime)
        call write_profile(x, nx, ux(:,jpro,kpro), uy(:,jpro,kpro), &
             uz(:,jpro,kpro), pp(:,jpro,kpro), "outputs/profil_x.dat")
        call write_profile(y, ny, ux(ipro,:,kpro), uy(ipro,:,kpro), &
             uz(ipro,:,kpro), pp(ipro,:,kpro), "outputs/profil_y.dat")
        call write_profile(z, nz, ux(ipro,jpro,:), uy(ipro,jpro,:), &
             uz(ipro,jpro,:), pp(ipro,jpro,:), "outputs/profil_z.dat")
        call statistics_calc(ux, uy, uz, nx, ny, nz, &
             dx, dy, dz, re, time)
     end if
     if (itime == nsve) then
        call save_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time)
     end if
     call CPU_TIME(end_time)
     call calcul_cpu_time(go_time, start_time, end_time, itime, &
          (itstop-itstart), sum_elapsed_time)
  end do
  stop
end program osinco3d
