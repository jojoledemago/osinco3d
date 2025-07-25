program osinco3d
  use IOfunctions
  use initialization
  use initial_conditions
  use integration
  use visualization
  use les_turbulence
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
  if ((iles+nscr) == 2) then
     call calc_visu_data_size(datasize, &
          nx, ny, nz, itstop, itstart, nfre, 8)
  else if ((iles+nscr) == 1) then
     call calc_visu_data_size(datasize, &
          nx, ny, nz, itstop, itstart, nfre, 7)
  else 
     call calc_visu_data_size(datasize, &
          nx, ny, nz, itstop, itstart, nfre, 6)
  end if
  call write_visu_data_size(datasize)
  ! Check out dir existences
  ! Set the simulation type
  if (typesim > 0) then
     call set_initialization_type(typesim)
     call init_condition(ux, uy, uz, pp, phi, x, y, z, &
          nx, ny, nz, l0, ratio, nscr)
     time0 = 0.d0
  else
     call get_filename(fields_name_file)
     call read_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time0, fields_name_file)
     time = time0
     phi = 0.d0
  end if
  if (ici == 2 .or. ici == 3) then
     call add_turbulent_init(ux, uy, uz, &
          nx, ny, nz, dy, u0, init_noise_x, init_noise_y, init_noise_z, typesim)
  end if
  if (ici == 1 .or. ici == 3) then
     call add_oscillations_init(ux, uy, uz, &
          nx, ny, nz, dy, u0, init_noise_x, init_noise_y, init_noise_z, typesim)
  end if
  if (sim2d == 1) then
     call apply_2dsim(uz)
  end if
  call print_velocity_values(ux, uy, uz)
  call rotational(rotx, roty, rotz, ux, uy, uz, dx, dy, dz, nx, ny, nz)
  call calculate_Q_criterion(q_criterion, ux, uy, uz, dx, dy, dz, nx, ny, nz)
  call divergence(divu, ux, uy, uz, &
       dx, dy, dz, nx, ny, nz, 1)
  divu_stats = function_stats(divu, nx, ny, nz)
  call print_divu_statistics(divu_stats, .false.)
  call compute_cfl(cflx, cfly, cflz, ux, uy, uz, dx, dy, dz, dt)
  call print_cfl(cflx, cfly, cflz)
  if (nscr == 1) then
     call visu(phi(:,:,kpro), x, y, nx, ny, num)
  else 
     call visu(rotz(:,:,kpro), x, y, nx, ny, num)
  end if
  call save_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time, 0)
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
  call write_all_data(ux, uy, uz, rotx, roty, rotz, &
       q_criterion, pp, phi, nu_t, numx, nscr, iles)
  call write_xdmf(nx, ny, nz, dx, dy, dz, x0, y0, z0, numx, nscr, iles, time)
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
     write(*,'(A12, I6, A1, I6)') "Iteration: ", itime, "/", itstop
     write(*,'(A8, F10.3, A1, F6.0)') "TIME = ", time, "/", time0 + itstop * dt
     write(*,*) "========================"
     call old_values(ux, uy, uz, old_ux, old_uy, old_uz, nx, ny, nz)
     call predict_velocity(ux_pred, uy_pred, uz_pred, ux, uy, uz, &
          fux, fuy, fuz, re, adt, bdt, cdt, itime, itscheme, &
          dx, dy, dz, nx, ny, nz, iles, cs, delta, nu_t)
     call correct_pression(pp, ux_pred, uy_pred, uz_pred, dx, dy, dz, &
          nx, ny, nz, dt, omega, eps, kmax, idyn, multigrid)
     call correct_velocity(ux, uy, uz, ux_pred, uy_pred, uz_pred, &
          pp, dt, dx, dy, dz, nx, ny, nz)
     if (nscr == 1) then
        call transeq(phi, ux, uy, uz, src, fphi, re, sc, adt, bdt, cdt, &
             itime, itscheme, dx, dy, dz, nx, ny, nz, iles, nu_t)
     end if
     call divergence(divu_pred, ux_pred, uy_pred, uz_pred, &
          dx, dy, dz, nx, ny, nz, 1)
     divu_stats = function_stats(divu_pred, nx, ny, nz)
     call print_divu_statistics(divu_stats, .true.)
     call divergence(divu, ux, uy, uz, &
          dx, dy, dz, nx, ny, nz, 1)
     divu_stats = function_stats(divu, nx, ny, nz)
     call print_divu_statistics(divu_stats, .false.)

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
             rotx(:,:,kpro), roty(:,:,kpro), rotz(:,:,kpro), &
             q_criterion(:,:,kpro), divu(:,:,kpro), phi(:,:,kpro), &
             time, "outputs/solution_xy.dat")
        call visualize_2d(y, z, ny, nz, &
             ux(ipro,:,:), uy(ipro,:,:), uz(ipro,:,:), pp(ipro,:,:), &
             rotx(ipro,:,:), roty(ipro,:,:), rotz(ipro,:,:), & 
             q_criterion(ipro,:,:), divu(ipro,:,:), phi(ipro,:,:), &
             time, "outputs/solution_yz.dat")
        call visualize_2d(x, z, nx, nz, &
             ux(:,jpro,:), uy(:,jpro,:), uz(:,jpro,:), pp(:,jpro,:), &
             rotx(:,jpro,:), roty(:,jpro,:), rotz(:,jpro,:), &
             q_criterion(:,jpro,:), divu(:,jpro,:), phi(:,jpro,:), &
             time, "outputs/solution_xz.dat")
        call write_all_data(ux, uy, uz, rotx, roty, rotz, &
             q_criterion, pp, phi, nu_t, numx, nscr, iles)
        call write_xdmf(nx, ny, nz, dx, dy, dz, x0, y0, z0, numx, nscr, iles, time)

        if (nscr == 1) then
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
        if (time > initstat) then
           call statistics_calc(ux, uy, uz, nx, ny, nz, &
                dx, dy, dz, re, time)
        end if
     end if
     if (modulo(itime, nsve) == 0) then
        call save_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time, itime)
     end if
     call CPU_TIME(end_time)
     call calcul_cpu_time(go_time, start_time, end_time, itime, &
          (itstop - itstart), sum_elapsed_time)
  end do
  stop
end program osinco3d
