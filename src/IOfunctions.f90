module IOfunctions
  use initialization
  implicit none

contains

  subroutine print_osinco3d_title()
    character(len=64) :: user_name, machine_name
    integer :: istat

    ! Initialize the character variables
    user_name = ''
    machine_name = ''

    ! Récupérer le nom d'utilisateur
    call get_environment_variable('USER', user_name)
    if (trim(user_name) == '') then
       call get_environment_variable('LOGNAME', user_name)  ! Alternative for some systems
    end if

    ! Récupérer le nom de la machine
    call get_machine_name(machine_name, istat)

    ! Clear the screen
    call execute_command_line('clear')

    ! Afficher le titre ASCII
    print *, ""
    print *, " OOO   SSSS IIIII N   N  CCCC  OOO   3333  DDDD "
    print *, "O   O S       I   NN  N C     O   O      3 D   D"
    print *, "O   O  SSS    I   N N N C     O   O  3333  D   D"
    print *, "O   O     S   I   N  NN C     O   O      3 D   D"
    print *, " OOO  SSSS  IIIII N   N  CCCC  OOO   3333  DDDD "
    print *, ""
    print *, "    User: ", trim(user_name)
    print *, " Machine: ", trim(machine_name)
    print *, ""
  end subroutine print_osinco3d_title

  subroutine print_parameters()
    use initialization, only: nx, ny, nz, xlx, yly, zlz, x0, y0, z0, &
         nbcx1, nbcxn, nbcy1, nbcyn, sim2d, &
         itscheme, cfl, itstart, itstop, &
         u0, l0, re, typesim, omega, eps, kmax, &
         nfre, xpro, ypro, zpro, &
         iin, inflow_noise, init_noise, ratio
    implicit none

    ! Formats
    character(len=20) :: int_format, float_format1, float_format2
    int_format = "(I10)"
    float_format1 = "(F10.2)" ! For values greater than 1
    float_format2 = "(E17.5)" ! For values less than or equal to 1

    ! Print domain parameters
    print *, "Domain Parameters:"
    write(*,*) "  Grid points (nx, ny, nz):"
    write(*,int_format, advance="no") nx 
    write(*,int_format, advance="no") ny 
    write(*,int_format) nz
    write(*,*) "  Domain dimensions (xlx, yly, zlz):"
    write(*,float_format1, advance="no") xlx
    write(*,float_format1, advance="no") yly
    write(*,float_format1) zlz
    write(*,*) "  Domain origin (x0, y0, z0):"
    write(*,float_format1, advance="no") x0
    write(*,float_format1, advance="no") y0
    write(*,float_format1) z0
    print *, ""

    ! Print boundary conditions
    print *, "Boundary Conditions:"
    write(*,*) "  X-axis: start (nbcx1), end (nbcxn):"
    write(*,int_format, advance="no") nbcx1
    write(*,int_format) nbcxn
    write(*,*) "  Y-axis: start (nbcy1), end (nbcyn):"
    write(*,int_format, advance="no") nbcy1
    write(*,int_format) nbcyn
    write(*,*) "  Z-axis: start (nbcz1), end (nbczn):"
    write(*,int_format, advance="no") nbcz1
    write(*,int_format) nbczn
    write(*,*) "  Simulation dimensionality (sim2d):"
    write(*,int_format) sim2d
    print *, ""

    ! Print time integration parameters
    print *, "Time Integration Parameters:"
    write(*,*) "  Time scheme (itscheme):"
    write(*,int_format) itscheme
    write(*,*) "  CFL number (cfl):"
    write(*,float_format1) cfl
    write(*,*) "  Time steps: start (itstart), stop (itstop):"
    write(*,int_format, advance="no") itstart
    write(*,int_format) itstop
    print *, ""

    ! Print flow parameters
    print *, "Flow Parameters:"
    write(*,*) "  Initial velocity (u0):"
    write(*,float_format1) u0
    write(*,*) "  Characteristic length (l0):"
    write(*,float_format1) l0
    write(*,*) "  Reynolds number (re):"
    write(*,float_format1) re
    write(*,*) "  Simulation type (typesim):"
    write(*,int_format) typesim
    print *, ""

    ! Print Poisson solver parameters
    print *, "Poisson Solver Parameters:"
    write(*,*) "  Relaxation parameter (omega):"
    write(*,float_format1) omega
    write(*,*) "  Convergence criterion (eps):"
    write(*,float_format2) eps
    write(*,*) "  Maximum iterations (kmax):"
    write(*,int_format) kmax
    print *, ""

    ! Print Scalar parameters
    write(*,*) "  Flag for scalar transport resolution"
    write(*,int_format) nscr
    write(*,*) "  Schmidt number (re):"
    write(*,float_format1) sc
    print *, ""


    ! Print visualization parameters
    print *, "Visualization Parameters:"
    write(*,*) "  Frequency (nfre):"
    write(*,int_format) nfre
    write(*,*) "  Profile coordinates (xpro, ypro, zpro):"
    write(*,float_format1, advance="no") xpro
    write(*,float_format1, advance="no") ypro
    write(*,float_format1) zpro
    print *, ""

    ! Print inflow initialization parameters
    print *, "Inflow Initialization Parameters:"
    write(*,*) "  Inflow condition (iin):"
    write(*,int_format) iin
    write(*,*) "  Inflow noise intensity (inflow_noise):"
    write(*,float_format1) inflow_noise
    write(*,*) "  Initialization noise intensity (init_noise):"
    write(*,float_format1) init_noise
    write(*,*) "  Ratio (ratio):"
    write(*,float_format1) ratio
    print *, ""

  end subroutine print_parameters

  subroutine print_variables()
    use initialization, only: nx, ny, nz, xlx, yly, zlz, x, y, z, u0, re, cfl, dt, &
         adt, bdt, cdt, t_ref, dx, dy, dz, itstart, itstop, ipro, jpro, kpro, itscheme
    implicit none
    real(kind=8) :: dmin, cnu

    ! Formats
    character(len=20) :: int_format, float_format1, float_format2
    int_format = "(I10)"
    float_format1 = "(F7.2)" ! For values greater than 1
    float_format2 = "(E14.5)" ! For values less than or equal to 1

    ! Calculate minimum grid spacing
    dmin = min(dx, min(dy, dz))

    ! Calculate cnu
    cnu = 1.d0 / re * (dt / (dmin * dmin))

    ! Adjust dt if cnu > 0.1
    if (cnu > 0.1d0) then
       dt = 0.1d0 * dmin * dmin * re 
       print *, "Predominant viscous terms, cnu =", cnu
    end if

    ! Print calculated variables
    print *, "Calculated Variables:"

    write(*,*) "  Grid spacing (dx, dy, dz):"
    write(*, float_format2) dx
    write(*, float_format2) dy
    write(*, float_format2) dz
    write(*,*) "  Time step (dt):"
    write(*, float_format2) dt
    write(*,*) "  Coefficients for time integration (adt, bdt, cdt):"
    write(*, float_format2) adt(itscheme)
    write(*, float_format2) bdt(itscheme)
    write(*, float_format2) cdt(itscheme)
    write(*,*) "  Simulation time (time):"
    write(*,*) "  t_ref"
    write(*, '(F9.1)') t_ref
    write(*,*) "  t_end"
    write(*, '(F9.1)') dt * (itstop - itstart)
    write(*,*) "  Profile index (ipro, jpro, kpro):"
    write(*,'(A5,I4,A4,F7.3)') "   x(", ipro,") = ", x(ipro)
    write(*,'(A5,I4,A4,F7.3)') "   y(", jpro,") = ", y(jpro)
    write(*,'(A5,I4,A4,F7.3)') "   z(", kpro,") = ", z(kpro)
    print *, ""
  end subroutine print_variables

  subroutine get_machine_name(machine_name, istat)
    character(len=64), intent(out) :: machine_name
    integer, intent(out) :: istat
    character(len=256) :: command
    character(len=64) :: line
    integer :: unit

    command = 'hostname > temp_hostname.txt'
    call execute_command_line(command, wait=.true., exitstat=istat)

    if (istat == 0) then
       open(newunit=unit, file='temp_hostname.txt', status='old', action='read')
       read(unit,'(A)') line
       close(unit)
       call execute_command_line('rm temp_hostname.txt')
       machine_name = trim(adjustl(line))
    else
       machine_name = 'Unknown'
    end if
  end subroutine get_machine_name

  subroutine print_residuals(res_u, res_v, res_w, aa, ia, ja, ka, bb, ib, jb, kb, cc, ic, jc, kc)
    implicit none

    real(kind=8), intent(in) :: res_u, res_v, res_w, aa, bb, cc
    integer, intent(in) :: ia, ja, ka, ib, jb, kb, ic, jc, kc
    write(*,'(a20,e12.6)') 'Residual L2 on u: ', res_u
    write(*,'(a20,e12.6)') 'Residual L2 on v: ', res_v
    write(*,'(a20,e12.6)') 'Residual L2 on w: ', res_w
    write(*,*) ''
    write(*,'(a22,e12.6,a21,i4,a2,i4,a2,i4)') 'Residual Linf on u: ', aa, ', at point (i,j,k) : ', ia, ', ', ja, ', ', ka
    write(*,'(a22,e12.6,a21,i4,a2,i4,a2,i4)') 'Residual Linf on v: ', bb, ', at point (i,j,k) : ', ib, ', ', jb, ', ', kb
    write(*,'(a22,e12.6,a21,i4,a2,i4,a2,i4)') 'Residual Linf on w: ', cc, ', at point (i,j,k) : ', ic, ', ', jc, ', ', kc
    write(*,*) ''

  end subroutine print_residuals

  subroutine save_residu(it, res_u, res_v, res_w)
    integer, intent(in) :: it
    real(kind=8), intent(in) :: res_u, res_v, res_w
    character(len=11), parameter :: filename = "residus.dat"
    integer, parameter :: iunit = 27
    integer :: ierr

    open(unit=iunit, file=filename, status='unknown', action='write', iostat=ierr, position='append')

    ! Vérifiez si l'ouverture du fichier a réussi
    if (ierr /= 0) then
       print *, 'Erreur lors de l''ouverture du fichier.'
       stop
    end if

    write(iunit, '(1(I7, 1x), 3(e12.6, 1x))') it, res_u, res_v, res_w

    close(iunit)

    return
  end subroutine save_residu

  subroutine print_cpu_time(elapsed_time, remaining_cpu_hm, &
       time_since_start_hm, mean_elapsed_time)
    real(kind=8), intent(in) :: elapsed_time, mean_elapsed_time
    character(len=7), intent(in) :: remaining_cpu_hm, time_since_start_hm

    print *, "* CPU time of the simulation"
    write(*, '(A, f07.4, A, f07.4, A)') "  Elapsed time for the loop: ", elapsed_time, " s (", mean_elapsed_time, " s)"
    write(*, '(A, 5A)') "  Remaining CPU Time (h:m): ", remaining_cpu_hm
    write(*, '(A, 5A)') "  Time since start (h:m): ", time_since_start_hm
    print *, "" 

  end subroutine print_cpu_time

  subroutine print_divu_statistics(divu_stats, star)
    real(kind=8), dimension(6), intent(in) :: divu_stats
    logical, intent(in) :: star

    ! Affichage des statistiques de div_u
    if (star) then
       print *, "* Statistics of div_u*:"
       write(*, '(A17, ES10.3, A9, ES10.3)') &
            "div(u*): max = ", max(abs(divu_stats(1)), divu_stats(2)), ", mean = ", divu_stats(3)
       write(*,'(A15, I5, A1, I5, A1, I5)') "Max value at:", int(divu_stats(4)), ',', int(divu_stats(5)), ',', int(divu_stats(6))
    else
       print *, "* Statistics of div_u:"
       write(*, '(A16, ES10.3, A9, ES10.3)') &
            "div(u): max = ", max(abs(divu_stats(1)), divu_stats(2)), ", mean = ", divu_stats(3)
       write(*,'(A15, I5, A1, I5, A1, I5)') "Max value at:", int(divu_stats(4)), ',', int(divu_stats(5)), ',', int(divu_stats(6))
    end if
    print *, ""
  end subroutine print_divu_statistics

  subroutine print_velocity_values(ux, uy, uz)
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)

    write(*, '(A26, 3E13.5)') "* Velocities U, V, W min:", minval(ux), minval(uy), minval(uz)
    write(*, '(A26, 3E13.5)') "* Velocities U, V, W max:", maxval(ux), maxval(uy), maxval(uz)
    print *, ""

    return
  end subroutine print_velocity_values

  subroutine  print_scalar_values(phi)
    real(kind=8), intent(in) :: phi(:,:,:)

    write(*, '(A14, 1E13.5)') "* Scalar min:", minval(phi)
    write(*, '(A14, 1E13.5)') "* Scalar max:", maxval(phi)
    print *, ""

    return
  end subroutine print_scalar_values

  subroutine print_in_outflow_rate(um_in,um_out)
    real(kind=8), intent(in) :: um_in, um_out

    write(*,'(A24, 2F8.5, E11.3)') "* Flow rate x I/O/O-I: ", &
         um_in, um_out, um_out-um_in
    print*, ""

    return
  end subroutine print_in_outflow_rate

  subroutine print_cfl(cflx, cfly, cflz)
    real(kind=8), intent(in) :: cflx, cfly, cflz
    write(*, '(A23, 3E13.5)') "* CFL U, CFL V, CFL W:", cflx, cfly, cflz
    print *, ""

    return
  end subroutine print_cfl

  subroutine save_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time)
    !> Save the flow state in a binary file to eventually restart the 
    !> simulation
    !> INPUT
    !> x       : array of coordinates in x-direction
    !> y       : array of coordinates in y-direction
    !> z       : array of coordinates in z-direction
    !> ux      : array of velocity components in x-direction
    !> uy      : array of velocity components in y-direction
    !> uz      : array of velocity components in z-direction
    !> pp      : array of pressure values
    !> nx      : number of grid points in x-direction
    !> ny      : number of grid points in y-direction
    !> nz      : number of grid points in z-direction
    !> time    : time from the beginning of the simulation
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: x(nx), y(ny), z(nz)
    real(kind=8), intent(in) :: ux(nx, ny, nz), uy(nx, ny, nz), uz(nx, ny, nz), pp(nx, ny, nz)
    real(kind=8), intent(in) :: time

    character(len=10), parameter :: filename = "fields.bin"
    integer, parameter :: iunit = 101
    integer :: ios

    print *, "* Save flow state in :", filename

    open(unit=iunit, file=filename, status='replace', access='stream', &
         form='unformatted', action='write', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file: ", trim(filename)
       return
    endif

    write(iunit) time
    write(iunit) nx, ny, nz
    write(iunit) x, y, z
    write(iunit) ux, uy, uz, pp

    close(iunit)
    return
  end subroutine save_fields

  subroutine read_fields(x, y, z, ux, uy, uz, pp, nx, ny, nz, time)
    !> Read the flow state from a binary file to restart the simulation
    !> INPUT
    !> nx      : number of grid points in x-direction
    !> ny      : number of grid points in y-direction
    !> nz      : number of grid points in z-direction
    !> x       : array of coordinates in x-direction
    !> y       : array of coordinates in y-direction
    !> z       : array of coordinates in z-direction
    !> ux      : array of velocity components in x-direction
    !> uy      : array of velocity components in y-direction
    !> uz      : array of velocity components in z-direction
    !> pp      : array of pressure values
    !> time    : time from the beginning of the simulation
    !>
    !> OUTPUT
    !> x       : array of coordinates in x-direction read from file
    !> y       : array of coordinates in y-direction read from file
    !> z       : array of coordinates in z-direction read from file
    !> ux      : array of velocity components in x-direction read from file
    !> uy      : array of velocity components in y-direction read from file
    !> uz      : array of velocity components in z-direction read from file
    !> pp      : array of pressure values read from file
    !> time    : simulation time read from file

    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(inout) :: x(nx), y(ny), z(nz)
    real(kind=8), intent(inout) :: ux(nx, ny, nz), uy(nx, ny, nz)
    real(kind=8), intent(inout) :: uz(nx, ny, nz), pp(nx, ny, nz)
    real(kind=8), intent(inout) :: time

    character(len=10), parameter :: filename = "fields.bin"
    integer, parameter :: iunit = 102
    integer :: ios, nx_r, ny_r, nz_r

    print *, "* Read flow state from :", filename

    open(unit=iunit, file=filename, status='old', access='stream', &
         form='unformatted', action='read', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file: ", trim(filename)
       return
    endif

    read(iunit) time
    read(iunit) nx_r, ny_r, nz_r
    if (nx /= nx_r .or. ny /= ny_r .or. nz /= nz_r) then
       print *, "* Problem: number of cells are different in parameters and fields.bin"
       write(*,*) "nx, nx read:", nx, nx_r
       write(*,*) "ny, ny read:", ny, ny_r
       write(*,*) "nz, nz read:", nz, nz_r
       close(iunit)
       return
    endif

    read(iunit) x, y, z
    read(iunit) ux, uy, uz, pp

    close(iunit)

    return
  end subroutine read_fields

  subroutine check_directories()
    implicit none
    character(len=100) :: directory_name
    character(len=100) :: command
    logical :: dir_exists

    ! Checking and creating 'images' directory
    directory_name = 'images/'
    inquire(file=trim(directory_name), exist=dir_exists)
    if (.not. dir_exists) then
       command = 'mkdir ' // trim(directory_name)
       call system(command)
       print *, 'Directory created: ', trim(directory_name)
    else
       print *, 'Directory exists: ', trim(directory_name)
    end if

    ! Checking and creating 'outputs' directory
    directory_name = 'outputs/'
    inquire(file=trim(directory_name), exist=dir_exists)
    if (.not. dir_exists) then
       command = 'mkdir ' // trim(directory_name)
       call system(command)
       print *, 'Directory created: ', trim(directory_name)
    else
       print *, 'Directory exists: ', trim(directory_name)
    end if
    print *, ""
    return

  end subroutine check_directories

end module IOfunctions

