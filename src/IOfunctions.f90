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
         nbcx1, nbcxn, nbcy1, nbcyn, nbcz1, nbczn, sim2d, &
         itscheme, dt, itstart, itstop, &
         u0, l0, re, typesim, omega, eps, kmax, &
         nfre, xpro, ypro, zpro, &
         init_noise_x, init_noise_y, init_noise_z, ratio, &
         iles, cs, idyn, nscr, sc
    implicit none

    ! Separator pattern
    character(len=40) :: sep_line
    sep_line = "========================================"

    ! Domain parameters
    print *, sep_line
    print *, " DOMAIN PARAMETERS "
    print *, sep_line
    write(*,'(A,I10,2X,I10,2X,I10)') " Grid points (nx, ny, nz): ", nx, ny, nz
    write(*,'(A,3F10.3)') " Domain dimensions (xlx, yly, zlz): ", xlx, yly, zlz
    write(*,'(A,3F10.3)') " Domain origin (x0, y0, z0): ", x0, y0, z0
    print *, ""

    ! Boundary conditions
    print *, sep_line
    print *, " BOUNDARY CONDITIONS "
    print *, sep_line
    write(*,'(A,I10,2X,I10)') " X-axis: start (nbcx1), end (nbcxn): ", nbcx1, nbcxn
    write(*,'(A,I10,2X,I10)') " Y-axis: start (nbcy1), end (nbcyn): ", nbcy1, nbcyn
    write(*,'(A,I10,2X,I10)') " Z-axis: start (nbcz1), end (nbczn): ", nbcz1, nbczn
    write(*,'(A,I10)') " Simulation dimensionality (sim2d): ", sim2d
    print *, ""

    ! Time integration parameters
    print *, sep_line
    print *, " TIME INTEGRATION PARAMETERS "
    print *, sep_line
    write(*,'(A,I10)') " Time scheme (itscheme): ", itscheme
    write(*,'(A,E17.5)') " Time step (dt): ", dt
    write(*,'(A,I10,2X,I10)') " Time steps: start (itstart), stop (itstop): ", itstart, itstop
    print *, ""

    ! Flow parameters
    print *, sep_line
    print *, " FLOW PARAMETERS "
    print *, sep_line
    write(*,'(A,F10.3)') " Initial velocity (u0): ", u0
    write(*,'(A,F10.3)') " Characteristic length (l0): ", l0
    write(*,'(A,F10.3)') " Reynolds number (re): ", re
    write(*,'(A,I10)') " Simulation type (typesim): ", typesim
    print *, ""

    ! Poisson solver parameters
    print *, sep_line
    print *, " POISSON SOLVER PARAMETERS "
    print *, sep_line
    write(*,'(A,F10.3)') " Relaxation parameter (omega): ", omega
    write(*,'(A,E17.5)') " Convergence criterion (eps): ", eps
    write(*,'(A,I10)') " Maximum iterations (kmax): ", kmax
    if (idyn == 1) then
       print *, "Dynamic procedure for omega: ON"
    else 
       print *, "Dynamic procedure for omega: OFF"
    end if
    print *, ""

    ! LES parameters (if relevant)
    if (iles == 1) then
       print *, sep_line
       print *, " LARGE EDDY SIMULATION (LES) PARAMETERS "
       print *, sep_line
       ! Ajoute ici les paramètres LES si nécessaire, par ex :
       write(*,'(A,F10.3)') " Smagorinsky constant (cs): ", cs
       print *, ""
    end if

    ! Scalar transport parameters
    print *, sep_line
    print *, "  SCALAR TRANSPORT PARAMETERS "
    print *, sep_line
    write(*,'(A,I10)') " Scalar transport resolution flag (nscr): ", nscr
    write(*,'(A,F10.3)') " Schmidt number (sc): ", sc
    print *, ""

    ! Visualization parameters
    print *, sep_line
    print *, " VISUALIZATION PARAMETERS "
    print *, sep_line
    write(*,'(A,I10)') " Frequency (nfre): ", nfre
    write(*,'(A)') " Profile coordinates (xpro, ypro, zpro): "
    write(*,'(3F10.3)') xpro, ypro, zpro
    print *, ""

    ! Inflow initialization parameters
    print *, sep_line
    print *, " INFLOW INITIALIZATION PARAMETERS "
    print *, sep_line
    write(*,'(A)') " Initialization noise intensity (init_noise_x, y, z): "
    write(*,'(3F10.3)') init_noise_x, init_noise_y, init_noise_z
    write(*,'(A,F10.3)') " Ratio: ", ratio
    print *, sep_line
    print *, ""

  end subroutine print_parameters

  subroutine print_variables()
    use initialization, only: nx, ny, nz, xlx, yly, zlz, x, y, z, u0, re, cfl, cnu, dt, &
         adt, bdt, cdt, t_ref, dx, dy, dz, itstart, itstop, ipro, jpro, kpro, itscheme, delta, iles, cs
    implicit none

    real(kind=8) :: dmin, cnux, cnuy, cnuz
    character(len=40) :: sep_line

    sep_line       = "========================================"

    ! Compute grid minimum spacing and CFL diffusivity
    dmin = min(dx, min(dy, dz))
    cnu = 1.d0 / re * (dt / (dmin * dmin))

    ! Adjust dt if necessary
    if (cnu > 0.1d0) then
       dt = 0.1d0 * dmin * dmin * re 
       print *, "Predominant viscous terms, adjusted dt"
       print *, "Updated cnu =", cnu
    end if

    cnux = 1.d0 / re * (dt / (dx * dx))
    cnuy = 1.d0 / re * (dt / (dy * dy))
    cnuz = 1.d0 / re * (dt / (dz * dz))

    ! ====== PRINTING START ======
    print *, sep_line
    print *, " CALCULATED VARIABLES"
    print *, sep_line

    ! Grid spacing
    write(*,'(A)') " Grid spacing (dx, dy, dz): " 
    write(*,'(E10.3, E10.3, E10.3)') dx, dy, dz
    print *, ""

    ! CFL numbers
    write(*,'(A,E10.3)') " CFL number (cfl): ", cfl
    write(*,'(A,E10.3)') " CFL diffusivity in x  (cnux): ", cnux
    write(*,'(A,E10.3)') " CFL diffusivity in y  (cnuy): ", cnuy
    write(*,'(A,E10.3)') " CFL diffusivity in z  (cnuz): ", cnuz

    ! Time integration coefficients
    print *, sep_line
    print *, " TIME INTEGRATION COEFFICIENTS"
    print *, sep_line
    write(*,'(A,E10.3)') " adt: ", adt(itscheme)
    write(*,'(A,E10.3)') " bdt: ", bdt(itscheme)
    write(*,'(A,E10.3)') " cdt: ", cdt(itscheme)

    ! Time reference and total simulated time
    print *, sep_line
    print *, " TIME INFORMATION"
    print *, sep_line
    write(*,'(A,F9.1)') " Reference time (t_ref): ", t_ref
    write(*,'(A,F9.1)') " Total simulated time (t_end): ", dt * (itstop - itstart)

    ! Profile point info
    print *, sep_line
    print *, " PROFILE POSITION"
    print *, sep_line
    write(*,'(A,I4,A,F7.3)') " x(", ipro, ") = ", x(ipro)
    write(*,'(A,I4,A,F7.3)') " y(", jpro, ") = ", y(jpro)
    write(*,'(A,I4,A,F7.3)') " z(", kpro, ") = ", z(kpro)

    ! LES info (if used)
    if (iles == 1) then
       print *, sep_line
       print *, " LES PARAMETERS"
       print *, sep_line
       write(*,'(A,E12.5)') " Filter size (delta): ", delta
       write(*,'(A,F5.2)') " Smagorinsky constant (cs): ", cs
    end if

    print *, sep_line
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
    character(len=13), intent(in) :: remaining_cpu_hm, time_since_start_hm

    print *, "* CPU time of the simulation"
    write(*, '(A, f07.4, A, f07.4, A)') "  Elapsed time for the loop: ", elapsed_time, " s (", mean_elapsed_time, " s)"
    write(*, '(A,13A)') "  Remaining CPU Time (h:m): ", remaining_cpu_hm
    write(*, '(A,13A)') "  Time since start (h:m): ", time_since_start_hm
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

  subroutine print_nu_t_statistics(nu_t_stats)
    real(kind=8), dimension(6), intent(in) :: nu_t_stats

    print *, " Statistics of nu_t:"
    write(*, '(A22, ES10.3, ES10.3, ES10.3)') &
         " nu_t min, max, mean:", nu_t_stats(1), nu_t_stats(2), nu_t_stats(3)
  end subroutine print_nu_t_statistics


  subroutine print_velocity_values(ux, uy, uz)
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)

    write(*, '(A26, 3F12.8)') "* Velocities U, V, W min:", minval(ux), minval(uy), minval(uz)
    write(*, '(A26, 3F12.8)') "* Velocities U, V, W max:", maxval(ux), maxval(uy), maxval(uz)
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

  subroutine save_fields(x, y, z, ux, uy, uz, pp, phi, nx, ny, nz, time, itime)
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
    !> phi     : array of passive scalar
    !> nx      : number of grid points in x-direction
    !> ny      : number of grid points in y-direction
    !> nz      : number of grid points in z-direction
    !> time    : time from the beginning of the simulation
    !> itiem   : index of the time step
    integer, intent(in) :: nx, ny, nz, itime
    real(kind=8), intent(in) :: x(nx), y(ny), z(nz)
    real(kind=8), intent(in) :: ux(nx, ny, nz), uy(nx, ny, nz), uz(nx, ny, nz), pp(nx, ny, nz), phi(nx, ny, nz)
    real(kind=8), intent(in) :: time

    character(len=30) :: filename
    integer, parameter :: iunit = 101
    integer :: ios

    write(filename, '(A,I0.6,A)') "fields_", itime, ".bin"
    print *, "* Save flow state in: ", filename

    open(unit=iunit, file=filename, status='replace', access='stream', &
         form='unformatted', action='write', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file: ", trim(filename)
       return
    endif

    write(iunit) time
    write(iunit) nx, ny, nz
    write(iunit) x, y, z
    write(iunit) ux, uy, uz, pp, phi

    close(iunit)
    return

  end subroutine save_fields

  subroutine read_fields(x, y, z, ux, uy, uz, pp, phi, nx, ny, nz, time, filename)
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
    !> phi     : array of passive scalar
    !> time    : time from the beginning of the simulation
    !> filename: name of binary fields file
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
    character(len=30), intent(in) :: filename
    real(kind=8), intent(inout) :: x(nx), y(ny), z(nz)
    real(kind=8), intent(inout) :: ux(nx, ny, nz), uy(nx, ny, nz)
    real(kind=8), intent(inout) :: uz(nx, ny, nz), pp(nx, ny, nz)
    real(kind=8), intent(inout) :: phi(nx, ny, nz)
    real(kind=8), intent(inout) :: time

    integer, parameter :: iunit = 102
    integer :: ios, nx_r, ny_r, nz_r

    print *, "* Read flow state from: ", filename

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
    read(iunit) ux, uy, uz, pp, phi

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

  subroutine write_statistics(t, ek, epst, eps, dzeta, usm, vsm, wsm, &
       duxdx_mean, duxdy_mean, duxdz_mean, &
       duydx_mean, duydy_mean, duydz_mean, &
       duzdx_mean, duzdy_mean, duzdz_mean)

    !> Write statistical quantities of the flow field to a file.
    !>
    !> INPUT:
    !>   t          - Time
    !>   ek         - Kinetic energy
    !>   epst       - Turbulent dissipation rate
    !>   eps        - Dissipation rate
    !>   dzeta      - Enstrophy
    !>   usm        - Mean velocity in x-direction
    !>   vsm        - Mean velocity in y-direction
    !>   wsm        - Mean velocity in z-direction
    !>   duxdx_mean - Mean derivative of ux with respect to x
    !>   duxdy_mean - Mean derivative of ux with respect to y
    !>   duxdz_mean - Mean derivative of ux with respect to z
    !>   duydx_mean - Mean derivative of uy with respect to x
    !>   duydy_mean - Mean derivative of uy with respect to y
    !>   duydz_mean - Mean derivative of uy with respect to z
    !>   duzdx_mean - Mean derivative of uz with respect to x
    !>   duzdy_mean - Mean derivative of uz with respect to y
    !>   duzdz_mean - Mean derivative of uz with respect to z
    !>
    !> OUTPUT:
    !>   Appends the statistical quantities to the file "outputs/stats.dat".

    real(kind=8), intent(in) :: t, ek, epst, eps, dzeta, usm, vsm, wsm
    real(kind=8), intent(in) :: duxdx_mean, duxdy_mean, duxdz_mean
    real(kind=8), intent(in) :: duydx_mean, duydy_mean, duydz_mean
    real(kind=8), intent(in) :: duzdx_mean, duzdy_mean, duzdz_mean

    ! Declare the variables
    integer :: iunit, ios
    parameter (iunit = 10) ! File unit number

    ! Open the file in append mode
    open(unit=iunit, file="outputs/stats.dat", status="unknown", position="append", iostat=ios)

    ! Check if the file was successfully opened
    if (ios /= 0) then
       print *, "Error: Unable to open the file 'outputs/stats.dat'"
       return
    end if

    ! Write the statistics to the file
    write(iunit, '(17es21.12)') t, ek, epst, eps, dzeta, usm, vsm, wsm, &
         duxdx_mean, duxdy_mean, duxdz_mean, &
         duydx_mean, duydy_mean, duydz_mean, &
         duzdx_mean, duzdy_mean, duzdz_mean

    ! Close the file
    close(iunit)

    return

  end subroutine write_statistics

  subroutine header_stats()
    integer :: i, ierr
    ! Define the header lines
    character(len=200) :: header(17)

    header(1)  = "# Column 1  : time t"
    header(2)  = "# Column 2  : kinetic energy E_k [=(u^2+v^2+w^2)/2]"
    header(3)  = "# Column 3  : dissipation epsilon_t [=-dE_k/dt]"
    header(4)  = "# Column 4  : dissipation epsilon [= nu ((du/dx)^2+(du/dy)^2+" // &
         "(du/dz)^2+(dv/dx)^2+(dv/dy)^2+(dv/dz)^2+" // &
         "(dw/dx)^2+(dw/dy)^2+(dw/dz)^2)]"
    header(5)  = "# Column 5  : enstrophy Dzeta [=2 nu epsilon]"
    header(6)  = "# Column 6  : mean square u^2"
    header(7)  = "# Column 7  : mean square v^2"
    header(8)  = "# Column 8  : mean square w^2"
    header(9)  = "# Column 9  : mean square (du/dx)^2"
    header(10) = "# Column 10 : mean square (du/dy)^2"
    header(11) = "# Column 11 : mean square (du/dz)^2"
    header(12) = "# Column 12 : mean square (dv/dx)^2"
    header(13) = "# Column 13 : mean square (dv/dy)^2"
    header(14) = "# Column 14 : mean square (dv/dz)^2"
    header(15) = "# Column 15 : mean square (dw/dx)^2"
    header(16) = "# Column 16 : mean square (dw/dy)^2"
    header(17) = "# Column 17 : mean square (dw/dz)^2"

    ! Open the file for writing
    open(unit=10, file='outputs/stats.dat', status='replace', action='write', &
         form='formatted', iostat=ierr)

    ! Check for errors in opening the file
    if (ierr /= 0) then
       print *, "Error opening file outputs/stats.dat"
       return
    endif

    ! Write the header to the file
    do i = 1, 17
       write(10, '(A)') trim(header(i))
    end do

    ! Close the file
    close(10)

  end subroutine header_stats

  subroutine write_visu_data_size(datasize)

    real(kind=8), intent(in) :: datasize

    print *,"Visu data size:"
    write(*,'(F8.3, 2A)') datasize, ' Go'
    print*, ""

  end subroutine write_visu_data_size

  subroutine write_velocity_diverged()

    print *, "Velocity has diverged"
    return
  end subroutine write_velocity_diverged

  subroutine print_noise_gene(direction)
    character(len=2), intent(in) :: direction

    write(*, '(A55, A2)') " * Generation of the noise for the velocity component: ", direction

    return
  end subroutine print_noise_gene

  subroutine get_filename(filename)
    ! Argument
    character(30), intent(inout) :: filename  ! Output filename

    ! Local variables
    integer :: num_args  ! Number of command-line arguments

    ! Get the number of command-line arguments
    num_args = COMMAND_ARGUMENT_COUNT()

    ! If an argument is provided, use it as the filename
    if (num_args >= 1) then
       call GET_COMMAND_ARGUMENT(1, filename)
    else
       ! Default filename if no argument is provided
       filename = "fields.bin"
    end if

  end subroutine get_filename

end module IOfunctions

