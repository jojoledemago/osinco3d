module visualization
  use utils
  implicit none

contains

  subroutine visualize_init(x1, x2, n1, n2, ux, uy, uz, pp, filename)
    !> Visualize the flow field by outputting data to a file.
    !> This subroutine outputs the flow velocities and pressure along with their
    !> respective coordinates to a specified file.
    !
    !> INPUT:
    !> x1(n1)     : First coordinate axis values
    !> x2(n2)     : Second coordinate axis values
    !> n1         : Number of points along the first coordinate axis
    !> n2         : Number of points along the second coordinate axis
    !> ux(n1, n2) : X-component of velocity at each grid point
    !> uy(n1, n2) : Y-component of velocity at each grid point
    !> uz(n1, n2) : Z-component of velocity at each grid point (could be zero in 2D cases)
    !> pp(n1, n2) : Pressure at each grid point
    !> filename   : Name of the file to write the output
    !
    !> OUTPUT:
    !> None (outputs to a file)

    real(kind=8), intent(in) :: ux(:,:), uy(:,:), uz(:,:), pp(:,:)
    real(kind=8), intent(in) :: x1(:), x2(:)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n1, n2
    integer :: i1, i2, io_stat
    integer, parameter :: iuf = 11  ! File unit

    open(unit=iuf, file=trim(filename), status='replace', iostat=io_stat)
    if (io_stat /= 0) then
       print *, "Error opening file: ", trim(filename)
       return
    endif
    !
    write(iuf, *) "# x1, x2, ux, uy, uz, pp, i1, i2"
    do i1 = 1, n1
       do i2 = 1, n2
          write(iuf, '(6(F12.6,1x), 2(I6,1x))') x1(i1), x2(i2), &
               ux(i1,i2), uy(i1,i2), uz(i1,i2), pp(i1,i2), i1, i2
       end do
       write(iuf, *)
    end do
    close(iuf)
    return

  end subroutine visualize_init

  subroutine visualize_2d(x1, x2, n1, n2, ux, uy, uz, pp, rotx, roty, rotz, &
       Q, divu, phi, time, filename)
    !> Visualize 2D flow field data by outputting it to a file in a plain text format.
    !> This subroutine outputs the flow velocities, pressure, vorticity, and Q-criterion
    !> at each grid point along with the grid coordinates and the current simulation time.
    !
    !> INPUT:
    !> x1(n1)     : Array of grid points along the first coordinate axis (e.g., x-axis).
    !> x2(n2)     : Array of grid points along the second coordinate axis (e.g., y-axis).
    !> n1         : Number of grid points along the first coordinate axis.
    !> n2         : Number of grid points along the second coordinate axis.
    !> ux(n1, n2) : X-component of velocity at each grid point.
    !> uy(n1, n2) : Y-component of velocity at each grid point.
    !> uz(n1, n2) : Z-component of velocity at each grid point (usually zero in 2D cases).
    !> pp(n1, n2) : Pressure at each grid point.
    !> rotx(n1, n2) : X-component of vorticity at each grid point.
    !> roty(n1, n2) : Y-component of vorticity at each grid point.
    !> rotz(n1, n2) : Z-component of vorticity at each grid point.
    !> Q(n1, n2)    : Q-criterion values at each grid point, used to identify vortices.
    !> divu(n1, n2) : Divergence of velocity values 
    !> time       : Simulation time associated with the data, used for output labeling.
    !> filename   : Name of the file to which the output will be written.
    !
    !> OUTPUT:
    !> Outputs data to a text file specified by 'filename'. The file format includes
    !> headers for each column and data formatted in rows corresponding to each grid point.
    !> The file includes the current simulation time at the beginning.

    real(kind=8), intent(in) :: ux(:,:), uy(:,:), uz(:,:), pp(:,:)
    real(kind=8), intent(in) :: rotx(:,:), roty(:,:), rotz(:,:), Q(:,:)
    real(kind=8), intent(in) :: divu(:,:), phi(:,:)
    real(kind=8), intent(in) :: x1(:), x2(:)
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n1, n2
    integer :: i1, i2, io_stat
    integer, parameter :: iuf = 11  ! File unit identifier

    print *, "* Save 2D solution in ", trim(filename)
    open(unit=iuf, file=trim(filename), status='replace', iostat=io_stat)
    if (io_stat /= 0) then
       print *, "Error opening file: ", trim(filename)
       return
    endif

    write(iuf, *) "# Time : ", time
    write(iuf, *) "# x1, x2, ux, uy, uz, pp, rotx, roty, rotz, Q, divu, phi, i1, i2"
    do i1 = 1, n1
       do i2 = 1, n2
          write(iuf, '(12(F12.6,1x), 2(I6,1x))') x1(i1), x2(i2), &
               ux(i1,i2), uy(i1,i2), uz(i1,i2), pp(i1,i2), &
               rotx(i1,i2), roty(i1,i2), rotz(i1,i2), &
               Q(i1,i2), divu(i1,i2), phi(i1, i2), i1, i2
       end do
       write(iuf, *)
    end do
    close(iuf)
    return

  end subroutine visualize_2d

  subroutine write_vtk(x, y, z, ux, uy, uz, pp, rotx, roty, rotz, Q, nx, ny, nz, timestep)
    !> Write simulation data to a VTK file for visualization in tools like Paraview.
    !> The output includes velocity vectors, vorticity vectors, and scalar pressure.
    !
    !> INPUT:
    !> x, y, z (nx, ny, nz): Coordinate arrays of the simulation domain.
    !> ux, uy, uz (nx, ny, nz): Velocity components at each grid point.
    !> pp (nx, ny, nz): Pressure field at each grid point.
    !> rotx, roty, rotz (nx, ny, nz): Vorticity components at each grid point.
    !> Q(nx, ny, nz) : Q-criterion at each grid point
    !> nx, ny, nz: Dimensions of the grid.
    !> timestep: Current time step number, used for file naming.
    !
    !> OUTPUT:
    !> A VTK file named 'output_<timestep>.vtk' containing the structured grid data.

    integer, intent(in) :: nx, ny, nz, timestep
    real(kind=8), dimension(:), intent(in) :: x, y, z
    real(kind=8), dimension(:,:,:), intent(in) :: ux, uy, uz, pp
    real(kind=8), dimension(:,:,:), intent(in) :: rotx, roty, rotz
    real(kind=8), dimension(:,:,:), intent(in) :: Q
    character(len=32) :: filename
    integer :: i, j, k
    real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax

    write(filename, '(A,I5.5,A)') 'vtk/output_', timestep, '.vtk'
    open(unit=10, file=trim(filename), form='formatted', status='replace')
    print *, "* Write vtk format data in ", filename

    write(10, '(A)') '# vtk DataFile Version 3.0'
    write(10, '(A)') 'CFD Simulation Data'
    write(10, '(A)') 'ASCII'  ! Change this to BINARY for binary files
    write(10, '(A)') 'DATASET STRUCTURED_POINTS'
    write(10, '(A,3I6)') 'DIMENSIONS ', nx, ny, nz
    xmin = minval(x); xmax = maxval(x)
    ymin = minval(y); ymax = maxval(y)
    zmin = minval(z); zmax = maxval(z)
    write(10, '(A,3F12.6)') 'ORIGIN ', xmin, ymin, zmin
    write(10, '(A,3F12.6)') 'SPACING ', (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1), (zmax-zmin)/(nz-1)

    write(10, '(A,I8)') 'POINT_DATA ', nx*ny*nz
    write(10, '(A)') 'VECTORS velocity float'
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             write(10, '(3F12.6)') ux(i,j,k), uy(i,j,k), uz(i,j,k)
          end do
       end do
    end do

    write(10, '(A)') 'VECTORS vorticity float'
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             write(10, '(3F12.6)') rotx(i,j,k), roty(i,j,k), rotz(i,j,k)
          end do
       end do
    end do

    write(10, '(A)') 'SCALARS pressure float 1'
    write(10, '(A)') 'LOOKUP_TABLE default'
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             write(10, '(F12.6)') pp(i,j,k)
          end do
       end do
    end do

    write(10, '(A)') 'SCALARS Q-criterion float 1'
    write(10, '(A)') 'LOOKUP_TABLE default'
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             write(10, '(F12.6)') Q(i,j,k)
          end do
       end do
    end do

    close(10)
  end subroutine write_vtk

  subroutine visu(f, x, y, nx, ny, num)
    integer, intent(in) :: nx, ny
    integer, intent(inout) :: num
    real(kind=8), intent(in) :: f(:,:), x(:), y(:)

    !> Local variables
    real(kind=8) :: max_value
    character(len=4) :: suffixe
    integer :: i, j

    print *, '* ENTREE VISU'
    open(10,file='tampon.dat',form='formatted')
    max_value = maxval(f)
    do j=2,ny-1
       do i=2,nx-1
          write(10,*) x(i), y(j), f(i,j) / max_value
       enddo
       write(10,*)
    enddo
    close(10)
    call system('gnuplot visu.plt')
    write(suffixe,'(i4)') num+1000
    call system('mv tampon.jpeg images/image'//suffixe(1:4)//'.jpeg')
    call system('rm tampon.dat')
    num=num+1

    return
  end subroutine visu

  subroutine write_binary(filename, array)
    implicit none
    character(len=*), intent(in) :: filename
    real(kind=8), intent(in) :: array(:,:,:)
    integer :: iunit, ios

    open(newunit=iunit, file=filename, status='replace', access='stream', form='unformatted', action='write', iostat=ios)
    if (ios /= 0) then
       print *, 'Error opening file: ', filename
       return
    endif

    write(iunit) array
    close(iunit)

    return
  end subroutine write_binary

  subroutine write_all_data(ux, uy, uz, rotx, roty, rotz, qcriterion, pp, phi, num, nscr)
    implicit none
    real(kind=8), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(kind=8), intent(in) :: rotx(:,:,:), roty(:,:,:), rotz(:,:,:)
    real(kind=8), intent(in) :: qcriterion(:,:,:), pp(:,:,:), phi(:,:,:)
    integer, intent(in) :: num, nscr
    character(len=100) :: binaryname


    write(binaryname, '(A,I0,A)') 'outputs/ux_', num, '.bin'
    call write_binary(binaryname, ux)
    write(binaryname, '(A,I0,A)') 'outputs/uy_', num, '.bin'
    call write_binary(binaryname, uy)
    write(binaryname, '(A,I0,A)') 'outputs/uz_', num, '.bin'
    call write_binary(binaryname, uz)
    write(binaryname, '(A,I0,A)') 'outputs/pp_', num, '.bin'
    call write_binary(binaryname, pp)
    write(binaryname, '(A,I0,A)') 'outputs/rotx_', num, '.bin'
    call write_binary(binaryname, rotx)
    write(binaryname, '(A,I0,A)') 'outputs/roty_', num, '.bin'
    call write_binary(binaryname, roty)
    write(binaryname, '(A,I0,A)') 'outputs/rotz_', num, '.bin'
    call write_binary(binaryname, rotz)
    write(binaryname, '(A,I0,A)') 'outputs/vort_', num, '.bin'
    call write_binary(binaryname, sqrt(rotx**2 + roty**2 + rotz**2))
    write(binaryname, '(A,I0,A)') 'outputs/qcrit_', num, '.bin'
    call write_binary(binaryname, qcriterion)
    if (nscr == 1) then
       write(binaryname, '(A,I0,A)') 'outputs/phi_', num, '.bin'
       call write_binary(binaryname, phi)
    end if

    return
  end subroutine write_all_data

  subroutine write_xdmf(nx, ny, nz, dx, dy, dz, x0, y0, z0, num, nscr)
    implicit none
    integer, intent(in) :: nx, ny, nz
    integer, intent(inout) ::  num, nscr
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), intent(in) :: x0, y0, z0
    integer :: iunit, ios
    character(len=100) :: binaryname, filename

    write(filename, '(A,I0,A)') 'outputs/output_', num, '.xdmf'
    print *, "* Write data formated in ", filename
    open(newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       print *, 'Error opening file: ', filename
       return
    endif

    write(iunit, '(A)') '<?xml version="1.0" ?>'
    write(iunit, '(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(iunit, '(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(iunit, '(A)') '<Domain>'
    write(iunit, '(A,I5,I5,I5,A)') '<Topology name="topo" TopologyType="3DRectMesh" Dimensions="', nz, ny, nx, '"/>'
    write(iunit, '(A)') '<Geometry name="geo" GeometryType="ORIGIN_DXDYDZ">'
    write(iunit, '(A)') '<DataItem Dimensions="3" Format="XML">'
    write(iunit, '(3F25.12)') x0, y0, z0
    write(iunit, '(A)') '</DataItem>'
    write(iunit, '(A)') '<DataItem Dimensions="3" Format="XML">'
    write(iunit, '(3F25.12)') dx, dy, dz
    write(iunit, '(A)') '</DataItem>'
    write(iunit, '(A)') '</Geometry>'
    write(iunit, '(A,I0,A)') '<Grid Name="', num, '" GridType="Uniform">'
    write(iunit, '(A)') '<Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(iunit, '(A)') '<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

    write(binaryname, '(A,I0,A)') 'ux_', num, '.bin'
    call write_data_item(iunit, binaryname, 'ux', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'uy_', num, '.bin'
    call write_data_item(iunit, binaryname, 'uy', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'uz_', num, '.bin'
    call write_data_item(iunit, binaryname, 'uz', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'pp_', num, '.bin'
    call write_data_item(iunit, binaryname, 'pp', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'rotx_', num, '.bin'
    call write_data_item(iunit, binaryname, 'rotx', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'roty_', num, '.bin'
    call write_data_item(iunit, binaryname, 'roty', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'rotz_', num, '.bin'
    call write_data_item(iunit, binaryname, 'rotz', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'vort_', num, '.bin'
    call write_data_item(iunit, binaryname, 'vort', nx, ny, nz)
    write(binaryname, '(A,I0,A)') 'qcrit_', num, '.bin'
    call write_data_item(iunit, binaryname, 'qcrit', nx, ny, nz)
    if (nscr == 1) then
       write(binaryname, '(A,I0,A)') 'phi_', num, '.bin'
       call write_data_item(iunit, binaryname, 'phi', nx, ny, nz)
    end if

    write(iunit, '(A)') '</Grid>'
    write(iunit, '(A)') '</Domain>'
    write(iunit, '(A)') '</Xdmf>'

    close(iunit)
    num = num+1

    return
  end subroutine write_xdmf

  subroutine write_data_item(iunit, filename, dataname, nx, ny, nz)
    implicit none
    integer, intent(in) :: iunit, nx, ny, nz
    character(len=*), intent(in) :: filename, dataname

    write(iunit, '(A,A,A)') '<Attribute Name="', dataname, '" Center="Node">'
    write(iunit, '(A,I5,I5,I5,A)') '<DataItem Format="Binary"&
         & DataType="Float" Precision="8" Endian="little" Seek="0" Dimensions="', nz, ny, nx, '">'
    write(iunit, '(A)') filename
    write(iunit, '(A)') '</DataItem>'
    write(iunit, '(A)') '</Attribute>'
  end subroutine write_data_item

  subroutine write_profile(x, n, u, v, w, p, filename)
    character(len=*) :: filename
    real(kind=8), dimension(:), intent(in) :: x, u, v, w, p
    integer, intent(in) :: n
    integer :: i
    integer, parameter :: iunit=16

    open(unit=iunit, file=filename, status='replace')
    do i = 1, n
       write(iunit, '(6(f12.6))') x(i), u(i), v(i), w(i), p(i)
    end do
    close(iunit)
    return
  end subroutine write_profile

end module visualization

