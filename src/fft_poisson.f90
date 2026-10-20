module fft_poisson
  !> Direct (FFT based) solver for the 3D Poisson equation
  !>
  !>     Laplacian(pp) = rhs
  !>
  !> on a UNIFORM node-centred grid, with the same boundary-condition
  !> convention and the same calling sequence as multigrid_poisson:
  !>
  !>   MG_PERIODIC (0) : periodic, nodes 1..n, period n*h          -> DFT (c2c)
  !>   MG_NEUMANN  (1) : free-slip wall, dp/dn = 0, nodes ON walls -> DCT-I
  !>                     (FFTW_REDFT00, logical period 2*(n-1))
  !>
  !> Any combination of the two BC in x, y, z is accepted, and any grid size
  !> (n=1 in a direction is accepted: 2D computations).
  !> Best speed: n = 2**a 3**b 5**c (periodic), n-1 = 2**a 3**b 5**c (Neumann).
  !>
  !> WHICH OPERATOR IS INVERTED  (this is the point of this module)
  !>
  !> composed = .true. (default), order = 2, 4 or 6:
  !>    The symbol of D.G is used, with D and G the centred first-derivative
  !>    stencils of the given order (the ones of the divergence and of the
  !>    pressure gradient of the projection step):
  !>        k'(theta) h = 3/2 sin(t) - 3/10 sin(2t) + 1/30 sin(3t)    (order 6)
  !>                    = 4/3 sin(t) - 1/6 sin(2t)                    (order 4)
  !>                    = sin(t)                                      (order 2)
  !>        lambda = k'x**2 + k'y**2 + k'z**2   (Laplacian = -lambda)
  !>    If the correction step uses G of the same order and the same ghost
  !>    convention as D (periodic wrap / mirror about the wall nodes), then
  !>    div(u^{n+1}) = 0 to machine precision (about 1e-13 in practice).
  !>    Symbol zeros (mean mode, Nyquist modes) are removed: the right-hand side
  !>    D.u* has no component there if D has the same order.
  !>
  !> composed = .false., order = 2 or 4:
  !>    The compact Laplacian of multigrid_poisson (3 or 5 points per
  !>    direction) is inverted exactly. This reproduces the multigrid solution
  !>    without the tolerance, but it does NOT cancel div(u^{n+1}) if D and G
  !>    are of different (wider) stencils.
  !>
  !> As in multigrid_poisson: the mean of rhs is removed (compatibility) and pp
  !> is returned with zero (trapezoid-weighted) mean: the DCT-I mode m = 0 is
  !> exactly the trapezoid weighted sum, so the two conventions coincide.
  !>
  !> eps, kmax and the initial value of pp are not used (direct solver).
  !> Plans and work arrays are built once and rebuilt only if the grid, BC,
  !> spacing or operator change. Not thread safe (module save variables).
  !>
  !> Requires: libfftw3 (and libfftw3_threads if poisson_fft_nthreads > 1).

  use, intrinsic :: iso_c_binding
  use multigrid_poisson, only: MG_PERIODIC, MG_NEUMANN
  implicit none
  private
  include 'fftw3.f03'

  public :: poisson_fft_solver, poisson_fft_free
  public :: poisson_fft_nthreads

  !> number of FFTW threads (1 = sequential). Set before the first call.
  integer, save :: poisson_fft_nthreads = 1

  ! Hierarchy / plans kept between calls
  logical,      save :: ready     = .false.
  integer,      save :: key_n(3)  = -1
  integer,      save :: key_bc(3) = -1
  integer,      save :: key_order = -1
  logical,      save :: key_comp  = .false.
  real(kind=8), save :: key_h(3)  = -1.d0

  logical,      save :: has_per = .false.
  integer,      save :: nthr_used = 0

  real(C_DOUBLE),            allocatable, save :: wr(:,:,:)   ! real work array
  complex(C_DOUBLE_COMPLEX), allocatable, save :: wc(:,:,:)   ! complex work array
  real(kind=8),              allocatable, save :: sym(:,:)    ! symbol per direction
  real(kind=8), save :: inv_norm = 1.d0
  real(kind=8), save :: lam_tol  = 0.d0

  type(C_PTR), save :: p_r2r(3) = C_NULL_PTR   ! DCT-I (Neumann directions)
  type(C_PTR), save :: p_fwd(3) = C_NULL_PTR   ! c2c forward (periodic directions)
  type(C_PTR), save :: p_bwd(3) = C_NULL_PTR   ! c2c backward

contains

  ! ---------------------------------------------------------------------------
  !  Public driver (same calling sequence as poisson_multigrid_solver,
  !  plus the optional "composed")
  ! ---------------------------------------------------------------------------
  subroutine poisson_fft_solver(pp, rhs, dx, dy, dz, nx, ny, nz, &
       bc, verbose, order, composed)
    real(kind=8), intent(inout) :: pp(:,:,:)
    real(kind=8), intent(in)    :: rhs(:,:,:)
    real(kind=8), intent(in)    :: dx, dy, dz
    integer,  intent(in)    :: nx, ny, nz
    integer,  intent(in)    :: bc(3)
    logical,  intent(in),  optional :: verbose
    integer,  intent(in),  optional :: order           ! default 6
    logical,  intent(in),  optional :: composed        ! default .true.

    integer  :: n(3), d, i, j, k, ord
    real(kind=8) :: h(3), lam, fac, sj, sk
    logical  :: talk, comp

    print *, "* Poisson FFT solver start"
    talk = .false.
    if (present(verbose)) talk = verbose
    ord = 6
    if (present(order)) ord = order
    comp = .true.
    if (present(composed)) comp = composed

    if (comp) then
       if (ord /= 2 .and. ord /= 4 .and. ord /= 6) then
          write(*,*) 'poisson_fft: order must be 2, 4 or 6 (composed D.G)'
          stop 1
       end if
    else
       if (ord /= 2 .and. ord /= 4) then
          write(*,*) 'poisson_fft: order must be 2 or 4 (compact Laplacian)'
          stop 1
       end if
    end if
    do d = 1, 3
       if (bc(d) /= MG_PERIODIC .and. bc(d) /= MG_NEUMANN) then
          write(*,*) 'poisson_fft: bc must be MG_PERIODIC or MG_NEUMANN'
          stop 1
       end if
    end do

    n = (/ nx, ny, nz /)
    h = (/ dx, dy, dz /)

    if (any(shape(pp) < n) .or. any(shape(rhs) < n)) then
       write(*,*) 'poisson_fft: pp / rhs smaller than (nx,ny,nz)'
       stop 1
    end if
    if (all(n == 1)) then
       pp(1,1,1) = 0.d0
       return
    end if

    ! (Re)build symbol, work arrays and plans if needed
    if (.not. ready .or. any(n /= key_n) .or. any(bc /= key_bc) &
         .or. any(abs(h - key_h) > 1.d-14*abs(h)) .or. ord /= key_order &
         .or. (comp .neqv. key_comp)) then
       call poisson_fft_free()
       call build_symbol(n, h, bc, ord, comp)
       call build_plans(n, bc)
       key_n = n;  key_bc = bc;  key_h = h;  key_order = ord;  key_comp = comp
       ready = .true.
       if (talk) then
          if (comp) then
             write(*,'(A,I2,A)') '   FFT solver: symbol of D.G, first derivative of order ', ord, &
                  ' (consistent projection)'
          else
             write(*,'(A,I2,A)') '   FFT solver: compact Laplacian of order ', ord, &
                  ' (same operator as multigrid_poisson)'
          end if
          write(*,'(A,3I2,A)') '   BC (0 = periodic, 1 = free-slip) :', bc, ' (x, y, z)'
       end if
    end if

    ! ---- Forward transforms ------------------------------------------------
    wr = rhs(1:nx,1:ny,1:nz)
    do d = 1, 3
       if (n(d) > 1 .and. bc(d) == MG_NEUMANN) call fftw_execute_r2r(p_r2r(d), wr, wr)
    end do

    ! ---- Division by the symbol, inverse transforms ------------------------
    if (has_per) then
       wc = cmplx(wr, 0.d0, kind=C_DOUBLE_COMPLEX)
       do d = 1, 3
          if (n(d) > 1 .and. bc(d) == MG_PERIODIC) call fftw_execute_dft(p_fwd(d), wc, wc)
       end do

       do k = 1, nz
          sk = sym(k,3)
          do j = 1, ny
             sj = sym(j,2) + sk
             do i = 1, nx
                lam = sym(i,1) + sj
                if (lam > lam_tol) then
                   fac = -inv_norm / lam         ! Laplacian = -lambda
                   wc(i,j,k) = wc(i,j,k) * fac
                else
                   wc(i,j,k) = (0.d0, 0.d0)      ! mean / blind modes
                end if
             end do
          end do
       end do

       do d = 1, 3
          if (n(d) > 1 .and. bc(d) == MG_PERIODIC) call fftw_execute_dft(p_bwd(d), wc, wc)
       end do
       wr = real(wc, kind=C_DOUBLE)
    else
       do k = 1, nz
          sk = sym(k,3)
          do j = 1, ny
             sj = sym(j,2) + sk
             do i = 1, nx
                lam = sym(i,1) + sj
                if (lam > lam_tol) then
                   fac = -inv_norm / lam
                   wr(i,j,k) = wr(i,j,k) * fac
                else
                   wr(i,j,k) = 0.d0
                end if
             end do
          end do
       end do
    end if

    do d = 1, 3
       if (n(d) > 1 .and. bc(d) == MG_NEUMANN) call fftw_execute_r2r(p_r2r(d), wr, wr)
    end do

    pp(1:nx,1:ny,1:nz) = wr

  end subroutine poisson_fft_solver

  ! ---------------------------------------------------------------------------
  !  Release plans and work arrays
  ! ---------------------------------------------------------------------------
  subroutine poisson_fft_free()
    integer :: d
    do d = 1, 3
       if (c_associated(p_r2r(d))) call fftw_destroy_plan(p_r2r(d))
       if (c_associated(p_fwd(d))) call fftw_destroy_plan(p_fwd(d))
       if (c_associated(p_bwd(d))) call fftw_destroy_plan(p_bwd(d))
       p_r2r(d) = C_NULL_PTR
       p_fwd(d) = C_NULL_PTR
       p_bwd(d) = C_NULL_PTR
    end do
    if (allocated(wr))  deallocate(wr)
    if (allocated(wc))  deallocate(wc)
    if (allocated(sym)) deallocate(sym)
    ready = .false.
  end subroutine poisson_fft_free

  ! ---------------------------------------------------------------------------
  !  Symbol of the operator (positive: Laplacian = -lambda)
  ! ---------------------------------------------------------------------------
  subroutine build_symbol(n, h, bc, order, composed)
    integer,  intent(in) :: n(3), bc(3), order
    real(kind=8), intent(in) :: h(3)
    logical,  intent(in) :: composed

    real(kind=8), parameter :: pi = acos(-1.d0)
    integer  :: d, m, nn
    real(kind=8) :: th, f, s, scale
    logical  :: nyq

    allocate(sym(maxval(n), 3))
    sym = 0.d0
    scale = 0.d0

    do d = 1, 3
       nn = n(d)
       if (nn <= 1) cycle
       scale = scale + 1.d0 / (h(d)*h(d))
       do m = 0, nn-1
          if (bc(d) == MG_PERIODIC) then
             th  = 2.d0 * pi * real(m, kind=8) / real(nn, kind=8)
             nyq = (mod(nn,2) == 0 .and. m == nn/2)
          else
             th  = pi * real(m, kind=8) / real(nn-1, kind=8)
             nyq = (m == nn-1)
          end if
          if (composed) then
             select case (order)
             case (2)
                f = sin(th)
             case (4)
                f = (4.d0/3.d0)*sin(th) - sin(2.d0*th)/6.d0
             case default
                f = 1.5d0*sin(th) - 0.3d0*sin(2.d0*th) + sin(3.d0*th)/30.d0
             end select
             s = (f/h(d))**2
             if (nyq .or. m == 0) s = 0.d0       ! exact zeros of the symbol
          else
             if (order == 2) then
                s = 2.d0*(1.d0 - cos(th)) / (h(d)*h(d))
             else
                s = (30.d0 - 32.d0*cos(th) + 2.d0*cos(2.d0*th)) / (12.d0*h(d)*h(d))
             end if
             if (m == 0) s = 0.d0
          end if
          sym(m+1, d) = s
       end do
    end do

    lam_tol = 1.d-12 * scale
  end subroutine build_symbol

  ! ---------------------------------------------------------------------------
  !  FFTW plans, one 1D transform per direction, batched over the two others
  !  with the guru interface (arbitrary strides, no loop over planes)
  ! ---------------------------------------------------------------------------
  subroutine build_plans(n, bc)
    integer, intent(in) :: n(3), bc(3)

    type(fftw_iodim) :: dims(1), howm(2)
    integer(C_FFTW_R2R_KIND) :: knd(1)
    integer :: d, e1, e2, stride(3), ierr

    stride(1) = 1
    stride(2) = n(1)
    stride(3) = n(1) * n(2)

    has_per = any(bc == MG_PERIODIC .and. n > 1)

    allocate(wr(n(1), n(2), n(3)))
    wr = 0.d0
    if (has_per) then
       allocate(wc(n(1), n(2), n(3)))
       wc = (0.d0, 0.d0)
    end if

    ! Normalisation of one forward + backward pass
    inv_norm = 1.d0
    do d = 1, 3
       if (n(d) <= 1) cycle
       if (bc(d) == MG_PERIODIC) then
          inv_norm = inv_norm / real(n(d), kind=8)
       else
          inv_norm = inv_norm / real(2*(n(d)-1), kind=8)
       end if
    end do

    if (poisson_fft_nthreads > 1 .and. nthr_used /= poisson_fft_nthreads) then
       ierr = fftw_init_threads()
       call fftw_plan_with_nthreads(int(poisson_fft_nthreads, C_INT))
       nthr_used = poisson_fft_nthreads
    end if

    do d = 1, 3
       if (n(d) <= 1) cycle
       e1 = mod(d, 3) + 1
       e2 = mod(d+1, 3) + 1

       dims(1)%n  = n(d)
       dims(1)%is = stride(d)
       dims(1)%os = stride(d)
       howm(1)%n  = n(e1)
       howm(1)%is = stride(e1)
       howm(1)%os = stride(e1)
       howm(2)%n  = n(e2)
       howm(2)%is = stride(e2)
       howm(2)%os = stride(e2)

       if (bc(d) == MG_NEUMANN) then
          knd(1) = FFTW_REDFT00
          p_r2r(d) = fftw_plan_guru_r2r(1, dims, 2, howm, wr, wr, knd, FFTW_ESTIMATE)
          if (.not. c_associated(p_r2r(d))) then
             write(*,*) 'poisson_fft: cannot create DCT-I plan, direction ', d
             stop 1
          end if
       else
          p_fwd(d) = fftw_plan_guru_dft(1, dims, 2, howm, wc, wc, FFTW_FORWARD,  FFTW_ESTIMATE)
          p_bwd(d) = fftw_plan_guru_dft(1, dims, 2, howm, wc, wc, FFTW_BACKWARD, FFTW_ESTIMATE)
          if (.not. c_associated(p_fwd(d)) .or. .not. c_associated(p_bwd(d))) then
             write(*,*) 'poisson_fft: cannot create DFT plan, direction ', d
             stop 1
          end if
       end if
    end do
  end subroutine build_plans

end module fft_poisson
