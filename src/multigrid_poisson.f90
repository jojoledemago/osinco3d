module multigrid_poisson
  !> Geometric multigrid (V-cycles) for the 3D Poisson equation
  !>
  !>     Laplacian(pp) = rhs
  !>
  !> on a uniform node-centred grid, with FOURTH-ORDER finite differences
  !> (default) or second-order ones.
  !>
  !> Discretisation (in each direction, h = grid spacing):
  !>   order 2 :  ( p(i-1) - 2 p(i) + p(i+1) ) / h**2                    3 points
  !>   order 4 :  ( -p(i-2) + 16 p(i-1) - 30 p(i) + 16 p(i+1) - p(i+2) )
  !>              / (12 h**2)                                            5 points
  !>
  !> Boundary conditions, chosen independently in each direction:
  !>   MG_PERIODIC (0) : periodic wrap-around (nodes 1..n, period n*h)
  !>   MG_NEUMANN  (1) : free-slip wall, dp/dn = 0, nodes ON the walls.
  !>                     Ghost values are mirror images:
  !>                       p(0) = p(2), p(-1) = p(3)   (and same at the other wall)
  !>                     This is accurate to 4th order because a free-slip wall is
  !>                     a symmetry plane: p is an even function of the wall-normal
  !>                     distance, so all its odd normal derivatives vanish.
  !>                     (NOT valid for a wall where dp/dn = 0 only is imposed
  !>                     but the solution is not symmetric: it falls back to
  !>                     2nd-order accuracy near the wall.)
  !>
  !> Singular problem (no Dirichlet condition): pp is defined up to a constant.
  !> The solver removes the (trapezoid-weighted) mean of rhs and returns pp
  !> with zero (trapezoid-weighted) mean.
  !>
  !> Multigrid strategy:
  !>   - finest level  : selected order (4 by default), residual AND smoothing
  !>   - coarser levels: 2nd-order operator (cheaper, and equally good as a
  !>                     coarse-grid correction for the smooth error)
  !>   - the converged solution satisfies the fine-grid 4th-order equation
  !>   - smoother : Gauss-Seidel by colours, (i+j+k) mod 3 on the 4th-order level
  !>                (valid colouring for the 5-point stencil), mod 2 (red-black)
  !>                on the 2nd-order levels
  !>   - transfers: full weighting (restriction), trilinear (prolongation)
  !>   - coarsest level: a fixed number of Gauss-Seidel sweeps
  !>
  !> Coarsening (factor 2, independently in each direction):
  !>   periodic : n even (n >= 4)    -> n/2
  !>   Neumann  : n-1 even (n >= 3)  -> (n-1)/2 + 1
  !> A direction that cannot be coarsened any more is kept as is (semi-
  !> coarsening). Best efficiency: n = 2**k (periodic), n = 2**k+1 (Neumann).
  !> n = 1 in a direction is accepted (2D computations).

  implicit none
  private

  public :: poisson_multigrid_solver
  public :: MG_PERIODIC, MG_NEUMANN

  integer, parameter :: MG_PERIODIC = 0
  integer, parameter :: MG_NEUMANN  = 1

  integer, parameter :: MAX_LEVELS = 32
  integer, parameter :: NU_PRE     = 2   ! pre-smoothing sweeps
  integer, parameter :: NU_POST    = 2   ! post-smoothing sweeps

  type :: mg_level
     integer  :: order = 2                ! 2 or 4
     integer  :: ncol  = 2                ! number of colours of the smoother
     integer  :: n(3)  = 1                ! number of nodes
     integer  :: nc(3) = 1                ! number of nodes of next coarser level
     real(kind=8) :: coef(3)  = 0.d0          ! coefficient of p(i+-1)  (0 if n = 1)
     real(kind=8) :: coef2(3) = 0.d0          ! coefficient of p(i+-2)  (0 for order 2)
     real(kind=8) :: diag     = 1.d0          ! centre coefficient of the stencil
     logical  :: coarsen(3) = .false.     ! is direction d coarsened -> next level
     ! neighbour indices (BC included) and trapezoid weights
     integer,  allocatable :: im(:,:), ip(:,:)       ! i-1, i+1
     integer,  allocatable :: im2(:,:), ip2(:,:)     ! i-2, i+2
     real(kind=8), allocatable :: wt(:,:)
     ! transfer operators to the next coarser level
     integer,  allocatable :: pc1(:,:), pc2(:,:)     ! prolongation: coarse indices
     real(kind=8), allocatable :: pw(:,:)                ! weight of pc2 (weight of pc1 = 1-pw)
     integer,  allocatable :: rn(:,:,:)              ! restriction: fine indices (-1:1, coarse i, d)
     real(kind=8) :: rw(-1:1,3) = 0.d0                   ! restriction weights
     ! fields
     real(kind=8), allocatable :: u(:,:,:), f(:,:,:), r(:,:,:)
  end type mg_level

  ! Hierarchy kept between calls (rebuilt only if grid / BC / order change)
  type(mg_level), allocatable, save :: lev(:)
  integer,  save :: nlev   = 0
  integer,  save :: key_n(3)  = -1
  integer,  save :: key_bc(3) = -1
  integer,  save :: key_order = -1
  real(kind=8), save :: key_h(3)  = -1.d0

contains

  ! ---------------------------------------------------------------------------
  !  Public driver
  ! ---------------------------------------------------------------------------
  subroutine poisson_multigrid_solver(pp, rhs, dx, dy, dz, nx, ny, nz, &
       bc, eps, kmax, iter_out, res_out, verbose, order)
    !> INPUT:
    !> pp   : initial guess (e.g. previous pressure)
    !> rhs  : right-hand side  (Laplacian(pp) = rhs)
    !> dx, dy, dz : grid spacings
    !> nx, ny, nz : number of nodes
    !> bc   : bc(1:3) = boundary condition in x, y, z (MG_PERIODIC / MG_NEUMANN)
    !> eps  : tolerance on the RELATIVE max-norm residual
    !>        max|rhs - L pp| <= eps * max|rhs - mean(rhs)|
    !> kmax : maximum number of V-cycles
    !> order (optional) : 4 (default) or 2, order of the finite differences
    !> OUTPUT:
    !> pp   : solution, with zero (weighted) mean
    !> iter_out (optional) : number of V-cycles performed
    !> res_out  (optional) : final relative residual
    !> verbose  (optional) : print the convergence history

    real(kind=8), intent(inout) :: pp(:,:,:)
    real(kind=8), intent(in)    :: rhs(:,:,:)
    real(kind=8), intent(in)    :: dx, dy, dz
    integer,  intent(in)    :: nx, ny, nz
    integer,  intent(in)    :: bc(3)
    real(kind=8), intent(in)    :: eps
    integer,  intent(in)    :: kmax
    integer,  intent(out), optional :: iter_out
    real(kind=8), intent(out), optional :: res_out
    logical,  intent(in),  optional :: verbose
    integer,  intent(in),  optional :: order

    integer  :: n(3), iter, l, ord
    real(kind=8) :: h(3), fnorm, res, res_old, tol
    logical  :: talk

    print *, "* Poisson multigrid solver start"
    talk = .false.
    if (present(verbose)) talk = verbose
    ord = 4
    if (present(order)) ord = order
    if (ord /= 2 .and. ord /= 4) then
       write(*,*) 'multigrid_poisson: order must be 2 or 4'
       stop 1
    end if

    n = (/ nx, ny, nz /)
    h = (/ dx, dy, dz /)

    if (talk) write(*,'(A14,I2,A,I5)') ' with order :', ord, ', kmax: ', kmax
    if (talk) write(*,'(A24,ES10.3)')   ' convergence criteria: ', eps

    if (any(shape(pp) < n) .or. any(shape(rhs) < n)) then
       write(*,*) 'multigrid_poisson: pp / rhs smaller than (nx,ny,nz)'
       stop 1
    end if
    if (all(n == 1)) then
       pp(1,1,1) = 0.d0
       if (present(iter_out)) iter_out = 0
       if (present(res_out))  res_out  = 0.d0
       return
    end if

    ! (Re)build the grid hierarchy if needed
    if (.not. allocated(lev) .or. any(n /= key_n) .or. any(bc /= key_bc) &
         .or. any(abs(h - key_h) > 1.d-14*abs(h)) .or. ord /= key_order) then
       call mg_setup(n, h, bc, ord)
       key_n = n;  key_bc = bc;  key_h = h;  key_order = ord
       if (talk) then
          write(*,'(A,I3,A,I2,A)') '   Multigrid: ', nlev, ' levels, finest level order ', &
               ord, ', coarser levels order 2'
          do l = 1, nlev
             write(*,'(A,I3,A,3I6)') '   level ', l, ' : n = ', lev(l)%n
          end do
       end if
    end if

    ! Finest level: initial guess and compatible right-hand side
    lev(1)%u = pp(1:nx,1:ny,1:nz)
    lev(1)%f = rhs(1:nx,1:ny,1:nz)
    lev(1)%f = lev(1)%f - wmean(lev(1), lev(1)%f)
    fnorm = maxval(abs(lev(1)%f))

    if (fnorm <= tiny(1.d0)) then
       ! Zero source: solution is a constant
       pp(1:nx,1:ny,1:nz) = 0.d0
       if (present(iter_out)) iter_out = 0
       if (present(res_out))  res_out  = 0.d0
       return
    end if

    tol = eps * fnorm
    call residual(lev(1))
    res = maxval(abs(lev(1)%r))
    res_old = res
    iter = 0
    if (talk) write(*,'(A8,I2,A13,ES12.4)') ' cycle ', 0, ', rel. res = ', res/fnorm

    do while (res > tol .and. iter < kmax)
       iter = iter + 1
       res_old = res
       call vcycle()
       call residual(lev(1))
       res = maxval(abs(lev(1)%r))
       if (talk) then
          if (mod(iter, 5) == 0) write(*,'(A8,I2,A13,ES12.4,A11,F7.4)') ' cycle ', iter, &
               ', rel. res = ', res/fnorm, ', factor = ', res/res_old
       end if
    end do

    if (res > tol) then
       write(*,'(A,I5,A,ES10.3)') ' multigrid_poisson: NOT converged after ', &
            iter, ' V-cycles, rel. res = ', res/fnorm
    else
       if (talk) write(*,'(A8,I2,A13,ES12.4,A11,F7.4)') ' cycle ', iter, &
            ', rel. res = ', res/fnorm, ', factor = ', res/res_old
    end if

    ! Fix the arbitrary constant: zero weighted mean
    lev(1)%u = lev(1)%u - wmean(lev(1), lev(1)%u)
    pp(1:nx,1:ny,1:nz) = lev(1)%u

    if (present(iter_out)) iter_out = iter
    if (present(res_out))  res_out  = res / fnorm
  end subroutine poisson_multigrid_solver

  ! ---------------------------------------------------------------------------
  !  One V-cycle
  ! ---------------------------------------------------------------------------
  subroutine vcycle()
    integer :: l

    do l = 1, nlev-1
       call smooth(lev(l), NU_PRE)
       call residual(lev(l))
       call restrict(lev(l), lev(l+1))
       lev(l+1)%u = 0.d0
    end do

    call coarse_solve(lev(nlev))

    do l = nlev-1, 1, -1
       call prolong_add(lev(l), lev(l+1))
       call smooth(lev(l), NU_POST)
    end do
  end subroutine vcycle

  ! ---------------------------------------------------------------------------
  !  Coarsest level: remove the mean of f (singular problem) and relax
  ! ---------------------------------------------------------------------------
  subroutine coarse_solve(L)
    type(mg_level), intent(inout) :: L
    integer :: nsweep

    L%f = L%f - wmean(L, L%f)
    nsweep = max(20, 10*maxval(L%n))
    call smooth(L, nsweep)
  end subroutine coarse_solve

  ! ---------------------------------------------------------------------------
  !  Gauss-Seidel by colours (colour = mod(i+j+k, ncol)), in place:
  !    u_new = ( f - sum_neighbours coef*u ) / diag
  !  ncol = 2 (red-black) for the 3-point stencil, 3 for the 5-point stencil.
  !  The terms with coef2 vanish identically on 2nd-order levels.
  ! ---------------------------------------------------------------------------
  subroutine smooth(L, nsweep)
    type(mg_level), intent(inout) :: L
    integer, intent(in) :: nsweep

    integer  :: s, color, ncol, i, j, k, i0
    integer  :: jm1, jp1, jm2, jp2, km1, kp1, km2, kp2
    real(kind=8) :: cx, cy, cz, cx2, cy2, cz2, invd

    cx  = L%coef(1);   cy  = L%coef(2);   cz  = L%coef(3)
    cx2 = L%coef2(1);  cy2 = L%coef2(2);  cz2 = L%coef2(3)
    invd = 1.d0 / L%diag
    ncol = L%ncol

    do s = 1, nsweep
       do color = 0, ncol-1
          do k = 1, L%n(3)
             km1 = L%im(k,3);  kp1 = L%ip(k,3)
             km2 = L%im2(k,3); kp2 = L%ip2(k,3)
             do j = 1, L%n(2)
                jm1 = L%im(j,2);  jp1 = L%ip(j,2)
                jm2 = L%im2(j,2); jp2 = L%ip2(j,2)
                i0 = 1 + modulo(color - 1 - j - k, ncol)   ! first i with mod(i+j+k,ncol)==color
                do i = i0, L%n(1), ncol
                   L%u(i,j,k) = ( L%f(i,j,k) &
                        - cx *(L%u(L%im(i,1),j,k)  + L%u(L%ip(i,1),j,k))  &
                        - cx2*(L%u(L%im2(i,1),j,k) + L%u(L%ip2(i,1),j,k)) &
                        - cy *(L%u(i,jm1,k) + L%u(i,jp1,k))               &
                        - cy2*(L%u(i,jm2,k) + L%u(i,jp2,k))               &
                        - cz *(L%u(i,j,km1) + L%u(i,j,kp1))               &
                        - cz2*(L%u(i,j,km2) + L%u(i,j,kp2)) ) * invd
                end do
             end do
          end do
       end do
    end do
  end subroutine smooth

  ! ---------------------------------------------------------------------------
  !  r = f - L u
  ! ---------------------------------------------------------------------------
  subroutine residual(L)
    type(mg_level), intent(inout) :: L

    integer  :: i, j, k
    integer  :: jm1, jp1, jm2, jp2, km1, kp1, km2, kp2
    real(kind=8) :: cx, cy, cz, cx2, cy2, cz2

    cx  = L%coef(1);   cy  = L%coef(2);   cz  = L%coef(3)
    cx2 = L%coef2(1);  cy2 = L%coef2(2);  cz2 = L%coef2(3)

    do k = 1, L%n(3)
       km1 = L%im(k,3);  kp1 = L%ip(k,3)
       km2 = L%im2(k,3); kp2 = L%ip2(k,3)
       do j = 1, L%n(2)
          jm1 = L%im(j,2);  jp1 = L%ip(j,2)
          jm2 = L%im2(j,2); jp2 = L%ip2(j,2)
          do i = 1, L%n(1)
             L%r(i,j,k) = L%f(i,j,k) - ( &
                  cx *(L%u(L%im(i,1),j,k)  + L%u(L%ip(i,1),j,k))  + &
                  cx2*(L%u(L%im2(i,1),j,k) + L%u(L%ip2(i,1),j,k)) + &
                  cy *(L%u(i,jm1,k) + L%u(i,jp1,k))               + &
                  cy2*(L%u(i,jm2,k) + L%u(i,jp2,k))               + &
                  cz *(L%u(i,j,km1) + L%u(i,j,kp1))               + &
                  cz2*(L%u(i,j,km2) + L%u(i,j,kp2))               + &
                  L%diag*L%u(i,j,k) )
          end do
       end do
    end do
  end subroutine residual

  ! ---------------------------------------------------------------------------
  !  Full-weighting restriction of L%r  ->  Lc%f
  ! ---------------------------------------------------------------------------
  subroutine restrict(L, Lc)
    type(mg_level), intent(in)    :: L
    type(mg_level), intent(inout) :: Lc

    integer  :: I, J, K, a, b, c, jj, kk
    real(kind=8) :: s, wc, wbc

    do K = 1, Lc%n(3)
       do J = 1, Lc%n(2)
          do I = 1, Lc%n(1)
             s = 0.d0
             do c = -1, 1
                wc = L%rw(c,3)
                if (wc < tiny(1.d0)) cycle
                kk = L%rn(c,K,3)
                do b = -1, 1
                   if (L%rw(b,2) < tiny(1.d0)) cycle
                   jj  = L%rn(b,J,2)
                   wbc = wc * L%rw(b,2)
                   do a = -1, 1
                      if (L%rw(a,1) < tiny(1.d0)) cycle
                      s = s + wbc * L%rw(a,1) * L%r(L%rn(a,I,1), jj, kk)
                   end do
                end do
             end do
             Lc%f(I,J,K) = s
          end do
       end do
    end do
  end subroutine restrict

  ! ---------------------------------------------------------------------------
  !  Trilinear prolongation of Lc%u, added to L%u
  ! ---------------------------------------------------------------------------
  subroutine prolong_add(L, Lc)
    type(mg_level), intent(inout) :: L
    type(mg_level), intent(in)    :: Lc

    integer  :: i, j, k, i1, i2, j1, j2, k1, k2
    real(kind=8) :: wi, wj, wk, v00, v10, v01, v11

    do k = 1, L%n(3)
       k1 = L%pc1(k,3);  k2 = L%pc2(k,3);  wk = L%pw(k,3)
       do j = 1, L%n(2)
          j1 = L%pc1(j,2);  j2 = L%pc2(j,2);  wj = L%pw(j,2)
          do i = 1, L%n(1)
             i1 = L%pc1(i,1);  i2 = L%pc2(i,1);  wi = L%pw(i,1)
             ! interpolation in x, then y, then z
             v00 = (1.d0-wi)*Lc%u(i1,j1,k1) + wi*Lc%u(i2,j1,k1)
             v10 = (1.d0-wi)*Lc%u(i1,j2,k1) + wi*Lc%u(i2,j2,k1)
             v01 = (1.d0-wi)*Lc%u(i1,j1,k2) + wi*Lc%u(i2,j1,k2)
             v11 = (1.d0-wi)*Lc%u(i1,j2,k2) + wi*Lc%u(i2,j2,k2)
             L%u(i,j,k) = L%u(i,j,k) + (1.d0-wk)*((1.d0-wj)*v00 + wj*v10) &
                  +        wk*((1.d0-wj)*v01 + wj*v11)
          end do
       end do
    end do
  end subroutine prolong_add

  ! ---------------------------------------------------------------------------
  !  Weighted mean (weights 1/2 on Neumann walls: the discrete operator is
  !  self-adjoint for this inner product, so this is the right compatibility
  !  condition and the right "zero mean" definition)
  ! ---------------------------------------------------------------------------
  function wmean(L, a) result(m)
    type(mg_level), intent(in) :: L
    real(kind=8), intent(in) :: a(:,:,:)
    real(kind=8) :: m
    integer  :: i, j, k
    real(kind=8) :: s, sw, w, wjk

    s = 0.d0;  sw = 0.d0
    do k = 1, L%n(3)
       do j = 1, L%n(2)
          wjk = L%wt(j,2) * L%wt(k,3)
          do i = 1, L%n(1)
             w  = L%wt(i,1) * wjk
             s  = s  + w * a(i,j,k)
             sw = sw + w
          end do
       end do
    end do
    m = s / sw
  end function wmean

  ! ---------------------------------------------------------------------------
  !  Hierarchy construction
  ! ---------------------------------------------------------------------------
  subroutine mg_setup(n_in, h_in, bc, ord)
    integer,  intent(in) :: n_in(3), bc(3), ord
    real(kind=8), intent(in) :: h_in(3)

    integer  :: n(3), d
    real(kind=8) :: h(3)

    if (allocated(lev)) deallocate(lev)
    allocate(lev(MAX_LEVELS))

    n = n_in
    h = h_in
    nlev = 0
    do
       nlev = nlev + 1
       if (nlev == 1) then
          call init_level(lev(nlev), n, h, bc, ord)   ! finest: requested order
       else
          call init_level(lev(nlev), n, h, bc, 2)     ! coarser: 2nd order
       end if
       do d = 1, 3
          lev(nlev)%coarsen(d) = can_coarsen(n(d), bc(d))
       end do
       if (.not. any(lev(nlev)%coarsen) .or. nlev == MAX_LEVELS) exit
       call build_transfer(lev(nlev), bc)
       do d = 1, 3
          if (lev(nlev)%coarsen(d)) then
             n(d) = lev(nlev)%nc(d)
             h(d) = 2.d0 * h(d)
          end if
       end do
    end do
    lev(nlev)%coarsen = .false.
  end subroutine mg_setup

  logical function can_coarsen(n, bc)
    integer, intent(in) :: n, bc
    if (bc == MG_PERIODIC) then
       can_coarsen = (n >= 4 .and. mod(n,2) == 0)
    else
       can_coarsen = (n >= 3 .and. mod(n-1,2) == 0)
    end if
  end function can_coarsen

  integer function wrap_index(j, n, b)
    !> Map a (1-based) index j, possibly outside 1..n, back into 1..n:
    !>   periodic : wrap-around
    !>   Neumann  : mirror about node 1 and node n (nodes ON the walls)
    integer, intent(in) :: j, n, b
    integer :: p, jj

    if (n == 1) then
       wrap_index = 1
    else if (b == MG_PERIODIC) then
       wrap_index = modulo(j-1, n) + 1
    else
       p  = 2*(n-1)
       jj = modulo(j-1, p)
       if (jj > n-1) jj = p - jj
       wrap_index = jj + 1
    end if
  end function wrap_index

  subroutine init_level(L, n, h, bc, order)
    type(mg_level), intent(inout) :: L
    integer,  intent(in) :: n(3), bc(3), order
    real(kind=8), intent(in) :: h(3)

    integer :: d, i, m

    m = maxval(n)
    L%n     = n
    L%order = order
    allocate(L%im(m,3), L%ip(m,3), L%im2(m,3), L%ip2(m,3), L%wt(m,3))
    L%im = 1;  L%ip = 1;  L%im2 = 1;  L%ip2 = 1;  L%wt = 1.d0

    do d = 1, 3
       if (n(d) > 1) then
          if (order == 4) then
             L%coef(d)  =  4.d0 / (3.d0*h(d)*h(d))
             L%coef2(d) = -1.d0 / (12.d0*h(d)*h(d))
          else
             L%coef(d)  = 1.d0 / (h(d)*h(d))
             L%coef2(d) = 0.d0
          end if
       else
          L%coef(d)  = 0.d0
          L%coef2(d) = 0.d0
       end if
       do i = 1, n(d)
          L%im(i,d)  = wrap_index(i-1, n(d), bc(d))
          L%ip(i,d)  = wrap_index(i+1, n(d), bc(d))
          L%im2(i,d) = wrap_index(i-2, n(d), bc(d))
          L%ip2(i,d) = wrap_index(i+2, n(d), bc(d))
       end do
       if (n(d) > 1 .and. bc(d) == MG_NEUMANN) then
          L%wt(1,d)    = 0.5d0
          L%wt(n(d),d) = 0.5d0
       end if
    end do
    L%diag = -2.d0 * (sum(L%coef) + sum(L%coef2))

    if (order == 4) then
       L%ncol = 3
    else
       L%ncol = 2
    end if

    allocate(L%u(n(1),n(2),n(3)), L%f(n(1),n(2),n(3)), L%r(n(1),n(2),n(3)))
    L%u = 0.d0;  L%f = 0.d0;  L%r = 0.d0
  end subroutine init_level

  subroutine build_transfer(L, bc)
    !> Restriction / prolongation tables from level L to the next coarser one.
    !> Coarse node I coincides with fine node 2*I-1 in a coarsened direction.
    type(mg_level), intent(inout) :: L
    integer, intent(in) :: bc(3)

    integer :: d, i, ic, nf, nc, m

    m = maxval(L%n)
    allocate(L%pc1(m,3), L%pc2(m,3), L%pw(m,3), L%rn(-1:1,m,3))
    L%pc1 = 1;  L%pc2 = 1;  L%pw = 0.d0;  L%rn = 1

    do d = 1, 3
       nf = L%n(d)
       if (L%coarsen(d)) then
          if (bc(d) == MG_PERIODIC) then
             nc = nf / 2
          else
             nc = (nf - 1) / 2 + 1
          end if
          L%rw(-1,d) = 0.25d0;  L%rw(0,d) = 0.5d0;  L%rw(1,d) = 0.25d0
          do ic = 1, nc
             i = 2*ic - 1
             L%rn(-1,ic,d) = L%im(i,d)
             L%rn( 0,ic,d) = i
             L%rn( 1,ic,d) = L%ip(i,d)
          end do
          do i = 1, nf
             if (mod(i,2) == 1) then
                L%pc1(i,d) = (i+1)/2
                L%pc2(i,d) = (i+1)/2
                L%pw(i,d)  = 0.d0
             else
                L%pc1(i,d) = i/2
                L%pc2(i,d) = i/2 + 1
                L%pw(i,d)  = 0.5d0
                if (L%pc2(i,d) > nc) L%pc2(i,d) = 1   ! periodic wrap
             end if
          end do
       else
          nc = nf
          L%rw(-1,d) = 0.d0;  L%rw(0,d) = 1.d0;  L%rw(1,d) = 0.d0
          do ic = 1, nc
             L%rn(:,ic,d) = ic
             L%pc1(ic,d)  = ic
             L%pc2(ic,d)  = ic
             L%pw(ic,d)   = 0.d0
          end do
       end if
       L%nc(d) = nc
    end do
  end subroutine build_transfer

end module multigrid_poisson
