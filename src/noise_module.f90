module noise_module
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
contains

  subroutine noise_generator_1D(n, k0, s, noise, seed_input)
    !> 1D random noise generator in physical space from a spectral distribution
    !> with energy envelope E(k) ∝ k² / (1 + (k/k₀)^|s|)
    !> Allows both positive and negative spectral slope s
    !> Based on inverse FFT using FFTW3 (complex to real transform)
    !>
    !> INPUT:  n          : number of spatial points (grid size)
    !>         k0         : wavenumber at which energy is peaked
    !>         s          : spectral slope (can be negative)
    !>         seed_input : optionnal seed for Initialize random seed
    !> OUTPUT: noise      : real-valued 1D noise in physical space (size n)

    integer, intent(in) :: n
    integer, intent(in), optional :: seed_input
    real(kind=8), intent(in) :: k0, s
    real(kind=8), intent(out) :: noise(n)

    complex(c_double_complex), allocatable :: spectrum(:)
    type(C_PTR) :: plan
    real(kind=8) :: rand_r, rand_i
    integer :: k, nk, seed_used
    real(kind=8) :: amplitude, envelope

    ! Variables for random seed initialization
    integer :: seed(8), i

    ! Allocate spectral array size for real FFT (half spectrum + 1)

    nk = n/2 + 1
    allocate(spectrum(nk))

    ! Initialize random seed using system clock
    if (present(seed_input)) then
       seed_used = seed_input
    else
       seed_used = 12345  ! default value
    end if

    call system_clock(count=seed(1))
    do i = 2, 8
       seed(i) = seed(i - 1) + i * seed_used
    end do
    call random_seed(put=seed)

    ! Generate spectrum with stable bounded energy envelope
    do k = 0, nk - 1
       call random_number(rand_r)
       call random_number(rand_i)
       rand_r = 2.0d0 * rand_r - 1.0d0
       rand_i = 2.0d0 * rand_i - 1.0d0

       if (k /= 0) then
          envelope = real(k*k, kind=8) / (1.0d0 + (real(k, kind=8)/k0)**abs(s))
          amplitude = sqrt(envelope)
       else
          amplitude = 0.0d0
       end if

       spectrum(k + 1) = cmplx(rand_r, rand_i, kind=c_double_complex) * amplitude
    end do

    ! Force Nyquist frequency (last mode) to be real for even-sized transforms
    if (mod(n, 2) == 0) spectrum(nk) = cmplx(real(spectrum(nk)), 0.0d0, kind=c_double_complex)

    ! Create and execute inverse FFT plan (complex-to-real)
    plan = fftw_plan_dft_c2r_1d(n, spectrum, noise, FFTW_ESTIMATE)
    call fftw_execute_dft_c2r(plan, spectrum, noise)

    ! Normalize output by grid size
    noise = noise / n

    ! Clean up FFT plan and memory
    call fftw_destroy_plan(plan)
    deallocate(spectrum)

  end subroutine noise_generator_1D

  !---------------------------------------------------------------------------
  !> Initialise le générateur aléatoire de façon reproductible.
  !---------------------------------------------------------------------------
  subroutine init_noise_seed(iseed)
    integer, intent(in) :: iseed
    integer :: n, i
    integer, allocatable :: seed(:)
    real(kind=8) :: r

    call random_seed(size=n)
    allocate(seed(n))
    do i = 1, n
       seed(i) = iseed + 37 * (i - 1)
    end do
    call random_seed(put=seed)
    deallocate(seed)
    ! Rodage du générateur (graines faiblement entropiques)
    do i = 1, 10000
       call random_number(r)
    end do
  end subroutine init_noise_seed

  !---------------------------------------------------------------------------
  !> Spectre d'énergie modèle (forme seule, la normalisation est fixée
  !> ensuite par urms) : Passot-Pouquet
  !>   E(k) ~ (k/kp)^4 exp(-2 (k/kp)^2)
  !> Remplacer ici pour changer de spectre (von Karman-Pao, etc.).
  !---------------------------------------------------------------------------
  pure function model_spectrum(k, kpeak) result(e)
    real(kind=8), intent(in) :: k, kpeak
    real(kind=8) :: e
    e = (k / kpeak)**4 * exp(-2.d0 * (k / kpeak)**2)
  end function model_spectrum

  !---------------------------------------------------------------------------
  !> Bruit blanc gaussien (Box-Muller), moyenne nulle, variance 1.
  !---------------------------------------------------------------------------
  subroutine gaussian_noise(f, n)
    integer(kind=8), intent(in) :: n
    real(kind=8), intent(out) :: f(*)
    integer(kind=8) :: m
    real(kind=8) :: r1, r2, twopi

    twopi = 2.d0 * acos(-1.d0)
    do m = 1, n
       call random_number(r1)
       call random_number(r2)
       r1 = max(r1, 1.d-300)
       f(m) = sqrt(-2.d0 * log(r1)) * cos(twopi * r2)
    end do
  end subroutine gaussian_noise

  !---------------------------------------------------------------------------
  !> Champ de vitesse solénoïdal (div u = 0 spectralement), isotrope,
  !> à spectre E(k) imposé, de moyenne nulle et de rms = urms.
  !>
  !> INPUT:
  !> nx, ny, nz   : nombre de points
  !> lx, ly, lz   : longueurs du domaine périodique (lx = nx*dx)
  !> kpeak        : nombre d'onde du pic du spectre
  !> kcut         : nombre d'onde de coupure (modes |k| > kcut annulés)
  !> urms         : vitesse rms cible, urms = sqrt(<u_i u_i>/3)
  !>
  !> OUTPUT:
  !> ux, uy, uz   : composantes de la vitesse
  !---------------------------------------------------------------------------
  subroutine generate_solenoidal_noise(ux, uy, uz, nx, ny, nz, &
       lx, ly, lz, kpeak, kcut, urms)
    real(kind=8), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer, intent(in) :: nx, ny, nz
    real(kind=8), intent(in) :: lx, ly, lz, kpeak, kcut, urms

    integer :: nxh, i, j, k, ic, jj, kk3
    integer(kind=8) :: ntot
    type(C_PTR) :: pr, pc, plan_fwd, plan_bwd
    real(C_DOUBLE), pointer :: rbuf(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: cbuf(:,:,:)
    complex(C_DOUBLE_COMPLEX), allocatable :: uh(:,:,:,:)
    complex(C_DOUBLE_COMPLEX) :: kdotu
    real(kind=8) :: twopi, kfx, kfy, kfz, kx, ky, kz, kmag, amp, nrm, fac, rms
    logical :: nyq

    twopi = 2.d0 * acos(-1.d0)
    kfx = twopi / lx
    kfy = twopi / ly
    kfz = twopi / lz
    nxh = nx / 2 + 1
    ntot = int(nx, 8) * int(ny, 8) * int(nz, 8)

    ! Tampons FFTW (FFTW est "row-major" : on inverse l'ordre des dimensions)
    pr = fftw_alloc_real(int(ntot, C_SIZE_T))
    call c_f_pointer(pr, rbuf, [nx, ny, nz])
    pc = fftw_alloc_complex(int(nxh, C_SIZE_T) * int(ny, C_SIZE_T) * int(nz, C_SIZE_T))
    call c_f_pointer(pc, cbuf, [nxh, ny, nz])
    plan_fwd = fftw_plan_dft_r2c_3d(nz, ny, nx, rbuf, cbuf, FFTW_ESTIMATE)
    plan_bwd = fftw_plan_dft_c2r_3d(nz, ny, nx, cbuf, rbuf, FFTW_ESTIMATE)

    allocate(uh(nxh, ny, nz, 3))

    ! 1-2. Bruit blanc -> espace spectral (une réalisation par composante)
    do ic = 1, 3
       call gaussian_noise(rbuf, ntot)
       call fftw_execute_dft_r2c(plan_fwd, rbuf, cbuf)
       uh(:,:,:,ic) = cbuf
    end do

    ! 3-4. Projection solénoïdale + amplitude imposée (phases aléatoires)
    do k = 1, nz
       kk3 = k - 1
       if (kk3 > nz / 2) kk3 = kk3 - nz
       kz = kk3 * kfz
       do j = 1, ny
          jj = j - 1
          if (jj > ny / 2) jj = jj - ny
          ky = jj * kfy
          do i = 1, nxh
             kx = (i - 1) * kfx
             kmag = sqrt(kx*kx + ky*ky + kz*kz)

             ! Modes de Nyquist (n pair) : annulés pour préserver div u = 0
             nyq = (mod(nx,2) == 0 .and. i - 1 == nx / 2) .or. &
                  (mod(ny,2) == 0 .and. j - 1 == ny / 2) .or. &
                  (mod(nz,2) == 0 .and. k - 1 == nz / 2)

             if (kmag < 1.d-12 .or. kmag > kcut .or. nyq) then
                uh(i,j,k,:) = (0.d0, 0.d0)
                cycle
             end if

             kdotu = (kx * uh(i,j,k,1) + ky * uh(i,j,k,2) + kz * uh(i,j,k,3)) &
                  / (kmag * kmag)
             uh(i,j,k,1) = uh(i,j,k,1) - kx * kdotu
             uh(i,j,k,2) = uh(i,j,k,2) - ky * kdotu
             uh(i,j,k,3) = uh(i,j,k,3) - kz * kdotu

             nrm = sqrt(abs(uh(i,j,k,1))**2 + abs(uh(i,j,k,2))**2 &
                  + abs(uh(i,j,k,3))**2)
             ! Énergie par mode ~ E(k) / (4 pi k^2)
             amp = sqrt(model_spectrum(kmag, kpeak) / (4.d0 * acos(-1.d0) * kmag**2))
             if (nrm > 1.d-30) then
                uh(i,j,k,:) = uh(i,j,k,:) * (amp / nrm)
             else
                uh(i,j,k,:) = (0.d0, 0.d0)
             end if
          end do
       end do
    end do

    ! 5. Retour dans l'espace physique
    do ic = 1, 3
       cbuf = uh(:,:,:,ic)
       call fftw_execute_dft_c2r(plan_bwd, cbuf, rbuf)
       select case (ic)
       case (1)
          ux(1:nx,1:ny,1:nz) = rbuf
       case (2)
          uy(1:nx,1:ny,1:nz) = rbuf
       case (3)
          uz(1:nx,1:ny,1:nz) = rbuf
       end select
    end do

    ! Renormalisation : urms = sqrt(<u_i u_i> / 3)
    rms = sqrt((sum(ux**2) + sum(uy**2) + sum(uz**2)) / (3.d0 * real(ntot, 8)))
    if (rms > 0.d0) then
       fac = urms / rms
       ux = ux * fac
       uy = uy * fac
       uz = uz * fac
    end if

    deallocate(uh)
    call fftw_destroy_plan(plan_fwd)
    call fftw_destroy_plan(plan_bwd)
    call fftw_free(pr)
    call fftw_free(pc)

  end subroutine generate_solenoidal_noise

end module noise_module

