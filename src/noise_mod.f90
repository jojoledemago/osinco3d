module noise_mod
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
contains

  subroutine noise_generator_1D(n, k0, s, noise)
    !> 1D random noise generator in physical space from a spectral distribution
    !> with energy envelope E(k) ∝ k² / (1 + (k/k₀)^|s|)
    !> Allows both positive and negative spectral slope s
    !> Based on inverse FFT using FFTW3 (complex to real transform)
    !>
    !> INPUT:  n      : number of spatial points (grid size)
    !>         k0     : wavenumber at which energy is peaked
    !>         s      : spectral slope (can be negative)
    !> OUTPUT: noise  : real-valued 1D noise in physical space (size n)

    integer, intent(in) :: n
    real(kind=8), intent(in) :: k0, s
    real(kind=8), intent(out) :: noise(n)

    complex(c_double_complex), allocatable :: spectrum(:)
    type(C_PTR) :: plan
    real(c_double) :: rand_r, rand_i
    integer :: k, nk
    real(kind=8) :: amplitude, envelope

    ! Variables for random seed initialization
    integer :: seed(8), i

    ! Allocate spectral array size for real FFT (half spectrum + 1)

    nk = n/2 + 1
    allocate(spectrum(nk))

    ! Initialize random seed using system clock
    call system_clock(count=seed(1))
    do i = 2, 8
       seed(i) = seed(i - 1) + i * 12345
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

end module noise_mod

