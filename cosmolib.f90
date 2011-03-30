! vim: set filetype=fortran et ts=2 sw=2 sts=2 :
module cosmolib
  ! class to calculate distances. This uses gauss-legendre integration
  ! extremely fast and accurate
  !
  ! For integration, 5 points is essentially exact and very fast.

  implicit none

  ! class variables
  integer*8, save, private :: has_been_init = 0
  integer*8, save :: npts
  real*8, private, save, dimension(:), allocatable :: xxi, wwi

  real*8, save :: H0
  real*8, save :: omega_m
  real*8, save :: omega_l

  ! The hubble distance c/H0
  real*8, save :: DH

  ! use in scinv for dlens in Mpc
  real*8, parameter :: four_pi_G_over_c_squared = 6.0150504541630152e-07

  ! for integral calculations
  real*8, private :: f1,f2,z,ezinv

contains

  ! you must initialize
  subroutine cosmo_init(H0_new, omega_m_new, npts_new)

    real*8, intent(in) :: H0_new, omega_m_new
    integer*8, intent(in) :: npts_new

    H0      = H0_new
    omega_m = omega_m_new
    omega_l = 1.0-omega_m
    npts    = npts_new

    DH = 2.99792458e5/100.0

    !call set_cosmo_weights(npts)

    has_been_init = 1

  end subroutine cosmo_init

end module cosmolib
