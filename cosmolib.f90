! vim: set filetype=fortran
module cosmolib
    ! class to calculate distances.
    !
    ! This uses gauss-legendre integration extremely fast and accurate
    !
    ! For 1/E(z) integration, 5 points is good to 1.e-8

    implicit none

    ! class variables
    integer*8, save, private :: has_been_init = 0
    integer*8, save :: npts
    integer*8, save :: vnpts
    real*8, private, save, dimension(:), allocatable :: xxi, wwi
    real*8, private, save, dimension(:), allocatable :: vxxi, vwwi

    real*8, save :: H0
    real*8, save :: omega_m
    real*8, save :: omega_l
    real*8, save :: omega_k
    logical, save :: flat

    real*8, save :: sqrt_omega_k_over_DH
    real*8, save :: sqrt_omega_k_over_DH_inv


    ! use in scinv for dlens in Mpc
    real*8, parameter :: four_pi_G_over_c_squared = 6.0150504541630152e-07_8
    ! The hubble distance c/H0
    real*8, parameter :: c = 2.99792458e5_8
    real*8, save :: DH

    ! for integral calculations
    real*8, private :: f1,f2,z,ezinv


    ! for integration
    real*8, parameter, public :: M_PI    = 3.141592653589793238462643383279502884197_8


contains

    ! you must initialize
    subroutine cosmo_init(flat_new, H0_new, omega_m_new, &
                          npts_new, vnpts_new, &
                          omega_k_new, omega_l_new)

        logical, intent(in) :: flat_new

        real*8, intent(in) :: H0_new, omega_m_new
        integer*8, intent(in) :: npts_new, vnpts_new

        ! we only need these if the universe is not flat
        real*8, intent(in), optional :: omega_l_new, omega_k_new

        flat    = flat_new
        H0      = H0_new
        omega_m = omega_m_new
        npts    = npts_new
        vnpts    = vnpts_new

        DH = c/H0

        if (flat) then
            omega_l = 1.0-omega_m
            sqrt_omega_k_over_DH = 0
            sqrt_omega_k_over_DH_inv = 0
        else
            ! it will be up to the python wrapper to ensure the right inputs
            if (present(omega_l_new) .neqv. .true.) then
                print '(a)','You must enter omega_k if not flat'
                call exit(45)
            endif
            if (present(omega_l_new) .neqv. .true.) then
                print '(a)','You must enter omega_l if not flat'
                call exit(45)
            endif
            omega_k=omega_k_new
            omega_l = omega_l_new

            if (omega_k > 0) then
                sqrt_omega_k_over_DH = sqrt(omega_k)/DH
            else
                sqrt_omega_k_over_DH = sqrt(-omega_k)/DH
            endif
            sqrt_omega_k_over_DH_inv = 1./sqrt_omega_k_over_DH 

        endif

        call set_cosmo_weights(npts, vnpts)

        has_been_init = 1

    end subroutine cosmo_init


    ! comoving distance
    ! variants for a vector of zmax and zmin,zmax both vectors
    real*8 function cdist(zmin, zmax)
        ! comoving distance
        real*8, intent(in) :: zmin, zmax
        cdist = DH*ez_inverse_integral(zmin, zmax)
    end function cdist

    subroutine cdist_vec(zmin, zmax, n, dc)
        real*8, intent(in) :: zmin
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: dc

        integer*8 i

        do i=1,n
            dc(i) = DH*ez_inverse_integral(zmin, zmax(i))
        end do
    end subroutine cdist_vec

    subroutine cdist_2vec(zmin, zmax, n, dc)
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: dc

        integer*8 i

        do i=1,n
            dc(i) = cdist(zmin(i), zmax(i))
        end do
    end subroutine cdist_2vec



    real*8 function tcdist(zmin, zmax) result(d)
        ! useful for calculating transverse comoving distance at zmax.  When
        ! zmin is not zero, useful in angular diameter distances

        real*8, intent(in) :: zmin, zmax

        d = cdist(zmin, zmax)

        if (flat) then
             !just use comoving distance already calculated
        else if (omega_k > 0) then
            d = sinh( d*sqrt_omega_k_over_DH )*sqrt_omega_k_over_DH_inv
        else
            d =  sin( d*sqrt_omega_k_over_DH )*sqrt_omega_k_over_DH_inv
        endif

    end function tcdist

    subroutine tcdist_vec(zmin, zmax, n, dm)
        integer*8, intent(in) :: n
        real*8, intent(in) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: dm

        integer*8 i

        do i=1,n
            dm(i) = tcdist(zmin, zmax(i))
        enddo

    end subroutine tcdist_vec

    subroutine tcdist_2vec(zmin, zmax, n, dm)
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: dm

        integer*8 i

        do i=1,n
            dm(i) = tcdist(zmin(i),zmax(i))
        enddo

    end subroutine tcdist_2vec




    real*8 function angdist(zmin, zmax) result(d)
        ! angular diameter distance
        real*8, intent(in) :: zmin, zmax
        if (flat) then
            ! this is just the comoving distance over 1+zmax
            d = DH*ez_inverse_integral(zmin, zmax)/(1.+zmax)
        else
            d = tcdist(zmin, zmax)/(1.+zmax)
        endif
    end function angdist

    subroutine angdist_vec(zmin, zmax, n, da)
        integer*8, intent(in) :: n
        real*8, intent(in) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: da

        integer*8 i

        do i=1,n
            da(i) = angdist(zmin, zmax(i))
        enddo

    end subroutine angdist_vec

    subroutine angdist_2vec(zmin, zmax, n, da)
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: da

        integer*8 i

        do i=1,n
            da(i) = angdist(zmin(i),zmax(i))
        enddo

    end subroutine angdist_2vec


    real*8 function lumdist(zmin, zmax) result(d)
        ! angular diameter distance
        real*8, intent(in) :: zmin, zmax
        d = angdist(zmin, zmax)*(1.+zmax)**2
    end function lumdist

    subroutine lumdist_vec(zmin, zmax, n, d)
        integer*8, intent(in) :: n
        real*8, intent(in) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: d

        integer*8 i

        do i=1,n
            d(i) = angdist(zmin, zmax(i))*(1.+zmax(i))**2
        enddo

    end subroutine lumdist_vec

    subroutine lumdist_2vec(zmin, zmax, n, d)
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: d

        integer*8 i

        do i=1,n
            d(i) = angdist(zmin(i),zmax(i))*(1.+zmax(i))**2
        enddo

    end subroutine lumdist_2vec


    real*8 function dv(z)
        ! comoving volume element at redshift z
        real*8, intent(in) :: z
        real*8 ezinv, da

        da = angdist(0.0_8, z)
        ezinv = ez_inverse(z)

        dv = DH*da**2*ezinv*(1.0+z)**2

    end function dv

    subroutine dv_vec(z, n, dvvec)
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: z
        real*8, intent(inout), dimension(n) :: dvvec
        integer*8 i

        do i=1,n
            dvvec(i) = dv( z(i) ) 
        enddo
    end subroutine dv_vec




    real*8 function scinv(zl, zs)
        ! inverse critical density
        real*8, intent(in) :: zl,zs

        real*8 dl, ds, dls

        if (zs <= zl) then 
            scinv=0.0
            return
        end if


        if (flat) then
            ! we can save some computation in the flat case
            dl = cdist(0.0_8, zl)
            ds = cdist(0.0_8, zs)
            scinv = dl/(1.+zl)*(ds-dl)/ds * four_pi_G_over_c_squared
        else
            dl  = angdist(0.0_8, zl)
            ds  = angdist(0.0_8, zs)
            dls = angdist(zl, zs)

            scinv = dls*dl/ds*four_pi_G_over_c_squared
        endif

    end function scinv

    subroutine scinv_vec(zmin, zmax, n, sc_inv)
        integer*8, intent(in) :: n
        real*8, intent(in) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: sc_inv

        integer*8 i

        do i=1,n
            sc_inv(i) = scinv(zmin, zmax(i))
        enddo

    end subroutine scinv_vec

    subroutine scinv_2vec(zmin, zmax, n, sc_inv)
        integer*8, intent(in) :: n
        real*8, intent(in), dimension(n) :: zmin
        real*8, intent(in), dimension(n) :: zmax
        real*8, intent(inout), dimension(n) :: sc_inv

        integer*8 i

        do i=1,n
            sc_inv(i) = scinv(zmin(i),zmax(i))
        enddo

    end subroutine scinv_2vec




    real*8 function ez_inverse(z)
        real*8, intent(in) :: z

        if (flat) then
            ez_inverse = omega_m*(1.+z)**3 + omega_l
        else
            ez_inverse = omega_m*(1.+z)**3 + omega_k*(1.+z)**2 + omega_l
        endif
        ez_inverse = sqrt(1.0/ez_inverse)
    end function ez_inverse

    subroutine ez_inverse_vec(z, n, ez)
        integer*8, intent(in) :: n
        real*8, dimension(n), intent(in) :: z
        real*8, dimension(n), intent(inout) :: ez
        integer*8 i

        do i=1,n
            ez(i) = ez_inverse( z(i) ) 
        enddo
    end subroutine ez_inverse_vec


    real*8 function ez_inverse_integral(zmin, zmax) result(val)
        real*8, intent(in) :: zmin, zmax
        integer*8 i


        f1 = (zmax-zmin)/2.
        f2 = (zmax+zmin)/2.

        val = 0.0

        do i=1,npts
            z = xxi(i)*f1 + f2
            ezinv = ez_inverse(z)

            val = val + f1*ezinv*wwi(i);
        end do

    end function ez_inverse_integral


    real*8 function volume(zmin, zmax)
        real*8, intent(in) :: zmin, zmax
        real*8 f1, f2, z, v

        integer*8 i


        f1 = (zmax-zmin)/2.
        f2 = (zmax+zmin)/2.

        volume = 0
        do i=1,vnpts
            z = vxxi(i)*f1 + f2
            volume = volume + f1*dv(z)*vwwi(i)
        enddo
    end function volume















    subroutine set_cosmo_weights(npts_new, vnpts_new)

        integer*8, intent(in) :: npts_new
        integer*8, intent(in) :: vnpts_new
        npts = npts_new
        vnpts = vnpts_new

        if (allocated(xxi)) then
            deallocate(xxi)
        endif
        allocate(xxi(npts)); xxi=0
        if (allocated(wwi)) then
            deallocate(wwi)
        endif
        allocate(wwi(npts)); wwi=0

        if (allocated(vxxi)) then
            deallocate(vxxi)
        endif
        allocate(vxxi(vnpts)); vxxi=0

        if (allocated(vwwi)) then
            deallocate(vwwi)
        endif
        allocate(vwwi(vnpts)); vwwi=0

        call gauleg(-1.0_8, 1.0_8, npts, xxi, wwi)
        call gauleg(-1.0_8, 1.0_8, vnpts, vxxi, vwwi)

    end subroutine set_cosmo_weights

    subroutine print_weights()
        integer*8 i

        do i=1,size(wwi)
            print '("xxi: ",F15.8,"  wwi: ",F15.8)',xxi(i),wwi(i)
        enddo
    end subroutine


    ! from numerical recipes
    subroutine gauleg(x1, x2, npts, x, w)

        integer*8, intent(in) :: npts
        real*8, intent(in) :: x1, x2
        real*8, intent(inout), dimension(npts) :: x, w
        

        integer*8 :: i, j, m
        real*8 :: xm, xl, z1, z, p1, p2, p3, pp, EPS, abszdiff


        pp = 0.0
        EPS = 4.e-11

        m = (npts + 1)/2

        xm = (x1 + x2)/2.0
        xl = (x2 - x1)/2.0
        z1 = 0.0

        do i=1,m

            z=cos( M_PI*(i-0.25)/(npts+.5) )

            abszdiff = abs(z-z1)

            do while (abszdiff > EPS) 

                p1 = 1.0
                p2 = 0.0
                do j=1,npts
                    p3 = p2
                    p2 = p1
                    p1 = ( (2.0*j - 1.0)*z*p2 - (j-1.0)*p3 )/j
                end do
                pp = npts*(z*p1 - p2)/(z*z -1.)
                z1=z
                z=z1 - p1/pp

                abszdiff = abs(z-z1)

            end do

            x(i) = xm - xl*z
            x(npts+1-i) = xm + xl*z
            w(i) = 2.0*xl/( (1.-z*z)*pp*pp )
            w(npts+1-i) = w(i)


        end do

    end subroutine gauleg


end module cosmolib
