from __future__ import print_function
from . import _cosmolib

import numpy
from numpy import isscalar, log10

class Cosmo(dict):
    def __init__(self, 
                 omega_m=0.3, 
                 omega_l=0.7,
                 omega_k=0.0,
                 H0=100.0,
                 flat=True,
                 npts=5,
                 vnpts=10):

        # If flat is specified, make sure omega_l = 1-omega_m
        # and omega_k=0
        omega_m, omega_l, omega_k = \
                self.extract_omegas(omega_m,omega_l,omega_k,flat)


        self['H0'] = H0 
        self['omega_m'] = omega_m
        self['omega_l'] = omega_l
        self['omega_k'] = omega_k
        self['flat'] = flat
        self['npts'] = npts
        self['vnpts'] = vnpts
        _cosmolib.cosmolib.cosmo_init(flat,H0, omega_m, npts, vnpts, omega_k, omega_l)
        self['DH'] = _cosmolib.cosmolib.dh

        self.Distmod = self.distmod

    def extract_omegas(self, omega_m, omega_l, omega_k, flat):
        """
        If flat is specified, make sure omega_l = 1-omega_m
        and omega_k=0
        """
        if flat or omega_k == 0:
            omega_l = 1.0-omega_m
            omega_k = 0.0
        return omega_m, omega_l, omega_k

    def DH(self):
        """
        Return the Hubble distance in Mpc/h.  You can also access this
        with
            c = cosmology.Cosmo()
            c['DH']
        """

        return self['DH']



    def Dc(self, zmin, zmax):
        """
        Calculate the comoving distance from zmin to zmax in units of Mpc.

        Parameters
        ----------
        zmin, zmax: scalars or arrays
            The following combinations are supported
                1) Two scalars
                2) zmin a scalar and zmax an array
                2) Both arrays of the same length.

        """

        if isscalar(zmin) and isscalar(zmax):
            # two scalars of any kind.
            dc = _cosmolib.cosmolib.cdist(zmin, zmax)

        elif isscalar(zmin) and not isscalar(zmax):
            # scalar for zmin, array for zmax
            dc=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.cdist_vec(zmin, zmax, dc)

        elif not isscalar(zmin) and not isscalar(zmax):
            # both arrays: must be same length
            if len(zmin) != len(zmax):
                raise ValueError("If zmin and zmax are arrays, they must be same length")
            dc=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.cdist_2vec(zmin, zmax, dc)
        else:
            raise ValueError("zmin,zmax should be two scalars, zmin scalar zmax array, or both arrays")

        return dc

    def Dm(self, zmin, zmax):
        """
        Calculate the transvers comoving distance from zmin to zmax in units of Mpc.


        Useful for calculating transverse comoving distance at zmax.  When zmin
        is not zero, useful in calculating angular diameter distances

        Parameters
        ----------
        zmin, zmax: scalars or arrays
            The following combinations are supported
                1) Two scalars
                2) zmin a scalar and zmax an array
                2) Both arrays of the same length.

        """


        if isscalar(zmin) and isscalar(zmax):
            # two scalars of any kind.
            dm = _cosmolib.cosmolib.tcdist(zmin, zmax)

        elif isscalar(zmin) and not isscalar(zmax):
            # scalar for zmin, array for zmax
            dm=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.tcdist_vec(zmin, zmax, dm)

        elif not isscalar(zmin) and not isscalar(zmax):
            # both arrays: must be same length
            if len(zmin) != len(zmax):
                raise ValueError("If zmin and zmax are arrays, they must be same length")
            dm=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.tcdist_2vec(zmin, zmax, dm)
        else:
            raise ValueError("zmin,zmax should be two scalars, zmin scalar zmax array, or both arrays")

        return dm

    def Da(self, zmin, zmax):
        """
        Calculate the angular diameter distance from zmin to zmax in units of Mpc.


        Parameters
        ----------
        zmin, zmax: scalars or arrays
            The following combinations are supported
                1) Two scalars
                2) zmin a scalar and zmax an array
                2) Both arrays of the same length.

        """


        if isscalar(zmin) and isscalar(zmax):
            # two scalars of any kind.
            da = _cosmolib.cosmolib.angdist(zmin, zmax)

        elif isscalar(zmin) and not isscalar(zmax):
            # scalar for zmin, array for zmax
            da=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.angdist_vec(zmin, zmax, da)

        elif not isscalar(zmin) and not isscalar(zmax):
            # both arrays: must be same length
            if len(zmin) != len(zmax):
                raise ValueError("If zmin and zmax are arrays, they must be same length")
            da=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.angdist_2vec(zmin, zmax, da)
        else:
            raise ValueError("zmin,zmax should be two scalars, zmin scalar zmax array, or both arrays")

        return da

    def Dl(self, zmin, zmax):
        """
        Calculate the luminosity distance from zmin to zmax in units of Mpc.


        Parameters
        ----------
        zmin, zmax: scalars or arrays
            The following combinations are supported
                1) Two scalars
                2) zmin a scalar and zmax an array
                2) Both arrays of the same length.

        """

        if isscalar(zmin) and isscalar(zmax):
            # two scalars of any kind.
            d = _cosmolib.cosmolib.lumdist(zmin, zmax)

        elif isscalar(zmin) and not isscalar(zmax):
            # scalar for zmin, array for zmax
            d=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.lumdist_vec(zmin, zmax, d)

        elif not isscalar(zmin) and not isscalar(zmax):
            # both arrays: must be same length
            if len(zmin) != len(zmax):
                raise ValueError("If zmin and zmax are arrays, they must be same length")
            d=numpy.zeros(len(zmax), dtype='f8')
            _cosmolib.cosmolib.lumdist_2vec(zmin, zmax, d)
        else:
            raise ValueError("zmin,zmax should be two scalars, zmin scalar zmax array, or both arrays")

        return d



    def distmod(self, z):
        """
        Calculate the distance modulus to the given redshift.

        Parameters
        ----------
        z: scalar or array
            The redshift 
        """

        dmpc = self.Dl(0.0, z)
        dpc = dmpc*1.e6
        dm = 5.0*log10(dpc/10.0)
        return dm      


    def sigmacritinv(self, zl, zs):
        """
        Calculate the inverse critical density for the lens and source redshifts


        Parameters
        ----------
        zl, zs: scalars or arrays
            The following combinations are supported
                1) Two scalars
                2) zl a scalar and zs an array
                2) Both arrays of the same length.

        """


        if isscalar(zl) and isscalar(zs):
            # two scalars of any kind.
            scinv = _cosmolib.cosmolib.scinv(zl, zs)

        elif isscalar(zl) and not isscalar(zs):
            # scalar for zl, array for zs
            scinv=numpy.zeros(len(zs), dtype='f8')
            _cosmolib.cosmolib.scinv_vec(zl, zs, scinv)

        elif not isscalar(zl) and not isscalar(zs):
            # both arrays: must be same length
            if len(zl) != len(zs):
                raise ValueError("If zl and zs are arrays, they must be same length")
            scinv=numpy.zeros(len(zs), dtype='f8')
            _cosmolib.cosmolib.scinv_2vec(zl, zs, scinv)
        else:
            raise ValueError("zl,zs should be two scalars, zl scalar zs array, or both arrays")

        return scinv


    def Ez_inverse(self, z):
        """
        Integrate kernel 1/E(z) from 0 to z.
        
        1/E(z) is used for distance calculations in FRW.

        Parameters
        ----------
        z: scalar or array
            The redshift 
        """

        if isscalar(z):
            ez = _cosmolib.cosmolib.ez_inverse(z)
        else:
            ez = numpy.zeros(len(z), dtype='f8')
            _cosmolib.cosmolib.ez_inverse_vec(z, ez)

        return ez

    def Ezinv_integral(self, zmin, zmax):
        """
        Integrate kernel 1/E(z) from zmin to zmax.

        1/E(z) is used for distance calculations in FRW.

        Parameters
        ----------
        zmin,zmax: scalars
            The redshifts
        """

        return _cosmolib.cosmolib.ez_inverse_integral(zmin, zmax)



def test_cosmo(omega_k=None):
    if omega_k is not None:
        flat=False
    else:
        flat=True
        omega_k=0.0

    c=Cosmo(flat=flat, omega_k=omega_k)

    print("  Testing Da,Dl")
    da=c.Da(0.1, 0.5)
    da=c.Da(0.1, [0.4, 0.5])
    da=c.Da([0.1,0.2], [0.4, 0.5])

    dl=c.Dl(0.1, 0.5)
    dl=c.Dl(0.1, [0.4, 0.5])
    dl=c.Dl([0.1,0.2], [0.4, 0.5])

    print("  Testing sicmacrit inverse")
    scinv=c.sigmacritinv(0.1, 0.5)
    scinv=c.sigmacritinv(0.1, [0.4, 0.5])
    scinv=c.sigmacritinv([0.1,0.2], [0.4, 0.5])


def test():
    print("Testing flat")
    test_cosmo()

    omega_k=0.05
    print("Testing non-flat, omega_k:",omega_k)
    test_cosmo(omega_k=omega_k)
