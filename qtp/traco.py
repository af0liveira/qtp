import warnings

import numpy as np

from scipy import integrate
from scipy.interpolate import UnivariateSpline

class TransCoeff(object):
    """Transmission coefficient for a particle through a 1D potential barrier.

    Parameters
    ----------
    peb : PEB instance
        The representation of the 1D potential energy barrier.

    pmass : float
        Mass in m_e units of the particle tunneling through the potential barrier.
        1 m_e = 9.10938356e-31 kg.
    """

    def __init__(self, peb, pmass, **kwargs):
        self.peb = peb
        self.pmass = pmass

    def _get_traco_type01(self, zval):
        """Calculate T(E) for a symmetric PEB with single maximum at z=0."""

        _zero = np.finfo(float).resolution # on a MacBook Pro 2017, this is ca. 1e-15
        _tiny= 1e-6

        def fun(x):
            U = self.peb(x)
            E = self.peb(zval)
            if E > U or abs(E-U) < _tiny:
                ret = 0
            else:
                ret = np.sqrt(U - E)
            return ret

        if abs(zval) < _tiny:
            ln_traco = 0
        else:
            zlims = self.peb.get_zfromUvalue(np.round(self.peb(zval), 8))
            if len(zlims) != 2 and len(zlims) > 0:
                z_ = np.max(np.abs(zlims))
                zlims = [-z_, z_]
            assert len(zlims) == 2, f"Unexpected number of points with same energy. {zlims, zval}"
            intgl = integrate.quad(fun, zlims[0], zlims[1], limit=500)
            ln_traco = -2*np.sqrt(2*self.pmass) * intgl[0]

        return ln_traco, np.exp(ln_traco)

    def _get_traco_type02(self, zval):
        """Calculate T(E) for a symmetric PEB with two maxima when E > U(0).
        
        In this case, T(E) is scaled by 1/2 to represent a 50% chance of
        backscattering after passing throught the first barrier.
        """

        _zero = np.finfo(float).resolution # on a MacBook Pro 2017, this is ca. 1e-15
        _tiny= 1e-6

        def fun(x):
            U = self.peb(x)
            E = self.peb(zval)
            if E > U or abs(E-U) < _tiny:
                ret = 0
            else:
                ret = np.sqrt(U - E)
            return ret

        z_Umax, Umax = self.peb.get_Umax()

        if abs(self.peb(zval) - Umax) < _tiny:
            ln_traco = np.log(0.5)
        else:
            zlims = [z for z in self.peb.get_zfromUvalue(np.round(self.peb(zval), 8)) if z <= 0]
            limit_case = len(zlims)==1 \
                            and abs(self.peb(0) - self.peb(zlims[0])) < _tiny \
                            and zlims[0] <= z_Umax
            if limit_case:
                zlims.append(0)
            assert len(zlims) == 2,\
                    f"Unexpected number of points with same energy for z <= 0. (zval: {zval})"
            intgl = integrate.quad(fun, zlims[0], zlims[1], limit=500)
            ln_traco = np.log(0.5) - 2*np.sqrt(2*self.pmass) * intgl[0]

        return ln_traco, np.exp(ln_traco)

    def __call__(self, zval):
        """Calculate the transmission coefficient for E = U(zval).

        Parameters
        ----------
        zval : float
            The value of z (in bohr) for which E = U(z).
            This will be used in conjunction with `peb` to determine the
            integration limits when calculating T(E).

        Return
        ------
        ln(T(E)), T(E)
        """
        _zero = np.finfo(float).resolution # on a MacBook Pro 2017, this is ca. 1e-15
        _tiny= 1e-6

        U0 = self.peb(0)    # potential energy at z=0
        E = self.peb(zval)  # particle's energy
        z_Umax, Umax = self.peb.get_Umax()

        # We do bold assumptions here:
        # - it's always assumed that U(z) = U(-z) and U(z) >= 0;
        # - if max. at U(z=0), assume no further maxima;
        # - if max. not at U(z=0), assume a single local minimum at z=0.
        if abs(z_Umax) < _tiny:
            peb_type = 'type01'
        else:
            if E < U0:
                peb_type = 'type01'
            else:
                peb_type = 'type02'
        
        assert peb_type in self.__class__._get_traco, \
                f"Unrecognized PEB type: '{peb_type}'"

        return self._get_traco[peb_type](self, zval)

    _get_traco = {'type01': _get_traco_type01,
                  'type02': _get_traco_type02}