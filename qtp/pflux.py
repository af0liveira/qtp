import numpy as np
from scipy import integrate

class PFlux(object):
    """Particle flux (current) through a 1D potential energy barrier.

    Parameters
    ----------
    traco : TransCoeff instance
        Calculator for the particle's transmission coefficient T(E).
    """

    def __init__(self, traco, **kwargs):
        self.traco = traco
        self.peb = traco.peb
        self.pmass = traco.pmass
        self.z_Umax, self.Umax = self.peb.get_Umax()

    def __call__(self, beta):
        """Calculate the classical and non-classical flux components.

        The classical flux is calculated as

        .. math::
            j_c(\\beta) = \\exp(-\\beta U_{max}) / \\sqrt{2 \\pi m \\beta}

        whereas the non-classical flux is

        .. math::
            j_q(\\beta) = \\sqrt{frac{\\beta}{2 \\pi m}} \\int_{\\infty}^{0}
            T(U(z)) \\exp(-\\beta U(z)) \frac{\\partial{U}}{\partial{z}} dz

        The total flux is, thus,

        .. math::
            j_{tot} = j_c + j_q

        Parameters
        ----------
        beta : float
            Inverse temperature in 1/K.

        Return
        ------
        j_c, j_q, j_c + j_q : tuple of floats
        """
        _zero = np.finfo(float).resolution # on a MacBook Pro 2017, this is ca. 1e-15
        _tiny= 1e-6

        def fun(x):
            traco = self.traco(x)[1]
            U = self.peb(x)
            dU = self.peb(x, der=1)
            return traco * np.exp(-beta * U) * dU

        zlims = self.peb.get_zfromUvalue(_zero)
        assert(zlims[0] < 0), "Problems finding integration limits."
        zlims = [zlims[0], self.z_Umax]

        intgl = integrate.quad(fun, zlims[0], zlims[1])

        j_q = np.sqrt(beta/(2*np.pi * self.pmass)) * intgl[0]
        j_c = 1/np.sqrt(2*np.pi * self.pmass * beta) * np.exp(-beta * self.Umax)

        return j_c, j_q, j_c + j_q