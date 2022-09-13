import numpy as np
from scipy.interpolate import UnivariateSpline

class Arrhenius(object):
    """Calculator of Arrhenius activation energies and coefficients.
    
    Parameters
    ----------
    betas : list
        Values of 1/T in atomic units.
    ln_ks : list
        Values of log(k) corresponding to `betas`.
    """

    def __init__(self, betas, ln_ks):
        
        # make sure betas in ascending order
        self.betas, self.ln_ks = zip(*sorted(zip(betas, ln_ks)))
        self.fun = UnivariateSpline(self.betas, self.ln_ks, s=0)

    def __call__(self, beta):
        """Return the activation energy and coefficient for beta."""
        k = np.exp(self.fun(beta))
        e_act = -self.fun(beta, nu=1)
        coeff = k*np.exp(beta*e_act)
        return e_act, coeff