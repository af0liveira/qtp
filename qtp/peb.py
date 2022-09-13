import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize

class PEB(object):
    """Class for representing 1D potential energy barriers
    
    Parameters
    ----------
    zvals : list
        List of coordinates in bohr
    Uvals : list
        List of potential energies U(z) in hartree; must correspond to
        `zvals`

    Notes
    -----
    The input data must be compatible to a symmetry potential barrier, such
    that $U(-z) = U(z)$.
    The input data will be shifted such that Umin = 0; in addition, the data
    will be symmetrized to ensure the symmetry of U(z).
    For z values outside the `zvals` range, U(z) := 0.
    """

    def __init__(self, zvals, Uvals):

        assert len(zvals) == len(Uvals), \
                "Uvals and zvals must have the same number of items!"

        def symmetrize_data(zvals, Uvals):
            _tiny = 1.0e-6
            pts = [(abs(z), U) for z, U in zip(zvals, Uvals)]
            right_pts = sorted(set(pts))
            tmp_ = []
            for x in set([z for z,u in right_pts]):
                sel = [(z, u) for z, u in right_pts if abs(x-z) < _tiny]
                zz, uu = zip(*sel)
                tmp_.append((np.mean(zz), np.mean(uu)))
            right_pts = tmp_
            left_pts = [(-z, U) for z, U in right_pts if abs(z) > _tiny]
            all_pts = sorted(left_pts + right_pts)
            znew, Unew = zip(*all_pts)
            Unew = list(np.array(Unew) - np.min(Unew))
            return znew, Unew

        self.zvals, self.Uvals = symmetrize_data(zvals, Uvals)
        self.pebspl = UnivariateSpline(self.zvals, self.Uvals,
                                       ext='zeros', k=3, s=0)

    def __call__(self, z, der=0, ext=1):
        """Evaluate PEB or its derivatives at point z.
    
        Parameters
        ----------
        z : array-like
            Point(s) to the evaluate the function at.

        Other parameters
        ----------------
        der : int, optional
            Order of the derivative of the spline.

        ext : int, optional
            Control the value returned for z values out of the interval
            defined in the knot sequence.
            0: extrapolate; 1: return zero; 2: raise ValueError;
            3: return boundary value.
        """
        return self.pebspl(z, nu=der, ext=ext)

    def get_zfromUvalue(self, U):
        """Return the value(s) of z that give the requested U value."""
        zvals = np.array(self.zvals)
        uvals = np.array(self.Uvals) - U
        spl = UnivariateSpline(zvals, uvals, ext='zeros', k=3, s=0)
        return spl.roots()

    def get_Umax(self):
        """Find point of maximum in the potential energy barrier.

        Return
        ------
        (z, U(z)) tuple for the point of maximum U(z) value.
        """

        _small = 1.0e-4

        # The `root()` method is only supported in order 3 splines; 
        # thus, we need to make sure that d1 is order 3.
        spl = UnivariateSpline(self.zvals, self.Uvals, ext='zeros', k=4, s=0)
        d1 = spl.derivative(n=1)
        d2 = spl.derivative(n=2)

        critical_pts = d1.roots()
        minmax_test = d2(critical_pts)

        max_pts = [x for x, t in zip(critical_pts, minmax_test) if t < -1*_small]
        max_vals = spl(max_pts)

        max_idx = max_vals.argmax()
        return max_pts[max_idx], max_vals[max_idx]


        



