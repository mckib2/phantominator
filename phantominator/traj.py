'''Sample k-space trajectories.'''

from functools import partial

import numpy as np

def radial(sx, spokes, golden=False):
    '''k-space coordinates for a radial trajectory.

    Parameters
    ----------
    sx : int
        Number of samples along spoke.
    spokes : int
        Number of spokes.
    golden : bool, optional
        Use Golden angle ratio sampling.

    Returns
    -------
    kx, ky : array_like
        k-space coordinates for sampling along a radial trajectory.

    Notes
    -----
    This is a very simple radial trajectory mainly for demonstration
    purposes.
    '''

    N = sx*spokes
    kx = np.zeros(N)
    ky = np.zeros(N)
    x = np.linspace(-1, 1, sx)*sx/2 # scale fac to match BART's traj
    y = np.zeros(sx)

    # How we get
    if golden:
        GA = np.pi*(3 - np.sqrt(5)) # calculate GA for use if needed
        gettheta = partial(_nextspoke_golden, GA=GA)
    else:
        gettheta = partial(_nextspoke, spokes=spokes)

    for ii in range(spokes):
        # theta = ii/spokes*np.pi
        theta = gettheta(ii)
        ct = np.cos(theta)
        st = np.sin(theta)
        kx[ii*sx:(ii+1)*sx] = x*ct - y*st
        ky[ii*sx:(ii+1)*sx] = x*st + y*ct
    return(kx, ky)

def _nextspoke(ii, spokes):
    return ii/spokes*np.pi

def _nextspoke_golden(ii, GA):
    return GA*ii
