'''Sample k-space trajectories.'''

import numpy as np

def radial(sx, spokes):
    '''k-space coordinates for a radial trajectory.

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
    for ii in range(spokes):
        theta = ii/spokes*np.pi
        ct = np.cos(theta)
        st = np.sin(theta)
        kx[ii*sx:(ii+1)*sx] = x*ct - y*st
        ky[ii*sx:(ii+1)*sx] = x*st + y*ct
    return(kx, ky)
