'''Simple dynamic numerical phantom.'''

import numpy as np

def dynamic(N, nt):
    '''Concentric circles with dynamic movement.

    Parameters
    ----------
    N : int
        In-plane dimension.
    nt : int
        Number of time points.

    Returns
    -------
    ph : array_like
        Phantom of size (N, N, nt), time is last dimension.
    '''

    ph = np.zeros((N, N, nt))
    X, Y = np.meshgrid(
        np.linspace(-1, 1, N),
        np.linspace(-1, 1, N))

    # Make an outer circle
    thickness = .25
    idx0 = X**2 + Y**2 <= 1
    idx1 = X**2 + Y**2 >= (1 - thickness)**2
    idx = np.logical_and(idx0, idx1)
    ph[idx, :] = 1

    # Make another outer circle
    idx0 = X**2 + Y**2 <= (1 - thickness)**2
    idx1 = X**2 + Y**2 >= (1 - 2*thickness)**2
    idx = np.logical_and(idx0, idx1)
    ph[idx, :] = .2

    # Make inner circle
    r = np.cos(np.arange(nt)*2*np.pi/nt) + 2
    r = r/np.max(r)*.4
    thickness = .15
    for ii in range(nt):
        idx0 = X**2 + Y**2 <= r[ii]**2
        idx1 = X**2 + Y**2 >= (r[ii] - thickness)**2
        idx = np.logical_and(idx0, idx1)
        ph[idx, ii] = .8

    return ph


if __name__ == '__main__':
    pass
