'''Analytical form of Shepp-Logan kspace.'''

import numpy as np
from scipy.special import j1 # pylint: disable=E0611

from phantominator.ct_shepp_logan import (
    _shepp_logan_params_2d, _modified_shepp_logan_params_2d)

def kspace_shepp_logan(kx, ky, modified=True, E=None):
    '''2D Shepp-Logan phantom kspace measurements at points (kx, ky).

    Parameters
    ----------
    kx, ky : array_like
        1D arrays corresponding to kspace coordinates.  Coordinate
        units are the same as those returned by BART's traj function.
    modified : bool, optional
        Use original grey-scale values as given or use modified
        values for better contrast.
    E : array_like, optional
        See ct_shepp_logan for details.

    Returns
    -------
    k : array_like
        1D array of kspace measurements corresponding to coordinates
        (kx, ky).

    References
    ----------
    .. [1] Van de Walle, Rik, et al. "Reconstruction of MR images
           from data acquired on a general nonregular grid by
           pseudoinverse calculation." IEEE transactions on medical
           imaging 19.12 (2000): 1160-1167.
    '''

    assert kx.size == ky.size, 'kx and ky must be the same size!'

    # Get the ellipse parameters the user asked for
    if E is None:
        if modified:
            E = _modified_shepp_logan_params_2d()
        else:
            E = _shepp_logan_params_2d()

    # Extract params and get dims right to vectorize everything
    grey = E[:, 0][None, :]
    major = E[:, 1][None, :]
    minor = E[:, 2][None, :]
    xs = E[:, 3][None, :]
    ys = E[:, 4][None, :]
    alphas = E[:, 5][None, :]

    # Same for kspace coordinates.  Factor of 2 to match up with
    # BART's traj function
    kx = kx[:, None]/2
    ky = ky[:, None]/2

    # Sum of ellipses
    return np.sum(_kspace_ellipse(
        kx, ky, xs, ys, grey, major, minor, alphas), axis=-1)

def _kspace_ellipse(kx, ky, xc, yc, rho, A, B, alpha):
    '''Generates Fourier transform of a general ellipse.

    Notes
    -----
    Implements equation [21] in [1]_.
    '''

    k = kx + 1j*ky
    theta = np.angle(k)
    k = np.abs(k)
    t = np.sqrt(xc**2 + yc**2)
    gamma = np.arctan2(yc, xc)
    athetak = _a(A, B, theta, alpha)*k
    return np.exp(-1j*2*np.pi*k*t*np.cos(gamma - theta))*rho*A*B*j1(
        2*np.pi*athetak)/athetak

def _a(A, B, theta, alpha):
    return np.sqrt(
        A**2*np.cos(theta - alpha)**2 + B**2*np.sin(theta - alpha)**2)

if __name__ == '__main__':

    # Example usage
    from phantominator.traj import radial
    sx, spokes = 128, 128
    kx, ky = radial(sx, spokes)

    from bart import bart # pylint: disable=E0401
    traj = np.concatenate((
        kx.reshape((1, sx, spokes)),
        ky.reshape((1, sx, spokes)),
        np.zeros((1, sx, spokes))), axis=0)
    nufft = lambda x0: bart(
        1, 'nufft -i -t -d %d:%d:1' % (sx, sx),
        traj, x0.reshape((1, sx, spokes, 1))).squeeze()

    k = kspace_shepp_logan(kx, ky)

    import matplotlib.pyplot as plt
    # plt.scatter(kx, ky, 1)
    # plt.show()
    plt.imshow(np.abs(nufft(k)))
    plt.show()
