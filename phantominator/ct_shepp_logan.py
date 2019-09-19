# -*- coding: utf-8 -*-
'''The canonical Shepp-Logan phantom used for CT simulations.'''

import numpy as np

def ct_shepp_logan(
        N, modified=True, E=None, ret_E=False, zlims=(-1, 1)):
    '''Generate a Shepp-Logan phantom of size (N, N).

    Parameters
    ----------
    N : int or array_like
        Matrix size, (N, N) or (M, N) or (L, M, N).
    modified : bool, optional
        Use original grey-scale values as given in [1]_.  Most
        implementations use modified values for better contrast (for
        example, see [3]_ and [4]_).
    E : array_like, optional
        For 2D: ex6 numeric matrix defining e ellipses.  The six
        columns of E are:

            - Length of the horizontal semiaxis of the ellipse
            - Length of the vertical semiaxis of the ellipse
            - x-coordinate of the center of the ellipse (in [-1, 1])
            - y-coordinate of the center of the ellipse (in [-1, 1])
            - Angle between the horizontal semiaxis of the ellipse
              and the x-axis of the image (in rad)

        For 3D: ex8 numeric matrix defining e ellipsoids.  The
        columns are:

            - Additive intensity value of the ellipse
            - x principal axis of the ellipsoid
            - y principal axis of the ellipsoid
            - z principal axis of the ellipsoid
            - x-coordinate of the center of the ellipsoid (in [-1, 1])
            - y-coordinate of the center of the ellipsoid (in [-1, 1])
            - z-coordinate of the center of the ellipsoid (in [-1, 1])
            - Angle of the ellipsoid (in rad)

    ret_E : bool, optional
        Return the matrix E used to generate the phantom, ph.
    zlims : tuple, optional
        Only for 3D.  Specify bounds along z.  Often we only want the
        middle portion of a 3D phantom, e.g., zlim=(-.5, .5).

    Returns
    -------
    ph : array_like
        The Shepp-Logan phantom.
    E : array_like, optional
        The ellipse parameters used to generate ph.

    Notes
    -----
    This much abused phantom is due to [1]_.  The tabulated values in
    the paper are reproduced in the Wikipedia entry [2]_.  The
    original values do not produce great contrast, so modified values
    are used by default (see Table B.1 in [5]_ or implementations
    [3]_ and [4]_).

    The 3D version follows values tabulated in [6]_.  The modified
    version uses the same adjustments to grayscale.

    References
    ----------
    .. [1] Shepp, Lawrence A., and Benjamin F. Logan. "The Fourier
           reconstruction of a head section." IEEE Transactions on
           nuclear science 21.3 (1974): 21-43.
    .. [2] https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom
    .. [3] https://sigpy.readthedocs.io/en/latest/_modules/sigpy/
           sim.html#shepp_logan
    .. [4] http://www.mathworks.com/matlabcentral/fileexchange/
           9416-3d-shepp-logan-phantom
    .. [5] Toft, Peter Aundal, and John Aasted Sørensen. "The Radon
           transform-theory and implementation." (1996).
    .. [6] Koay, Cheng Guan, Joelle E. Sarlls, and Evren Özarslan.
           "Three‐dimensional analytical magnetic resonance imaging
           phantom in the Fourier domain." Magnetic Resonance in
           Medicine: An Official Journal of the International Society
           for Magnetic Resonance in Medicine 58.2 (2007): 430-436.
    '''

    # Get size of phantom
    if np.isscalar(N):
        M, N = N, N
        is2D = True
    else:
        if len(N) == 2:
            M, N = N[:]
            is2D = True
        elif len(N) == 3:
            L, M, N = N[:]
            is2D = False
        else:
            raise ValueError('Dimension must be scalar, 2D, or 3D!')

    # Give back either a 2D or 3D phantom
    if is2D:
        return ct_shepp_logan_2d(M, N, modified, E, ret_E)
    return ct_shepp_logan_3d(L, M, N, modified, E, ret_E, zlims)

def ct_shepp_logan_2d(M, N, modified, E, ret_E):
    '''Make a 2D phantom.'''

    # Get the ellipse parameters the user asked for
    if E is None:
        if modified:
            E = _modified_shepp_logan_params_2d()
        else:
            E = _shepp_logan_params_2d()

    # Extract params
    grey = E[:, 0]
    major = E[:, 1]
    minor = E[:, 2]
    xs = E[:, 3]
    ys = E[:, 4]
    theta = E[:, 5]

    # 2x2 square => FOV = (-1, 1)
    X, Y = np.meshgrid( # meshgrid needs linspace in opposite order
        np.linspace(-1, 1, N),
        np.linspace(-1, 1, M))
    ph = np.zeros((M, N))
    ct = np.cos(theta)
    st = np.sin(theta)

    for ii in range(E.shape[0]):
        xc, yc = xs[ii], ys[ii]
        a, b = major[ii], minor[ii]
        ct0, st0 = ct[ii], st[ii]

        # Find indices falling inside the ellipse
        idx = (
            ((X - xc)*ct0 + (Y - yc)*st0)**2/a**2 +
            ((X - xc)*st0 - (Y - yc)*ct0)**2/b**2 <= 1)

        # Sum of ellipses
        ph[idx] += grey[ii]

    if ret_E:
        return(ph, E)
    return ph

def ct_shepp_logan_3d(L, M, N, modified, E, ret_E, zlims):
    '''Make a 3D phantom.'''

    # Make sure zlims are appropriate
    assert len(zlims) == 2, (
        'zlims must be a tuple with 2 entries: upper and lower '
        'bounds!')
    assert zlims[0] <= zlims[1], (
        'zlims: lower bound must be first entry!')

    # Get parameters from paper if None provided
    if E is None:
        if modified:
            E = _modified_shepp_logan_params_3d()
        else:
            E = _shepp_logan_params_3d()

    # Extract some parameters so we can use them
    gray = E[:, 0]
    xaxis = E[:, 1]
    yaxis = E[:, 2]
    zaxis = E[:, 3]
    xs = E[:, 4]
    ys = E[:, 5]
    zs = E[:, 6]
    theta = E[:, 7]

    # Initialize array
    X, Y, Z = np.meshgrid( # meshgrid does X, Y backwards
        np.linspace(-1, 1, M),
        np.linspace(-1, 1, L),
        np.linspace(zlims[0], zlims[1], N))
    ct = np.cos(theta)
    st = np.sin(theta)
    ph = np.zeros((L, M, N))

    # Put ellipses where they need to be
    for ii in range(E.shape[0]):
        xc, yc, zc = xs[ii], ys[ii], zs[ii]
        a, b, c = xaxis[ii], yaxis[ii], zaxis[ii]
        ct0, st0 = ct[ii], st[ii]

        # Find indices falling inside the ellipsoid, ellipses only
        # rotated in xy plane
        idx = (
            ((X - xc)*ct0 + (Y - yc)*st0)**2/a**2 +
            ((X - xc)*st0 - (Y - yc)*ct0)**2/b**2 +
            (Z - zc)**2/c**2 <= 1)

        # Add ellipses together
        ph[idx] += gray[ii]

    if ret_E:
        return(ph, E)
    return ph

def _shepp_logan_params_2d():
    '''Return parameters for original Shepp-Logan phantom.

    Returns
    -------
    E : array_like, shape (10, 6)
        Parameters for the 10 ellipses used to construct the phantom.
    '''

    E = np.zeros((10, 6)) # (10, [A, a, b, xc, yc, theta])
    E[:, 0] = [2, -.98, -.02, -.02, .01, .01, .01, .01, .01, .01]
    E[:, 1] = [
        .69, .6624, .11, .16, .21, .046, .046, .046, .023, .023]
    E[:, 2] = [.92, .874, .31, .41, .25, .046, .046, .023, .023, .046]
    E[:, 3] = [0, 0, .22, -.22, 0, 0, 0, -.08, 0, .06]
    E[:, 4] = [0, -.0184, 0, 0, .35, .1, -.1, -.605, -.605, -.605]
    E[:, 5] = np.deg2rad([0, 0, -18, 18, 0, 0, 0, 0, 0, 0])
    return E

def _modified_shepp_logan_params_2d():
    '''Return parameters for modified Shepp-Logan phantom.

    Returns
    -------
    E : array_like, shape (10, 6)
        Parameters for the 10 ellipses used to construct the phantom.
    '''
    E = _shepp_logan_params_2d()
    E[:, 0] = [1, -0.8, -0.2, -0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    return E

def _shepp_logan_params_3d():
    '''Ellipsoid parameters for 3D Shepp-Logan.

    Returns
    -------
    E : array_like, shape (10, 8)
        Parameters for the 10 ellipsoids used to construct phantom.
    '''

    E = np.zeros((10, 8))
    # [ gray, x, y, z, a, b, c, theta ]
    E[0, :] = [2, 0.69, 0.92, 0.9, 0, 0, 0, 0]
    E[1, :] = [-.8, 0.6624, 0.874, 0.88, 0, 0, 0, 0]
    E[2, :] = [-.2, 0.41, 0.16, 0.21, -0.22, 0, -0.25, 3*np.pi/5]
    E[3, :] = [-.2, 0.31, 0.11, 0.22, 0.22, 0, -0.25, 2*np.pi/5]
    E[4, :] = [.2, 0.21, 0.25, 0.5, 0, 0.35, -0.25, 0]
    E[5, :] = [.2, 0.046, 0.046, 0.046, 0, 0.1, -0.25, 0]
    E[6, :] = [.1, 0.046, 0.023, 0.02, -0.08, -0.65, -0.25, 0]
    E[7, :] = [.1, 0.046, 0.023, 0.02, 0.06, -0.65, -0.25, np.pi/2]
    E[8, :] = [.2, 0.056, 0.04, 0.1, 0.06, -0.105, 0.625, np.pi/2]
    E[9, :] = [-.2, 0.056, 0.056, 0.1, 0, 0.1, 0.625, 0]
    return E

def _modified_shepp_logan_params_3d():
    '''Return parameters for modified Shepp-Logan phantom.

    Returns
    -------
    E : array_like, shape (10, 8)
        Parameters for the 10 ellipsoids used to construct phantom.
    '''
    E = _shepp_logan_params_3d()
    E[:, 0] = [1, -0.8, -0.2, -0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    return E


if __name__ == '__main__':
    pass
