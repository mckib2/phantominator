'''The canonical Shepp-Logan phantom.'''

import numpy as np

def shepp_logan(N, modified=True, E=None, ret_E=False):
    '''Generate a Shepp-Logan phantom of size (N, N).

    Parameters
    ----------
    N : int
        Matrix size, (N, N).
    modified : bool, optional
        Use original grey-scale values as given in [1]_.  Most
        implementations use modified values for better contrast (for
        example, see [3]_ and [4]_).
    E : array_like
        ex6 numeric matrix defining e ellipses.  The six columns of
        E are:

            - Additive intensity value of the ellipse
            - Length of the horizontal semiaxis of the ellipse
            - Length of the vertical semiaxis of the ellipse
            - x-coordinate of the center of the ellipse (in [-1, 1])
            - y-coordinate of the center of the ellipse (in [-1, 1])
            - Angle between the horizontal semiaxis of the ellipse
              and the x-axis of the image (in rad)

    ret_E : bool, optional
        Return the matrix E used to generate the phantom, ph.

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
    '''

    # Get the ellipse parameters the user asked for
    if E is None:
        if modified:
            E = modified_shepp_logan_params()
        else:
            E = shepp_logan_params()

    # Extract params
    grey = E[:, 0]
    major = E[:, 1]
    minor = E[:, 2]
    xs = E[:, 3]
    ys = E[:, 4]
    theta = E[:, 5]

    # 2x2 square => FOV = (-1, 1)
    xx = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(xx, xx)
    ph = np.zeros((N, N))
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

def shepp_logan_params():
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
    E[:, 2] = [.92, .874, .31, .41, .25, .046, .046, .023, .023, .023]
    E[:, 3] = [0, 0, .22, -.22, 0, 0, 0, -.08, 0, .06]
    E[:, 4] = [0, -.0184, 0, 0, .35, .1, -.1, -.605, -.605, -.605]
    E[:, 5] = np.deg2rad([0, 0, -18, 18, 0, 0, 0, 0, 0, 0])
    return E

def modified_shepp_logan_params():
    '''Return parameters for  modified Shepp-Logan phantom.

    Returns
    -------
    E : array_like, shape (10, 6)
        Parameters for the 10 ellipses used to construct the phantom.
    '''
    E = shepp_logan_params()
    E[:, 0] = [1, -0.8, -0.2, -0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    return E

if __name__ == '__main__':
    pass
