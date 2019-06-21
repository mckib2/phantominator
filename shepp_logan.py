'''The canonical Shepp-Logan phantom.'''

import numpy as np

def shepp_logan(N, orig_vals=False):
    '''Generate a Shepp-Logan phantom of size (N, N).

    Parameters
    ----------
    N : int
        Matrix size, (N, N).
    orig_vals : bool, optional
        Use original grey-scale values as given in [1]_.  Most
        implementations use modified values for better contrast (for
        example, see [3]_ and [4]_)

    Returns
    -------
    ph : array_like
        The Shepp-Logan phantom.

    Notes
    -----
    This much abused phantom is due to [1]_.  The tabulated values in
    the paper are reproduced in the Wikipedia entry [2]_.  The
    original values do not produce great contrast, so modified values
    are used by default.

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
    '''

    # Centers of ellipses
    ctrs = [
        (0, 0),
        (0, -0.0184),
        (0.22, 0),
        (-0.22, 0),
        (0, 0.35),
        (0, 0.1),
        (0, -0.1),
        (-0.08, -0.605),
        (0, -0.605),
        (0.06, -0.605)]

    # Major-axes
    major = [
        0.69, 0.6624, 0.11, 0.16, 0.21, 0.046, 0.046, 0.046, 0.023,
        0.023]

    # Minor-axes
    minor = [
        0.92, 0.874, 0.31, 0.41, 0.25, 0.046, 0.046, 0.023, 0.023,
        0.023]

    # Angle of rotation (in deg)
    theta = np.deg2rad([0, 0, -18, 18, 0, 0, 0, 0, 0, 0])

    # Grey level
    if orig_vals:
        grey = [
            2, -.98, -.02, -.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    else:
        # Better contrast:
        grey = [1, -0.8, -0.2, -0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    # 2x2 square => FOV = (-1, 1)
    xx = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(xx, xx)
    ph = np.zeros((N, N))
    ct = np.cos(theta)
    st = np.sin(theta)

    for ii, ctr in enumerate(ctrs):
        xc, yc = ctr
        a, b = major[ii], minor[ii]
        ct0, st0 = ct[ii], st[ii]

        # Find indices falling inside the ellipse
        idx = (
            ((X - xc)*ct0 + (Y - yc)*st0)**2/a**2 +
            ((X - xc)*st0 - (Y - yc)*ct0)**2/b**2 <= 1)

        # Sum of ellipses
        ph[idx] += grey[ii]

    return ph

if __name__ == '__main__':
    pass
