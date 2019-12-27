# -*- coding: utf-8 -*-
'''Shepp-Logan phantom for use with MR simulations.'''

import numpy as np

def mr_shepp_logan(N, E=None, B0=3, T2star=False, zlims=(-1, 1)):
    '''Shepp-Logan phantom with MR tissue parameters.

    Parameters
    ----------
    N : int or array_like
        Matrix size, (N, N, N), or (L, M, N).
    E : array_like, optional
        ex13 numeric matrix defining e ellipses.  The columns of E
        are:

            - x-coordinate of the center of the ellipsoid (in [-1, 1])
            - y-coordinate of the center of the ellipsoid (in [-1, 1])
            - z-coordinate of the center of the ellipsoid (in [-1, 1])
            - x principal axis of the ellipsoid
            - y principal axis of the ellipsoid
            - z principal axis of the ellipsoid
            - Angle of the ellipsoid (in rad)
            - spin density, M0
            - Parameter A for T1 determination
            - Parameter C for T1 determination
            - Explicit T1 value (in sec, or np.nan if model is used)
            - T2 value (in sec)
            - chi (change in magnetic susceptibility)

        If spin density is negative, M0, T1, and T2 will be subtracted
        instead of cummulatively added.
    B0 : float, optimal
        Field strength (in Tesla).
    T2star : bool, optional
        Use magnetic susceptibility values to return T2star values
        instead of T2. Gyromagnetic ratio is assumed to be that of
        hydrogen.
    zlims : tuple, optional
        Only for 3D.  Specify bounds along z.  Often we only want the
        middle portion of a 3D phantom, e.g., zlim=(-.5, .5).

    Returns
    -------
    M0 : array_like
        The proton density.
    T1 : array_like
        The T1 values.
    T2 : array_like
        The T2 values. If T2star is True, then these will be T2 star
        values.

    Notes
    -----
    Implements the phantoms described in [1]_.

    If parameters A, C are given and T1 is None, T1 is determined
    according to the equation:

        T1 = A*B0^C

    References
    ----------
    .. [1] Gach, H. Michael, Costin Tanase, and Fernando Boada.
           "2D & 3D Shepp-Logan phantom standards for MRI." 2008 19th
           International Conference on Systems Engineering. IEEE,
           2008.
    '''

    # Determine size of phantom
    if np.isscalar(N):
        L, M, N = N, N, N
    else:
        L, M, N = N[:]

    # Make sure zlims are appropriate
    assert len(zlims) == 2, (
        'zlims must be a tuple with 2 entries: upper and lower '
        'bounds!')
    assert zlims[0] <= zlims[1], (
        'zlims: lower bound must be first entry!')

    # Get parameters from paper if None provided
    if E is None:
        E = _ellipsoid_parameters()

    # Extract some parameters so we can use them
    xs = E[:, 0]
    ys = E[:, 1]
    zs = E[:, 2]
    xaxis = E[:, 3]
    yaxis = E[:, 4]
    zaxis = E[:, 5]
    theta = E[:, 6]
    M0 = E[:, 7]
    As = E[:, 8]
    Cs = E[:, 9]
    T1 = E[:, 10]
    T2 = E[:, 11]
    chis = E[:, 12]

    # Initialize array
    X, Y, Z = np.meshgrid( # meshgrid does X, Y backwards
        np.linspace(-1, 1, M),
        np.linspace(-1, 1, L),
        np.linspace(zlims[0], zlims[1], N))
    ct = np.cos(theta)
    st = np.sin(theta)
    sgn = np.sign(M0)
    T1s = np.zeros((L, M, N))
    T2s = np.zeros((L, M, N))
    M0s = np.zeros((L, M, N))

    # We'll need the gyromagnetic ratio if returning T2star values
    if T2star:
        # see https://en.wikipedia.org/wiki/Gyromagnetic_ratio:
        gamma0 = 267.52219 # 10^6 rad⋅s−1⋅T⋅−1

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

        # Add ellipses together -- subtract of M0 is negative
        M0s[idx] += M0[ii]

        # Use T2star values if user asked for them
        if T2star:
            T2s[idx] += sgn[ii]/(1/T2[ii] + gamma0*np.abs(
                B0*chis[ii]))
        else:
            T2s[idx] += sgn[ii]*T2[ii]

        # Use T1 model if not given explicit T1 value
        if np.isnan(T1[ii]):
            T1s[idx] += sgn[ii]*As[ii]*(B0**Cs[ii])
        else:
            T1s[idx] += sgn[ii]*T1[ii]

    return(M0s, T1s, T2s)


def _ellipsoid_parameters():
    '''Returns parameters of ellipsoids.

    Returns
    -------
    E : array_like
        Parameters for the ellipsoids used to construct the phantom.
    '''
    params = _mr_relaxation_parameters()

    E = np.zeros((15, 13))
    # [:, [x, y, z, a, b, c, theta, m0, A, C, (t1), t2, chi]]
    E[0, :] = [
        0, 0, 0, .72, .95, .93, 0, .8,
        *params['scalp']]
    E[1, :] = [
        0, 0, 0, .69, .92, .9, 0, .12,
        *params['marrow']]
    E[2, :] = [
        0, -.0184, 0, .6624, .874, .88, 0, .98,
        *params['csf']]
    E[3, :] = [
        0, -.0184, 0, .6524, .864, .87, 0, .745,
        *params['gray-matter']]
    E[4, :] = [
        -.22, 0, -.25, .41, .16, .21, np.deg2rad(-72), .98,
        *params['csf']]
    E[5, :] = [
        .22, 0, -.25, .31, .11, .22, np.deg2rad(72), .98,
        *params['csf']]
    E[6, :] = [
        0, .35, -.25, .21, .25, .35, 0, .617,
        *params['white-matter']]
    E[7, :] = [
        0, .1, -.25, .046, .046, .046, 0, .95,
        *params['tumor']]
    E[8, :] = [
        -.08, -.605, -.25, .046, .023, .02, 0, .95,
        *params['tumor']]
    E[9, :] = [
        .06, -.605, -.25, .046, .023, .02, np.deg2rad(-90), .95,
        *params['tumor']]
    E[10, :] = [
        0, -.1, -.25, .046, .046, .046, 0, .95,
        *params['tumor']]
    E[11, :] = [
        0, -.605, -.25, .023, .023, .023, 0, .95,
        *params['tumor']]
    E[12, :] = [
        .06, -.105, .0625, .056, .04, .1, np.deg2rad(-90), .93,
        *params['tumor']]
    E[13, :] = [
        0, .1, .625, .056, .056, .1, 0, .98,
        *params['csf']]
    E[14, :] = [
        .56, -.4, -.25, .2, .03, .1, np.deg2rad(70), .85,
        *params['blood-clot']]

    # Need to subtract some ellipses here...
    Eneg = np.zeros(E.shape)
    for ii in range(E.shape[0]):

        # Ellipsoid geometry
        Eneg[ii, :7] = E[ii, :7]

        # Tissue property differs after 4th subtracted ellipsoid
        if ii > 3:
            Eneg[ii, 7:] = E[3, 7:]
        else:
            Eneg[ii, 7:] = E[ii-1, 7:]

    # Throw out first as we skip this one in the paper's table
    Eneg = Eneg[1:, :]

    # Spin density is negative for subtraction
    Eneg[:, 7] *= -1

    # Paper doesn't use last blood-clot ellipsoid
    E = E[:-1, :]
    Eneg = Eneg[:-1, :]

    # Put both ellipsoid groups together
    E = np.concatenate((E, Eneg), axis=0)

    return E

def _mr_relaxation_parameters():
    '''Returns MR relaxation parameters for certain tissues.

    Returns
    -------
    params : dict
        Gives entries as [A, C, (t1), t2, chi]

    Notes
    -----
    If t1 is None, the model T1 = A*B0^C will be used.  If t1 is not
    np.nan, then specified t1 will be used.
    '''

    # params['tissue-name'] = [A, C, (t1 value if explicit), t2, chi]
    params = dict()
    params['scalp'] = [.324, .137, np.nan, .07, -7.5e-6]
    params['marrow'] = [.533, .088, np.nan, .05, -8.85e-6]
    params['csf'] = [np.nan, np.nan, 4.2, 1.99, -9e-6]
    params['blood-clot'] = [1.35, .34, np.nan, .2, -9e-6]
    params['gray-matter'] = [.857, .376, np.nan, .1, -9e-6]
    params['white-matter'] = [.583, .382, np.nan, .08, -9e-6]
    params['tumor'] = [.926, .217, np.nan, .1, -9e-6]
    return params
