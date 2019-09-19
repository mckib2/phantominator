'''Analytical form of Shepp-Logan kspace.'''

from time import time

import numpy as np
from scipy.special import j1 # pylint: disable=E0611

from phantominator.ct_shepp_logan import (
    _shepp_logan_params_2d, _modified_shepp_logan_params_2d)
from phantominator.sens_coeffs import _sens_coeffs, _sens_info

def kspace_shepp_logan(kx, ky, modified=True, E=None, ncoil=None):
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
    grey = E[:, 0]
    major = E[:, 1]
    minor = E[:, 2]
    xs = E[:, 3]
    ys = E[:, 4]
    alphas = E[:, 5]

    # We want sensitivity maps!
    if ncoil is not None:
        MAX_COIL, NUM_COEFF = _sens_info()
        assert ncoil <= MAX_COIL, (
            'Only %d coils possible to simulate!' % MAX_COIL)
        t0 = time()

        # Build up the coefficient matrix, we'll do all coils for
        # each ellipse for coefficiency
        coeffs = np.zeros((ncoil, NUM_COEFF), dtype='complex')
        for cc in range(ncoil):
            coeffs[cc, :] = _sens_coeffs(cc)

        # Add up all the ellipse kspaces with all coils
        val = np.zeros((kx.size, ncoil), dtype='complex')
        for ii in range(E.shape[0]):
            # Have to get screwy with the center coordinates to make
            # it work.  Not sure what's different between us and
            # MATLAB script...
            val += _kspace_ellipse_sens(
                kx, ky, ys[ii]/2, -xs[ii]/2, grey[ii], major[ii],
                minor[ii], alphas[ii], coeffs).T

        print('Took %g seconds to simulate %d coils' % (
            time() - t0, ncoil))
        return val

    # Same for kspace coordinates.  Factor of 2 to match up with
    # BART's traj function
    kx = kx[:, None]/2
    ky = ky[:, None]/2

    # Sum of ellipses
    return np.sum(_kspace_ellipse(
        kx, ky, xs, ys, grey[None, :], major[None, :],
        minor[None, :], alphas[None, :]), axis=-1)

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

def _kspace_ellipse_sens(kx, ky, xc, yc, rho, A, B, theta, coeffs):
    '''Fourier transform of ellipse with polynomial sensitivity map.
    '''

    width = np.array([B, A])
    ct, st = np.cos(theta), np.sin(theta)
    Rmat = np.array([[ct, -st], [st, ct]])
    Dmat = np.diag(width/2)

    return rho*MRDataEllipseSinusoidal(
        kx, ky, Dmat, Rmat, xc, yc, coeffs)

def MRDataEllipseSinusoidal(kx, ky, Dmat, Rmat, xc, yc, coeffs):
    '''Sinusoidal model'''
    N = coeffs.shape[1]
    L = np.floor(np.sqrt(N))

    k = np.concatenate((kx[None, :], ky[None, :]), axis=0)
    Nk = k.shape[1]

    ky, kx = np.meshgrid(
        np.linspace(-np.floor(L/2), np.floor((L-1)/2), L),
        np.linspace(-np.floor(L/2), np.floor((L-1)/2), L))
    kx = kx.flatten('C')/2
    ky = ky.flatten('C')/2
    kx = np.tile(kx[:, None], (1, Nk)) + np.tile(k[0, :], (N, 1))
    ky = np.tile(ky[:, None], (1, Nk)) + np.tile(k[1, :], (N, 1))
    kxy = np.concatenate(
        (kx.flatten('C')[None, :], ky.flatten('C')[None, :]), axis=0)

    wu = -2*np.pi*Dmat @ Rmat.conj().T @ kxy
    # wux = np.reshape(wu[0, :], (N, Nk), 'C')
    # wuy = np.reshape(wu[1, :], (N, Nk), 'C')
    # modw = np.sqrt(wux**2 + wuy**2)
    modw = np.abs(wu[0, :] + 1j*wu[1, :]).reshape((N, Nk))

    # The robust Bessel function
    Gval = np.zeros(modw.shape)
    ind_big = np.abs(modw) >= 1e-10
    Gval[ind_big] = j1(modw[ind_big])/modw[ind_big]
    Gval[~ind_big] = .5*(1 - (modw[~ind_big]/2)**2/2)

    return 2*np.pi*np.linalg.det(Dmat)*coeffs @ (
        np.exp(2*np.pi*1j*(xc*kx + yc*ky))*Gval)

if __name__ == '__main__':
    pass
    # # Example usage
    # from phantominator.traj import radial
    # sx, spokes, ncoil = 128, 128, 8
    # kx, ky = radial(sx, spokes)
    # k = kspace_shepp_logan(kx, ky, ncoil=ncoil)
    #
    # import matplotlib.pyplot as plt
    # sos = lambda x0: np.sqrt(np.sum(np.abs(x0)**2, axis=-1))
    # coil_ims = gridder(kx, ky, k, sx, sx)
    # plt.imshow(sos(coil_ims))
    # # for ii in range(ncoil):
    # #     plt.subplot(1, ncoil, ii+1)
    # #     plt.imshow(sos(coil_ims[..., ii][..., None]))
    # plt.show()
