'''Analytical form of Shepp-Logan kspace.'''

from time import time
from math import tau

import numpy as np
from scipy.special import j1 # pylint: disable=E0611

from phantominator.ct_shepp_logan import (
    _shepp_logan_params_2d, _modified_shepp_logan_params_2d)
from phantominator.sens_coeffs import _sens_coeffs, _sens_info


def kspace_shepp_logan(
        kx, ky=None, modified=True, E=None, ncoil=None):
    '''2D Shepp-Logan phantom kspace measurements at points (kx, ky).

    Parameters
    ----------
    kx, ky : array_like
        1D arrays corresponding to kspace coordinates.  Coordinate
        units are the same as those returned by BART's traj function.
        If ky=None, then kx should be a complex-valued array
        corresponding to kx + 1j*ky.
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

    if ky is None:
        k = np.asarray(kx)
    else:
        kx, ky = np.asarray(kx), np.asarray(ky)
        assert kx.shape == ky.shape, (
            'kx and ky must be the same size!')
        k = kx + 1j*ky

    # Get the ellipse parameters the user asked for
    if E is None:
        if modified:
            E = _modified_shepp_logan_params_2d()
        else:
            E = _shepp_logan_params_2d()

    # We want sensitivity maps!
    if ncoil is not None:
        # Extract params and get dims right to vectorize everything
        grey = E[:, 0]
        major = E[:, 1]
        minor = E[:, 2]
        xs = E[:, 3]
        ys = E[:, 4]
        alphas = E[:, 5]

        MAX_COIL, NUM_COEFF = _sens_info()
        assert ncoil <= MAX_COIL, (
            f'Only {MAX_COIL} coils possible to simulate!')
        t0 = time()

        # Build up the coefficient matrix, we'll do all coils for
        # each ellipse at the same time for efficiency
        coeffs = np.zeros((ncoil, NUM_COEFF), dtype=np.complex)
        for cc in range(ncoil):
            coeffs[cc, :] = _sens_coeffs(cc)

        # Add up all the ellipse kspaces with all coils
        val = np.zeros((k.size, ncoil), dtype=np.complex)
        for ii in range(E.shape[0]):
            # Have to get screwy with the center coordinates to make
            # it work.  Not sure what's different between us and
            # MATLAB script...
            val += _kspace_ellipse_sens(
                k, ys[ii]/2, -xs[ii]/2, grey[ii], major[ii],
                minor[ii], alphas[ii], coeffs).T

        print('Took %g seconds to simulate %d coils' % (
            time() - t0, ncoil))
        return val

    # Same for kspace coordinates.  Factor of 2 to match up with
    # BART's traj function
    # FIXME: Scaling removed so that this matches Matlab mriphantom.
    #        Which was is correct?
    #k /= 2

    # Sum of ellipses
    return np.sum(_kspace_ellipse(k, E), axis=-1)

def _kspace_ellipse(k, E):
    '''Generates Fourier transform of (an array of) a general ellipse.

    Notes
    -----
    Implements equation [21] in [1]_.
    '''
    k = np.asarray(k)
    rho, A, B, xc, yc, alpha = np.asarray(E).T
    E = rho*A*B
    Ec = xc + 1j*yc
    ret = np.empty(shape=np.shape(k) + np.shape(E), dtype=np.complex)
    zero = np.isclose(k, 0)
    ret[zero] = .5*tau*E  # lim a->0: j1(tau * a) / a = .5 * tau
    k = k[~zero, None]
    if k.size:
        k, theta = abs(k), np.angle(k)
        t, gamma = abs(Ec), np.angle(Ec)
        athetak = _a(A, B, theta, alpha)*k
        rotation = np.exp(-1j*tau*k*t*np.cos(gamma - theta))
        ret[~zero] = rotation*E*j1(tau*athetak)/athetak
    return ret

def _a(A, B, theta, alpha):
    return np.sqrt(
        A**2*np.cos(theta - alpha)**2 + B**2*np.sin(theta - alpha)**2)

def _kspace_ellipse_sens(k, xc, yc, rho, A, B, theta, coeffs):
    '''Fourier transform of ellipse with polynomial sensitivity map.
    '''

    width = np.array([B, A])
    ct, st = np.cos(theta), np.sin(theta)
    Rmat = np.array([[ct, -st], [st, ct]])
    Dmat = np.diag(width/2)

    return rho*MRDataEllipseSinusoidal(k, Dmat, Rmat, xc, yc, coeffs)

def MRDataEllipseSinusoidal(k, Dmat, Rmat, xc, yc, coeffs):
    '''Sinusoidal model'''
    N = coeffs.shape[1]
    L = np.floor(np.sqrt(N))

    k = np.concatenate((k.real[None, :], k.imag[None, :]), axis=0)
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

    wu = -tau * Dmat @ Rmat.conj().T @ kxy
    modw = np.abs(wu[0, :] + 1j*wu[1, :]).reshape((N, Nk))

    # The robust Bessel function
    Gval = np.zeros(modw.shape)
    ind_big = np.abs(modw) >= 1e-10
    Gval[ind_big] = j1(modw[ind_big])/modw[ind_big]
    Gval[~ind_big] = .5*(1 - (modw[~ind_big]/2)**2/2)

    return tau*np.linalg.det(Dmat)*coeffs @ (
        np.exp(1j*tau*(xc*kx + yc*ky))*Gval)

if __name__ == '__main__':
    pass
    # # Example usage (requires pygrappa package to be installed!)
    # from phantominator.traj import radial
    # from scipy.cluster.vq import whiten # pylint: disable=C0412
    # import matplotlib.pyplot as plt
    # from pygrappa import radialgrappaop, grog # pylint: disable=E0611
    #
    # sx, spokes, ncoil = 288, 72, 8
    # kx, ky = radial(sx, spokes)
    # kx = np.reshape(kx, (sx, spokes), 'F').flatten()
    # ky = np.reshape(ky, (sx, spokes), 'F').flatten()
    # k = kspace_shepp_logan(kx, ky, ncoil=ncoil)
    # k = whiten(k)
    #
    # # Grid and check out the results:
    # Gx, Gy = radialgrappaop(
    #     np.reshape(kx, (sx, spokes)),
    #     np.reshape(ky, (sx, spokes)),
    #     np.reshape(k, (sx, spokes, ncoil)))
    # coil_ims = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(
    #     grog(kx, ky, k, sx, sx, Gx, Gy),
    #     axes=(0, 1)), axes=(0, 1)), axes=(0, 1))
    # sos = lambda x0: np.sqrt(np.sum(np.abs(x0)**2, axis=-1))
    # # plt.imshow(sos(coil_ims))
    # for ii in range(ncoil):
    #     plt.subplot(1, ncoil, ii+1)
    #     plt.imshow(sos(coil_ims[..., ii][..., None]))
    # plt.show()
