"""Replicate paper figures from [1]_.

References
----------
.. [1] Gach, H. Michael, Costin Tanase, and Fernando Boada.
       "2D & 3D Shepp-Logan phantom standards for MRI." 2008 19th
       International Conference on Systems Engineering. IEEE,
       2008.
"""

from time import time

import numpy as np
import matplotlib.pyplot as plt
from ssfp import spoiled_gre

from phantominator import shepp_logan


def se90(T1, T2, TR, TE, M0=1):
    """Spin echo simulation assuming 90 deg flip angle

    Parameters
    ----------
    T1 : array_like
        longitudinal relaxation time.
    T2 : array_like
        transverse relaxation time.
    TR : float
        repetition time.
    TE : float
        echo time
    M0 : array_like, optional
        Proton density.

    Returns
    -------
    S : array_like
        Simulated magnitude spin echo image.
    """

    # Make sure we don't divide by zero
    idx1 = np.nonzero(T1)
    idx2 = np.nonzero(T2)
    E1 = np.zeros(T1.shape)
    E1[idx1] = np.exp(-TR/T1[idx1])
    E2 = np.zeros(T2.shape)
    E2[idx2] = np.exp(-TE/T2[idx2])

    S = M0*(1 - E1)*E2
    return S


if __name__ == '__main__':

    # Generate a single slice of 3D Shepp-Logan phantom, z=-.25
    # corresponds to the familiar 2D CT phantom
    t0 = time()
    N = 256
    M0, T1, T2 = shepp_logan((N, N, 1), MR=True, zlims=(-.25, .25))
    print(f"Took {time() - t0} seconds")

    # This is a 3D simulation, but only 1 slice, so remove singleton
    # dimension at the end for simplicity:
    M0, T1, T2 = np.squeeze(M0), np.squeeze(T1), np.squeeze(T2)

    # Rotate to match orientation shown in paper
    M0, T1, T2 = (
        np.rot90(M0, k=2), np.rot90(T1, k=2), np.rot90(T2, k=2))

    # T1 weighted GRE sequence simulation
    TR = 0.035
    TE = 0.01
    FA = np.deg2rad(40)
    T1w = spoiled_gre(T1, T2, TR, TE, alpha=FA, M0=M0)

    # T2 weighted SE sequence simulation
    TR = 4
    TE = 0.1
    T2w = se90(T1, T2, TR, TE, M0)

    # Take a look
    nx, ny = 2, 2
    plt_opts = {
        'cmap': 'gray',
    }
    plt.subplot(nx, ny, 1)
    plt.imshow(M0, **plt_opts)
    plt.axis('off')
    plt.title('Proton Density')

    plt.subplot(nx, ny, 2)
    plt.imshow(T1w, **plt_opts)
    plt.axis('off')
    plt.title('T1 weighted')

    plt.subplot(nx, ny, 3)
    plt.imshow(T2w, **plt_opts)
    plt.axis('off')
    plt.title('T2 weighted')
    plt.show()
