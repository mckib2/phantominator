"""Show how to generate radial kspace with coil sensitivities."""

import numpy as np
from scipy.cluster.vq import whiten # pylint: disable=C0412
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pygrappa import radialgrappaop, grog  # pylint: disable=E0611

from phantominator import kspace_shepp_logan
from phantominator.traj import radial


if __name__ == '__main__':

    # Example usage (requires pygrappa package to be installed!)
    sx, spokes, ncoil = 288, 72, 8
    kx, ky = radial(sx, spokes)
    kx = np.reshape(kx, (sx, spokes), 'F').flatten()
    ky = np.reshape(ky, (sx, spokes), 'F').flatten()
    k = kspace_shepp_logan(kx, ky, ncoil=ncoil)
    k = whiten(k)

    # Grid via GROG and check out the results:
    Gx, Gy = radialgrappaop(
        np.reshape(kx, (sx, spokes)),
        np.reshape(ky, (sx, spokes)),
        np.reshape(k, (sx, spokes, ncoil)))
    coil_ims = np.abs(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(
        grog(kx, ky, k, sx, sx, Gx, Gy),
        axes=(0, 1)), axes=(0, 1)), axes=(0, 1)))

    # Some code to look at the animation
    fig = plt.figure()
    ax = plt.imshow(coil_ims[..., 0], cmap='gray')

    def init():
        """Initialize ax data."""
        ax.set_array(coil_ims[..., 0])
        return(ax,)

    def animate(frame):
        """Update frame."""
        ax.set_array(coil_ims[..., frame])
        return(ax,)

    anim = FuncAnimation(
        fig, animate, init_func=init, frames=ncoil, interval=50,
        blit=True)
    plt.show()
