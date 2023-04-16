"""Basic usage of dynamic phantom."""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from phantominator import dynamic


if __name__ == '__main__':
    # Create a phantom with nt time points
    N = 128
    nt = 30
    ph = dynamic(N, nt)

    # Some code to look at the animation
    fig = plt.figure()
    ax = plt.imshow(ph[..., 0], cmap='gray')
    plt.title("2+1D simulation")

    def init():
        """Initialize ax data."""
        ax.set_array(ph[..., 0])
        return(ax,)

    def animate(frame):
        """Update frame."""
        ax.set_array(ph[..., frame])
        return(ax,)

    anim = FuncAnimation(
        fig, animate, init_func=init, frames=nt, interval=50,
        blit=True)
    plt.show()
