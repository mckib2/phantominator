"""Example demonstrating how to make a Shepp-Logan phantom."""

import numpy as np
import matplotlib.pyplot as plt

from phantominator import shepp_logan


if __name__ == '__main__':

    # The original Shepp-Logan
    ph = shepp_logan(128, modified=False)
    plt.subplot(1, 2, 1)
    plt.title('Shepp-Logan')
    plt.imshow(ph, cmap='gray')

    # Modified Shepp-Logan for better contrast
    ph = shepp_logan(128, modified=True)
    plt.subplot(1, 2, 2)
    plt.title('Modified Shepp-Logan')
    plt.imshow(ph, cmap='gray')
    plt.show()

    # Generate phantoms with different sizes
    ph = shepp_logan((128, 256))
    plt.title('Strange sizes')
    plt.imshow(ph, cmap='gray')
    plt.show()

    # Get a 3D phantom
    ph = shepp_logan((128, 128, 20), zlims=(-.5, .5))

    # Fancy dancing to nicely show all slices on same plot
    nx = int(np.ceil(np.sqrt(ph.shape[-1])))
    fig = plt.figure(figsize=(nx, nx))
    for ii in range(ph.shape[-1]):
        ax = fig.add_subplot(nx, nx, ii+1)
        plt.imshow(ph[..., ii], cmap='gray')
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_aspect('equal')
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()
