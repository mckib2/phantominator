"""Demonstrate how to use MR Shepp-Logan."""

import matplotlib.pyplot as plt

from phantominator import shepp_logan


if __name__ == '__main__':

    # Get proton density, T1, and T2 maps for (L, M, N) sized Shepp-
    # Logan phantom.
    L, M, N = 256, 255, 20
    M0, T1, T2 = shepp_logan((L, M, N), MR=True, zlims=(-.25, .25))

    # Take a look
    nx, ny = 2, 2
    plt.subplot(nx, ny, 1)
    plt.imshow(M0[..., 0])
    plt.title('Proton Density')
    plt.axis('off')

    plt.subplot(nx, ny, 2)
    plt.imshow(T1[..., 0])
    plt.title('T1')
    plt.axis('off')

    plt.subplot(nx, ny, 3)
    plt.imshow(T2[..., 0])
    plt.title('T2')
    plt.axis('off')

    plt.show()
