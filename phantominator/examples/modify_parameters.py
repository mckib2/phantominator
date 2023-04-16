"""Example demonstrating how to modify parameters."""

import numpy as np
import matplotlib.pyplot as plt

from phantominator import shepp_logan, mr_ellipsoid_parameters


if __name__ == '__main__':

    # Get stock ellipsoid parameters
    E = mr_ellipsoid_parameters()

    # 8th column is spin density, there are 13 ellipsoids.  Let's
    # change spin density of the first 5 ellipsoids
    E[:5, 7] = np.linspace(.1, .9, 5)

    # Make the Shepp-Logan phantom with modified spin densities
    N = 256
    M0, T1, T2 = shepp_logan(
        (N, N, 1), MR=True, E=E, zlims=(-.25, -.25))
    M0, T1, T2 = M0[..., 0], T1[..., 0], T2[..., 0]

    # Take a gander
    plt.imshow(M0)
    plt.title('Even-more-modified Shepp-Logan')
    plt.axis('off')
    plt.colorbar()
    plt.show()
