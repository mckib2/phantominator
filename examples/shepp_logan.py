'''Example demonstrating how to make a Shepp-Logan phantom.'''

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
