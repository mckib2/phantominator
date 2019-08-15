'''Demonstrate how to use MR Shepp-Logan.'''

from phantominator import shepp_logan

if __name__ == '__main__':

    # Get proton density, T1, and T2 maps for (L, M, N) sized Shepp-
    # Logan phantom.
    L, M, N = 256, 255, 20
    M0, T1, T2 = shepp_logan((L, M, N), MR=True)
