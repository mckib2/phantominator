'''Demonstrate how to use MR Shepp-Logan.'''

from phantominator import mr_shepp_logan

if __name__ == '__main__':

    N = 128
    M0, T1, T2 = mr_shepp_logan(N, T2star=True)

    print(M0.shape)

    from mr_utils import view
    view(M0, montage_axis=-1)
