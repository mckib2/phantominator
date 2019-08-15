'''Driving function.

Notes
-----
Calls to CT, MR, 2D, and 3D variants of Shepp-Logan will be handled
by this handler function.  It serves to point to the correct
implementation of Shepp-Logan to use given a set of parameters.
'''

from phantominator import mr_shepp_logan, ct_shepp_logan

def shepp_logan(*args, **kwargs):
    '''Shepp-Logan phantom.

    Notes
    -----
    See phantominator.mr_shepp_logan() and
    phantominator.ct_shepp_logan() for docstrings explaining usage.
    '''

    MR = kwargs.get('MR', False)
    if 'MR' in kwargs:
        del kwargs['MR']

    if MR:
        return mr_shepp_logan(*args, **kwargs)
    return ct_shepp_logan(*args, **kwargs)
