Installation
============

.. code-block:: bash

    pip install phantominator

The goal is to have easy installation and usage for everyone.  If
something doesn't work, please open an issue and/or submit a pull
request.  We'll get it figured out.

About
=====

Python package for easy generation of numerical phantoms.  I often
need a simple image to try something out on.  In MATLAB, I would use
the `phantom` command to quickly get something to work with.  In
Python, it's not quite so easy, so I made this package that's quick
to install and use so there's as little friction as possible.  There
are other implementations of Shepp-Logan available from other
projects, but they are usually not as easy to install or include other
things that I don't want for this project.

This package offers both CT and MR versions.

Going forward, this module will no longer support Python 2.  Please do
the world a favor and move on to Python 3.

Usage
=====

Also see the `examples` module and docstrings.  The interface for CT
phantom generation is similar to MATLAB's `phantom` function.

Basic usage:

.. code-block:: python

    # CT phantom
    from phantominator import shepp_logan
    ph = shepp_logan(128)

    # MR phantom (returns proton density, T1, and T2 maps)
    M0, T1, T2 = shepp_logan((128, 128, 20), MR=True)

The Shepp-Logan 3D phantom has ellipsoids in [-1, 1] along the z-axis.
The 2D Shepp-Logan exists at z=-0.25, so if we want just a subset
along the z-axis with the first slice being the traditional 2D
phantom, we can use the `zlims` option:

.. code-block:: python

    from phantominator import shepp_logan
    M0, T1, T2 = shepp_logan((64, 64, 5), MR=True, zlims=(-.25, .25))

We can also generate simple oscillating concentric circles:

.. code-block:: python

    # Dynamic (concentric circles), 20 time frames
    from phantominator import dynamic
    ph = dynamic(128, 20)

If we want to modify ellipse/ellipsoid parameters or we just want to
see what they are.  For example, we can get the MR ellipsoid
parameters like this:

.. code-block:: python

    from phantominator import mr_ellipsoid_parameters
    E = mr_ellipsoid_parameters()

See docstrings for `ct_shepp_logan`, and `mr_shepp_logan` for how
the array `E` is structured.  It follows conventions from MATLAB's
phantom function.

Arbitrary k-space sampling is supported for the single coil 2D
Shepp-Logan phantom:

.. code-block:: python

    # Given k-space coordinates (kx, ky), where kx and ky are 1D
    # arrays using the same unit conventions as BART's traj command,
    # we can find the corresponding k-space measurements:
    from phantominator import kspace_shepp_logan
    k = kspace_shepp_logan(kx, ky)
