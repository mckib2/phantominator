Installation
============

.. code-block:: bash

    pip install phantominator

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

Usage
=====

Basic usage:

.. code-block:: python

    # CT phantom
    from phantominator import shepp_logan
    ph = shepp_logan(128)

    # MR phantom (returns proton density, T1, and T2 maps)
    M0, T1, T2 = shepp_logan((128, 128, 20), MR=True)

    # Dynamic (concentric circles), 20 time frames
    from phantominator import dynamic
    ph = dynamic(128, 20)

Also see the `examples` module and docstrings.  The interface for CT
phantom generation is similar to MATLAB's `phantom` function.
