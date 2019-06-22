'''Setup.py'''

from distutils.core import setup
from setuptools import find_packages

setup(
    name='phantominator',
    version='0.0.6',
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/phantominator',
    license='GPLv3',
    description='Generate numerical phantoms.',
    long_description=open('README.rst').read(),
    install_requires=[
        "numpy>=1.16.2",
        "matplotlib>=3.0.3"
    ],
    python_requires='>=3.6',
)
