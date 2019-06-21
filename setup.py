'''Setup.py'''

from distutils.core import setup
from setuptools import find_packages

setup(
    name='phantominator',
    version='0.0.1',
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/phantominator',
    license='',
    description='Generate numerical phantoms.',
    long_description='', #open('README.rst').read(),
    install_requires=[
        "numpy>=1.16.2",
    ],
    python_requires='>=3.6',
)
