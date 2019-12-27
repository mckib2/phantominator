'''Setup.py'''

from distutils.core import setup
from setuptools import find_packages

setup(
    name='phantominator',
    version='0.5.0',
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/phantominator',
    license='GPLv3',
    description='Generate numerical phantoms.',
    long_description=open('README.rst').read(),
    install_requires=[
        "numpy>=1.17.4",
        "scipy>=1.3.3",
        "matplotlib>=2.1.1"
    ],
    python_requires='>=2.7', # 2.7.15+, to be exact...
)
