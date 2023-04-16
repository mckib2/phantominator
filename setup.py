"""Setup.py"""

from distutils.core import setup
from setuptools import find_packages


setup(
    name="phantominator",
    version="0.7.0",
    author="Nicholas McKibben",
    author_email="nicholas.bgp@gmail.com",
    packages=find_packages(),
    scripts=[],
    url="https://github.com/mckib2/phantominator",
    license="MIT",
    description="Generate numerical phantoms.",
    long_description=open("README.rst").read(),
    install_requires=[
        "numpy>=1.24.2",
        "scipy>=1.9.1",
        "matplotlib>=3.7.1",
        "ssfp>=1.0.0",
    ],
    python_requires=">=3.8",
)
