#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-


try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension


class get_numpy_include_dirs(object):
    """Returns a list of include directories for Numpy after lazy loading;
    Ensure that Numpy is installed before referencing it."""
    
    def __str__(self):
        import numpy as np
        return np.get_include()


arc_ext = Extension(
    "arc.arc_c_extensions",
    sources=["arc/arc_c_extensions.c"],
    extra_compile_args=["-Wall", "-O3"],
    include_dirs=[get_numpy_include_dirs()],
)


setup(
    name="ARC-Alkali-Rydberg-Calculator",
    version="3.2.8",
    description="Alkali Rydberg Calculator",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="BSD3",
    keywords=[
        "rydberg",
        "physics",
        "stark maps",
        "atom interactions",
        "quantum optics",
        "van der Waals",
        "scientific",
        "atomic sensors",
        "quantum simulator",
        "alkali atoms",
        "alkaline atoms",
        "divalent atoms",
        "quantum computing",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "License :: OSI Approved :: BSD License",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Development Status :: 5 - Production/Stable",
    ],
    setup_requires=["oldest-supported-numpy"],
    url="https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator",
    download_url="https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/archive/refs/tags/v3.2.8.tar.gz",
    author="Nikola Sibalic,  Elizabeth J. Robertson, Jonathan D. Pritchard, Robert M. Potvliege, Matthew P. A. Jones, Charles S. Adams, Kevin J. Weatherill, and contributors",
    author_email="nikolasibalic@physics.org",
    packages=["arc", "arc.advanced"],
    package_data={
        "arc": ["data/*", "data/refractive_index_data/*", "arc_c_extensions.c"]
    },
    install_requires=[
        "scipy>=0.18.1",
        "numpy>=1.16.0",
        "matplotlib>=1.5.3",
        "sympy>=1.1.1",
        "lmfit>=0.9.0",
    ],
    zip_safe=False,
    ext_modules=[arc_ext],
)
