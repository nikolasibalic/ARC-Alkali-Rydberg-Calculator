#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-


import sys
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
from numpy.distutils.misc_util import get_numpy_include_dirs


arc_ext = Extension(
            'arc.arc_c_extensions',
            sources = ['arc/arc_c_extensions.c'],
            extra_compile_args = ['-Wall','-O3'],
            include_dirs=get_numpy_include_dirs(),
        )


setup(
    name="ARC-Alkali-Rydberg-Calculator",
    version="2.0.11",
    description="Alkali Rydberg Calculator",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    license="BSD3",
    keywords=["rydberg", "physics", "stark maps", "atom interactions",
              "quantum optics", "van der Waals", "scientific","atomic-sensors",
              "quantum-simulator","alkali-atoms"],
    url="https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator",
    download_url="https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/archive/2.0.11.tar.gz",
    author = 'Nikola Sibalic, Jonathan D. Pritchard, Charles S. Adams, Kevin J. Weatherill',
    author_email = 'nikolasibalic@physics.org',

    packages=['arc'],

    package_data={'arc': ['data/*', 'arc_c_extensions.c']},

    zip_safe=False,
    ext_modules=[arc_ext],

)
