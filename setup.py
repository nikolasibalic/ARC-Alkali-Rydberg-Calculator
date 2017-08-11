#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-


import sys
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
    name="arc",
    version="1.2.1",
    description="Alkali Rydberg Calculator",
    license="BSD3",
    keywords=["rydberg", "physics"],
    url="http://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/",

    packages=['arc'],

    package_data={'arc': ['data/*', 'arc_c_extensions.c']},

    include_package_data=True,
    zip_safe=False,
    ext_modules=[arc_ext],
)

