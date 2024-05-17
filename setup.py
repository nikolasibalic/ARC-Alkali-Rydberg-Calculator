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
    ext_modules=[arc_ext],
)
