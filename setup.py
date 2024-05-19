#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-


from setuptools import setup, Extension
import platform


class get_numpy_include_dirs(object):
    """Returns a list of include directories for Numpy after lazy loading;
    Ensure that Numpy is installed before referencing it."""

    def __str__(self):
        import numpy as np

        return np.get_include()


if platform.system() == "Windows":
    c_args = ["/Wall", "/O2"]
else:
    c_args = ["-Wall", "-O3"]

arc_ext = Extension(
    "arc.arc_c_extensions",
    sources=["arc/arc_c_extensions.c"],
    extra_compile_args=c_args,
    include_dirs=[get_numpy_include_dirs()],
)


setup(
    ext_modules=[arc_ext],
)
