Installation instructions
=========================
Prerequisite: Python
--------------------

Install Python and packages for scientific computing in Python (scipy, numpy, matplotlib). The package is tested and works with **both Python 2.7 and Python 3.5**.  We recommend installing Python distributions that comes with Numpy that is connected to the optimized numeric libraries like ATLAS. One such distribution is `Anaconda <https://www.continuum.io/downloads>`_, that provides `ATLAS <https://anaconda.org/anaconda/atlas>`_ support and optimized math kernel.


Download the ARC library/package
--------------------------------

`Download latest release for your operating system <https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/releases>`_, unzip the archive and set the folder somewhere within the Python package search path or directly in your project directory. Simply import and use the module::

    from arc import *
    # write your code that uses ARC then.

It is important that package is stored somewhere where user has write permissions, so that it can update the databases with atomic properties. **This is the end of the standard installation for majority of the users.**


Installation of the package globally with setup.py
--------------------------------------------------

**This is tested on Linux so far**

Do this only if you have Linux/UNIX (as it is tested on it) and you are sure that you don't want to change underlying ARC code. 
Make sure you have C compiler and `python development headers <[https://anaconda.org/StatisKit/python-dev>`_ installed. To compile and install for local user ARC call from terminal::

    python setup.py build
    python setup.py install

Databases that have to be changed with new values will be locally copied from package data location to `~.arc-data` folder when arc is used for the first time.

Compiling C extension
----------------------

If you do need to compile C extension yourself, this is how to do it without 
installing globally the package (as in the previos section 
"Installation of the package globally with setup.py").
Optimized version of the Numerov is provided as the C code `arc_c_extensions.c`.
**You don't need to perform this step** of manual compilation of that code if you
followed recommended installation instruction by downloading **precompiled
binary distribution** for the latest `release <https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/releases>`_ .
Note that path to arc directory **should not contain spaces** in order
to setupc.py script to work.

**For Windows users**

If precompiled binaries don't work, please contact developers. Compiling Numpy C
extensions on Windows is a bit complicated due to use of C89 standard (instead of C99). Procedure is the following.
One needs to use `MSVC compiler <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_
in order to compile Numpy extension for Python 2.7 under Windows. For other
Python versions (3.5) find correct compiler `here <https://www.scipy.org/scipylib/building/windows.html#microsoft-visual-c-msvc>`_ .
After installation of the compiler, find in Start menu "Visual C++ 2008 32-bit Command Prompt"
(for 32-bit Python) or "Visual C++ 2008 64-bit Command Prompt" (for 64-bit Python).
Set the following variables set in the command prompt environment::

  SET DISTUTILS_USE_SDK=1
  SET MSSdk=1
  python setupc.py build_ext --inplace

This should build C Numpy extension (implementing Numerov integration)
under Windows. We recommend, however, using
pre-build binaries available on the `release page <https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/releases>`_ .

**For Linux users**

Download and install GNU C compiler. Then with terminal open, navigate to arc folder where `setupc.py` file is located execute::

    python setupc.py build_ext --inplace


**For Mac users**

Download and install GNU C compiler. Then with terminal open, navigate to arc folder where `setupc.py` file is located execute::

    python setupc.py build_ext --inplace

**Slow alternative: Numerov implemented in pure Python**

Alternative solution, if you don't want to compile anything, is to use pure Python implementation of the Numerov, provided in the package. This is done by passing `cpp_numerov = False` flag whenever atoms are initialized, e.g::

    atom = Rubidium(cpp_numerov=False)

This is not recommended option for complex calculations, since it will run much more slowly then optimized C version, but is fine if you need just a few numbers.

**Finally...**

That is all, enjoy using ARC package. Check :ref:`get-started-page` to see some ideas where to start.
