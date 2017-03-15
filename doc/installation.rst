Installation instructions
=========================
Prerequisite: Python
--------------------

Install Python and packages for scientific computing in Python (scipy, numpy, matplotlib). The package is tested with Python 2.* version, and currently is not supported for Python 3.* version.  We recommend installing Python distributions that comes with Numpy that is connected to the optimized numeric libraries like ATLAS. One such distribution is `Anaconda <https://www.continuum.io/downloads>`_, that provides `ATLAS <https://anaconda.org/anaconda/atlas>`_ support and optimized math kernel.


Download the ARC library/package
--------------------------------

`Download the module from the main site <https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator>`_, unzip the archive and set the folder somewhere within the Python package search path or directly in your project directory. Simply import and use the module::

    >>> from arc import *
    >>> # write your code that uses ARC then.

It is important that package is stored somewhere where user has write permissions, so that it can update the databases with atomic properties.



Compiling C extension
---------------------

Optimized version of the Numerov is provided as the C code `arc_c_extensions.c`.
Note that path to arc directory **should not contain spaces** in order
to setupc.py script to work.

**For Windows users**

Download and install `MinGW <http://www.mingw.org/>`_ compiler, or some other distribution of GNU C compiler. In the command prompt, navigate to the arc folder where `setupc.py` file is located and execute::

    python setupc.py build_ext --inplace

**For Linux users**

Download and install GNU C compiler. Then with terminal open, navigate to arc folder where `setupc.py` file is located execute::

    python setupc.py build_ext --inplace


**For MAC users**

Download and install GNU C compiler. Then with terminal open, navigate to arc folder where `setupc.py` file is located execute::

    python setupc.py build_ext --inplace

**For people who don't want to compile anything**

Alternative solution, if you don't want to compile anything, is to use pure Python implementation of the Numerov, provided in the package. This is done by passing `cpp_numerov = False` flag whenever atoms are initialized, e.g::

    atom = Rubidium(cpp_numerov=False)

This is not recommended option for complex calculations, since it will run much more slowly then optimized C version, but is fine if you need just a few numbers.

**Finally...**

That is all, enjoy using ARC package. Check :ref:`get-started-page` to see some ideas where to start.
