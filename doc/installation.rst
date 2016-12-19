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

That is all, enjoy using ARC package. Check :ref:`get-started-page` to see some ideas where to start. If you are **not** using 64-bit Windows/Mac/Linux, please continue reading below.



Precompiled modules on unsupported operating systems 
----------------------------------------------------

Optimized of the Numerov is provided as the C++ code `nvwcpp.cpp`. If you are using 64-bit Windows, Mac or Linux operating system, package will recognize the system and use appropriate precompuled binary provided in the package - i.e. you should not have to do anything additional manually, this is out-of-the box feature.

If you are using some other operating system, or 32-bit version of the operating systems above, to run the optimized precompiled version (recommended) do the following:

**For Windows users**

Downlaod and install `MinGW <http://www.mingw.org/>`_ compiler, or some other distribution of GNU C++ compiler. In the command propt, navigate to the arc folder where `nvwcpp.cpp` file is located and execute::

    g++ -O3 nvwcpp.cpp -o nvwcpp_win

**For Linux users**

Download and install GNU C++ compiler. Then with terminal open, navigate to arc folder where `nvwcpp.cpp` file is located execute::

    g++ -O3 nvwcpp.cpp -o nvwcpp_linux


**For MAC users**

Download and install GNU C++ compiler. Then with terminal open, navigate to arc folder where `nvwcpp.cpp` file is located execute::

    g++ -O3 nvwcpp.cpp -o nvwcpp_mac
    
** For people who don't want to compile anythig **
    
Alternative solution, if you don't want to compile anything, is to use pure Python implementation of the Numerov, provided in the package. This is done by passing `cpp_numerov = False` flag whenever atoms are initialized, e.g::

    atom = Rubidium(cpp_numerov=False)

This is not recommended option for complex calculations, since it will run much more slowly then optimized C++ version, but is fine if you need just a few numbers.
