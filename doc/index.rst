.. Atom calculator documentation master file, created by
   sphinx-quickstart on Fri Jul 29 12:19:10 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*****************
ARC documentation
*****************

.. figure:: ./logo.png
    :width: 94px
    :align: right
    :height: 106px
    :alt: alternate text
    :figclass: align-right

ARC (Alkali.ne Rydberg Calculator)  is package of routines written in Python,
using object-oriented programming (OOP) to make modular, reusable and extendable
collection of routines and data for performing useful calculations of
single atom and two-atom properties, like level diagrams, interactions and
transition strengths for alkali and alkaline earth metal atoms.

Contents
========

.. toctree::
   :maxdepth: 2

   installation
   getting_started
   detailed_doc
   contribute

.. note::
    Support for Alkaline Earth atoms and some of the functions
    (`Wavefunction`, `AtomSurfaceVdW`,
    `OpticalLattice1D`, `DynamicPolarizability`, and optical materials
    properties)in this modules will be
    released in the forthcoming ARC 3.0.0 version.
    To used them now as a beta feature do::
     from arc.beta import *
    Publication describing this upgrade is in preparation (check this place
    soon). For now cite as "E. J. Robertson, N. Šibalić, R. M. Potvliege and
    M. P. A. Jones, *in preparation* (2020)".

Package structure
=================

.. figure:: ./overview_of_modules3.png
    :width: 800px
    :align: center
    :height: 480px
    :alt: module overview image
    :figclass: align-center

    Overview of modules and interdependencies in the :obj:`arc` package. Click on image to enlarge.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Credits
=======

:Authors:
    Nikola Šibalić,
    Jonathan D. Pritchard,
    Charles S. Adams,
    Kevin J. Weatherill

:Licence: BSD 3-Clause

:Version: 3.0.0 of 2020/01/17
