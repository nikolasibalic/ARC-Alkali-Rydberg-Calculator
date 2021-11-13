Getting started with ARC
=========================

 .. _get-started-page:

IPython notebook with examples
------------------------------

`Rydberg atoms - a primer`_ introduces Rydberg atoms and ARC package, and is a
good **starting point** to learn how to use ARC to get relevant information about
**alkali atoms** and **Rydberg states** in general.

`An introduction to ARC 3.0: Alkali.ne Rydberg Calculator`_ introduces features
added in ARC 3.0 version: support for **divalent atoms**, inter-species
calculations, atom-surface interactions, dynamic polarizability calculations
(AC Stark Shift), wave function plotting, and methods for work with optical
lattices.
   
`ARC 3.1 update: support for hyperfine structure for alkali atoms`_ expands
support for **alkali metalk**: hyperfine structure is added,
and functions for dealing with Raman transitions and level structures in strong
magnetic fields (Breit-Rabi diagrams).

.. _`Rydberg atoms - a primer`: ./Rydberg_atoms_a_primer_notebook.ipynb

.. _`An introduction to ARC 3.0: Alkali.ne Rydberg Calculator`: ./ARC_3_0_introduction.ipynb

.. _`ARC 3.1 update: support for hyperfine structure for alkali atoms`: ./ARC_3_1_additions.ipynb

.. _`Jupyter`: https://jupyter.org/

Click on the corresponding topic below to open static (HTML) version of the
notebooks. If you want directly .ipynb format, open directly files form ARC
repository.

.. toctree::
   
    Rydberg_atoms_a_primer_notebook
    ARC_3_0_introduction
    ARC_3_1_additions

On demand examples from online Atom calculator
----------------------------------------------

You can try using the package without installing anything on your computer. Simply point your web browser from your computer, tablet or phone to  `atomcalc.jqc.org.uk <https://atomcalc.jqc.org.uk>`_ and use ARC online.

Online version also generates the correct code necessary for answering the questions you ask, which can be downladed and used as a starting point for running the package locally on your computer.

Frequently asked questions (FAQ)
--------------------------------

If you have a question how to do a common calculation, we recommend checking above mentioned `Rydberg atoms - a primer`_ IPython notebook. For general questions about the package usage check here:

**1. How to save calculation (or matrix) for later use?**

Calculations of pair-state interactions :obj:`PairStateInteractions` and Stark maps :obj:`StarkMap` can be easily saved at any point by calling :obj:`alkali_atom_functions.saveCalculation` . This can be loaded later by using :obj:`alkali_atom_functions.loadSavedCalculation` and calculation can be continued from that point.


**2. How to export results?**

If you want to export results e.g. for analysis and plotting in other programs, you can use :obj:`calculations_atom_pairstate.PairStateInteractions.exportData` and :obj:`calculations_atom_single.StarkMap.exportData` to export results of Stark map and Pair-state interaction calculations in **.csv** format. See documentation of corresponding functions for more details.

**3. Calculation is not outputting anything? How long does it take for calculation to finish?**

Most of the functions have `progressOutput` and `debugOutput` as an optional parameter (by default set to False) - check documentation of individual functions for details. We recommend setting at least `progressOutput=True` so that you have minimum output about the status of calculations. This often displays percentage of the current calculation that is finished, that you can use to estimate total time. Setting `debugOutput=True` outputs even more verbose output, like states in the selected basis, and individual coupling strengths etc.
