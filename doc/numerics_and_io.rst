Numerics and IO functions
=========================

Numerical routines, input and output of calculations (saving and loading
calculations), printing state strings etc.

.. currentmodule:: arc.alkali_atom_functions

Numerics
--------

.. autosummary::
    :toctree: generated/

    NumerovBack

Input, output, formating
------------------------

.. autosummary::
    :toctree: generated/

    saveCalculation
    loadSavedCalculation
    printState
    printStateString
    printStateStringLatex
    printStateLetter
    formatNumberSI

Citations
---------

.. currentmodule:: arc._database

.. autosummary::
    :toctree: generated/

    getCitationForARC


HTML output
-----------

Developed for Atom calculator (ARC web application access).

.. currentmodule:: arc.web_functionality

.. autosummary::
    :toctree: generated/

    htmlLiteratureOutput
    rabiFrequencyWidget
    printValueString
    plotStarkMap
    plotInteractionLevels
    webPlot