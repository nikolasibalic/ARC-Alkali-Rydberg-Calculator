
PairStateInteractions
=====================

.. currentmodule:: arc.calculations_atom_pairstate

.. autoclass:: PairStateInteractions
    :members: __init__
    :exclude-members: __init__

.. autosummary::
    :toctree: calculations_atom_pairstate.PairStateInteractions/

Calculate
---------

Full numerical calculation

.. autosummary::
    :toctree: generated/

    ~PairStateInteractions.defineBasis
    ~PairStateInteractions.getC6perturbatively
    ~PairStateInteractions.getLeRoyRadius
    ~PairStateInteractions.diagonalise
cd
Analysis of dominant angular channels and perturbative calculation

.. autosummary::
    :toctree: generated/

    ~PairStateInteractions.getC6perturbativelyAngularChannel
    ~PairStateInteractions.calculateAngularChannelData
    ~PairStateInteractions.loadAngularChannelData
    ~PairStateInteractions.getC6perturbatively_anglePairs

Analyse
-------


.. autosummary::
    :toctree: generated/

    ~PairStateInteractions.getC6fromLevelDiagram
    ~PairStateInteractions.getC3fromLevelDiagram
    ~PairStateInteractions.getVdwFromLevelDiagram
    ~PairStateInteractions.exportData


Visualise
---------

.. autosummary::
    :toctree: generated/

    ~PairStateInteractions.plotLevelDiagram
    ~PairStateInteractions.showPlot
    ~PairStateInteractions.savePlot


Attributes
----------

Internal variables of the class. This is for low-level access to intermediate
results (low level API).

.. autosummary::
    :toctree: generated/

    ~PairStateInteractions.atom1
    ~PairStateInteractions.s1
    ~PairStateInteractions.interactionsUpTo
    ~PairStateInteractions.basisStates
    ~PairStateInteractions.originalPairStateIndex
    ~PairStateInteractions.channel
    ~PairStateInteractions.coupling
    ~PairStateInteractions.matrixElement
    ~PairStateInteractions.matDiagonal
    ~PairStateInteractions.matR
    
    
    
    
    
