
AC Stark Maps
=============

.. currentmodule:: arc.calculations_atom_single

The following classes provide methods for calculating AC Stark Maps for Rydberg manifolds
using Floquet or other approximate methods.

Basis Generation
----------------

.. autoclass:: StarkBasisGenerator
    :members: __init__
    :exclude-members: __init__

.. autosummary::
    :toctree: calculations_atom_single.StarkBasisGenerator/

.. autosummary::
    :toctree: generated/

    ~StarkBasisGenerator.defineBasis

Attributes
++++++++++

.. autosummary::
    :toctree: generated/

    ~StarkBasisGenerator.basisStates
    ~StarkBasisGenerator.H
    ~StarkBasisGenerator.V
    ~StarkBasisGenerator.bareEnergies

Shirley's Method
----------------

.. autoclass:: ShirleyMethod
    :members: __init__
    :exclude-members: __init__

.. autosummary::
    :toctree: calculations_atom_single.ShirleyMethod/

Calculate
+++++++++

.. autosummary::
    :toctree: generated/

    ~ShirleyMethod.defineBasis
    ~ShirleyMethod.defineShirleyHamiltonian
    ~ShirleyMethod.diagonalise

Attributes
++++++++++

.. autosummary::
    :toctree: generated/

    ~StarkBasisGenerator.basisStates
    ~ShirleyMethod.H0
    ~ShirleyMethod.B
    ~ShirleyMethod.dT
    ~ShirleyMethod.eFields
    ~ShirleyMethod.freqs
    ~ShirleyMethod.targetShifts
    ~ShirleyMethod.eigs
    ~ShirleyMethod.eigVectors
    ~ShirleyMethod.transProbs

RWA Approximation
-----------------

.. autoclass:: RWAStarkShift
    :members: __init__
    :exclude-members: __init__

.. autosummary::
    :toctree: calculations_atom_single.RWAStarkShift/

Calculate
+++++++++

.. autosummary::
    :toctree: generated/

    ~RWAStarkShift.defineBasis
    ~RWAStarkShift.findDipoleCoupledStates
    ~RWAStarkShift.makeRWA

Attributes
++++++++++

.. autosummary::
    :toctree: generated/
    
    ~StarkBasisGenerator.basisStates
    ~RWAStarkShift.dipoleCoupledStates
    ~RWAStarkShift.dipoleCoupledFreqs
    ~RWAStarkShift.starkShifts