
AC Stark Maps
=============

.. currentmodule:: arc.calculations_atom_single

The following classes provide methods for calculating AC Stark Maps for Rydberg manifolds
using Floquet or other approximate methods.

Basis Generation
----------------

.. autoclass:: BasisGenerator
    :members: __init__
    :exclude-members: __init__

.. autosummary::
    :toctree: calculations_atom_single.BasisGenerator/

.. autosummary::
    :toctree: generated/

    ~BasisGenerator.defineBasis

Attributes
++++++++++

.. autosummary::
    :toctree: generated/

    ~BasisGenerator.basisStates
    ~BasisGenerator.H
    ~BasisGenerator.V

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

    ~ShirleyMethod.defineShirleyHamiltonian
    ~ShirleyMethod.diagonalise

Attributes
++++++++++

.. autosummary::
    :toctree: generated/

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
+++++++++++++++++

.. autoclass:: RWAModel
    :members: __init__
    :exclude-members: __init__

.. autosummary::
    :toctree: calculations_atom_single.RWAModel/

Calculate
+++++++++

.. autosummary::
    :toctree: generated/

    ~RWAModel.findDipoleCoupledStates
    ~RWAModel.RWAModel

Attributes
++++++++++

.. autosummary::
    :toctree: generated/

    ~RWAModel.dipoleCoupledStates
    ~RWAModel.dipoleCoupledFreqs
    ~RWAModel.starkShifts