Single atom calculations
========================

.. note::
  Some of the functions (`Wavefunction`, `AtomSurfaceVdW`,
  `OpticalLattice1D`, `DynamicPolarizability`, and optical materials
  properties) in this module are part of new ARC 3.0.0 version.
  See more at `E. J. Robertson, N. Šibalić, R. M. Potvliege and
  M. P. A. Jones, *Computer Physics Communications* **261**, 107814 (2021)
  `https://doi.org/10.1016/j.cpc.2020.107814 <https://doi.org/10.1016/j.cpc.2020.107814>`_ .

Overview
--------

.. currentmodule:: arc.calculations_atom_single

.. rubric:: Wavefunction Methods

.. autosummary::

    Wavefunction.getRtimesPsiSpherical
    Wavefunction.getRtimesPsi
    Wavefunction.getPsi
    Wavefunction.getRtimesPsiSquaredInPlane
    Wavefunction.plot2D
    Wavefunction.plot3D

.. rubric:: StarkMap Methods

.. autosummary::

    StarkMap.defineBasis
    StarkMap.diagonalise
    StarkMap.plotLevelDiagram
    StarkMap.showPlot
    StarkMap.savePlot
    StarkMap.exportData
    StarkMap.getPolarizability
    StarkMap.getState

.. rubric:: LevelPlot Methods

LevelPlot is also called Grotrian diagram, or term diagram.

.. autosummary::

    LevelPlot.makeLevels
    LevelPlot.drawLevels
    LevelPlot.showPlot

.. rubric:: AtomSurfaceVdW Methods

.. autosummary::

    AtomSurfaceVdW.getC3contribution
    AtomSurfaceVdW.getStateC3

.. rubric:: OpticalLattice1D Methods

.. autosummary::

    OpticalLattice1D.getRecoilEnergy
    OpticalLattice1D.getTrappingFrequency
    OpticalLattice1D.defineBasis
    OpticalLattice1D.diagonalise
    OpticalLattice1D.plotLevelDiagram
    OpticalLattice1D.BlochWavefunction
    OpticalLattice1D.getWannierFunction

.. rubric:: DynamicPolarizability Methods

.. autosummary::

    DynamicPolarizability.defineBasis
    DynamicPolarizability.getPolarizability
    DynamicPolarizability.plotPolarizability

.. currentmodule:: arc.materials

.. rubric:: Optical material properties

.. autosummary::
    OpticalMaterial
    Air
    Sapphire

Detailed documentation
----------------------

.. automodule:: arc.calculations_atom_single
   :members:

.. automodule:: arc.materials
   :show-inheritance:
   :members:
