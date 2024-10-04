# Google-style docstrings are used for documenting code
# for more on those, see http://www.sphinx-doc.org/en/stable/ext/example_google.html#example-google
# this is automatically parsed by Aphinx module in Python, using napoleon addon
from __future__ import division, print_function, absolute_import

__version__ = "3.6.0"

__all__ = [
    "AlkaliAtom",
    "printState",
    "printStateString",
    "printStateStringLatex",
    "printStateLetter",
    "formatNumberSI",
    "Hydrogen",
    "Caesium",
    "Cesium",
    "Rubidium85",
    "Rubidium",
    "Rubidium87",
    "Lithium6",
    "Lithium7",
    "Sodium",
    "Potassium",
    "Potassium39",
    "Potassium40",
    "Potassium41",
    "Strontium88",
    "Calcium40",
    "Ytterbium174",
    "Ylm",
    "Wavefunction",
    "StarkMap",
    "LevelPlot",
    "AtomSurfaceVdW",
    "OpticalLattice1D",
    "DynamicPolarizability",
    "StarkBasisGenerator",
    "ShirleyMethod",
    "RWAStarkShift",
    "PairStateInteractions",
    "StarkMapResonances",
    "Wigner3j",
    "Wigner6j",
    "TriaCoeff",
    "CG",
    "WignerDmatrix",
    "C_k",
    "C_c",
    "C_h",
    "C_e",
    "C_m_e",
    "getCitationForARC",
]

from arc.alkali_atom_functions import (
    AlkaliAtom,
    printState,
    printStateString,
    printStateStringLatex,
    printStateLetter,
    formatNumberSI,
)
from arc.alkali_atom_data import (
    Hydrogen,
    Caesium,
    Cesium,
    Rubidium85,
    Rubidium,
    Rubidium87,
    Lithium6,
    Lithium7,
    Sodium,
    Potassium,
    Potassium39,
    Potassium40,
    Potassium41,
)
from arc.divalent_atom_data import Strontium88, Calcium40, Ytterbium174
from arc.calculations_atom_single import (
    Ylm,
    Wavefunction,
    StarkMap,
    LevelPlot,
    AtomSurfaceVdW,
    OpticalLattice1D,
    DynamicPolarizability,
    StarkBasisGenerator,
    ShirleyMethod,
    RWAStarkShift,
)
from arc.calculations_atom_pairstate import (
    PairStateInteractions,
    StarkMapResonances,
)
from arc.wigner import Wigner3j, Wigner6j, TriaCoeff, CG, WignerDmatrix
from arc._database import getCitationForARC
from scipy.constants import k as C_k
from scipy.constants import c as C_c
from scipy.constants import h as C_h
from scipy.constants import e as C_e
from scipy.constants import m_e as C_m_e
