import sqlite3
import numpy as np

sqlite3.register_adapter(np.float64, float)
sqlite3.register_adapter(np.float32, float)
sqlite3.register_adapter(np.int64, int)
sqlite3.register_adapter(np.int32, int)


class UsedModulesARC(object):
    alkali_atoms = False
    divalent_atoms = False
    arc3_0_methods = False
    hyperfine = False
    advanced_getPopulationTime = False
    ac_stark = False


def getCitationForARC():
    """Returns recommended citation for ARC based on used ARC modules."""

    citation = ""

    if UsedModulesARC.divalent_atoms or UsedModulesARC.arc3_0_methods:
        citation += """- E. J. Robertson, N. Šibalić, R. M. Potvliege, M. P. A. Jones, ARC 3.0: An expanded Python toolbox for atomic physics calculations, Computer Physics Communications 261, 107814 (2021) https://doi.org/10.1016/j.cpc.2020.107814\n"""
        if UsedModulesARC.hyperfine:
            citation += """- N. Šibalić, J. D. Pritchard, K. J. Weatherill, C. S. Adams, ARC: An open-source library for calculating properties of alkali Rydberg atoms, Computer Physics Communications 220, 319 (2017), https://doi.org/10.1016/j.cpc.2017.06.015\n"""
    elif UsedModulesARC.alkali_atoms:
        citation += """- N. Šibalić, J. D. Pritchard, K. J. Weatherill, C. S. Adams, ARC: An open-source library for calculating properties of alkali Rydberg atoms, Computer Physics Communications 220, 319 (2017), https://doi.org/10.1016/j.cpc.2017.06.015\n"""

    if UsedModulesARC.advanced_getPopulationTime:
        citation += """- M. Archimi, C. Simonelli, L. Di Virgilio, A. Greco, M. Ceccanti, E. Arimondo, D. Ciampini, I. I. Ryabtsev, I. I. Beterov, and O. Morsch, Phys. Rev. A 100, 030501(R) (2019) https://doi.org/10.1103/PhysRevA.100.030501\n"""

    if UsedModulesARC.ac_stark:
        citation += """- D. H. Meyer, Z. A. Castillo, K. C. Cox, P. D. Kunz, J. Phys. B: At. Mol. Opt. Phys, 53, 034001 (2020) https://doi.org/10.1088/1361-6455/ab6051"""

    return citation
