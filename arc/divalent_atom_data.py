# -*- coding: utf-8 -*-

"""
Data sources
-------------

.. [#c1] J. A. Armstrong, J. J. Wynne and P. Esherick,
        "Bound, odd-parity J = 1 spectra of the alkaline earths: Ca, Sr,
        and Ba",
        *J. Opt. Soc. Am.* **69**, 211-230 (1979)


.. [#c2] R.Beigang, K.Lücke, A.Timmermann, P.J.West and D.Frölich,
        Determination of absolute level energies of 5sns1S0 and 5snd1D2
        Rydberg series of Sr,
        *Opt. Commun.* **42**, 19 (1982).

.. [#c3] J. E. Sansonetti and G Nave,
        Wavelengths, Transition Probabilities, and Energy Levels for the
        Spectrum of Neutral Strontium (Sr I),
        Journal of Physical and Chemical Reference Data **39**, 033103 (2010).

.. [#c4] Baig M, Yaseen M, Nadeem A, Ali R. and Bhatti S.
        Three-photon excitation of strontium Rydberg levels,
        *Optics Communications* **156**, 279 (1998)

.. [#c5] P. Esherick, J. J. Wynne and J A Armstrong,
        Spectroscopy of 3P0 states of alkaline earths,
        *Optics Letters* **1**, 19 (1977).

.. [#c6] P Esherick,
        Bound, even-parity J = 0 and J = 2 spectra of Sr,
        *Physical Review A* **15**, 1920 (1977).

.. [#c7] R. Beigang and D. Schmidt,
        Two-Channel MQDT Analysis of Bound 5snd 3D1,3 Rydberg States of
        Strontium,
        *Physica Scripta* **27**, 172 (1983).

.. [#c8] J R. Rubbmark and S. A. Borgstr¨om,
         Rydberg Series in Strontium Found in Absorption by Selectively,
         Laser-Excited Atoms.
         *Physica Scripta* **18**, 196 (1978)

.. [#c9] Beigang R, Lucke K, Schmidt D, Timmermann A. and West P. J,
        One-Photon Laser Spectroscopy of Rydberg Series from Metastable Levels
        in Calcium and Strontium,
        *Phys. Scr.* **26**, 183 (1982)

.. [#c10] L. Couturier, I. Nosske, F. Hu, C. Tan, C. Qiao, Y. H. Jiang, P. Chen
        and M. Weidemüller.
        Measurement of the strontium triplet Rydberg series by depletion
        spectroscopy of ultracold atoms
        http://arxiv.org/abs/1810.07611

.. [#yb1] H. Maeda, Y. Matsuo, M. Takami and A. Suzuki,
        Optical-microwave double-resonance spectroscopy of highly excited
        Rydberg states of ytterbium,
        *Physical Review A* **45**, 1732 (1992)

.. [#yb2] M. Aymar, A. Debarre and O. Robaux,
        Highly excited levels of neutral ytterbium. II. Multichannel quantum
        defect analysis of odd- and even-parity spectra,
        *Journal of Physics B: Atomic and Molecular Physics* **13**,
        1089 (1980)
        https://doi.org/10.1088/0022-3700/13/6/016

.. [#yb3] H. Lehec, A. Zuliani, W. Maineult, E. Luc-Koenig, P. Pillet,
        P. Cheinet, F. Niyaz and T. F. Gallagher,
        Laser and microwave spectroscopy of even-parity Rydberg states of
        neutral ytterbium and multichannel-quantum-defect-theory analysis,
        *Physical Review A* **98**, 062506 (2018)

.. [#MT78] W. F. Meggers and J. L. Tech, *J. Res. Natl. Bur. Stand.* (U.S.)
        **83**, 13 (1978).

.. [#ca1] Thomas R. Gentile, Barbara J. Hughey, Daniel Kleppner
        and Theodore W. Ducas,
        Microwave spectroscopy of calcium Rydberg states,
        *Physical Review A* **42**, 440 (1990)

.. [#ca2] Masabumi Miyabe, Christopher Geppert, Masaaki Kato, Masaki Oba,
        Ikuo Wakaida, Kazuo Watanabe and Klaus D. A. Wendt,
        Determination of Ionization Potential of Calcium by High-Resolution
        Resonance Ionization Spectroscopy,
        *Journal of the Physical Society of Japan* **75**, 034302 (2006)
        https://doi.org/10.1143/JPSJ.75.034302

.. [#ca3] Meija, Juris et al,
        "Atomic weights of the elements 2013 (IUPAC Technical Report)",
        *Pure and Applied Chemistry* **88**,265 (2016)
        https://doi.org/10.1515/pac-2015-0305.

.. [#ca4] B.B. Zelener, S. A. Saakyan,V. A. Sautenkov, E. V. Vilshanskaya,
        B. V. Zelener and V. E. Fortov, *JETP Letters* **110**, 761 (2019)
        https://doi.org/10.1134/S0021364019240093

.. [#ca5] J. Sugar and C. Corliss, *J. Phys. Chem. Ref. Data* **14**,
        Suppl. 2 (1985)

.. [#pr] C.B.Alcock, V.P.Itkin, M.K.Horrigan,\
        *Canadian Metallurgical Quarterly*, **23**, 309 (1984)
        http://dx.doi.org/10.1179/cmq.1984.23.3.309

.. [#nist] NIST Standard reference database, https://dx.doi.org/10.18434/T4FW23
"""

from .divalent_atom_functions import *


class Strontium88(DivalentAtom):
    """
    Properties of Strontium 88 atoms
    """

    alphaC = 15

    ionisationEnergy = 1377012721e6 * C_h / C_e  #: (eV)  Ref. [#c10]_

    Z = 38
    I = 0.0

    #: Ref. [#c10]_
    scaledRydbergConstant = (
        109736.631
        * 1.0e2
        * physical_constants["inverse meter-electron volt relationship"][0]
    )

    quantumDefect = [
        [
            [3.269123, -0.177769, 3.4619, 0.0, 0.0, 0.0],
            [2.72415, -3.390, -220.0, 0.0, 0.0, 0.0],
            [2.384667, -42.03053, -619.0, 0.0, 0.0, 0.0],
            [0.090886, -2.4425, 61.896, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [3.3707725, 0.41979, -0.421377, 0.0, 0.0, 0.0],
            [2.88673, 0.433745, -1.800, 0.0, 0.0, 0.0],
            [2.675236, -13.23217, -4418.0, 0.0, 0.0, 0.0],
            [0.120588, -2.1847, 102.98, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [3.3707725, 0.41979, -0.421377, 0.0, 0.0, 0.0],
            [2.88265, 0.39398, -1.1199, 0.0, 0.0, 0.0],
            [2.661488, -16.8524, -6629.26, 0.0, 0.0, 0.0],
            [0.11899, -2.0446, 103.26, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [3.3707725, 0.41979, -0.421377, 0.0, 0.0, 0.0],
            [2.88163, -2.462, 145.18, 0.0, 0.0, 0.0],
            [2.655, -65.317, -13576.7, 0.0, 0.0, 0.0],
            [0.12000, -2.37716, 118.97, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    ]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{1},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{1},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""

    groundStateN = 5

    # levels that are for smaller n than ground level, but are above in energy
    # due to angular part
    extraLevels = [
        [4, 2, 3, 1],
        [4, 2, 1, 1],
        [4, 3, 3, 0],
        [4, 3, 4, 1],
        [4, 3, 3, 1],
        [4, 3, 2, 1],
        [4, 2, 2, 0],
    ]

    #: Sources Refs. [#c1]_, [#c2]_, [#c3]_, [#c4]_, [#c5]_, [#c6]_, [#c7]_,
    #: [#c8]_ , [#c10]_
    levelDataFromNIST = "sr_level_data.csv"

    precalculatedDB = "sr88_precalculated.db"
    dipoleMatrixElementFile = "sr_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "sr_quadrupole_matrix_elements.npy"

    literatureDMEfilename = "strontium_literature_dme.csv"

    elementName = "Sr88"
    meltingPoint = 777 + 273.15  #: in K

    #: Ref. [#nist]_
    mass = 87.905619 * physical_constants["atomic mass constant"][0]

    #: Quantum defect principal quantum number fitting ranges for different
    #: series
    defectFittingRange = {
        "1S0": [14, 34],
        "3S1": [15, 50],
        "1P1": [10, 29],
        "3P2": [19, 41],
        "3P1": [8, 21],
        "3P0": [8, 15],
        "1D2": [20, 50],
        "3D3": [20, 37],
        "3D2": [28, 50],
        "3D1": [28, 50],
        "1F3": [10, 28],
        "3F4": [10, 28],
        "3F3": [10, 24],
        "3F2": [10, 24],
    }

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Calculates pressure based on Ref. [#pr]_ (accuracy +- 5%).
        """
        if temperature < 298:
            print("WARNING: Sr vapour pressure below 298 K is unknown (small)")
            return 0
        if temperature < self.meltingPoint:
            return 10 ** (
                5.006
                + 9.226
                - 8572 / temperature
                - 1.1926 * log(temperature) / log(10.0)
            )
        else:
            raise ValueError(
                "ERROR: Sr vapour pressure above %.0f C is unknown" % self.meltingPoint
            )


class Calcium40(DivalentAtom):
    """
    Properties of Calcium 40 atoms
    """

    #: eV Ref. [#ca4]_
    ionisationEnergy = (
        49305.91966
        * 1e2
        * physical_constants["inverse meter-electron volt relationship"][0]
    )

    Z = 20
    I = 0

    #: eV Ref. [#ca2]_
    scaledRydbergConstant = (
        109735.81037
        * 1e2
        * physical_constants["inverse meter-electron volt relationship"][0]
    )

    quantumDefect = [
        [
            [2.33793, -0.1142, 0.0, 0.0, 0.0, 0.0],
            [1.885584, -0.3240, -23.8, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.09864, -1.29, 36, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [2.440956, 0.35, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.8833, -0.02, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [2.440956, 0.35, 0.0, 0.0, 0.0, 0.0],
            [1.964709, 0.228, 0.0, 0.0, 0.0, 0.0],
            [0.8859, 0.13, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [2.440956, 0.35, 0.0, 0.0, 0.0, 0.0],
            [1.9549, 2.5, -1.60e2, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    ]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{1},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{1},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""

    groundStateN = 4
    extraLevels = []  #: TODO unkown if such exist at time of writing

    #: Sources Refs. [#c1]_, [#c5]_, [#c9]_, [#ca1]_, [#ca5]_
    levelDataFromNIST = "ca_level_data.csv"

    precalculatedDB = "ca40_precalculated.db"

    dipoleMatrixElementFile = "ca_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "ca_quadrupole_matrix_elements.npy"

    literatureDMEfilename = "calcium_literature_dme.csv"

    elementName = "Ca40"
    meltingPoint = 842 + 273.15  #: in K

    #: Ref. [#nist]_
    mass = 39.962591 * physical_constants["atomic mass constant"][0]

    #: Quantum defect principal quantum number fitting ranges for different
    #: series
    defectFittingRange = {
        "1S0": [22, 55],
        "3S1": [22, 55],
        "1P1": [22, 55],
        "3P2": [8, 18],
        "3P1": [22, 55],
        "3D2": [22, 55],
        "3D1": [22, 55],
        "1F3": [20, 150],
    }

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Calculates pressure based on Ref. [#pr]_ (accuracy +- 5%).
        """
        if temperature < 298:
            print("WARNING: Ca vapour pressure below 298 K is unknown (small)")
            return 0
        if temperature < self.meltingPoint:
            return 10 ** (
                5.006
                + 10.127
                - 9517 / temperature
                - 1.4030 * log(temperature) / log(10.0)
            )
        else:
            raise ValueError(
                "ERROR: Ca vapour pressure above %.0f C is unknown" % self.meltingPoint
            )


class Ytterbium174(DivalentAtom):
    """
    Properties of Ytterbium 174 atoms
    """

    ionisationEnergy = (
        50443.07041
        * 1e2
        * physical_constants["inverse meter-electron volt relationship"][0]
    )
    #: eV Ref. [#yb3]_

    Z = 70
    I = 0

    #: eV Ref. [#yb3]_
    scaledRydbergConstant = (
        109736.96959
        * 1e2
        * physical_constants["inverse meter-electron volt relationship"][0]
    )

    quantumDefect = [
        [
            [4.278367, -5.60943, -258.5, 0.0, 0.0, 0.0],
            [3.953434, -10.58286, 728.100, 0.0, 0.0, 0.0],
            [2.7130117, -0.929878, -636.4, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [2.7485996, 0.0137, -106.55, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    ]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{1},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{1},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""

    groundStateN = 6
    extraLevels = []  #: TODO unkown if such exist at time of writing

    #: Sources Refs. [#yb1]_, [#yb2]_, [#yb3]_, [#MT78]_
    levelDataFromNIST = "yb_level_data.csv"

    precalculatedDB = "yb174_precalculated.db"
    dipoleMatrixElementFile = "yb_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "yb_quadrupole_matrix_elements.npy"

    literatureDMEfilename = "ytterbium_literature_dme.csv"

    elementName = "Yb174"
    meltingPoint = 819 + 273.15  #: in K

    #: Ref. [#nist]_
    mass = 173.9388664 * physical_constants["atomic mass constant"][0]

    #: Quantum defect principal quantum number fitting ranges for different
    #: series
    defectFittingRange = {
        "1S0": [34, 80],
        "1P1": [35, 54],
        "1D2": [40, 80],
        "3D2": [35, 80],
    }

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Calculates pressure based on Ref. [#pr]_ (accuracy +- 5%).
        """
        if temperature < 298:
            print("WARNING: Yb vapour pressure below 298 K is unknown (small)")
            return 0
        if temperature < 900:
            return 10 ** (
                5.006
                + 9.111
                - 8111 / temperature
                - 1.0849 * log(temperature) / log(10.0)
            )
        else:
            raise ValueError("ERROR: Yb vapour pressure above 900 K is unknown")
