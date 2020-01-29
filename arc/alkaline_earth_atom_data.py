# -*- coding: utf-8 -*-

"""
Data sources
-------------

.. [#c1] J. A. Armstrong, J. J. Wynne, and P. Esherick,
        "Bound, odd-parity J = 1 spectra of the alkaline earths: Ca, Sr, and Ba,"
        J. Opt. Soc. Am. 69, 211-230 (1979)


.. [#c2]R.Beigang K.Lücke A.Timmermann P.J.West D.Frölich
        Determination of absolute level energies of 5sns1S0 and 5snd1D2 Rydberg series of Sr
        Opt. Commun. 42, 19 1982.

.. [#c3] J E Sansonetti and G Nave. Wavelengths,
        Transition Probabilities, and Energy Levels for the Spectrum of Neutral Strontium (Sr I).
        Journal of Physical and Chemical Reference Data, 39:033103, 2010.

.. [#c4]Baig M Yaseen M Nadeem A Ali R Bhatti S
        Three-photon excitation of strontium Rydberg levels
        Optics Communications, vol: 156 (4-6) pp: 279-284, 1998

.. [#c5] P Esherick, J J Wynne, and J A Armstrong.
        Spectroscopy of 3P0 states of alkaline earths.
        Optics Letters, 1:19, 1977.

.. [#c6] P Esherick.
        Bound, even-parity J = 0 and J = 2 spectra of Sr.
        PhysicalReview A, 15:1920, 1977.

.. [#c7] R Beigang and D Schmidt.
        Two-Channel MQDT Analysis of Bound 5snd 3D1,3 Rydberg States of Strontium.
        Physica Scripta, 27:172, 1983.

.. [#c8]J R Rubbmark and S A Borgstr¨om.
        Rydberg Series in Strontium Found in Absorption by Selectively Laser-Excited Atoms.
        Physica Scripta, 18:196,1978

.. [#c9] Beigang R, Lucke K, Schmidt D, Timmermann A and West P J ¨
        One-Photon Laser Spectroscopy of Rydberg Series from Metastable Levels in Calcium and Strontium
        Phys. Scr. 26 183, 1982

.. [c10] L. Couturier, I. Nosske, F. Hu, C. Tan, C. Qiao, Y. H. Jiang, P. Chen, and M. Weidemüller.
        Measurement of the strontium triplet Rydberg series by depletion spectroscopy of ultracold atoms
        http://arxiv.org/abs/1810.07611

.. [#yb1] Optical-microwave double-resonance spectroscopy of highly excited Rydberg states of ytterbium
        H. Maeda Y. Matsuo M. Takami A. Suzuki
        Physical Review Avol. 45issue 3(1992)pp: 1732-1741Published by American Physical Society

.. [#yb2] Three-step laser spectroscopy and multichannel quantum defect analysis of odd-parity Rydberg states of neutral ytterbium
        M Aymar R J Champeau C Delsart O Robaux
        Journal of Physics B: Atomic and Molecular Physicsvol. 17issue 18(1984)pp: 3645-3661

.. [#yb3] Laser and microwave spectroscopy of even-parity Rydberg states of neutral ytterbium and multichannel-quantum-defect-theory analysis
        H. Lehec A. Zuliani W. Maineult E. Luc-Koenig P. Pillet P. Cheinet F. Niyaz T. F. Gallagher
        Physical Review A vol. 98 issue 6 (2018) pp: 062506 Published by American Physical Society

.. [#ca1] Microwave spectroscopy of calcium Rydberg states
        Thomas R. Gentile Barbara J. Hughey Daniel Kleppner Theodore W. Ducas
        Physical Review Avol. 42issue 1(1990)pp: 440-451Published by American Physical Society

.. [#ca2] Determination of Ionization Potential of Calcium by High-Resolution Resonance Ionization Spectroscopy
        Masabumi Miyabe, Christopher Geppert, Masaaki Kato, Masaki Oba, Ikuo Wakaida, Kazuo Watanabe, Klaus D. A. Wendt
        Journal of the Physical Society of Japan, 75, 034302 (2006) 10.1143/JPSJ.75.034302

.. [#ca3] Meija, Juris; et al. (2016). "Atomic weights of the elements 2013 (IUPAC Technical Report)". Pure and Applied Chemistry. 88 (3): 265–91. doi:10.1515/pac-2015-0305.


.. [#pr] C.B.Alcock, V.P.Itkin, M.K.Horrigan,\
        *Canadian Metallurgical Quarterly*, **23**, 309 (1984)
        http://dx.doi.org/10.1179/cmq.1984.23.3.309
"""

from .alkaline_earth_atom_functions import *


class Strontium88(AlkalineEarthAtom):
    """
    Properties of Strontium 88 atoms
    """

    alphaC = 15

    ionisationEnergy = 5.69486740       #: (eV)  Ref. [#c3]_

    Z = 38
    I = 0.

    #: TODO source
    scaledRydbergConstant = 109736.627 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    quantumDefect = [[[3.26923346261, -0.252029996277, 12.6529707842,
                       0.0, 0.0, 0.0],
                      [2.73329407388, -5.97060805042,
                       - 40.2119216814, 0.0, 0.0, 0.0],
                      [2.38878451407, -48.8061795134,
                          0.128122744619, 0.0, 0.0, 0.0],
                      [0.0921617687317, -2.89264811181,
                          98.7654059257, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[3.26923346261, -0.252029996277, 12.6529707842,
                       0.0, 0.0, 0.0],
                      [2.88651200565, 0.442418088474,
                       - 1.78011356853, 0.0, 0.0, 0.0],
                      [2.66410028047, -0.209248799245,
                       - 8204.50063566, 0.0, 0.0, 0.0],
                      [0.121637481305, -2.57038002901,
                          133.391164866, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[3.26923346261, -0.252029996277, 12.6529707842,
                       0.0, 0.0, 0.0],
                      [2.88243476924, 0.390633779457,
                       - 0.456452019507, 0.0, 0.0, 0.0],
                      [2.6617115361, -15.7900189481,
                       - 7520.34023263, 0.0, 0.0, 0.0],
                      [0.121231620172, -2.84141643149, 169.2044499,
                       0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[3.37044793739, 0.534941918165, -11.1375820384,
                       0.0, 0.0, 0.0],
                      [2.87169218101, 0.451857877887,
                       - 1.64370227677, 0.0, 0.0, 0.0],
                      [2.60700545163, -18.3685877987,
                       - 24643.6669608, 0.0, 0.0, 0.0],
                      [0.121422935185, -2.86416051794,
                          157.683861744, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{0},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{0},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""

    groundStateN = 5

    # levels that are for smaller n than ground level, but are above in energy
    # due to angular part
    extraLevels = [[4, 2, 3, 1], [4, 2, 1, 1],
                   [4, 3, 3, 0],
                   [4, 3, 4, 1], [4, 3, 3, 1], [4, 3, 2, 1],
                   [4, 2, 2, 0]
                   ]

    levelDataFromNIST = "sr_level_data.csv"

    precalculatedDB = "sr_precalculated.db"
    dipoleMatrixElementFile = "sr_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "sr_quadrupole_matrix_elements.npy"

    literatureDMEfilename = 'strontium_literature_dme.csv'

    elementName = 'Sr88'
    meltingPoint = 777 + 273.15  #: in K

    #: NIST Standard reference database, https://dx.doi.org/10.18434/T4FW23
    mass = 87.905619 * physical_constants["atomic mass constant"][0]

    #: TODO what is docstring here? (fitting ranges have been taken from Pauls fits and christophe)
    defectFittingRange = {"1S0": [14, 34], "3S1": [13, 45], "1P1": [14, 28],
                          "3P2": [8, 18], "3P1": [8, 22], "3P0": [8, 15],
                          "1D2": [36, 66], "3D3": [20, 45], "3D2": [22, 37],
                          "3D1": [20, 32], "1F3": [10, 25], "3F4": [10, 24],
                          "3F3": [10, 24], "3F2": [10, 24]}

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Calculates pressure based on Ref. [#pr]_ (accuracy +- 5%).
        """
        if temperature < 298:
            print("WARNING: Sr vapour pressure below 298 K is unknown (small)")
            return 0
        if temperature < self.meltingPoint:
            return 10**(5.006 + 9.226 - 8572 / temperature
                        - 1.1926 * log(temperature) / log(10.))
        else:
            raise ValueError("ERROR: Sr vapour pressure above %.0f C is unknown"
                             % self.meltingPoint)


class Calcium40(AlkalineEarthAtom):
    """
    Properties of Calcium 40 atoms
    """

    ionisationEnergy = 49305.924 / 8065.544  #: eV ref. [#ca2]_

    Z = 20
    I = 0

    #: TODO source
    scaledRydbergConstant = 109735.81 * 1e2 * \
        physical_constants["inverse meter-electron volt relationship"][0]

    quantumDefect = [[[2.33793, -3.96, 0.0, 0.0, 0.0, 0.0],
                      [1.885584, -0.114, -23.8, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.089, -2, 30, 0.0, 0.0, 0.0, 0.0]],
                     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.8833, -0.02, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [1.964709, 0.228, 0.0, 0.0, 0.0, 0.0],
                      [0.8859, 0.13, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[2.440956, 0.35, 0.0, 0.0, 0.0, 0.0],
                      [1.9549, 2.5, -1.60e2, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{0},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{0},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""

    groundStateN = 4
    extraLevels = []  #: TODO unkown if such exist at time of writing

    levelDataFromNIST = "ca_level_data.csv"

    precalculatedDB = "ca_precalculated.db"

    dipoleMatrixElementFile = "ca_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "ca_quadrupole_matrix_elements.npy"

    literatureDMEfilename = 'calcium_literature_dme.csv'

    elementName = 'Ca40'
    meltingPoint = 842 + 273.15  #: in K

    #: NIST Standard Reference Database 108, https://dx.doi.org/10.18434/T4FW23
    mass = 39.962591 * physical_constants["atomic mass constant"][0]

    #: TODO - what is here docstring? (fitting ranges have been taken from Pauls fits and christophe)
    defectFittingRange = {"1S0": [14, 34], "3S1": [13, 45], "1P1": [14, 28],
                          "3P2": [8, 18], "3P1": [8, 22], "3P0": [8, 15],
                          "1D2": [36, 66], "3D3": [20, 45], "3D2": [22, 37],
                          "3D1": [20, 32], "1F3": [10, 25], "3F4": [10, 24],
                          "3F3": [10, 24], "3F2": [10, 24]}

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Calculates pressure based on Ref. [#pr]_ (accuracy +- 5%).
        """
        if temperature < 298:
            print("WARNING: Ca vapour pressure below 298 K is unknown (small)")
            return 0
        if temperature < self.meltingPoint:
            return 10**(5.006 + 10.127 - 9517 / temperature
                        - 1.4030 * log(temperature) / log(10.))
        else:
            raise ValueError("ERROR: Ca vapour pressure above %.0f C is unknown"
                             % self.meltingPoint)


class Ytterbium174(AlkalineEarthAtom):
    """
    Properties of Ytterbium 174 atoms
    """

    ionisationEnergycm = 50443.08  # cm-1  ref. [#yb3]
    ionisationEnergy = ionisationEnergycm / 8065.544  # eV ref.

    Z = 70
    I = 0

    #: TODO source
    scaledRydbergConstant = 109736.627 * 1e2 * \
        physical_constants["inverse meter-electron volt relationship"][0]

    quantumDefect = [[[4.27914, -7.06, 565, 0.0, 0.0, 0.0],
                      [3.95433, -12.33, 1729, 0.0, 0.0, 0.0],
                      [2.71363, -2.01, 0.00E+00, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{0},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{0},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""

    groundStateN = 6
    extraLevels = []  #: TODO unkown if such exist at time of writing

    # : Store list of filenames to read in no data
    levelDataFromNIST = "yb_level_data.csv"

    precalculatedDB = "yb_precalculated.db"
    dipoleMatrixElementFile = "yb_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "yb_quadrupole_matrix_elements.npy"

    literatureDMEfilename = 'ytterbium_literature_dme.csv'

    elementName = 'Yb174'
    meltingPoint = 819 + 273.15  #: in K

    #: NIST Standard reference database, https://dx.doi.org/10.18434/T4FW23
    mass = 173.9388664 * \
        physical_constants["atomic mass constant"][0]

    #: TODO what is docstring here? (fitting ranges have been taken from Pauls fits and christophe)
    defectFittingRange = {"1S0": [15, 43], "1P1": [
        35, 53], "1D2": [28, 75], "3D2": [10, 52]}

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Calculates pressure based on Ref. [#pr]_ (accuracy +- 5%).
        """
        if temperature < 298:
            print("WARNING: Yb vapour pressure below 298 K is unknown (small)")
            return 0
        if temperature < 900:
            return 10**(5.006 + 9.111 - 8111 / temperature
                        - 1.0849 * log(temperature) / log(10.))
        else:
            raise ValueError(
                "ERROR: Yb vapour pressure above 900 K is unknown")
