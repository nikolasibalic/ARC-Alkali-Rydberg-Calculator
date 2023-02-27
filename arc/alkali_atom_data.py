# -*- coding: utf-8 -*-
"""
    This module specifies properties of individual alkali metals.

    If you want to change e.g. coefficients used for model potential, \
    quantum defects, or other numerical values, this is the place to look at.

    **How to delete precalculated dipole/quadrupole matrix elements values
    and/or start a new database?** To delete precalculated values, simply
    delete files, whose names are stated in `dipoleMatrixElementFile`,
    `quadrupoleMatrixElementFile` and `precalculatedDB` variables for the
    corresponding atom type, from data/ folder. Alternatively, if you
    want to keep old values, but want to also start completely new
    calculation of dipole matrix elements (e.g. because you changed
    parameters of energy levels significantly or model potential parameters),
    simply set new values for `dipoleMatrixElementFile`,
    `quadrupoleMatrixElementFile` and `precalculatedDB` variables.

    Note that by default isotopes of Rubidium and Potassium
    are sharing precalculated dipole and quadrupole matrix elements. This is
    because the small energy level differences typically don't change this
    matrix elements within a typical accuracy.

    Data sources
    -------------

    .. [#c1] M. Marinescu, H. R. Sadeghpour, and A. Dalgarno, \
            *Phys.Rev.A* **49**, 982 (1994)
            https://doi.org/10.1103/PhysRevA.49.982

    .. [#Weber1987] K.-H. Weber and Craig J. Sansonetti,
            *Phys.Rev.A* **35**, 4650 (1987)

    .. [#c3] C.B.Alcock, V.P.Itkin, M.K.Horrigan,\
            *Canadian Metallurgical Quarterly*, **23**, 309 (1984)
            http://dx.doi.org/10.1179/cmq.1984.23.3.309

    .. [#c4] Wenhui Li, I. Mourachko, M. W. Noel, and T. F. Gallagher, \
            *Phys. Rev. A* **67**, 052502 (2003)
            https://doi.org/10.1103/PhysRevA.67.052502

    .. [#c5] Jianing Han, Yasir Jamil, D. V. L. Norum, Paul J. Tanner, \
            and T. F. Gallagher, *Phys. Rev. A* **74**, 054502 (2006)
            https://doi.org/10.1103/PhysRevA.74.054502

    .. [#Mack2011] Markus Mack, Florian Karlewski, Helge Hattermann, Simone Hockh,
            Florian Jessen, Daniel Cano, and Jozsef Fortagh, *Phys. Rev. A* **83**,
            052515 (2011), https://doi.org/10.1103/PhysRevA.83.052515

    .. [#Afrousheh2006a] K. Afrousheh, P. Bohlouli-Zanjani, J. A. Petrus, and
            J. D. D. Martin, *Phys. Rev. A* **74**, 062712 (2006)
            https://doi.org/10.1103/PhysRevA.74.062712

    .. [#Raithel2020] K. Moore, A. Duspayev, R. Cardman, and G. Raithel,
             *Phys. Rev. A* **102**, 062817 (2020),
             https://doi.org/10.1103/PhysRevA.102.062817

    .. [#c6] P. Goy, J. Liang, M. Gross, and S. Haroche,\
            *Phys. Rev. A* **34**, 2889 (1986)
            https://doi.org/10.1103/PhysRevA.34.2889

    .. [#jd2016] Johannes Deiglmayr, Holger Herburger, Heiner Sassmannshausen,
        Paul Jansen, Hansjurg Schmutz, Frederic Merkt, *Phys. Rev. A* **93**, 013424 (2016)
        https://doi.org/10.1103/PhysRevA.93.013424

    .. [#Lorenzen1984] C. -J. Lorenzen, and K. Niemax, *Z. Phys. A* **315**,
        127 (1984) https://doi.org/10.1007/BF01419370

    .. [#c7] C. -J. Lorenzen, and K. Niemax, *Physica Scripta* **27**, 300 (1983)

    .. [#c8] NIST, P. Mohr and S. Kotochigova, unpublished calculations (2000).
        The wavelengths for the Balmer-alpha and Balmer-beta transitions at 6563
        and 4861 :math:`\\unicode{xC5}` include only the stronger components of
        more extensive fine structures.

    .. [#c11] R. L. Kelly, *J. Phys. Chem. Ref. Data* **16**, Suppl. 1 (1987).

    .. [#c14] J. S. Coursey, D. J. Schwab, J. J. Tsai, and R. A. Dragoset,
        (2015), Atomic Weights and Isotopic Compositions (version 4.1).
        Online Available: http://physics.nist.gov/Comp (2017, March, 14).
        National Institute of Standards and Technology, Gaithersburg, MD.

    .. [#Sanguinetti2009] B. Sanguinetti, H. O. Majeed, M. L. Jones and
        B. T. H. Varcoe, *J. Phys. B* **42**, 165004 (2009)
        http://iopscience.iop.org/article/10.1088/0953-4075/42/16/165004/meta

    .. [#SteckCs] D. Steck, Cesium D Line Data v2.2.1
        https://steck.us/alkalidata/cesiumnumbers.1.6.pdf

    .. [#Steck85Rb] D. Steck, Rubidium 85 D Line Data v2.2.1
        https://steck.us/alkalidata/rubidium85numbers.pdf

    .. [#Steck87Rb] D. Steck, Rubidium 87 D Line Data v2.2.1
        https://steck.us/alkalidata/rubidium87numbers.pdf

    .. [#Feiertag72] D. Feiertag, A. Sahm, and G. zu Putlitz, Z. Phys. 255, 93 (1972)

    .. [#Arrimondo77] E. Arrimondo, M. Inguscio and P. Violino, Rev. Mod. Phys. 49, 31 (1977)

    Module
    ------
"""

from __future__ import division, print_function, absolute_import

from .alkali_atom_functions import *

from scipy.constants import Rydberg as C_Rydberg
from scipy.constants import m_e as C_m_e


class Hydrogen(AlkaliAtom):
    """
        Properties of hydrogen atoms
    """
    ionisationEnergy = 13.598433  #: (eV), Ref. [#c8]_.
    Z = 1
    scaledRydbergConstant = 109677.5834 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    # NOTE: below data is from NIST Atomic Spectra Database (ASD, ver. 5.5.6)
    # Even more precise Hydrogen level data is available if needed on
    # U.D. Jentschura, S. Kotochigova, E.O. LeBigot, P.J. Mohr, and B.N. Taylor (2005),
    # The Energy Levels of Hydrogen and Deuterium (version 2.1). [Online]
    # Available: http://physics.nist.gov/HDEL [2018, May 3].
    # National Institute of Standards and Technology, Gaithersburg, MD.
    levelDataFromNIST = "h_NIST_level_data.ascii"
    NISTdataLevels = 12

    precalculatedDB = "h_precalculated.db"
    dipoleMatrixElementFile = "h_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "h_quadrupole_matrix_elements.npy"

    groundStateN = 1
    minQuantumDefectN = 8

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 1.00782503223 * physical_constants["atomic mass constant"][0]

    a1 = [0.0, 0.0, 0.0, 0.0]
    a2 = [0.0, 0.0, 0.0, 0.0]
    a3 = [0.0, 0.0, 0.0, 0.0]
    a4 = [0.0, 0.0, 0.0, 0.0]
    rc = [0.0, 0.0, 0.0, 0.0]

    def potential(self, l, s, j, r):
        # Returns total potential that electron feels = core potential + Spin-Orbit interaction
        return -self.Z / r + self.alpha**2 / (2.0 * r**3) * (j * (j + 1.0) - l * (l + 1.0) - s * (s + 1)) / 2.0

    def stateQuantumDefect(self, n, l, j):
        defect = 0.
        return defect


class Caesium(AlkaliAtom):
    """
        Properties of caesium atoms
    """

    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 15.6440
    """
        model potential parameters from [#c1]_

    """
    #
    a1 = [3.49546309, 4.69366096, 4.32466196, 3.01048361]
    """
        model potential parameters from [#c1]_

    """

    a2 = [1.47533800, 1.71398344, 1.61365288, 1.40000001]
    """
        model potential parameters from [#c1]_
    """
    a3 = [-9.72143084, -24.65624280, -6.70128850, -3.20036138]
    """
        model potential parameters from [#c1]_
    """
    a4 = [0.02629242, -0.09543125, -0.74095193, 0.00034538]
    """
        model potential parameters from [#c1]_
    """
    rc = [1.92046930, 2.13383095, 0.93007296, 1.99969677]
    """
        model potential parameters from [#c1]_
    """
    Z = 55

    I = 3.5  # 7/2

    #: (eV), Ref. [#jd2016]_.
    ionisationEnergy = 31406.4677325 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    NISTdataLevels = 25

    # first index [0]:  j-1/2    [1]: j+1/2
    # second index [0..4] : s,p,d,f,g
    # third index [delta0,delta2...]

    quantumDefect = [[[4.04935665, 0.2377037, 0.255401, 0.00378, 0.25486, 0.0],
                      [3.59158950, 0.360926, 0.41905, 0.64388, 1.45035, 0.0],
                      [2.4754562, 0.009320, -0.43498, -0.76358, -18.0061, 0.0],
                      [0.03341424, -0.198674, 0.28953, -0.2601, 0.0, 0.0],
                      [0.00703865, -0.049252, 0.01291, 0.0, 0.0, 0.0]],
                     [[4.04935665, 0.2377037, 0.255401, 0.00378, 0.25486, 0.0],
                      [3.5589599, 0.392469, -0.67431, 22.3531, -92.289, 0.0],
                      [2.46631524, 0.013577, -0.37457, -2.1867, -1.5532, -56.6739],
                      [0.03341424, -0.198674, 0.28953, -0.2601, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects for :math:`S_{1/2}`, :math:`nP_{1/2}`, :math:`D_{5/2}`,
        :math:`F_{5/2}` and :math:`G_{7/2}` are from [#Weber1987]_, while
        quantum defects for :math:`nP_{3/2}`,:math:`D_{3/2}` are from [#Lorenzen1984]_,

        Note:
            f_7/2 quantum defects are PUT TO BE EXACTLY the same as f_5/2 (~10MHz difference?!)
    """

    minQuantumDefectN = 9

    levelDataFromNIST = "cs_NIST_level_data.ascii"

    precalculatedDB = "cs_precalculated.db"

    dipoleMatrixElementFile = "cs_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "cs_quadrupole_matrix_elements.npy"

    literatureDMEfilename = 'caesium_literature_dme.csv'

    #: levels that are for smaller n than ground level, but are above in energy due to angular part
    extraLevels = [[5, 2, 2 + 0.5], [5, 2, 2 - 0.5], [5, 3, 3 + 0.5], [5, 3, 3 - 0.5],
                   [5, 4, 4 + 0.5], [5, 4, 4 - 0.5], [4, 3, 3 + 0.5], [4, 3, 3 - 0.5]]

    groundStateN = 6

    mass = 132.9054519610 * physical_constants["atomic mass constant"][0]
    abundance = 1.000

    #: in eV
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "Cs133"
    meltingPoint = 28.44 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "cs_hfs_data.csv"

    gL = 0.99999587      #: Electron orbital g-factor [#SteckCs]_
    gI = -0.00039885395  #: Nuclear g-factor  [#SteckCs]_


    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Cs in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Cs in liquid phase.

        """

        # returns pressure in Pa for temperature in K

        if temperature < self.meltingPoint:
            # Cs is in solid phase (from table 2. for recommended equations / +-5%)
            return 10.0**(2.881 + 4.711 - 3999. / temperature) * 133.322368

        elif temperature < 550. + 273.15:
            # Cs is in liquid phase (from table 3. of the cited reference "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.232 - 4062. / temperature -
                          1.3359 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: Cs vapour pressure above 550 C is unknown \
                    (limits of experimental interpolation)")
            return 0

    def getPressureOld(self, temperature):
        # returns pressure in Pa for temperature in K

        # from A.N.Nesmeyanov, Vapor Pressure of the Chemical Elements (Elsevier, Amsterdam, 1963). English edition
        # edited by Robert Gary
        # as was found in Steck Alkali metal data, revision 1.6, 14 October 2003

        print("WARNING: getPressureOld is provided just for reference for \
                the old versions of the programme")
        print("New programmes should use getPressure function instead !")

        if temperature < 28.44 + 273.15:
            # Cs is in solid phase
            return 10.0**(-219.482 + 1088.676 / temperature - 0.08336185 * temperature +
                          94.88752 * log(temperature) / log(10.0)) * 133.322368

        elif temperature < 671 + 273.15:
            # Cs is in liquid phase
            return 10.0**(8.22127 - 4006.048 / temperature - 0.00060194 * temperature -
                          0.19623 * log(temperature) / log(10.0)) * 133.322368
        else:
            print("ERROR: Cs vapour pressure above 671 C is unknown \
                    (limits of experimental interpolation)")
            return 0

class Cesium(Caesium):
    """
        support for American English spelling
    """
    pass

class Rubidium85(AlkaliAtom):
    """
        Properites of rubidium 85 atoms
    """

    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 9.0760
    """
        model potential parameters from [#c1]_

    """

    a1 = [3.69628474, 4.44088978, 3.78717363, 2.39848933]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.64915255, 1.92828831, 1.57027864, 1.76810544]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-9.86069196, -16.79597770, -11.65588970, -12.07106780]
    """
        model potential parameters from [#c1]_

    """
    a4 = [0.19579987, -0.8163314, 0.52942835, 0.77256589]
    """
        model potential parameters from [#c1]_

    """
    rc = [1.66242117, 1.50195124, 4.86851938, 4.79831327]
    """
        model potential parameters from [#c1]_

    """
    Z = 37

    I = 2.5  # 5/2

    NISTdataLevels = 77

    #: (eV) Ref. [#Sanguinetti2009]_
    ionisationEnergy = (1010.024700e12) / C_c \
        * physical_constants["inverse meter-electron volt relationship"][0]

    quantumDefect = [[[3.1311804, 0.1784, 0.0, 0.0, 0.0, 0.0],
                      [2.6548849, 0.2900, 0.0, 0.0, 0.0, 0.0],
                      [1.34809171, -0.60286, 0.0, 0.0, 0.0, 0.0],
                      [0.0165192, -0.085, 0.0, 0.0, 0.0, 0.0],
                      [0.0039990, -0.0202, 0.0, 0.0, 0.0, 0.0]],
                     [[3.1311804, 0.1784, 0.0, 0.0, 0.0, 0.0],
                      [2.6416737, 0.2950, 0.0, 0.0, 0.0, 0.0],
                      [1.34646572, -0.59600, 0.0, 0.0, 0.0, 0.0],
                      [0.0165437, -0.086, 0.0, 0.0, 0.0, 0.0],
                      [0.0039990, -0.0202, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects for :math:`nF` states are
        from [#c5]_. Quantum defects for :math:`nG` states are
        from [#Raithel2020]_. All other quantum defects are from from [#c4]_

    """

    levelDataFromNIST = "rb_NIST_level_data.ascii"
    dipoleMatrixElementFile = "rb_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "rb_quadrupole_matrix_elements.npy"

    minQuantumDefectN = 8

    precalculatedDB = "rb85_precalculated.db"

    literatureDMEfilename = 'rubidium_literature_dme.csv'

    #: levels that are for smaller n than ground level, but are above in energy due to angular part
    extraLevels = [[4, 2, 2 + 0.5], [4, 2, 2 - 0.5],
                   [4, 3, 3 + 0.5], [4, 3, 3 - 0.5]]

    groundStateN = 5

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 84.9117897379 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.7217

    #:  in eV
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "Rb85"
    meltingPoint = 39.31 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "rb85_hfs_data.csv"

    gL = 0.99999354      #: Electron orbital g-factor [#Steck85Rb]_
    gI = -0.00029364000  #: Nuclear g-factor  [#Steck85Rb]_

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Rb in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Rb in liquid phase.

        """

        if temperature < self.meltingPoint:
            # Rb is in solid phase (from table 2. for recommended equations / +-5%)
            return 10.0**(2.881 + 4.857 - 4215. / temperature) * 133.322368

        elif temperature < 550. + 273.15:
            # Rb is in liquid phase (from table 3. of the cited reference "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.316 - 4275. / temperature -
                          1.3102 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: Rb vapour pressure above 550 C is unknown \
                    (limits of experimental interpolation)")
            return 0


class Rubidium(Rubidium85):
    """
        backward compatibility:
        before there was only one Rubidium class, and that one corresponded
        to Rubidium85
    """
    pass


class Rubidium87(AlkaliAtom):
    """
        Properites of rubidium 87 atoms
    """

    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 9.0760
    """
        model potential parameters from [#c1]_

    """

    a1 = [3.69628474, 4.44088978, 3.78717363, 2.39848933]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.64915255, 1.92828831, 1.57027864, 1.76810544]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-9.86069196, -16.79597770, -11.65588970, -12.07106780]
    """
        model potential parameters from [#c1]_

    """
    a4 = [0.19579987, -0.8163314, 0.52942835, 0.77256589]
    """
        model potential parameters from [#c1]_

    """
    rc = [1.66242117, 1.50195124, 4.86851938, 4.79831327]
    """
        model potential parameters from [#c1]_

    """
    Z = 37

    I = 1.5  # 3/2

    NISTdataLevels = 77

    #: (eV) Ref. [#Mack2011]_
    ionisationEnergy = (1010.0291646e12 - 4.271676631815181e9) / C_c \
        * physical_constants["inverse meter-electron volt relationship"][0]

    quantumDefect = [[[3.1311804, 0.1784, 0.0, 0.0, 0.0, 0.0],
                      [2.6548849, 0.2900, 0.0, 0.0, 0.0, 0.0],
                      [1.34809171, -0.60286, 0.0, 0.0, 0.0, 0.0],
                      [0.0165192, -0.085, 0.0, 0.0, 0.0, 0.0],
                      [0.00405, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[3.1311804, 0.1784, 0.0, 0.0, 0.0, 0.0],
                      [2.6416737, 0.2950, 0.0, 0.0, 0.0, 0.0],
                      [1.34646572, -0.59600, 0.0, 0.0, 0.0, 0.0],
                      [0.0165437, -0.086, 0.0, 0.0, 0.0, 0.0],
                      [0.00405, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects for :math:`nF` states are
        from [#c5]_. Quantum defects for :math:`nG` states are
        from [#Afrousheh2006a]_. All other quantum defects are from from [#c4]_

    """

    levelDataFromNIST = "rb_NIST_level_data.ascii"
    dipoleMatrixElementFile = "rb_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "rb_quadrupole_matrix_elements.npy"

    minQuantumDefectN = 8

    precalculatedDB = "rb87_precalculated.db"

    literatureDMEfilename = 'rubidium_literature_dme.csv'

    #: levels that are for smaller n than ground level, but are above in energy due to angular part
    extraLevels = [[4, 2, 2 + 0.5], [4, 2, 2 - 0.5],
                   [4, 3, 3 + 0.5], [4, 3, 3 - 0.5]]

    groundStateN = 5

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 86.9091805310 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.2783

    #:  in eV  (M_ion core = m_atomic - m_electron)
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]
    #109736.605*1.e2 \
#        *physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "Rb87"
    meltingPoint = 39.31 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "rb87_hfs_data.csv"

    gL = 0.99999369      #: Electron orbital g-factor [#Steck87Rb]_
    gI = -0.0009951414   #: Nuclear g-factor  [#Steck87Rb]_

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Rb in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Rb in liquid phase.

        """

        if temperature < self.meltingPoint:
            # Rb is in solid phase (from table 2. for recommended equations / +-5%)
            return 10.0**(2.881 + 4.857 - 4215. / temperature) * 133.322368

        elif temperature < 550. + 273.15:
            # Rb is in liquid phase (from table 3. of the cited reference "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.316 - 4275. / temperature -
                          1.3102 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: Rb vapour pressure above 550 C is unknown \
                    (limits of experimental interpolation)")
            return 0


class Lithium6(AlkaliAtom):  # Li
    """
        Properties of lithium 6 atoms
    """

    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 0.1923
    """
        model potential parameters from [#c1]_

    """
    # model potential parameters from Marinescu et.al, PRA 49:982 (1994)
    a1 = [2.47718079, 3.45414648, 2.51909839, 2.51909839]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.84150932, 2.55151080, 2.43712450, 2.43712450]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-0.02169712, -0.21646561, 0.32505524, 0.32505524]
    """
        model potential parameters from [#c1]_

    """
    a4 = [-0.11988362, -0.06990078, 0.10602430, 0.10602430]
    """
        model potential parameters from [#c1]_

    """
    rc = [0.61340824, 0.61566441, 2.34126273, 2.34126273]
    """
        model potential parameters from [#c1]_

    """

    Z = 3

    I = 1.

    NISTdataLevels = 42

    # (eV) from Ref. [#c7]_
    ionisationEnergy = 43487.15 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    # PRA 34, 2889 (1986); and (for D_J and F_J) from Physica Scripta 27:300-305 (1983)
    quantumDefect = [[[0.3995101, 0.0290, 0.0, 0.0, 0.0, 0.0],
                      [0.0471835, -0.024, 0.0, 0.0, 0.0, 0.0],
                      [0.002129, -0.01491, 0.1759, -0.8507, 0.0, 0.0],
                      [-0.000077, 0.021856, -0.4211, 2.3891, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.3995101, 0.0290, 0.0, 0.0, 0.0, 0.0],
                      [0.0471720, -0.024, 0.0, 0.0, 0.0, 0.0],
                      [0.002129, -0.01491, 0.1759, -0.8507, 0.0, 0.0],
                      [-0.000077, 0.021856, -0.4211, 2.3891, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects for :math:`nS` and :math:`nP` are from Ref. [#c6]_ .
        Quantum defects for :math:`D_j` and :math:`F_j` are from Ref. [#c7]_
        (note that this defects in Ref. [#c7]_ are for Li7, differences
        are expected not be too big).

    """

    levelDataFromNIST = "li_NIST_level_data.ascii"
    dipoleMatrixElementFile = "li6_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "li6_quadrupole_matrix_elements.npy"

    minQuantumDefectN = 4

    precalculatedDB = "li6_precalculated.db"

    # levels that are for smaller n than ground level, but are above in energy
    # due to angular part
    extraLevels = []

    groundStateN = 2

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 6.0151228874 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.0759

    gL = 1 - physical_constants["electron mass"][0]/mass

    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "Li6"
    meltingPoint = 180.54 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "li6_hfs_data.csv"

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 3.
            (accuracy +-1 %) are used both for liquid and solid phase of Li.

        """

        if temperature < self.meltingPoint:
            # Li is in solid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 7.790 - 8423. / temperature -
                          0.7074 * log(temperature) / log(10.)) * 133.322368

        elif temperature < 1000. + 273.15:
            # Li is in liquid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.409 - 8320. / temperature -
                          1.0255 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: Li vapour pressure above 1000 C is unknown \
                    (limits of experimental interpolation)")
            return 0


class Lithium7(AlkaliAtom):  # Li
    """
        Properties of lithium 7 atoms
    """
    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    # model potential parameters from Marinescu et.al, PRA 49:982 (1994)
    alphaC = 0.1923
    """
        model potential parameters from [#c1]_

    """
    a1 = [2.47718079, 3.45414648, 2.51909839, 2.51909839]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.84150932, 2.55151080, 2.43712450, 2.43712450]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-0.02169712, -0.21646561, 0.32505524, 0.32505524]
    """
        model potential parameters from [#c1]_

    """
    a4 = [-0.11988362, -0.06990078, 0.10602430, 0.10602430]
    """
        model potential parameters from [#c1]_

    """
    rc = [0.61340824, 0.61566441, 2.34126273, 2.34126273]
    """
        model potential parameters from [#c1]_

    """

    Z = 3

    I = 1.5  # 3/2

    NISTdataLevels = 42
    ionisationEnergy = 5.391719  #: (eV) NIST Ref. [#c11]_.

    quantumDefect = [[[0.3995101, 0.0290, 0.0, 0.0, 0.0, 0.0],
                      [0.0471780, -0.024, 0.0, 0.0, 0.0, 0.0],
                      [0.002129, -0.01491, 0.1759, -0.8507, 0.0, 0.0],
                      [-0.000077, 0.021856, -0.4211, 2.3891, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.3995101, 0.0290, 0.0, 0.0, 0.0, 0.0],
                      [0.0471665, -0.024, 0.0, 0.0, 0.0, 0.0],
                      [0.002129, -0.01491, 0.1759, -0.8507, 0.0, 0.0],
                      [-0.000077, 0.021856, -0.4211, 2.3891, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects for :math:`nS` and :math:`nP` states are
        from Ref. [#c6]_. Quantum defects for :math:`D_j` and :math:`F_j`
        states are from [#c7]_.

    """

    levelDataFromNIST = "li_NIST_level_data.ascii"
    dipoleMatrixElementFile = "li7_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "li7_quadrupole_matrix_elements.npy"

    minQuantumDefectN = 4

    precalculatedDB = "li7_precalculated.db"

    # levels that are for smaller n than ground level,
    # but are above in energy due to angular part
    extraLevels = []

    groundStateN = 2

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 7.0160034366 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.9241

    gL = 1 - physical_constants["electron mass"][0]/mass

    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "Li7"
    meltingPoint = 180.54 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "li7_hfs_data.csv"

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature (in K).

            Uses equation and values from [#c3]_. Values from table 3.
            (accuracy +-1 %) are used for both liquid and solid phase of Li.

        """

        if temperature < self.meltingPoint:
            # Li is in solid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 7.790 - 8423. / temperature -
                          0.7074 * log(temperature) / log(10.)) * 133.322368

        elif temperature < 1000. + 273.15:
            # Li is in liquid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.409 - 8320. / temperature -
                          1.0255 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: Li vapour pressure above 1000 C is unknown \
                    (limits of experimental interpolation)")
            return 0


class Sodium(AlkaliAtom):  # Na23
    """
        Properties of sodium 23 atoms
    """

    #: ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 0.9448
    """
        model potential parameters from [#c1]_

    """
    a1 = [4.82223117, 5.08382502, 3.53324124, 1.11056646]
    """
        model potential parameters from [#c1]_

    """
    a2 = [2.45449865, 2.18226881, 2.48697936, 1.05458759]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-1.12255048, -1.19534623, -0.75688448, 1.73203428]
    """
        model potential parameters from [#c1]_

    """
    a4 = [-1.42631393, -1.03142861, -1.27852357, -0.09265696]
    """
        model potential parameters from [#c1]_

    """
    rc = [0.45489422, 0.45798739, 0.71875312, 28.6735059]
    """
        model potential parameters from [#c1]_

    """

    Z = 11

    I = 1.5  # 3/2

    NISTdataLevels = 20

    #: (eV) from Ref. [#c7]_
    ionisationEnergy = 41449.44 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    quantumDefect = [[[1.347964, 0.060673, 0.0233, -0.0085, 0.0, 0.0],
                      [0.855380, 0.11363, 0.0384, 0.1412, 0.0, 0.0],
                      [0.015543, -0.08535, 0.7958, -4.0513, 0.0, 0.0],
                      [0.001453, 0.017312, -0.7809, 7.021, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[1.347964, 0.060673, 0.0233, -0.0085, 0.0, 0.0],
                      [0.854565, 0.114195, 0.0352, 0.1533, 0.0, 0.0],
                      [0.015543, -0.08535, 0.7958, -4.0513, 0.0, 0.0],
                      [0.001453, 0.017312, -0.7809, 7.021, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        Quantum defects are from Ref. [#c7]_. Note that we are using modified
        Rydberg-Ritz formula. In literature both modified and non-modified
        coefficients appear. For more details about the two equations see
        page 301. of Ref. [#c7]_.
    """

    levelDataFromNIST = "na_NIST_level_data.ascii"
    dipoleMatrixElementFile = "na23_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "na23_quadrupole_matrix_elements.npy"

    precalculatedDB = "na23_precalculated.db"

    literatureDMEfilename = 'sodium_literature_dme.csv'

    # levels that are for smaller n than ground level, but are above in
    # energy due to angular part
    extraLevels = []

    groundStateN = 3

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 22.9897692820 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 1.00

    gL = 1 - physical_constants["electron mass"][0]/mass

    #: (eV)
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "Na23"
    meltingPoint = 97.794 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "na23_hfs_data.csv"

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Na in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Na in liquid phase.

        """

        if temperature < self.meltingPoint:
            # Na is in solid phase (from table 2. of the cited reference  / +- 5%)
            return 10.0**(2.881 + 5.298 - 5603. / temperature) * 133.322368

        elif temperature < 700. + 273.15:
            # Na is in liquid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.400 - 5634. / temperature -
                          1.1748 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: Na vapour pressure above 700 C is unknown \
                    (limits of experimental interpolation)")
            return 0


class Potassium39(AlkaliAtom):
    """
        Properties of potassium 39 atoms
    """
    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 5.3310
    """
        model potential parameters from [#c1]_

    """
    a1 = [3.56079437, 3.65670429, 4.12713694, 1.42310446]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.83909642, 1.67520788, 1.79837462, 1.27861156]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-1.74701102, -2.07416615, -1.69935174, 4.77441476]
    """
        model potential parameters from [#c1]_

    """
    a4 = [-1.03237313, -0.89030421, -0.98913582, -0.94829262]
    """
        model potential parameters from [#c1]_

    """
    rc = [0.83167545, 0.85235381, 0.83216907, 6.50294371]
    """
        model potential parameters from [#c1]_

    """

    Z = 19

    I = 1.5  # 3/2

    NISTdataLevels = 46

    #: (eV), weighted average of values in Ref. [#c7]_.
    ionisationEnergy = 35009.8139375 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    # quantum defects from Physica Scripta 27:300 (1983)
    quantumDefect = [[[2.1801985, 0.13558, 0.0759, 0.117, -0.206, 0.0],
                      [1.713892, 0.233294, 0.16137, 0.5345, -0.234, 0.0],
                      [0.27697, -1.024911, -0.709174, 11.839, -26.689, 0.0],
                      [0.010098, -0.100224, 1.56334, -12.6851, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[2.1801985, 0.13558, 0.0759, 0.117, -0.206, 0.0],
                      [1.710848, 0.235437, 0.11551, 1.1015, -2.0356, 0.0],
                      [0.2771580, -1.025635, -0.59201, 10.0053, -19.0244, 0.0],
                      [0.010098, -0.100224, 1.56334, -12.6851, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects from Ref. [#c7]_.
    """

    levelDataFromNIST = "k_NIST_level_data.ascii"
    dipoleMatrixElementFile = "k_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "k_quadrupole_matrix_elements.npy"

    precalculatedDB = "k39_precalculated.db"

    literatureDMEfilename = 'potassium_literature_dme.csv'

    #: levels that are for smaller n than ground level, but are above in energy due to angular part
    extraLevels = [[3, 2, 2 + 0.5], [3, 2, 2 - 0.5]]

    groundStateN = 4

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 38.9637064864 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.932581

    gL = 1 - physical_constants["electron mass"][0]/mass

    # in eV
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "K39"
    meltingPoint = 63.5 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "k39_hfs_data.csv"

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Na in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Na in liquid phase.

        """

        if temperature < self.meltingPoint:
            # K is in solid phase (from table 2. of the cited reference  / +- 5%)
            return 10.0**(2.881 + 4.961 - 4646. / temperature) * 133.322368

        elif temperature < 600. + 273.15:
            # K is in liquid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.233 - 4693. / temperature -
                          1.2403 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: K vapour pressure above 600 C is unknown \
                (limits of experimental interpolation)")
            return 0


class Potassium(Potassium39):
    """
      backward compatibility:
      before only one class for Potassium existed and
      it corresponded to Potassium 39
    """
    pass


class Potassium40(AlkaliAtom):
    """
        Properties of potassium 40 atoms
    """
    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 5.3310
    """
        model potential parameters from [#c1]_

    """
    a1 = [3.56079437, 3.65670429, 4.12713694, 1.42310446]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.83909642, 1.67520788, 1.79837462, 1.27861156]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-1.74701102, -2.07416615, -1.69935174, 4.77441476]
    """
        model potential parameters from [#c1]_

    """
    a4 = [-1.03237313, -0.89030421, -0.98913582, -0.94829262]
    """
        model potential parameters from [#c1]_

    """
    rc = [0.83167545, 0.85235381, 0.83216907, 6.50294371]
    """
        model potential parameters from [#c1]_

    """

    Z = 19

    I = 4

    NISTdataLevels = 46

    #: (eV), weighted average of values in Ref. [#c7]_.
    ionisationEnergy = 35009.8139375 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    # quantum defects from Physica Scripta 27:300 (1983)
    quantumDefect = [[[2.1801985, 0.13558, 0.0759, 0.117, -0.206, 0.0],
                      [1.713892, 0.233294, 0.16137, 0.5345, -0.234, 0.0],
                      [0.27697, -1.024911, -0.709174, 11.839, -26.689, 0.0],
                      [0.010098, -0.100224, 1.56334, -12.6851, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[2.1801985, 0.13558, 0.0759, 0.117, -0.206, 0.0],
                      [1.710848, 0.235437, 0.11551, 1.1015, -2.0356, 0.0],
                      [0.2771580, -1.025635, -0.59201, 10.0053, -19.0244, 0.0],
                      [0.010098, -0.100224, 1.56334, -12.6851, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects from Ref. [#c7]_.
    """

    levelDataFromNIST = "k_NIST_level_data.ascii"
    dipoleMatrixElementFile = "k_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "k_quadrupole_matrix_elements.npy"

    precalculatedDB = "k40_precalculated.db"

    literatureDMEfilename = 'potassium_literature_dme.csv'

    #: levels that are for smaller n than ground level, but are above in energy due to angular part
    extraLevels = [[3, 2, 2 + 0.5], [3, 2, 2 - 0.5]]

    groundStateN = 4

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 39.963998166 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.000117

    gL = 1 - physical_constants["electron mass"][0]/mass

    #: in eV
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "K40"
    meltingPoint = 63.5 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "k40_hfs_data.csv"

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Na in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Na in liquid phase.

        """

        if temperature < self.meltingPoint:
            # K is in solid phase (from table 2. of the cited reference  / +- 5%)
            return 10.0**(2.881 + 4.961 - 4646. / temperature) * 133.322368

        elif temperature < 600. + 273.15:
            # K is in liquid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.233 - 4693. / temperature -
                          1.2403 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: K vapour pressure above 600 C is unknown \
                (limits of experimental interpolation)")
            return 0


class Potassium41(AlkaliAtom):
    """
        Properties of potassium 41 atoms
    """
    # ALL PARAMETERES ARE IN ATOMIC UNITS (HATREE)
    alphaC = 5.3310
    """
        model potential parameters from [#c1]_

    """
    a1 = [3.56079437, 3.65670429, 4.12713694, 1.42310446]
    """
        model potential parameters from [#c1]_

    """
    a2 = [1.83909642, 1.67520788, 1.79837462, 1.27861156]
    """
        model potential parameters from [#c1]_

    """
    a3 = [-1.74701102, -2.07416615, -1.69935174, 4.77441476]
    """
        model potential parameters from [#c1]_

    """
    a4 = [-1.03237313, -0.89030421, -0.98913582, -0.94829262]
    """
        model potential parameters from [#c1]_

    """
    rc = [0.83167545, 0.85235381, 0.83216907, 6.50294371]
    """
        model potential parameters from [#c1]_

    """

    Z = 19

    I = 1.5  # 3/2

    NISTdataLevels = 46

    #: (eV), weighted average of values in Ref. [#c7]_.
    ionisationEnergy = 35009.8139375 * 1.e2\
        * physical_constants["inverse meter-electron volt relationship"][0]

    # quantum defects from Physica Scripta 27:300 (1983)
    quantumDefect = [[[2.1801985, 0.13558, 0.0759, 0.117, -0.206, 0.0],
                      [1.713892, 0.233294, 0.16137, 0.5345, -0.234, 0.0],
                      [0.27697, -1.024911, -0.709174, 11.839, -26.689, 0.0],
                      [0.010098, -0.100224, 1.56334, -12.6851, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[2.1801985, 0.13558, 0.0759, 0.117, -0.206, 0.0],
                      [1.710848, 0.235437, 0.11551, 1.1015, -2.0356, 0.0],
                      [0.2771580, -1.025635, -0.59201, 10.0053, -19.0244, 0.0],
                      [0.010098, -0.100224, 1.56334, -12.6851, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]
    """
        quantum defects from Ref. [#c7]_.
    """

    levelDataFromNIST = "k_NIST_level_data.ascii"
    dipoleMatrixElementFile = "k_dipole_matrix_elements.npy"
    quadrupoleMatrixElementFile = "k_quadrupole_matrix_elements.npy"

    precalculatedDB = "k41_precalculated.db"

    literatureDMEfilename = 'potassium_literature_dme.csv'

    #: levels that are for smaller n than ground level, but are above in energy due to angular part
    extraLevels = [[3, 2, 2 + 0.5], [3, 2, 2 - 0.5]]

    groundStateN = 4

    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    mass = 40.9618252579 * physical_constants["atomic mass constant"][0]
    #: source NIST, Atomic Weights and Isotopic Compositions [#c14]_
    abundance = 0.067302

    gL = 1 - physical_constants["electron mass"][0]/mass

    #: in eV
    scaledRydbergConstant = (mass - C_m_e) / (mass) * C_Rydberg\
        * physical_constants["inverse meter-electron volt relationship"][0]

    elementName = "K41"
    meltingPoint = 63.5 + 273.15  #: in K

    #: source of HFS magnetic dipole and quadrupole constants
    hyperfineStructureData = "k41_hfs_data.csv"

    def getPressure(self, temperature):
        """
            Pressure of atomic vapour at given temperature.

            Uses equation and values from [#c3]_. Values from table 2.
            (accuracy +- 5%) are used for Na in solid phase. Values from table 3.
            (accuracy +-1 %) are used for Na in liquid phase.

        """

        if temperature < self.meltingPoint:
            # K is in solid phase (from table 2. of the cited reference  / +- 5%)
            return 10.0**(2.881 + 4.961 - 4646. / temperature) * 133.322368

        elif temperature < 600. + 273.15:
            # K is in liquid phase (from table 3. of the cited reference
            # "precisely fitted equations / +- 1%)
            return 10.0**(2.881 + 8.233 - 4693. / temperature -
                          1.2403 * log(temperature) / log(10.)) * 133.322368
        else:
            print("ERROR: K vapour pressure above 600 C is unknown \
                (limits of experimental interpolation)")
            return 0
