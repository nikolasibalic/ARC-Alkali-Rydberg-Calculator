# -*- coding: utf-8 -*-
"""
Implements general single-atom calculations

This module calculates single (isolated) atom properties of all alkali metals in
general. For example, it calculates dipole matrix elements, quandrupole matrix
elements, etc.  Also, some helpful general functions are here, e.g. for saving
and loading calculations (single-atom and pair-state based), printing state
labels etc.


"""

from __future__ import division, print_function, absolute_import

import sqlite3
import csv
import gzip
from math import log, exp, sqrt
from mpmath import angerj
# for web-server execution, uncomment the following two lines
# import matplotlib
# matplotlib.use("Agg")
import numpy as np
import re
import shutil

from numpy.linalg import eigh

from .wigner import Wigner6j, Wigner3j, CG, WignerDmatrix
from scipy.constants import physical_constants, pi, epsilon_0, hbar
from scipy.constants import k as C_k
from scipy.constants import c as C_c
from scipy.constants import h as C_h
from scipy.constants import e as C_e
from scipy.constants import m_e as C_m_e

# for matrices
from scipy.sparse import csr_matrix
from numpy import floor

import sys
import os
if sys.version_info > (2,):
    xrange = range

import pickle
sqlite3.register_adapter(np.float64, float)
sqlite3.register_adapter(np.float32, float)
sqlite3.register_adapter(np.int64, int)
sqlite3.register_adapter(np.int32, int)

DPATH = os.path.join(os.path.expanduser('~'), '.arc-data')
__arc_data_version__ = 10


def setup_data_folder():
    """ Setup the data folder in the users home directory.

    """
    if not os.path.exists(DPATH):
        os.makedirs(DPATH)

    # check what is the local version of data
    copyDataLocally = True
    versionFile = os.path.join(DPATH, "version.txt")
    if os.path.exists(versionFile):
        with open(versionFile, "r") as f:
            version = int(f.readline())
        if (version == __arc_data_version__):
            copyDataLocally = False

    if copyDataLocally:
        dataFolder = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), "data")
        for fn in os.listdir(dataFolder):
            if os.path.isfile(os.path.join(dataFolder, fn)):
                shutil.copy(os.path.join(dataFolder, fn), DPATH)

        dataFolder = os.path.join(dataFolder, "refractive_index_data")
        refractiveIndexData = os.path.join(DPATH, "refractive_index_data")

        if not os.path.exists(refractiveIndexData):
            os.makedirs(refractiveIndexData)

        for fn in os.listdir(dataFolder):
            if os.path.isfile(os.path.join(dataFolder, fn)):
                shutil.copy(os.path.join(dataFolder, fn), refractiveIndexData)

        with open(versionFile, "w") as f:
            f.write("%d" % __arc_data_version__)


class AlkaliAtom(object):
    """
        Implements general calculations for alkali atoms.

        This abstract class implements general calculations methods.

        Args:
            preferQuantumDefects (bool):
                Use quantum defects for energy level calculations. If False,
                uses NIST ASD values where available. If True, uses quantum
                defects for energy calculations for principal quantum numbers
                equal or above :obj:`minQuantumDefectN` which is specified for
                each element separately. For principal quantum numbers below
                this value, NIST ASD values are used, since quantum defects
                don't reproduce well low-lying states. Default is True.
            cpp_numerov (bool):
                should the wavefunction be calculated with Numerov algorithm
                implemented in C++; if False, it uses pure Python
                implementation that is much slower. Default is True.

    """

    gS = 2.0023193043737  # : Electron Spin g-factor [Steck]
    gL = 1.0          #: Electron Orbital g-factor
    gI = 0.0          #: Nuclear g-factor

    # ALL PARAMETERS ARE IN ATOMIC UNITS (Hatree)
    alpha = physical_constants["fine-structure constant"][0]

    #: Model potential parameters fitted from experimental observations for
    #: different l (electron angular momentum)
    a1, a2, a3, a4, rc = [0], [0], [0], [0], [0]

    alphaC = 0.0    #: Core polarizability
    Z = 0.0       #: Atomic number
    I = 0.0       #: Nuclear spin

    #: state energies from NIST values
    #: sEnergy [n,l] = state energy for n, l, j = l-1/2
    #: sEnergy [l,n] = state energy for j = l+1/2
    sEnergy = 0
    NISTdataLevels = 0
    scaledRydbergConstant = 0  # : in eV

    #: Contains list of modified Rydberg-Ritz coefficients for calculating
    #: quantum defects for [[ :math:`S_{1/2},P_{1/2},D_{3/2},F_{5/2}`],
    #: [ :math:`S_{1/2},P_{3/2},D_{5/2},F_{7/2}`]]."""
    quantumDefect = [[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]

    #: location of stored NIST values of measured energy levels in eV
    levelDataFromNIST = ""
    #: location of hard-disk stored dipole matrix elements
    dipoleMatrixElementFile = ""
    #: location of hard-disk stored dipole matrix elements
    quadrupoleMatrixElementFile = ""

    dataFolder = DPATH

    # now additional literature sources of dipole matrix elements

    #: Filename of the additional literature source values of dipole matrix
    #: elements.
    #: These additional values should be saved as reduced dipole matrix
    #: elements in J basis.
    literatureDMEfilename = ""

    #: levels that are for smaller principal quantum number (n) than ground
    #: level, but are above in energy due to angular part
    extraLevels = []

    #: principal quantum number for the ground state
    groundStateN = 0

    #: swich - should the wavefunction be calculated with Numerov algorithm
    #: implemented in C++
    cpp_numerov = True

    mass = 0.  #: atomic mass in kg
    abundance = 1.0  #: relative isotope abundance

    elementName = "elementName"  #: Human-readable element name
    meltingPoint = 0  #: melting point of the element at standard conditions

    preferQuantumDefects = False

    #: minimal quantum number for which quantum defects can be used;
    #: uses measured energy levels otherwise
    minQuantumDefectN = 0

    #: file cotaining data on hyperfine structure (magnetic dipole A and
    #: magnetic quadrupole B constnats).
    hyperfineStructureData = ""

    def __init__(self, preferQuantumDefects=True, cpp_numerov=True):
        # should the wavefunction be calculated with Numerov algorithm
        # implemented in C; if false, it uses Python implementation
        # that is much slower

        self.cpp_numerov = cpp_numerov
        self.preferQuantumDefects = preferQuantumDefects

        self._databaseInit()
        c = self.conn.cursor()

        if self.cpp_numerov:
            from .arc_c_extensions import NumerovWavefunction
            self.NumerovWavefunction = NumerovWavefunction

        # load dipole matrix elements previously calculated
        data = []
        if (self.dipoleMatrixElementFile != ""):
            if preferQuantumDefects is False:
                self.dipoleMatrixElementFile = \
                    "NIST_" + self.dipoleMatrixElementFile
            try:
                data = np.load(os.path.join(self.dataFolder,
                                            self.dipoleMatrixElementFile),
                               encoding='latin1', allow_pickle=True)
            except IOError as e:
                print("Error reading dipoleMatrixElement File "
                      + os.path.join(self.dataFolder, self.dipoleMatrixElementFile))
                print(e)
        # save to SQLite database
        try:
            c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='dipoleME';''')
            if (c.fetchone()[0] == 0):
                # create table
                c.execute('''CREATE TABLE IF NOT EXISTS dipoleME
                 (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED,
                 j1_x2 TINYINT UNSIGNED,
                 n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED,
                 j2_x2 TINYINT UNSIGNED,
                 dme DOUBLE,
                 PRIMARY KEY (n1,l1,j1_x2,n2,l2,j2_x2)
                ) ''')
                if (len(data) > 0):
                    c.executemany(
                        'INSERT INTO dipoleME VALUES (?,?,?,?,?,?,?)', data)
                self.conn.commit()
        except sqlite3.Error as e:
            print("Error while loading precalculated values into the database")
            print(e)
            exit()

        # load quadrupole matrix elements previously calculated
        data = []
        if (self.quadrupoleMatrixElementFile != ""):
            if preferQuantumDefects is False:
                self.quadrupoleMatrixElementFile = \
                    "NIST_" + self.quadrupoleMatrixElementFile
            try:
                data = np.load(os.path.join(self.dataFolder,
                                            self.quadrupoleMatrixElementFile),
                               encoding='latin1', allow_pickle=True)

            except IOError as e:
                print("Error reading quadrupoleMatrixElementFile File "
                      + os.path.join(self.dataFolder,
                                     self.quadrupoleMatrixElementFile))
                print(e)
        # save to SQLite database
        try:
            c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='quadrupoleME';''')
            if (c.fetchone()[0] == 0):
                # create table
                c.execute('''CREATE TABLE IF NOT EXISTS quadrupoleME
                 (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED,
                 j1_x2 TINYINT UNSIGNED,
                 n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED,
                 j2_x2 TINYINT UNSIGNED,
                 qme DOUBLE,
                 PRIMARY KEY (n1,l1,j1_x2,n2,l2,j2_x2)
                ) ''')
                if (len(data) > 0):
                    c.executemany(
                        'INSERT INTO quadrupoleME VALUES (?,?,?,?,?,?,?)',
                        data)
                self.conn.commit()
        except sqlite3.Error as e:
            print("Error while loading precalculated values into the database")
            print(e)
            exit()

        self.sEnergy = np.array([[0.0] * (self.NISTdataLevels + 1)]
                                * (self.NISTdataLevels + 1))

        # Always load NIST data on measured energy levels;
        # Even when user wants to use quantum defects, qunatum defects for
        # lowest lying state are not always so accurate, so below the
        # minQuantumDefectN cut-off (defined for each element separately)
        # getEnergy(...) will always return measured,
        # not calculated energy levels
        if (self.levelDataFromNIST == ""):
            print(
                "NIST level data file not specified."
                + "Only quantum defects will be used.")
        else:
            levels = self._parseLevelsFromNIST(
                os.path.join(self.dataFolder, self.levelDataFromNIST))
            br = 0

            while br < len(levels):
                self._addEnergy(levels[br][0], levels[br]
                                [1], levels[br][2], levels[br][3])
                br = br + 1

        # read Literature values for dipole matrix elements
        self._readLiteratureValues()
        self._readHFSdata()

        return

    def _databaseInit(self):
        # SQL connection and cursor
        self.conn = sqlite3.connect(os.path.join(self.dataFolder,
                                                 self.precalculatedDB))

    def getPressure(self, temperature):
        """ Vapour pressure (in Pa) at given temperature

            Args:
                temperature (float): temperature in K
            Returns:
                float: vapour pressure in Pa
        """
        print("Error: getPressure to-be implement in child class (otherwise "
              + "this call is invalid for the specified atom")
        exit()

    def getNumberDensity(self, temperature):
        """ Atom number density at given temperature

            See `calculation of basic properties example snippet`_.

            .. _`calculation of basic properties example snippet`:
                ./Rydberg_atoms_a_primer.html#General-atomic-properties

            Args:
                temperature (float): temperature in K
            Returns:
                float: atom concentration in :math:`1/m^3`
        """
        return self.getPressure(temperature) / (C_k * temperature)

    def getAverageInteratomicSpacing(self, temperature):
        """
            Returns average interatomic spacing in atomic vapour

            See `calculation of basic properties example snippet`_.

            .. _`calculation of basic properties example snippet`:
                ./Rydberg_atoms_a_primer.html#General-atomic-properties

            Args:
                temperature (float): temperature of the atomic vapour

            Returns:
                float: average interatomic spacing in m
        """
        return (5. / 9.) * self.getNumberDensity(temperature)**(-1. / 3.)

    def corePotential(self, l, r):
        """ core potential felt by valence electron

            For more details about derivation of model potential see
            Ref. [#marinescu]_.

            Args:
                l (int): orbital angular momentum
                r (float): distance from the nucleus (in a.u.)
            Returns:
                float: core potential felt by valence electron (in a.u. ???)

            References:

                .. [#marinescu] M. Marinescu, H. R. Sadeghpour, and A. Dalgarno
                    PRA **49**, 982 (1994),
                    https://doi.org/10.1103/PhysRevA.49.982
        """

        return -self.effectiveCharge(l, r) / r - self.alphaC / (2 * r**4) * \
            (1 - exp(-(r / self.rc[l])**6))

    def effectiveCharge(self, l, r):
        """ effective charge of the core felt by valence electron

            For more details about derivation of model potential see
            Ref. [#marinescu]_.

            Args:
                l (int): orbital angular momentum
                r (float): distance from the nucleus (in a.u.)
            Returns:
                float: effective charge (in a.u.)
         """
        return 1.0 + (self.Z - 1) * exp(-self.a1[l] * r) - \
            r * (self.a3[l] + self.a4[l] * r) * exp(-self.a2[l] * r)

    def potential(self, l, s, j, r):
        """ returns total potential that electron feels

            Total potential = core potential + Spin-Orbit interaction

            Args:
                l (int): orbital angular momentum
                s (float): spin angular momentum
                j (float): total angular momentum
                r (float): distance from the nucleus (in a.u.)
            Returns:
                float: potential (in a.u.)
        """
        if l < 4:
            return self.corePotential(l, r) + self.alpha**2 / (2.0 * r**3) * \
                (j * (j + 1.0) - l * (l + 1.0) - s * (s + 1)) / 2.0
        else:
            # act as if it is a Hydrogen atom
            return -1. / r + self.alpha**2 / (2.0 * r**3) * \
                (j * (j + 1.0) - l * (l + 1.0) - s * (s + 1)) / 2.0

    def radialWavefunction(self, l, s, j, stateEnergy,
                           innerLimit, outerLimit, step):
        """
        Radial part of electron wavefunction

        Calculates radial function with Numerov (from outside towards the
        core). Note that wavefunction might not be calculated all the way to
        the requested `innerLimit` if the divergence occurs before. In that
        case third returned argument gives nonzero value, corresponding to the
        first index in the array for which wavefunction was calculated. For
        quick example see `Rydberg wavefunction calculation snippet`_.

        .. _`Rydberg wavefunction calculation snippet`:
            ./Rydberg_atoms_a_primer.html#Rydberg-atom-wavefunctions



        Args:
            l (int): orbital angular momentum
            s (float): spin angular momentum
            j (float): total angular momentum
            stateEnergy (float): state energy, relative to ionization
                threshold, should be given in atomic units (Hatree)
            innerLimit (float): inner limit at which wavefunction is requested
            outerLimit (float): outer limit at which wavefunction is requested
            step (flaot): radial step for integration mesh (a.u.)
        Returns:
            List[float], List[flaot], int:
                :math:`r`

                :math:`R(r)\\cdot r`

        .. note::
            Radial wavefunction is not scaled to unity! This normalization
            condition means that we are using spherical harmonics which are
            normalized such that
            :math:`\\int \\mathrm{d}\\theta~\\mathrm{d}\\psi~Y(l,m_l)^* \
            \\times Y(l',m_{l'})  =  \\delta (l,l') ~\\delta (m_l, m_{l'})`.

        Note:
            Alternative calculation methods can be added here (potenatial
            package expansion).

        """
        innerLimit = max(
            4. * step, innerLimit)  # prevent divergence due to hitting 0
        if self.cpp_numerov:
            # efficiant implementation in C
            if (l < 4):
                d = self.NumerovWavefunction(
                    innerLimit, outerLimit,
                    step, 0.01, 0.01,
                    l, s, j, stateEnergy, self.alphaC, self.alpha,
                    self.Z,
                    self.a1[l], self.a2[l], self.a3[l], self.a4[l],
                    self.rc[l],
                    (self.mass - C_m_e) / self.mass)
            else:
                d = self.NumerovWavefunction(
                    innerLimit, outerLimit,
                    step, 0.01, 0.01,
                    l, s, j, stateEnergy, self.alphaC, self.alpha,
                    self.Z, 0., 0., 0., 0., 0.,
                    (self.mass - C_m_e) / self.mass)

            psi_r = d[0]
            r = d[1]
            suma = np.trapz(psi_r**2, x=r)
            psi_r = psi_r / (sqrt(suma))
        else:
            # full implementation in Python
            mu = (self.mass - C_m_e) / self.mass

            def potential(x):
                r = x * x
                return -3. / (4. * r) + 4. * r * (
                    2. * mu * (stateEnergy - self.potential(l, s, j, r))
                    - l * (l + 1) / (r**2)
                )

            r, psi_r = NumerovBack(innerLimit, outerLimit, potential,
                                   step, 0.01, 0.01)

            suma = np.trapz(psi_r**2, x=r)
            psi_r = psi_r / (sqrt(suma))

        return r, psi_r

    def _parseLevelsFromNIST(self, fileData):
        """
            Parses the level energies from file listing the NIST ASD data

            Args:
                fileData (str): path to the file containing NIST ASD data for
                    the element
        """
        f = open(fileData, "r")
        l = 0
        n = 0
        levels = []
        for line in f:

            line = re.sub('[\[\]]', '', line)
            pattern = "\.\d*[spdfgh]"
            pattern2 = "\|\s+\d*/"
            pattern3 = "/\d* \|"
            pattern4 = "\| *\d*\.\d* *\|"
            match = re.search(pattern, line)
            if match is not None:
                n = int(line[match.start() + 1:match.end() - 1])
            if match is not None:
                ch = line[match.end() - 1:match.end()]
                if ch == "s":
                    l = 0
                elif ch == "p":
                    l = 1
                elif ch == "d":
                    l = 2
                elif ch == "f":
                    l = 3
                elif ch == "g":
                    l = 4
                elif ch == "h":
                    l = 5
                else:
                    print("Unidentified character in line:\n", line)
                    exit()

            match = re.search(pattern2, line)
            if match is not None:
                br1 = float(line[match.start() + 2:match.end() - 1])
                match = re.search(pattern3, line)
                br2 = float(line[match.start() + 1:match.end() - 2])
                match = re.search(pattern4, line)
                energyValue = float(line[match.start() + 1:match.end() - 1])
                levels.append([n, l, br1 / br2, energyValue])
        f.close()
        return levels

    def _addEnergy(self, n, l, j, energyNIST):
        """
            Adding energy levels

            Accepts energy level relative to **ground state**, and
            saves energy levels, relative to the **ionization treshold**.

            Args:
                energyNIST (float): energy relative to
                    the nonexcited level (= 0 eV)
        """
        #
        if abs(j - (l - 0.5)) < 0.001:
            # j =l-1/2
            self.sEnergy[n, l] = energyNIST - self.ionisationEnergy
        else:
            # j = l+1/2
            self.sEnergy[l, n] = energyNIST - self.ionisationEnergy

    def getTransitionWavelength(self, n1, l1, j1, n2, l2, j2, s=0.5, s2=None):
        """
            Calculated transition wavelength (in vacuum) in m.

            Returned values is given relative to the centre of gravity of the
            hyperfine-split states.

            Args:
                n1 (int): principal quantum number of the state
                    **from** which we are going
                l1 (int): orbital angular momentum of the state
                    **from** which we are going
                j1 (float): total angular momentum of the state
                    **from** which we are going
                n2 (int): principal quantum number of the state
                    **to** which we are going
                l2 (int): orbital angular momentum of the state
                    **to** which we are going
                j2 (float): total angular momentum of the state
                    **to** which we are going
                s (float): optional, spin of the intial state
                    (for Alkali this is fixed to 0.5)
                s2 (float): optional, spin of the final state.
                    If not set, defaults to same value as :obj:`s`

            Returns:
                float:
                    transition wavelength (in m). If the returned value is
                    negative, level from which we are going is **above**
                    the level to which we are going.
        """
        if s2 is None:
            s2 = s
        return (C_h * C_c) / ((self.getEnergy(n2, l2, j2, s=s2)
                               - self.getEnergy(n1, l1, j1, s=s)) * C_e)

    def getTransitionFrequency(self, n1, l1, j1, n2, l2, j2, s=0.5, s2=None):
        """
            Calculated transition frequency in Hz

            Returned values is given relative to the centre of gravity of the
            hyperfine-split states.

            Args:
                n1 (int): principal quantum number of the state
                    **from** which we are going
                l1 (int): orbital angular momentum of the state
                    **from** which we are going
                j1 (float): total angular momentum of the state
                    **from** which we are going
                n2 (int): principal quantum number of the state
                    **to** which we are going
                l2 (int): orbital angular momentum of the state
                    **to** which we are going
                j2 (float): total angular momentum of the state
                    **to** which we are going
                s (float): optional, spin of the intial state
                    (for Alkali this is fixed to 0.5)
                s2 (float): optional, spin of the final state
                    If not set, defaults to the same value as :obj:`s`

            Returns:
                float:
                    transition frequency (in Hz). If the returned value is
                    negative, level from which we are going is **above**
                    the level to which we are going.
        """
        if s2 is None:
            s2 = s
        return (self.getEnergy(n2, l2, j2, s=s2)
                - self.getEnergy(n1, l1, j1, s=s)) * C_e / C_h

    def getEnergy(self, n, l, j, s=0.5):
        """
            Energy of the level relative to the ionisation level (in eV)

            Returned energies are with respect to the center of gravity of the
            hyperfine-split states.
            If `preferQuantumDefects` =False (set during initialization)
            program will try use NIST energy value, if such exists,
            falling back to energy calculation with quantum defects if
            the measured value doesn't exist. For `preferQuantumDefects` =True,
            program will calculate energies from quantum defects
            (useful for comparing quantum defect calculations with measured
            energy level values) if the principal quantum number of the
            requested state is larger than the minimal quantum principal quantum
            number `self.minQuantumDefectN` which sets minimal quantum number
            for which quantum defects still give good estimate of state energy
            (below this value saved energies will be used if existing).

            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                s (float): optional, total spin angular momentum. Default value
                    of 0.5 is correct for Alkali atoms, and has to be specified
                    explicitly for divalent atoms.

            Returns:
                float: state energy (eV)
        """
        if l >= n:
            raise ValueError(
                "Requested energy for state l=%d >= n=%d !" % (l, n))

        # use NIST data ?
        if (not self.preferQuantumDefects
            or n < self.minQuantumDefectN) and (n <= self.NISTdataLevels) and \
                (abs(self._getSavedEnergy(n, l, j, s=s)) > 1e-8):
            return self._getSavedEnergy(n, l, j, s=s)

        # else, use quantum defects
        defect = self.getQuantumDefect(n, l, j, s=s)
        return -self.scaledRydbergConstant / ((n - defect)**2)

    def _getSavedEnergy(self, n, l, j, s=0.5):
        if abs(j - (l - 0.5)) < 0.001:
            # j = l-1/2
            return self.sEnergy[n, l]
        elif abs(j - (l + 0.5)) < 0.001:
            # j =l+1/2
            return self.sEnergy[l, n]
        else:
            raise ValueError("j (=%.1f) is not equal to l+1/2 nor l-1/2 (l=%d)"
                             % (j, l))

    def getQuantumDefect(self, n, l, j, s=0.5):
        """
            Quantum defect of the level.

            For an example, see `Rydberg energy levels example snippet`_.

            .. _`Rydberg energy levels example snippet`:
                ./Rydberg_atoms_a_primer.html#Rydberg-Atom-Energy-Levels

            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                s (float): (optional). Total spin angular momentum.
                    Default value of 0.5 correct for Alkali atoms. For divalent
                    atoms it has to be explicitly defined.

            Returns:
                float: quantum defect
        """
        defect = 0.0
        if (l < 5):
            # find correct part in table of quantum defects
            modifiedRRcoef = self.quantumDefect[int(floor(s) + s + j - l)][l]
            if (l < 3 and abs(modifiedRRcoef[0]) < 1e-9
                    and self.Z != 1):
                # it's not Hydrogen but for l in {s,p,d} quantum defect is 0
                raise ValueError("Quantum defects for requested state "
                                 + ("(n = %d, l = %d, j = %.1f, s=%.1f) are"
                                    % (n, l, j, s)) +
                                 " uknown. Aborting calculation.")
            defect = modifiedRRcoef[0] + \
                modifiedRRcoef[1] / ((n - modifiedRRcoef[0])**2) + \
                modifiedRRcoef[2] / ((n - modifiedRRcoef[0])**4) + \
                modifiedRRcoef[3] / ((n - modifiedRRcoef[0])**6) + \
                modifiedRRcoef[4] / ((n - modifiedRRcoef[0])**8) + \
                modifiedRRcoef[5] / ((n - modifiedRRcoef[0])**10)
        else:
            # use \delta_\ell = \delta_g * (4/\ell)**5
            # from https://journals.aps.org/pra/abstract/10.1103/PhysRevA.74.062712
            defect = self.quantumDefect[0][4][0] * (4 / l) ** 5
        return defect

    def getRadialMatrixElement(self,
                               n1, l1, j1,
                               n2, l2, j2,
                               s=0.5,
                               useLiterature=True):
        """
            Radial part of the dipole matrix element

            Calculates :math:`\\int \\mathbf{d}r~R_{n_1,l_1,j_1}(r)\\cdot \
                R_{n_1,l_1,j_1}(r) \\cdot r^3`.

            Args:
                n1 (int): principal quantum number of state 1
                l1 (int): orbital angular momentum of state 1
                j1 (float): total angular momentum of state 1
                n2 (int): principal quantum number of state 2
                l2 (int): orbital angular momentum of state 2
                j2 (float): total angular momentum of state 2
                s (float): optional, total spin angular momentum of state 1.
                    By default 0.5 for Alkali atoms.
                useLiterature (bool): optional, should literature values for
                    dipole matrix element be used if existing? If true,
                    compiled values stored in `literatureDMEfilename` variable
                    for a given atom (file is stored locally at ~/.arc-data/),
                    will be checked, and if the value is found, selects the
                    value with smallest error estimate (if there are multiple
                    entries). If no value is found, it will default to numerical
                    integration of wavefunctions. By default True.

            Returns:
                float: dipole matrix element (:math:`a_0 e`).
        """
        dl = abs(l1 - l2)
        dj = abs(j1 - j2)
        if not(dl == 1 and (dj < 1.1)):
            return 0

        if (self.getEnergy(n1, l1, j1, s=s)
                > self.getEnergy(n2, l2, j2, s=s)):
            temp = n1
            n1 = n2
            n2 = temp
            temp = l1
            l1 = l2
            l2 = temp
            temp = j1
            j1 = j2
            j2 = temp

        n1 = int(n1)
        n2 = int(n2)
        l1 = int(l1)
        l2 = int(l2)
        j1_x2 = int(round(2 * j1))
        j2_x2 = int(round(2 * j2))

        c = self.conn.cursor()

        if useLiterature:
            # is there literature value for this DME? If there is,
            # use the best one (smalles error)
            c.execute('''SELECT dme FROM literatureDME WHERE
             n1= ? AND l1 = ? AND j1_x2 = ? AND
             n2 = ? AND l2 = ? AND j2_x2 = ?
             ORDER BY errorEstimate ASC''', (n1, l1, j1_x2, n2, l2, j2_x2))
            answer = c.fetchone()
            if (answer):
                # we did found literature value
                return answer[0]

        # was this calculated before? If it was, retrieve from memory
        c.execute('''SELECT dme FROM dipoleME WHERE
         n1= ? AND l1 = ? AND j1_x2 = ? AND
         n2 = ? AND l2 = ? AND j2_x2 = ?''', (n1, l1, j1_x2, n2, l2, j2_x2))
        dme = c.fetchone()
        if (dme):
            return dme[0]

        step = 0.001
        r1, psi1_r1 = self.radialWavefunction(l1, 0.5, j1,
                                              self.getEnergy(
                                                  n1, l1, j1) / 27.211,
                                              self.alphaC**(1 / 3.0),
                                              2.0 * n1 * (n1 + 15.0), step)
        r2, psi2_r2 = self.radialWavefunction(l2, 0.5, j2,
                                              self.getEnergy(
                                                  n2, l2, j2) / 27.211,
                                              self.alphaC**(1 / 3.0),
                                              2.0 * n2 * (n2 + 15.0), step)

        upTo = min(len(r1), len(r2))

        # note that r1 and r2 change in same staps,
        # starting from the same value
        dipoleElement = np.trapz(
            np.multiply(np.multiply(psi1_r1[0:upTo], psi2_r2[0:upTo]),
                        r1[0:upTo]),
            x=r1[0:upTo]
        )

        c.execute(''' INSERT INTO dipoleME VALUES (?,?,?, ?,?,?, ?)''',
                       [n1, l1, j1_x2, n2, l2, j2_x2, dipoleElement])

        self.conn.commit()

        return dipoleElement

    def getQuadrupoleMatrixElement(self, n1, l1, j1, n2, l2, j2,
                                   s=0.5):
        """
            Radial part of the quadrupole matrix element

            Calculates :math:`\\int \\mathbf{d}r~R_{n_1,l_1,j_1}(r)\\cdot \
            R_{n_1,l_1,j_1}(r) \\cdot r^4`.
            See `Quadrupole calculation example snippet`_  .

            .. _`Quadrupole calculation example snippet`:
                ./Rydberg_atoms_a_primer.html#Quadrupole-matrix-elements

            Args:
                n1 (int): principal quantum number of state 1
                l1 (int): orbital angular momentum of state 1
                j1 (float): total angular momentum of state 1
                n2 (int): principal quantum number of state 2
                l2 (int): orbital angular momentum of state 2
                j2 (float): total angular momentum of state 2
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: quadrupole matrix element (:math:`a_0^2 e`).
        """

        dl = abs(l1 - l2)
        dj = abs(j1 - j2)
        if not ((dl == 0 or dl == 2 or dl == 1)and (dj < 2.1)):
            return 0

        if (self.getEnergy(n1, l1, j1, s=s) > self.getEnergy(n2, l2, j2, s=s)):
            temp = n1
            n1 = n2
            n2 = temp
            temp = l1
            l1 = l2
            l2 = temp
            temp = j1
            j1 = j2
            j2 = temp

        n1 = int(n1)
        n2 = int(n2)
        l1 = int(l1)
        l2 = int(l2)
        j1_x2 = int(round(2 * j1))
        j2_x2 = int(round(2 * j2))

        c = self.conn.cursor()
        # was this calculated before? If yes, retrieve from memory.
        c.execute('''SELECT qme FROM quadrupoleME WHERE
         n1= ? AND l1 = ? AND j1_x2 = ? AND
         n2 = ? AND l2 = ? AND j2_x2 = ?''', (n1, l1, j1_x2, n2, l2, j2_x2))
        qme = c.fetchone()
        if (qme):
            return qme[0]

        # if it wasn't, calculate now

        step = 0.001
        r1, psi1_r1 = self.radialWavefunction(l1, 0.5, j1,
                                              self.getEnergy(
                                                  n1, l1, j1) / 27.211,
                                              self.alphaC**(1 / 3.0),
                                              2.0 * n1 * (n1 + 15.0), step)
        r2, psi2_r2 = self.radialWavefunction(l2, 0.5, j2,
                                              self.getEnergy(
                                                  n2, l2, j2) / 27.211,
                                              self.alphaC**(1 / 3.0),
                                              2.0 * n2 * (n2 + 15.0), step)

        upTo = min(len(r1), len(r2))

        # note that r1 and r2 change in same staps,
        # starting from the same value
        quadrupoleElement = np.trapz(
            np.multiply(np.multiply(psi1_r1[0:upTo], psi2_r2[0:upTo]),
                        np.multiply(r1[0:upTo], r1[0:upTo])
                        ),
            x=r1[0:upTo]
        )

        c.execute(''' INSERT INTO quadrupoleME VALUES (?,?,?,?,?,?, ?)''',
                       [n1, l1, j1_x2, n2, l2, j2_x2, quadrupoleElement])

        self.conn.commit()

        return quadrupoleElement

    def getReducedMatrixElementJ_asymmetric(self, n1, l1, j1, n2, l2, j2,
                                            s=0.5):
        """
            Reduced matrix element in :math:`J` basis, defined in asymmetric
            notation.

            Note that notation for symmetric and asymmetricly defined
            reduced matrix element is not consistent in the literature. For
            example, notation is used e.g. in Steck [1]_ is precisely
            the oposite.

            Note:
                Note that this notation is asymmetric: :math:`( j||e r \
                ||j' ) \\neq ( j'||e r ||j )`.
                Relation between the two notation is :math:`\\langle j||er||j'\
                \\rangle=\\sqrt{2j+1} ( j ||er ||j')`.
                This function always returns value for transition from
                lower to higher energy state, independent of the order of
                states entered in the function call.

            Args:
                n1 (int): principal quantum number of state 1
                l1 (int): orbital angular momentum of state 1
                j1 (float): total angular momentum of state 1
                n2 (int): principal quantum number of state 2
                l2 (int): orbital angular momentum of state 2
                j2 (float): total angular momentum of state 2
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:
                    reduced dipole matrix element in Steck notation
                    :math:`( j || er || j' )` (:math:`a_0 e`).

            .. [1] Daniel A. Steck, "Cesium D Line Data," (revision 2.0.1,
                2 May 2008). http://steck.us/alkalidata
        """
        #
        if (self.getTransitionFrequency(n1, l1, j1,
                                        n2, l2, j2,
                                        s=s, s2=s) < 0):
            temp = n2
            n2 = n1
            n1 = temp
            temp = l1
            l1 = l2
            l2 = temp
            temp = j1
            j1 = j2
            j2 = temp
        return (-1)**(int((l2 + l1 + 3.) / 2. + s + j2)) *\
            sqrt((2.0 * j2 + 1.0) * (2.0 * l1 + 1.0)) *\
            Wigner6j(l1, l2, 1, j2, j1, s) *\
            sqrt(float(max(l1, l2)) / (2.0 * l1 + 1.0)) *\
            self.getRadialMatrixElement(n1, l1, j1, n2, l2, j2, s=s)

    def getReducedMatrixElementL(self, n1, l1, j1, n2, l2, j2, s=0.5):
        """
            Reduced matrix element in :math:`L` basis (symmetric notation)

            Args:
                n1 (int): principal quantum number of state 1
                l1 (int): orbital angular momentum of state 1
                j1 (float): total angular momentum of state 1
                n2 (int): principal quantum number of state 2
                l2 (int): orbital angular momentum of state 2
                j2 (float): total angular momentum of state 2

            Returns:
                float:
                    reduced dipole matrix element in :math:`L` basis
                    :math:`\\langle l || er || l' \\rangle` (:math:`a_0 e`).
        """

        return (-1)**l1 * sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0)) *\
            Wigner3j(l1, 1, l2, 0, 0, 0) *\
            self.getRadialMatrixElement(n1, l1, j1, n2, l2, j2, s=s)

    def getReducedMatrixElementJ(self, n1, l1, j1, n2, l2, j2,
                                 s=0.5):
        """
            Reduced matrix element in :math:`J` basis (symmetric notation)

            Args:
                n1 (int): principal quantum number of state 1
                l1 (int): orbital angular momentum of state 1
                j1 (float): total angular momentum of state 1
                n2 (int): principal quantum number of state 2
                l2 (int): orbital angular momentum of state 2
                j2 (float): total angular momentum of state 2
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:
                    reduced dipole matrix element in :math:`J` basis
                    :math:`\\langle j || er || j' \\rangle` (:math:`a_0 e`).
        """

        return (-1)**(int(l1 + s + j2 + 1.)) * sqrt((2. * j1 + 1.)
                                                    * (2. * j2 + 1.)) *\
            Wigner6j(j1, 1., j2, l2, s, l1) *\
            self.getReducedMatrixElementL(n1, l1, j1, n2, l2, j2, s=s)

    def getDipoleMatrixElement(self, n1, l1, j1, mj1, n2, l2, j2, mj2, q,
                               s=0.5):
        r"""
            Dipole matrix element
            :math:`\langle n_1 l_1 j_1 m_{j_1} |e\mathbf{r}|\
            n_2 l_2 j_2 m_{j_2}\rangle`
            in units of :math:`a_0 e`

            Args:
                n1. l1, j1, mj1: principal, orbital, total angular momentum,
                    and projection of total angular momentum for state 1
                n2. l2, j2, mj2: principal, orbital, total angular momentum,
                    and projection of total angular momentum for state 2
                q (int): specifies transition that the driving field couples to,
                    +1, 0 or -1 corresponding to driving :math:`\sigma^+`,
                    :math:`\pi` and :math:`\sigma^-` transitions respectively.
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: dipole matrix element( :math:`a_0 e`)

            Example:

                For example, calculation of :math:`5 S_{1/2}m_j=-\frac{1}{2}\
                \rightarrow  5 P_{3/2}m_j=-\frac{3}{2}`
                transition dipole matrix element for laser driving
                :math:`\sigma^-` transition::

                    from arc import *
                    atom = Rubidium()
                    # transition 5 S_{1/2} m_j=-0.5 -> 5 P_{3/2} m_j=-1.5
                    # for laser driving sigma- transition
                    print(atom.getDipoleMatrixElement(5,0,0.5,-0.5,5,1,1.5,-1.5,-1))


        """
        if abs(q) > 1.1:
            return 0
#        return (-1)**(int(j1 - mj1)) *\
#            Wigner3j(j1, 1, j2, -mj1, -q, mj2) *\
#            self.getReducedMatrixElementJ(n1, l1, j1, n2, l2, j2, s=s)
        return self.getSphericalDipoleMatrixElement(j1, mj1, j2, mj2, q) * \
            self.getReducedMatrixElementJ(n1, l1, j1, n2, l2, j2, s=s)

    def getDipoleMatrixElementHFS(self,
                                  n1, l1, j1, f1, mf1,
                                  n2, l2, j2, f2, mf2,
                                  q,
                                  s=0.5):
        r"""
        Dipole matrix element for hyperfine structure resolved transitions
        :math:`\langle n_1 l_1 j_1 f_1 m_{f_1} |e\mathbf{r}|\
        n_2 l_2 j_2 f_2 m_{f_2}\rangle`
        in units of :math:`a_0 e`

        For hyperfine resolved transitions, the dipole matrix element is
        :math:`\langle n_1,\ell_1,j_1,f_1,m_{f1} |  \
        \mathbf{\hat{r}}\cdot \mathbf{\varepsilon}_q  \
        | n_2,\ell_2,j_2,f_2,m_{f2} \rangle = (-1)^{f_1-m_{f1}} \
        \left( \
        \begin{matrix} \
        f_1 & 1 & f_2 \\ \
        -m_{f1} & q & m_{f2} \
        \end{matrix}\right) \
        \langle n_1 \ell_1 j_1 f_1|| r || n_2 \ell_2 j_2 f_2 \rangle,` where
        :math:`\langle n_1 \ell_1 j_1 f_1 ||r|| n_2 \ell_2 j_2 f_2 \rangle \
        = (-1)^{j_1+I+F_2+1}\sqrt{(2f_1+1)(2f_2+1)} ~ \
        \left\{ \begin{matrix}\
        F_1 & 1 & F_2 \\ \
        j_2 & I & j_1 \
        \end{matrix}\right\}~ \
        \langle n_1 \ell_1 j_1||r || n_2 \ell_2 j_2 \rangle.`


        Args:
            n1. l1, j1, f1, mf1: principal, orbital, total orbital,
                fine basis (total atomic) angular momentum,
                and projection of total angular momentum for state 1
            n2. l2, j2, f2, mf2: principal, orbital, total orbital,
                fine basis (total atomic) angular momentum,
                and projection of total angular momentum for state 2
            q (int): specifies transition that the driving field couples to,
                +1, 0 or -1 corresponding to driving :math:`\sigma^+`,
                :math:`\pi` and :math:`\sigma^-` transitions respectively.
            s (float): optional, total spin angular momentum of state.
                By default 0.5 for Alkali atoms.

        Returns:
            float: dipole matrix element( :math:`a_0 e`)


        """
#        dme = (- 1)**(f1 - mf1) * Wigner3j(f1, 1, f2, - mf1, -q, mf2)
#        dme *= (- 1)**(j1 + self.I + f2 + 1) * ((2. * f1 + 1)
#                                                * (2 * f2 + 1))**0.5
#        dme *= Wigner6j(f1, 1, f2, j2, self.I, j1)

        dme = self.getSphericalDipoleMatrixElement(f1, mf1, f2, mf2, q)
        dme *= self._reducedMatrixElementFJ(j1, f1, j2, f2)
        dme *= self.getReducedMatrixElementJ(n1, l1, j1, n2, l2, j2, s=s)
        return dme

    def getRabiFrequency(self,
                         n1, l1, j1, mj1,
                         n2, l2, j2, q,
                         laserPower, laserWaist,
                         s=0.5):
        """
            Returns a Rabi frequency for resonantly driven atom in a
            center of TEM00 mode of a driving field

            Args:
                n1,l1,j1,mj1 : state from which we are driving transition
                n2,l2,j2 : state to which we are driving transition
                q : laser polarization (-1,0,1 correspond to :math:`\\sigma^-`,
                    :math:`\\pi` and :math:`\\sigma^+` respectively)
                laserPower : laser power in units of W
                laserWaist : laser :math:`1/e^2` waist (radius) in units of m
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:
                    Frequency in rad :math:`^{-1}`. If you want frequency
                    in Hz, divide by returned value by :math:`2\\pi`
        """
        maxIntensity = 2 * laserPower / (pi * laserWaist**2)
        electricField = sqrt(2. * maxIntensity / (C_c * epsilon_0))
        return self.getRabiFrequency2(n1, l1, j1, mj1,
                                      n2, l2, j2, q,
                                      electricField,
                                      s=s)

    def getRabiFrequency2(self,
                          n1, l1, j1, mj1,
                          n2, l2, j2, q,
                          electricFieldAmplitude,
                          s=0.5):
        """
            Returns a Rabi frequency for resonant excitation with a given
            electric field amplitude

            Args:
                n1,l1,j1,mj1 : state from which we are driving transition
                n2,l2,j2 : state to which we are driving transition
                q : laser polarization (-1,0,1 correspond to :math:`\\sigma^-`,
                    :math:`\\pi` and :math:`\\sigma^+` respectively)
                electricFieldAmplitude : amplitude of electric field
                    driving (V/m)
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:
                    Frequency in rad :math:`^{-1}`. If you want frequency
                    in Hz, divide by returned value by :math:`2\\pi`
        """
        mj2 = mj1 + q
        if abs(mj2) - 0.1 > j2:
            return 0
        dipole = self.getDipoleMatrixElement(n1, l1, j1, mj1,
                                             n2, l2, j2, mj2, q,
                                             s=s) *\
            C_e * physical_constants["Bohr radius"][0]
        freq = electricFieldAmplitude * abs(dipole) / hbar
        return freq

    def getC6term(self, n, l, j, n1, l1, j1, n2, l2, j2, s=0.5):
        """
            C6 interaction term for the given two pair-states

            Calculates :math:`C_6` intaraction term for :math:`|n,l,j,n,l,j\
            \\rangle \\leftrightarrow |n_1,l_1,j_1,n_2,l_2,j_2\\rangle`.
            For details of calculation see Ref. [#c6r1]_.

            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                n1 (int): principal quantum number
                l1 (int): orbital angular momentum
                j1 (float): total angular momentum
                n2 (int): principal quantum number
                l2 (int): orbital angular momentum
                j2 (float): total angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:  :math:`C_6 = \\frac{1}{4\\pi\\varepsilon_0} \
                    \\frac{|\\langle n,l,j |er|n_1,l_1,j_1\\rangle|^2|\
                    \\langle n,l,j |er|n_2,l_2,j_2\\rangle|^2}\
                    {E(n_1,l_1,j_2,n_2,j_2,j_2)-E(n,l,j,n,l,j)}`
                (:math:`h` Hz m :math:`{}^6`).

            Example:
                We can reproduce values from Ref. [#c6r1]_ for C3 coupling
                to particular channels. Taking for example channels described
                by the Eq. (50a-c) we can get the values::

                    from arc import *

                    channels = [[70,0,0.5, 70, 1,1.5, 69,1, 1.5],\\
                                [70,0,0.5, 70, 1,1.5, 69,1, 0.5],\\
                                [70,0,0.5, 69, 1,1.5, 70,1, 0.5],\\
                                [70,0,0.5, 70, 1,0.5, 69,1, 0.5]]

                    print(" = = = Caesium = = = ")
                    atom = Caesium()
                    for channel in channels:
                        print("%.0f  GHz (mu m)^6" % ( atom.getC6term(*channel)
                                                      / C_h * 1.e27 ))

                    print("\\n = = = Rubidium  = = =")
                    atom = Rubidium()
                    for channel in channels:
                        print("%.0f  GHz (mu m)^6" % ( atom.getC6term(*channel)
                                                      / C_h * 1.e27 ))

                Returns::

                     = = = Caesium = = =
                    722  GHz (mu m)^6
                    316  GHz (mu m)^6
                    383  GHz (mu m)^6
                    228  GHz (mu m)^6

                     = = = Rubidium  = = =
                    799  GHz (mu m)^6
                    543  GHz (mu m)^6
                    589  GHz (mu m)^6
                    437  GHz (mu m)^6

                which is in good agreement with the values cited in the
                Ref. [#c6r1]_. Small discrepancies for Caesium originate from
                slightly different quantum defects used in calculations.


            References:
                .. [#c6r1] T. G. Walker, M. Saffman, PRA **77**, 032723 (2008)
                    https://doi.org/10.1103/PhysRevA.77.032723

        """
        d1 = self.getRadialMatrixElement(n, l, j, n1, l1, j1, s=s)
        d2 = self.getRadialMatrixElement(n, l, j, n2, l2, j2, s=s)
        d1d2 = 1 / (4.0 * pi * epsilon_0) * d1 * d2 * C_e**2 *\
            (physical_constants["Bohr radius"][0])**2
        return -d1d2**2 / (C_e * (self.getEnergy(n1, l1, j1, s=s)
                                  + self.getEnergy(n2, l2, j2, s=s)
                                  - 2 * self.getEnergy(n, l, j, s=s)))

    def getC3term(self, n, l, j, n1, l1, j1, n2, l2, j2, s=0.5):
        """
            C3 interaction term for the given two pair-states

            Calculates :math:`C_3` intaraction term for
                :math:`|n,l,j,n,l,j\\rangle \
                 \\leftrightarrow |n_1,l_1,j_1,n_2,l_2,j_2\\rangle`

            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                n1 (int): principal quantum number
                l1 (int): orbital angular momentum
                j1 (float): total angular momentum
                n2 (int): principal quantum number
                l2 (int): orbital angular momentum
                j2 (float): total angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:  :math:`C_3 = \\frac{\\langle n,l,j |er\
                |n_1,l_1,j_1\\rangle \
                \\langle n,l,j |er|n_2,l_2,j_2\\rangle}{4\\pi\\varepsilon_0}`
                (:math:`h` Hz m :math:`{}^3`).
        """
        d1 = self.getRadialMatrixElement(n, l, j, n1, l1, j1, s=s)
        d2 = self.getRadialMatrixElement(n, l, j, n2, l2, j2, s=s)
        d1d2 = 1 / (4.0 * pi * epsilon_0) * d1 * d2 * C_e**2 *\
            (physical_constants["Bohr radius"][0])**2
        return d1d2

    def getEnergyDefect(self, n, l, j, n1, l1, j1, n2, l2, j2,
                        s=0.5):
        """
            Energy defect for the given two pair-states (one of the state has
            two atoms in the same state)

            Energy difference between the states
            :math:`E(n_1,l_1,j_1,n_2,l_2,j_2) - E(n,l,j,n,l,j)`

            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                n1 (int): principal quantum number
                l1 (int): orbital angular momentum
                j1 (float): total angular momentum
                n2 (int): principal quantum number
                l2 (int): orbital angular momentum
                j2 (float): total angular momentum
                s (float): optional. Spin angular momentum
                    (default 0.5 for Alkali)

            Returns:
                float:  energy defect (SI units: J)
        """
        return C_e * (self.getEnergy(n1, l1, j1, s=s)
                      + self.getEnergy(n2, l2, j2, s=s)
                      - 2 * self.getEnergy(n, l, j, s=s))

    def getEnergyDefect2(self, n, l, j, nn, ll, jj, n1, l1, j1, n2, l2, j2,
                         s=0.5):
        """
            Energy defect for the given two pair-states

            Energy difference between the states
            :math:`E(n_1,l_1,j_1,n_2,l_2,j_2) - E(n,l,j,nn,ll,jj)`

            See `pair-state energy defects example snippet`_.

            .. _`pair-state energy defects example snippet`:
                ./Rydberg_atoms_a_primer.html#Rydberg-atom-interactions


            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                nn (int): principal quantum number
                ll (int): orbital angular momentum
                jj (float): total angular momentum
                n1 (int): principal quantum number
                l1 (int): orbital angular momentum
                j1 (float): total angular momentum
                n2 (int): principal quantum number
                l2 (int): orbital angular momentum
                j2 (float): total angular momentum
                s (float): optional. Spin angular momentum
                    (default 0.5 for Alkali)

            Returns:
                float:  energy defect (SI units: J)
        """
        return C_e * (self.getEnergy(n1, l1, j1, s=s)
                      + self.getEnergy(n2, l2, j2, s=s)
                      - self.getEnergy(n, l, j, s=s)
                      - self.getEnergy(nn, ll, jj, s=s))

    def updateDipoleMatrixElementsFile(self):
        """
            Updates the file with pre-calculated dipole matrix elements.

            This function will add the the file all the elements that have been
            calculated in the previous run, allowing quick access to them in
            the future calculations.
        """
        # obtain dipole matrix elements from the database

        dipoleMatrixElement = []
        c = self.conn.cursor()
        c.execute('''SELECT * FROM dipoleME ''')
        for v in c.fetchall():
            dipoleMatrixElement.append(v)

        # obtain quadrupole matrix elements from the database

        quadrupoleMatrixElement = []
        c.execute('''SELECT * FROM quadrupoleME ''')
        for v in c.fetchall():
            quadrupoleMatrixElement.append(v)

        # save dipole elements
        try:
            np.save(os.path.join(self.dataFolder,
                                 self.dipoleMatrixElementFile),
                    dipoleMatrixElement)
        except IOError as e:
            print("Error while updating dipoleMatrixElements File "
                  + self.dipoleMatrixElementFile)
            print(e)
        # save quadrupole elements
        try:
            np.save(os.path.join(self.dataFolder,
                                 self.quadrupoleMatrixElementFile),
                    quadrupoleMatrixElement)
        except IOError as e:
            print("Error while updating quadrupoleMatrixElements File "
                  + self.quadrupoleMatrixElementFile)
            print(e)

    def getTransitionRate(self, n1, l1, j1, n2, l2, j2, temperature=0.,
                          s=0.5):
        """
            Transition rate due to coupling to vacuum modes
            (black body included)

            Calculates transition rate from the first given state to the second
            given state :math:`|n_1,l_1,j_1\\rangle \\rightarrow \
            |n_2,j_2,j_2\\rangle` at given temperature due to interaction with
            the vacuum field. For zero temperature this returns Einstein A
            coefficient. For details of calculation see Ref. [#lf1]_ and
            Ref. [#lf2]_.
            See `Black-body induced population transfer example snippet`_.

            .. _`Black-body induced population transfer example snippet`:
                ./Rydberg_atoms_a_primer.html#Rydberg-Atom-Lifetimes

            Args:
                n1 (int): principal quantum number
                l1 (int): orbital angular momentum
                j1 (float): total angular momentum
                n2 (int): principal quantum number
                l2 (int): orbital angular momentum
                j2 (float): total angular momentum
                [temperature] (float): temperature in K
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:  transition rate in s :math:`{}^{-1}` (SI)

            References:
                .. [#lf1] C. E. Theodosiou, PRA **30**, 2881 (1984)
                    https://doi.org/10.1103/PhysRevA.30.2881

                .. [#lf2] I. I. Beterov, I. I. Ryabtsev, D. B. Tretyakov,\
                    and V. M. Entin, PRA **79**, 052504 (2009)
                    https://doi.org/10.1103/PhysRevA.79.052504
        """

        degeneracyTerm = 1.

        # find dipoleRadialPart
        if (self.getTransitionFrequency(n1, l1, j1, n2, l2, j2,
                                        s=s, s2=s) > 0):
            dipoleRadialPart = self.getReducedMatrixElementJ_asymmetric(
                n1, l1, j1,
                n2, l2, j2,
                s=s) *\
                C_e * (physical_constants["Bohr radius"][0])

        else:
            dipoleRadialPart = self.getReducedMatrixElementJ_asymmetric(
                n2, l2, j2,
                n1, l1, j1,
                s=s) *\
                C_e * (physical_constants["Bohr radius"][0])
            degeneracyTerm = (2. * j2 + 1.0) / (2. * j1 + 1.)

        omega = abs(
            2.0 * pi * self.getTransitionFrequency(n1, l1, j1, n2, l2, j2,
                                                   s=s, s2=s))

        modeOccupationTerm = 0.
        if (self.getTransitionFrequency(n1, l1, j1, n2, l2, j2,
                                        s=s, s2=s) < 0):
            modeOccupationTerm = 1.

        # only possible by absorbing thermal photons ?
        if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
            modeOccupationTerm += 1. / \
                (exp(hbar * omega / (C_k * temperature)) - 1.)

        return omega**3 * dipoleRadialPart**2 /\
            (3 * pi * epsilon_0 * hbar * C_c**3)\
            * degeneracyTerm * modeOccupationTerm

    def getStateLifetime(self, n, l, j, temperature=0, includeLevelsUpTo=0,
                         s=0.5):
        """
            Returns the lifetime of the state (in s)

            For non-zero temperatures, user must specify up to which principal
            quantum number levels, that is **above** the initial state, should
            be included in order to account for black-body induced transitions
            to higher lying states. See `Rydberg lifetimes example snippet`_.

            .. _`Rydberg lifetimes example snippet`:
                ./Rydberg_atoms_a_primer.html#Rydberg-Atom-Lifetimes

            Args:
                n, l, j (int,int,float): specifies state whose lifetime we are
                    calculating
                temperature : optional. Temperature at which the atom
                    environment is, measured in K. If this parameter
                    is non-zero, user has to specify transitions up to
                    which state (due to black-body decay) should be included
                    in calculation.
                includeLevelsUpTo (int): optional and not needed for atom
                    lifetimes calculated at zero temperature. At non zero
                    temperatures, this specify maximum principal quantum number
                    of the state to which black-body induced transitions will
                    be included. Minimal value of the parameter in that case is
                    :math:`n+1`
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.


            Returns:
                float:
                    State lifetime in units of s (seconds)

            See also:
                :obj:`getTransitionRate` for calculating rates of individual
                transition rates between the two states

        """
        if (temperature > 0.1 and includeLevelsUpTo <= n):
            raise ValueError(
                "For non-zero temperatures, user has to specify "
                + "principal quantum number of the maximum state *above* the "
                + "state for which we are calculating the lifetime. This is "
                + "in order to include black-body induced transitions to "
                + " higher lying up in energy levels.")
        elif (temperature < 0.1):
            includeLevelsUpTo = max(n, self.groundStateN)

        transitionRate = 0.

        for nto in xrange(max(self.groundStateN, l), includeLevelsUpTo + 1):

            # sum over all l-1
            if l > 0:
                lto = l - 1
                if lto > j - 0.5 - 0.1:
                    jto = j
                    transitionRate += self.getTransitionRate(n, l, j,
                                                             nto, lto, jto,
                                                             temperature,
                                                             s=s)
                jto = j - 1.
                if jto > 0:
                    transitionRate += self.getTransitionRate(n, l, j,
                                                             nto, lto, jto,
                                                             temperature,
                                                             s=s)

        for nto in xrange(max(self.groundStateN, l + 2),
                          includeLevelsUpTo + 1):
            # sum over all l+1
            lto = l + 1
            if lto - 0.5 - 0.1 < j:
                jto = j
                transitionRate += self.getTransitionRate(n, l, j,
                                                         nto, lto, jto,
                                                         temperature,
                                                         s=s)
            jto = j + 1
            transitionRate += self.getTransitionRate(n, l, j,
                                                     nto, lto, jto,
                                                     temperature,
                                                     s=s)
        # sum over additional states
        for state in self.extraLevels:
            if (abs(j - state[2]) < 1.1) and \
                    (abs(state[1] - l) < 1.1) and (abs(state[1] - l) > 0.9):
                transitionRate += self.getTransitionRate(
                    n, l, j,
                    state[0], state[1], state[2],
                    temperature,
                    s=s
                )

        # add something small decay (1e-50) rate to prevent division by zero
        return 1. / (transitionRate + 1e-50)

    def getRadialCoupling(self, n, l, j, n1, l1, j1, s=0.5):
        """
            Returns radial part of the coupling between two states (dipole and
            quadrupole interactions only)

            Args:
                n1 (int): principal quantum number
                l1 (int): orbital angular momentum
                j1 (float): total angular momentum
                n2 (int): principal quantum number
                l2 (int): orbital angular momentum
                j2 (float): total angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float:  radial coupling strength (in a.u.), or zero for
                forbidden transitions in dipole and quadrupole approximation.

        """
        dl = abs(l - l1)
        if (dl == 1 and abs(j - j1) < 1.1):
            return self.getRadialMatrixElement(n, l, j, n1, l1, j1,
                                               s=s)
        elif (dl == 0 or dl == 1 or dl == 2) and(abs(j - j1) < 2.1):
            # quadrupole coupling
            # return 0.
            return self.getQuadrupoleMatrixElement(n, l, j, n1, l1, j1,
                                                   s=s)
        else:
            # neglect octopole coupling and higher
            return 0

    def getAverageSpeed(self, temperature):
        """
            Average (mean) speed at a given temperature

            Args:
                temperature (float): temperature (K)

            Returns:
                float: mean speed (m/s)
        """
        return sqrt(8. * C_k * temperature / (pi * self.mass))

    def _readHFSdata(self):
        c = self.conn.cursor()
        c.execute('''DROP TABLE IF EXISTS hfsDataAB''')
        c.execute('''SELECT COUNT(*) FROM sqlite_master
                 WHERE type='table' AND name='hfsDataAB';''')
        if (c.fetchone()[0] == 0):
            # create table
            c.execute('''CREATE TABLE IF NOT EXISTS hfsDataAB
                (n TINYINT UNSIGNED, l TINYINT UNSIGNED, j_x2 TINYINT UNSIGNED,
                hfsA DOUBLE, hfsB DOUBLE,
                errorA DOUBLE, errorB DOUBLE,
                typeOfSource TINYINT,
                comment TINYTEXT,
                ref TINYTEXT,
                refdoi TINYTEXT
                );''')
            c.execute('''CREATE INDEX compositeIndexHFS
                ON hfsDataAB (n,l,j_x2);''')
        self.conn.commit()

        if (self.hyperfineStructureData == ""):
            return 0  # no file specified for literature values

        try:
            fn = open(os.path.join(self.dataFolder,
                                   self.hyperfineStructureData), 'r')
            dialect = csv.Sniffer().sniff(fn.read(2024), delimiters=";,\t")
            fn.seek(0)
            data = csv.reader(fn, dialect, quotechar='"')

            literatureHFS = []
            count = 0
            for row in data:
                if count != 0:
                    # if not header
                    n = int(row[0])
                    l = int(row[1])
                    j = float(row[2])
                    A = float(row[3])
                    B = float(row[4])
                    errorA = float(row[5])
                    errorB = float(row[6])
                    typeOfSource = row[7]
                    comment = row[8]
                    ref = row[9]
                    refdoi = row[10]

                    literatureHFS.append([n, l, j * 2, A, B,
                                          errorA, errorB,
                                          typeOfSource,
                                          comment,
                                          ref, refdoi])
                count += 1

            fn.close()
            try:
                if count > 1:
                    c.executemany('''INSERT INTO hfsDataAB
                                        VALUES (?,?,?,?,?,
                                                ?, ?,
                                                ?,?,?,?)''',
                                  literatureHFS)
                    self.conn.commit()
            except sqlite3.Error as e:
                if count > 0:
                    print("Error while loading precalculated values "
                          "into the database")
                    print(e)
                    return

        except IOError as e:
            print("Error reading literature values File "
                      + self.hyperfineStructureData)
            print(e)



    def _readLiteratureValues(self):
        # clear previously saved results, since literature file
        # might have been updated in the meantime
        c = self.conn.cursor()
        c.execute('''DROP TABLE IF EXISTS literatureDME''')
        c.execute('''SELECT COUNT(*) FROM sqlite_master
                        WHERE type='table' AND name='literatureDME';''')
        if (c.fetchone()[0] == 0):
            # create table
            c.execute('''CREATE TABLE IF NOT EXISTS literatureDME
             (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED, j1_x2 TINYINT UNSIGNED,
             n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED, j2_x2 TINYINT UNSIGNED,
             dme DOUBLE,
             typeOfSource TINYINT,
             errorEstimate DOUBLE,
             comment TINYTEXT,
             ref TINYTEXT,
             refdoi TINYTEXT
            );''')
            c.execute('''CREATE INDEX compositeIndex
            ON literatureDME (n1,l1,j1_x2,n2,l2,j2_x2); ''')
        self.conn.commit()

        if (self.literatureDMEfilename == ""):
            return 0  # no file specified for literature values

        try:
            fn = open(os.path.join(self.dataFolder,
                                   self.literatureDMEfilename), 'r')
            dialect = csv.Sniffer().sniff(fn.read(2024), delimiters=";,\t")
            fn.seek(0)
            data = csv.reader(fn, dialect, quotechar='"')

            literatureDME = []

            # i=0 is header
            i = 0
            for row in data:
                if i != 0:
                    n1 = int(row[0])
                    l1 = int(row[1])
                    j1 = float(row[2])
                    n2 = int(row[3])
                    l2 = int(row[4])
                    j2 = float(row[5])
                    if (
                        self.getEnergy(n1, l1, j1) > self.getEnergy(n2, l2, j2)
                    ):
                        temp = n1
                        n1 = n2
                        n2 = temp
                        temp = l1
                        l1 = l2
                        l2 = temp
                        temp = j1
                        j1 = j2
                        j2 = temp

                    # convered from reduced DME in J basis (symmetric notation)
                    # to radial part of dme as it is saved for calculated
                    # values
                    dme = float(row[6]) / (
                        (-1)**(int(l1 + 0.5 + j2 + 1.))
                        * sqrt((2. * j1 + 1.) * (2. * j2 + 1.))
                        * Wigner6j(j1, 1., j2, l2, 0.5, l1)
                        * (-1)**l1 * sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0))
                        * Wigner3j(l1, 1, l2, 0, 0, 0))

                    comment = row[7]
                    typeOfSource = int(row[8])  # 0 = experiment; 1 = theory
                    errorEstimate = float(row[9])
                    ref = row[10]
                    refdoi = row[11]

                    literatureDME.append([n1, l1, j1 * 2, n2, l2, j2 * 2, dme,
                                          typeOfSource, errorEstimate,
                                          comment, ref,
                                          refdoi])
                i += 1
            fn.close()

            try:
                if i > 1:
                    c.executemany('''INSERT INTO literatureDME
                                        VALUES (?,?,?,?,?,?,?,
                                                ?,?,?,?,?)''',
                                  literatureDME)
                    self.conn.commit()
            except sqlite3.Error as e:
                if i > 0:
                    print("Error while loading precalculated values "
                          "into the database")
                    print(e)
                    exit()

        except IOError as e:
            print("Error reading literature values File "
                  + self.literatureDMEfilename)
            print(e)

    def getLiteratureDME(self, n1, l1, j1, n2, l2, j2, s=0.5):
        """
            Returns literature information on requested transition.

            Args:
                n1,l1,j1: one of the states we are coupling
                n2,l2,j2: the other state to which we are coupling

            Returns:
                bool, float, [int,float,string,string,string]:
                    hasLiteratureValue?, dme, referenceInformation

                    **If Boolean value is True**, a literature value for
                    dipole matrix element was found and reduced DME in J basis
                    is returned as the number. The third returned argument
                    (array) contains additional information about the
                    literature value in the following order [ typeOfSource,
                    errorEstimate , comment , reference, reference DOI]
                    upon success to find a literature value for dipole matrix
                    element:

                    * typeOfSource=1 if the value is theoretical \
                        calculation; otherwise, if it is experimentally \
                        obtained value typeOfSource=0
                    * comment details where within the publication the \
                        value can be found
                    * errorEstimate is absolute error estimate
                    * reference is human-readable formatted reference
                    * reference DOI provides link to the publication.

                    **Boolean value is False**, followed by zero and an empty
                    array if no literature value for dipole matrix element is
                    found.

            Note:
                The literature values are stored in /data folder in
                <element name>_literature_dme.csv files as a ; separated
                values. Each row in the file consists of one literature entry,
                that has information in the following order:

                * n1
                * l1
                * j1
                * n2
                * l2
                * j2
                * dipole matrix element reduced l basis (a.u.)
                * comment (e.g. where in the paper value appears?)
                * value origin: 1 for theoretical; 0 for experimental values
                * accuracy
                * source (human readable formatted citation)
                * doi number (e.g. 10.1103/RevModPhys.82.2313 )

                If there are several values for a given transition, program
                outputs the value that has smallest error (under column
                accuracy). The list of values can be expanded - every time
                program runs this file is read and the list is parsed again
                for use in calculations.

        """

        if (self.getEnergy(n1, l1, j1) > self.getEnergy(n2, l2, j2)):
            temp = n1
            n1 = n2
            n2 = temp
            temp = l1
            l1 = l2
            l2 = temp
            temp = j1
            j1 = j2
            j2 = temp

        # is there literature value for this DME? If there is,
        # use the best one (wit the smallest error)

        j1_x2 = 2 * j1
        j2_x2 = 2 * j2
        c = self.conn.cursor()
        c.execute('''SELECT dme, typeOfSource,
                     errorEstimate ,
                     comment ,
                     ref,
                     refdoi FROM literatureDME WHERE
                     n1= ? AND l1 = ? AND j1_x2 = ? AND
                     n2 = ? AND l2 = ? AND j2_x2 = ?
                     ORDER BY errorEstimate ASC''',
                       (n1, l1, j1_x2, n2, l2, j2_x2))

        answer = c.fetchone()
        if (answer):
            # we did found literature value
            return True, answer[0], [answer[1], answer[2], answer[3],
                                     answer[4], answer[5]]

        # if we are here, we were unsucessfull in literature search
        # for this value
        return False, 0, []

    def getZeemanEnergyShift(self, l, j, mj, magneticFieldBz, s=0.5):
        r"""
            Retuns linear (paramagnetic) Zeeman shift.

            :math:`\mathcal{H}_P=\frac{\mu_B B_z}{\hbar}(\hat{L}_{\rm z}+\
            g_{\rm S}S_{\rm z})`

            Args:
                l (int): orbital angular momentum
                j (float): total angular momentum
                mj (float): projection of total angular momentum alon z-axis
                magneticFieldBz (float): applied magnetic field (alon z-axis
                    only) in units of T (Tesla)
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: energy offset of the state (in J)
        """
        prefactor = physical_constants["Bohr magneton"][0] * magneticFieldBz
        gs = - physical_constants["electron g factor"][0]
        sumOverMl = 0

        for ml in np.linspace(mj - s, mj + s, round(2 * s + 1)):
            if abs(ml) <= l + 0.1:
                ms = mj - ml
                sumOverMl += (ml + gs * ms) * \
                    abs(CG(l, ml, s, ms, j, mj))**2
        return prefactor * sumOverMl

    def _getRadialDipoleSemiClassical(self, n1, l1, j1, n2, l2, j2,
                                      s=0.5):

        # get the effective principal number of both states
        nu = np.sqrt( - self.scaledRydbergConstant
                     / self.getEnergy(n1, l1, j1, s=s))
        nu1 = np.sqrt( - self.scaledRydbergConstant
                      / self.getEnergy(n2, l2, j2, s=s))

        # get the parameters required to calculate the sum
        l_c = (l1 + l2 + 1.) / 2.
        nu_c = sqrt(nu * nu1)

        delta_nu = nu - nu1
        delta_l = l2 - l1

        # I am not sure if this correct

        gamma = (delta_l * l_c) / nu_c

        if delta_nu == 0:
            g0 = 1
            g1 = 0
            g2 = 0
            g3 = 0
        else:

            g0 = (1. / (3. * delta_nu)) * (
                angerj(delta_nu - 1., - delta_nu)
                - angerj(delta_nu + 1, - delta_nu))
            g1 = -(1. / (3. * delta_nu)) * (
                angerj(delta_nu - 1., - delta_nu)
                + angerj(delta_nu + 1, -delta_nu))
            g2 = g0 - np.sin(np.pi * delta_nu) / (np.pi * delta_nu)
            g3 = (delta_nu / 2.) * g0 + g1

        radial_ME = (3 / 2) * nu_c**2 * (1 - (l_c / nu_c)**(2))**0.5 * \
            (g0 + gamma * g1 + gamma**2 * g2 + gamma**3 * g3)
        return float(radial_ME)

    def _getRadialQuadrupoleSemiClassical(self, n1, l1, j1, n2, l2, j2,
                                          s=0.5):

        dl = abs(l2 - l1)

        nu = n1 - self.getQuantumDefect(n1, l1, j1, s=s)
        nu1 = n2 - self.getQuantumDefect(n2, l2, j2, s=s)

        # get the parameters required to calculate the sum
        l_c = (l1 + l2 + 1.) / 2.
        nu_c = np.sqrt(nu * nu1)

        delta_nu = nu - nu1
        delta_l = l2 - l1

        gamma = (delta_l * l_c) / nu_c

        if delta_nu == 0:
            q = np.array([1, 0, 0, 0])
        else:

            g0 = (1. / (3. * delta_nu)) * (
                angerj(delta_nu - 1., - delta_nu)
                - angerj(delta_nu + 1, -delta_nu))
            g1 = -(1. / (3. * delta_nu)) * (
                angerj(delta_nu - 1., - delta_nu)
                + angerj(delta_nu + 1, -delta_nu))

            q = np.zeros((4,))
            q[0] = -(6. / (5. * delta_nu)) * g1
            q[1] = -(6. / (5. * delta_nu)) * g0 + (6. / 5.) * \
                np.sin(np.pi * delta_nu) / (np.pi * delta_nu**2)
            q[2] = -(3. / 4.) * (6. / (5. * delta_nu) * g1 + g0)
            q[3] = 0.5 * (delta_nu * 0.5 * q[0] + q[1])

        sm = 0

        if dl == 0:
            quadrupoleElement = (5. / 2.) * nu_c**4 * \
                (1. - (3. * l_c**2) / (5 * nu_c**2))
            for p in range(0, 2, 1):
                sm += gamma**(2 * p) * q[2 * p]
            return quadrupoleElement * sm

        elif dl == 2:
            quadrupoleElement = (5. / 2.) * nu_c**4 * (
                1 - (l_c + 1) ** 2 / (nu_c**2))**0.5 * (1 - (l_c + 2)**2
                                                        / (nu_c**2))**0.5
            for p in range(0, 4):
                sm += gamma**(p) * q[p]
            return quadrupoleElement * sm
        else:
            return 0

    # Additional AMO Functions
    def getHFSCoefficients(self, n, l, j, s=None):
        """
            Returns hyperfine splitting coefficients for state :math:`n`,
            :math:`l`, :math:`j`.

            Args:
                n (int): principal quantum number
                l (int): orbital angular momentum
                j (float): total angular momentum
                s (float): (optional) total spin momentum

            Returns:
                float: A,B hyperfine splitting constants (in Hz)
        """
        c = self.conn.cursor()

        c.execute('''SELECT hfsA, hfsB FROM hfsDataAB WHERE
            n= ? AND l = ? AND j_x2 = ?''', (n, l, j*2))
        answer = c.fetchone()
        if (answer):
            # we did found literature value  (A and B respectively)
            return answer[0], answer[1]
        else:
            raise ValueError("There is no data available on HFS structure"
                             " of %s state" % printStateString(n,l,j,s=s))

    def _reducedMatrixElementFJ(self, j1, f1, j2, f2):

        sph = 0.0
        if((abs(f2 - f1) < 2) & (int(abs(j2 - j1)) < 2)):
            # Reduced Matrix Element <f||er||f'> in units of reduced matrix element <j||er||j'>
            sph = (-1.0)**(j1 + self.I + f2 + 1.0) * ((2. * f1 + 1) * (2 * f2 + 1)) ** \
                0.5 * Wigner6j(f1, 1, f2, j2, self.I, j1)

        return sph

    def getSphericalDipoleMatrixElement(self, j1, mj1, j2, mj2, q):
        # Spherical Component of Angular Matrix Element in units of reduced matrix element <j||er||j'>
        return (- 1)**(j1 - mj1) * Wigner3j(j1, 1, j2, -mj1, -q, mj2)

    def getSphericalMatrixElementHFStoFS(self, j1, f1, mf1, j2, mj2, q):
        r"""
             Spherical matrix element for transition from hyperfine  resolved state
             to unresolved fine-structure state
             :math:`\langle f,m_f \vert\mu_q\vert j',m_j'\rangle`
             in units of :math:`\langle j\vert\vert\mu\vert\vert j'\rangle`


             Args:
                 j1, f1, mf1: total orbital,
                     fine basis (total atomic) angular momentum,
                     and projection of total angular momentum for state 1
                 j2, mj2: total orbital,
                     fine basis (total atomic) angular momentum,
                     and projection of total orbital angular momentum for state 2
                 q (int): specifies transition that the driving field couples to,
                     +1, 0 or -1 corresponding to driving :math:`\sigma^+`,
                     :math:`\pi` and :math:`\sigma^-` transitions respectively.
                 s (float): optional, total spin angular momentum of state.
                     By default 0.5 for Alkali atoms.

             Returns:
                 float: spherical dipole matrix element( :math:`\langle j\vert\vert\mu\vert\vert j'\rangle`)
         """
        mf2 = mf1 + q
        mI = mf2 - mj2
        sph = 0.0
        if(abs(mI) <= self.I):
            for f2 in np.arange(max(self.I - j2, abs(mf2), f1 - 1), 1 + min(self.I + j2, f1 + 1)):
            	#Enforce Triangle Rule
            	if abs(j2-self.I)<= f2: 
	                # CG multiplied by <j1 f1 mf1|er_q|j2 f2 mf2> in units of <j1 || er || j2 >
	                sph += CG(j2, mj2, self.I, mI, f2, mf2) \
    	                * self.getSphericalDipoleMatrixElement(f1, mf1, f2, mf2, q) \
        	            * self._reducedMatrixElementFJ(j1, f1, j2, f2)

        return sph

    def getDipoleMatrixElementHFStoFS(self, n1, l1, j1, f1, mf1, n2, l2, j2, mj2, q, s=0.5):
        r"""
            Dipole matrix element for transition from hyperfine  resolved state
            to unresolved fine-structure state
            :math:`\langle n_1 l_1 j_1 f_1 m_{f_1} |e\mathbf{r}|\
            n_2 l_2 j_2 m_{j_2}\rangle`
            in units of :math:`a_0 e`

            For hyperfine resolved transitions, the dipole matrix element is
            :math:`\langle n_1,\ell_1,j_1,f_1,m_{f1} |  \
            \mathbf{\hat{r}}\cdot \mathbf{\varepsilon}_q  \
            | n_2,\ell_2,j_2,f_2,m_{f2} \rangle = (-1)^{f_1-m_{f1}} \
            \left( \
            \begin{matrix} \
            f_1 & 1 & f_2 \\ \
            -m_{f1} & q & m_{f2} \
            \end{matrix}\right) \
            \langle n_1 \ell_1 j_1 f_1|| r || n_2 \ell_2 j_2 f_2 \rangle,` where
            :math:`\langle n_1 \ell_1 j_1 f_1 ||r|| n_2 \ell_2 j_2 f_2 \rangle \
            = (-1)^{j_1+I+F_2+1}\sqrt{(2f_1+1)(2f_2+1)} ~ \
            \left\{ \begin{matrix}\
            F_1 & 1 & F_2 \\ \
            j_2 & I & j_1 \
            \end{matrix}\right\}~ \
            \langle n_1 \ell_1 j_1||r || n_2 \ell_2 j_2 \rangle.`


            Args:
                n1. l1, j1, f1, mf1: principal, orbital, total orbital,
                    fine basis (total atomic) angular momentum,
                    and projection of total angular momentum for state 1
                n2. l2, j2, mj2: principal, orbital, total orbital,
                    fine basis (total atomic) angular momentum,
                    and projection of total orbital angular momentum for state 2
                q (int): specifies transition that the driving field couples to,
                    +1, 0 or -1 corresponding to driving :math:`\sigma^+`,
                    :math:`\pi` and :math:`\sigma^-` transitions respectively.
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: dipole matrix element( :math:`a_0 e`)
        """
        return self.getSphericalMatrixElementHFStoFS(j1, f1, mf1, j2, mj2, q) \
            * self.getReducedMatrixElementJ(n1, l1, j1, n2, l2, j2, s=s)

    def getMagneticDipoleMatrixElementHFS(self, l, j, f1, mf1, f2, mf2, q, s=0.5):
        r"""

          Magnetic dipole matrix element :math:`\langle f_1,m_{f_1} \vert \mu_q \vert f_2,m_{f_2}\rangle` \for transitions from :math:`\vert f_1,m_{f_1}\rangle\rightarrow\vert f_2,m_{f_2}\rangle` within the same :math:`n,\ell,j` state in units of :math:`\mu_B B_q`.

            The magnetic dipole matrix element is given by
            :math:`\langle f_1,m_{f_1}\vert \mu_q \vert f_2,m_{f_2}\rangle = g_J \mu_B B_q (-1)^{f_2+j+I+1+f_1-m_{f_1}} \sqrt{(2f_1+1)(2f_2+1)j(j+1)(2j+1)} \begin{pmatrix}f_1&1&f_2\\-m_{f_1} & -q & m_{f_2}\end{pmatrix}                \begin{Bmatrix}f_1&1&f_2\\j & I & j\end{Bmatrix}`



            Args:
                l, j, f1, mf1: orbital, total orbital,
                    fine basis (total atomic) angular momentum,total anuglar momentum
                    and projection of total angular momentum for state 1
                f2,mf2: principal, orbital, total orbital,
                    fine basis (total atomic) angular momentum,
                    and projection of total orbital angular momentum for state 2
                q (int): specifies transition that the driving field couples to,
                    +1, 0 or -1 corresponding to driving :math:`\sigma^+`,
                    :math:`\pi` and :math:`\sigma^-` transitions respectively.
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: magnetic dipole matrix element (in units of :math:`\mu_BB_q`)
        """
        return self.getLandegj(l, j, s) * (-1)**(f2 + j + self.I + 1)\
            * np.sqrt((2 * f1 + 1) * (2 * f2 + 1) * j * (j + 1) * (2 * j + 1)) \
            * self.getSphericalDipoleMatrixElement(f1, mf1, f2, mf2, q) \
            * Wigner6j(f1, 1, f2, j, self.I, j)

    def getLandegj(self, l, j, s=0.5):
        r"""
            Lande g-factor :math:`g_J\simeq 1+\frac{j(j+1)+s(s+1)-l(l+1)}{2j(j+1)}`

            Args:
                l (float): orbital angular momentum
                j (float): total orbital angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: Lande g-factor ( :math:`g_J`)
        """
        return 1.0 + (j * (j + 1.0) + s * (s + 1.0) - l * (l + 1.0)) / (2.0 * j * (j + 1.0))

    def getLandegjExact(self, l, j, s=0.5):
        r"""
            Lande g-factor :math:`g_J=g_L\frac{j(j+1)-s(s+1)+l(l+1)}{2j(j+1)}+g_S\frac{j(j+1)+s(s+1)-l(l+1)}{2j(j+1)}`

            Args:
                l (float): orbital angular momentum
                j (float): total orbital angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: Lande g-factor ( :math:`g_J`)
        """
        return self.gL * (j * (j + 1.0) - s * (s + 1.0) + l * (l + 1.0)) / (2.0 * j * (j + 1.0)) \
            + self.gS * (j * (j + 1.0) + s * (s + 1.0) - l *
                         (l + 1.0)) / (2.0 * j * (j + 1.0))

    def getLandegf(self, l, j, f, s=0.5):
        r"""
            Lande g-factor :math:`g_F\simeq g_J\frac{f(f+1)-I(I+1)+j(j+1)}{2f(f+1)}`

            Args:
                l (float): orbital angular momentum
                j (float): total orbital angular momentum
                f (float): total atomic angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: Lande g-factor ( :math:`g_F`)
        """
        gf = self.getLandegj(l, j, s) * (f * (f + 1.0) - self.I *
                                         (self.I + 1.0) + j * (j + 1.0)) / (2.0 * f * (f + 1.0))
        return gf

    def getLandegfExact(self, l, j, f, s=0.5):
        r"""
            Lande g-factor :math:`g_F`
            :math:`g_F=g_J\frac{f(f+1)-I(I+1)+j(j+1)}{2f(f+1)}+g_I\frac{f(f+1)+I(I+1)-j(j+1)}{2f(f+1)}`

            Args:
                l (float): orbital angular momentum
                j (float): total orbital angular momentum
                f (float): total atomic angular momentum
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: Lande g-factor ( :math:`g_F`)
        """
        gf = self.getLandegjExact(l, j, s) * (f * (f + 1) - self.I * (self.I + 1) + j * (j + 1.0)) / (2 * f * (f + 1.0)) \
            + self.gI * (f * (f + 1.0) + self.I * (self.I + 1.0) -
                         j * (j + 1.0)) / (2.0 * f * (f + 1.0))
        return gf

    def getHFSEnergyShift(self, j, f, A, B=0, s=0.5):
        r"""
             Energy shift of HFS from centre of mass :math:`\Delta E_\mathrm{hfs}`

            :math:`\Delta E_\mathrm{hfs} = \frac{A}{2}K+B\frac{\frac{3}{2}K(K+1)-2I(I+1)J(J+1)}{2I(2I-1)2J(2J-1)}`

            where :math:`K=F(F+1)-I(I+1)-J(J+1)`

            Args:
                j (float): total orbital angular momentum
                f (float): total atomic angular momentum
                A (float): HFS magnetic dipole constant
                B (float): HFS magnetic quadrupole constant
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: Energy shift ( :math:`\Delta E_\mathrm{hfs}`)
        """
        K = f * (f + 1.0) - self.I * (self.I + 1.0) - j * (j + 1.0)
        Ehfs = A / 2.0 * K
        if abs(B) > 0:
            Ehfs += B * (3.0 / 2.0 * K * (K + 1) - 2.0 * self.I * (self.I + 1.0) * j *
                         (j + 1.0)) / (2.0 * self.I * (2.0 * self.I - 1.0) * 2.0 * j * (2.0 * j - 1))

        return Ehfs

    def getBranchingRatio(self, jg, fg, mfg, je, fe, mfe, s=0.5):
        r"""
             Branching ratio for decay from :math:`\vert j_e,f_e,m_{f_e} \rangle \rightarrow \vert j_g,f_g,m_{f_g}\rangle`

                :math:`b = \displaystyle\sum_q (2j_e+1)\left(\begin{matrix}f_1 & 1 & f_2 \\-m_{f1} & q & m_{f2}\end{matrix}\right)^2\vert \langle j_e,f_e\vert \vert er \vert\vert j_g,f_g\rangle\vert^2/|\langle j_e || er || j_g \rangle |^2`

            Args:
                jg, fg, mfg: total orbital, fine basis (total atomic) angular momentum,
                    and projection of total angular momentum for ground state
                je, fe, mfe: total orbital, fine basis (total atomic) angular momentum,
                and projection of total angular momentum for excited state
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: branching ratio
        """
        b = 0.0
        for q in np.arange(-1, 2):
            b += self.getSphericalDipoleMatrixElement(
                fg, mfg, fe, mfe, q)**2 * self._reducedMatrixElementFJ(jg, fg, je, fe)**2

        # Rescale
        return b * (2.0 * je + 1.0)

    def getSaturationIntensity(self, ng, lg, jg, fg, mfg, ne, le, je, fe, mfe, s=0.5):
        r"""
             Saturation Intensity :math:`I_\mathrm{sat}` for transition :math:`\vert j_g,f_g,m_{f_g}\rangle\rightarrow\vert j_e,f_e,m_{f_e}\rangle` in units of :math:`\mathrm{W}/\mathrm{m}^2`.

                :math:`I_\mathrm{sat} = \frac{c\epsilon_0\Gamma^2\hbar^2}{4\vert \epsilon_q\cdot\mathrm{d}\vert^2}`

            Args:
                ng, lg, jg, fg, mfg: total orbital, fine basis (total atomic) angular momentum,
                    and projection of total angular momentum for ground state
                ne, le, je, fe, mfe: total orbital, fine basis (total atomic) angular momentum,
                and projection of total angular momentum for excited state
                s (float): optional, total spin angular momentum of state.
                    By default 0.5 for Alkali atoms.

            Returns:
                float: Saturation Intensity in units of :math:`\mathrm{W}/\mathrm{m}^2`
        """
        q = mfe - mfg
        if abs(q) <= 1:
            d = self.getDipoleMatrixElementHFS(
                ng, lg, jg, fg, mfg, ne, le, je, fe, mfe, q) * C_e * physical_constants["Bohr radius"][0]
            Gamma = 1. / self.getStateLifetime(ne, le, je)
            Is = C_c * epsilon_0 * Gamma**2 * hbar**2 / (4.0 * d**2)
        else:
            raise ValueError("States not coupled")

        return Is

    def getSaturationIntensityIsotropic(self, ng, lg, jg, fg, ne, le, je, fe):
        r"""
             Isotropic Saturation Intensity :math:`I_\mathrm{sat}` for transition :math:`f_g\rightarrow f_e` averaged over all polarisations in units of :math:`\mathrm{W}/\mathrm{m}^2`.

                :math:`I_\mathrm{sat} = \frac{c\epsilon_0\Gamma^2\hbar^2}{4\vert \epsilon_q\cdot\mathrm{d}\vert^2}`

            Args:
                ng, lg, jg, fg, mfg: total orbital, fine basis (total atomic) angular momentum,
                    and projection of total angular momentum for ground state
                ne, le, je, fe, mfe: total orbital, fine basis (total atomic) angular momentum,
                and projection of total angular momentum for excited state

            Returns:
                float: Saturation Intensity in units of :math:`\mathrm{W}/\mathrm{m}^2`
        """
        d_iso_sq = 0.0
        for q in range(-1, 2):
            for mfg in range(-fg, fg + 1):
                d_iso_sq += self.getDipoleMatrixElementHFS(
                    ng, lg, jg, fg, mfg, ne, le, je, fe, mfg + q, q)**2

        # Avergage over (2fg+1) levels and 3 polarisationsand rescale
        d_iso_sq = d_iso_sq / 3.0 / \
            (2 * fg + 1) * (C_e * physical_constants["Bohr radius"][0])**2

        Gamma = 1. / self.getStateLifetime(ne, le, je)
        Is = C_c * epsilon_0 * Gamma**2 * hbar**2 / (4.0 * d_iso_sq)

        return Is

    def groundStateRamanTransition(self, Pa, wa, qa, Pb, wb, qb, Delta, f0, mf0, f1, mf1, ne, le, je):
        r"""
            Returns two-photon Rabi frequency :math:`\Omega_R`, differential AC Stark shift :math:`\Delta_\mathrm{AC}` and probability to scatter a photon during a :math:`\pi`-pulse :math:`P_\mathrm{sc}` for two-photon ground-state Raman transitions from :math:`\vert f_g,m_{f_g}\rangle\rightarrow\vert nL_{j_r} j_r,m_{j_r}\rangle` via an intermediate excited state :math:`n_e,\ell_e,j_e`.

                :math:`\Omega_R=\displaystyle\sum_{f_e,m_{f_e}}\frac{\Omega^a_{0\rightarrow f_e}\Omega^b_{1\rightarrow f_e}}{2(\Delta-\Delta_{f_e})},`

                :math:`\Delta_{\mathrm{AC}} = \displaystyle\sum_{f_e,m_{f_e}}\left[\frac{\vert\Omega^a_{0\rightarrow f_e}\vert^2-\vert\Omega^b_{1\rightarrow f_e}\vert^2}{4(\Delta-\Delta_{f_e})}+\frac{\vert\Omega^a_{1\rightarrow f_e}\vert^2}{4(\Delta+\omega_{01}-\Delta_{f_e})}-\frac{\vert\Omega^b_{0\rightarrow f_e}\vert^2}{4(\Delta-\omega_{01}-\Delta_{f_e})}\right],`

                :math:`P_\mathrm{sc} =\frac{\Gamma_e t_\pi}{2}\displaystyle\sum_{f_e,m_{f_e}}\left[\frac{\vert\Omega^a_{0\rightarrow f_e}\vert^2}{2(\Delta-\Delta_{f_e})^2}+\frac{\vert\Omega^b_{1\rightarrow f_e}\vert^2}{2(\Delta-\Delta_{f_e})^2}+\frac{\vert\Omega^a_{1\rightarrow f_e}\vert^2}{4(\Delta+\omega_{01}-\Delta_{f_e})^2}+\frac{\vert\Omega^b_{0\rightarrow f_e}\vert^2}{4(\Delta-\omega_{01}-\Delta_{f_e})^2}\right]`

                where :math:`\tau_\pi=\pi/\Omega_R`.

            .. figure:: ./GroundStateRaman.png
                :width: 250 px
                :alt: Schema of |0>-> -> |e> -> |1> transition
                :align: right

            Args:
                Pa:
                    power (W), of laser a :math:`\vert 0 \rangle\rightarrow\vert e\rangle`
                wa:
                    beam waist (m) of laser a :math:`\vert 0 \rangle\rightarrow\vert e\rangle`
                qa:
                    polarisation (+1, 0 or -1 corresponding to driving :math:`\sigma^+`, :math:`\pi` and :math:`\sigma^-`)
                    of laser a :math:`\vert 0 \rangle\rightarrow\vert e\rangle`
                Pb: power (W)  of laser b :math:`\vert 1 \rangle\rightarrow\vert e\rangle`
                wb: beam waist (m)  of laser b :math:`\vert 1 \rangle\rightarrow\vert e\rangle`
                qb: polarisation (+1, 0 or -1 corresponding to driving :math:`\sigma^+`, :math:`\pi` and :math:`\sigma^-`) of laser b :math:`\vert 1 \rangle\rightarrow\vert e\rangle`
                Delta : Detuning from excited state centre of mass (rad :math:`\mathrm{s}^{-1}`)
                f0,mf0: Lower hyperfine level
                f1,mf1: Upper hyperfine level
                ne, le, je: principal, orbital, total orbital quantum numbers of excited state

            Returns:
                float: Two-Photon Rabi frequency :math:`\Omega_R` (units :math:`\mathrm{rads}^{-1}`), differential AC Stark shift :math:`\Delta_\mathrm{AC}` (units :math:`\mathrm{rads}^{-1}`) and probability to scatter a photon during a :math:`\pi`-pulse :math:`P_\mathrm{sc}`
        """

        # Intensity/beam (W/m^2)
        Ia = 2.0 * Pa / (pi * wa**2)
        Ib = 2.0 * Pb / (pi * wb**2)
        # Electric field (V/m)
        Ea = np.sqrt(2.0 * Ia / (epsilon_0 * C_c))
        Eb = np.sqrt(2.0 * Ib / (epsilon_0 * C_c))
        # Reduced Matrix Element (au)
        ng = self.groundStateN
        lg = 0
        jg = 0.5
        rme_j = self.getReducedMatrixElementJ(ng, lg, jg, ne, le, je)
        # Rescale to (Cm)
        rme_j *= C_e * physical_constants["Bohr radius"][0]

        # Qubit level energy separation (rad s-1)
        [A, B] = self.getHFSCoefficients(ng, lg, jg)
        omega01 = (jg + self.I) * A * 2.0 * pi

        # Excited State Properties

        # Hyperfine Coefficients (Hz)
        [A, B] = self.getHFSCoefficients(ne, le, je)
        # Linewidth (rad s-1)
        Gamma = 1.0 / self.getStateLifetime(ne, le, je)

        # Initialise Output Variables
        OmegaR = np.zeros(np.shape(Delta))
        AC1 = np.zeros(np.shape(Delta))
        AC0 = np.zeros(np.shape(Delta))
        Pe = np.zeros(np.shape(Delta))

        # Loop over excited state energylevels
        for fe in range(int(abs(je - self.I)), int(1.0 + (je + self.I))):
            # Hyperfine energy shift (rad s-1)
            Ehfs = 2.0 * np.pi * self.getHFSEnergyShift(je, fe, A, B)
            for mfe in range(max(-fe, min(mf1, mf0) - 1), 1 + min(fe, max(mf1, mf0) + 1)):

                # Rabi frequency of each laser from each transition (rad s-1)
                Omaf0 = Ea * rme_j / hbar * self.getSphericalDipoleMatrixElement(f0, mf0, fe, mfe, qa) \
                    * self._reducedMatrixElementFJ(jg, f0, je, fe)

                Omaf1 = Ea * rme_j / hbar * self.getSphericalDipoleMatrixElement(f1, mf1, fe, mfe, qa) \
                    * self._reducedMatrixElementFJ(jg, f1, je, fe)

                Ombf0 = Eb * rme_j / hbar * self.getSphericalDipoleMatrixElement(f0, mf0, fe, mfe, qb) \
                    * self._reducedMatrixElementFJ(jg, f0, je, fe)

                Ombf1 = Eb * rme_j / hbar * self.getSphericalDipoleMatrixElement(f1, mf1, fe, mfe, qb) \
                    * self._reducedMatrixElementFJ(jg, f1, je, fe)

                # AC Stark shift on qubit states
                AC1 += Ombf1**2 / (4 * (Delta - Ehfs)) + \
                    Omaf1**2 / (4 * (Delta + omega01 - Ehfs))
                AC0 += Omaf0**2 / (4 * (Delta - Ehfs)) + \
                    Ombf0**2 / (4 * (Delta - omega01 - Ehfs))

                # Two-Photon Rabi Frequency
                OmegaR += Omaf0 * Ombf1 / (2 * (Delta - Ehfs))

                # Excitated state population Pe
                Pe += 0.5 * Omaf0**2 / (2 * (Delta - Ehfs)**2) + 0.5 * Ombf1**2 / (2 * (Delta - Ehfs)**2) \
                    + 0.5 * Omaf1**2 / (2 * (Delta + omega01 - Ehfs)**2) + \
                    0.5 * Ombf0**2 / (2 * (Delta - omega01 - Ehfs)**2)

        # Total Differential Shift
        AC = AC0 - AC1

        # Pi-rotation time (s)
        tau_pi = pi / abs(OmegaR)
        # Spontaneous Emission Probability
        Psc = Gamma * tau_pi * Pe

        return OmegaR, AC, Psc

    def twoPhotonRydbergExcitation(self, Pp, wp, qp, Pc, wc, qc, Delta, fg, mfg, ne, le, je, nr, lr, jr, mjr):
        r"""
                Returns two-photon Rabi frequency :math:`\Omega_R`, ground AC Stark shift :math:`\Delta_{\mathrm{AC}_g}`, Rydberg state AC Stark shift :math:`\Delta_{\mathrm{AC}_r}` and probability to scatter a photon during a :math:`\pi`-pulse :math:`P_\mathrm{sc}` for two-photon  excitation from :math:`\vert f_h,m_{f_g}\rangle\rightarrow \vert j_r,m_{j_r}\rangle` via intermediate excited state

                    :math:`\Omega_R=\displaystyle\sum_{f_e,m_{f_e}}\frac{\Omega_p^{g\rightarrow f_e}\Omega_c^{f_e\rightarrow r}}{2(\Delta-\Delta_{f_e})}`

                    :math:`\Delta_{\mathrm{AC}_g} = \displaystyle\sum_{f_e,m_{f_e}}\frac{\vert\Omega_p^{g\rightarrow f_e}\vert^2}{4(\Delta-\Delta_{f_e})}`

                    :math:`\Delta_{\mathrm{AC}_r} = \displaystyle\sum_{f_e,m_{f_e}}\frac{\vert\Omega_p^{g\rightarrow f_e}\vert^2}{4(\Delta-\Delta_{f_e})}``

                    :math:`P_\mathrm{sc} = \frac{\Gamma_et_\pi}{2}\displaystyle\sum_{f_e,m_{f_e}}\left[\frac{\vert\Omega_p^{g\rightarrow f_e}\vert^2}{2(\Delta-\Delta_{f_e})^2}+\frac{\vert\Omega_c^{f_e\rightarrow r}\vert^2}{2(\Delta-\Delta_{f_e})^2}\right]`

                where :math:`\tau_\pi=\pi/\Omega_R`.

                .. figure:: ./twophotonexcitation.png
                    :width: 150 px
                    :alt: Schema of |g-> -> |e> -> |r> transition
                    :align: right

                Args:
                    Pp: power (W) of probe laser :math:`\vert g \rangle\rightarrow\vert e\rangle`
                    wp: beam waist (m) of probe laser :math:`\vert g \rangle\rightarrow\vert e\rangle`
                    qp: polarisation (+1, 0 or -1 corresponding to driving :math:`\sigma^+`,:math:`\pi` and :math:`\sigma^-`) of probe laser :math:`\vert g \rangle\rightarrow\vert e\rangle`
                    Pb: power (W) of coupling laser :math:`\vert e\rangle\rightarrow\vert r\rangle`
                    wb: beam waist (m) of coupling laser :math:`\vert e\rangle\rightarrow\vert r\rangle`
                    qb: polarisation (+1, 0 or -1 corresponding to driving :math:`\sigma^+`,:math:`\pi` and :math:`\sigma^-`) of coupling laser :math:`\vert e\rangle\rightarrow\vert r\rangle`
                    Delta : Detuning from excited state centre of mass (rad s:math:`^{-1}`)
                    fg: ground state hyperfine state
                    mfg: projection of ground state hyperfine state
                    f1,mf1: upper hyperfine state
                    ne: principal quantum numbers of excited state 
                    le: orbital angular momentum of excited state 
                    je: total angular momentum of excited state
                    nr: principal quantum number of target Rydberg state
                    lr: orbital angular momentum of target Rydberg state
                    jr: total angular momentum of target Rydberg state
                    mjr: projection of total angular momenutm of target Rydberg state

                Returns:
                    float: Two-Photon Rabi frequency :math:`\Omega_R` (units :math:`\mathrm{rads}^{-1}`),
                    ground-state AC Stark shift :math:`\Delta_{\mathrm{AC}_g}` (units :math:`\mathrm{rads}^{-1}`) Rydberg-state AC Stark shift :math:`\Delta_{\mathrm{AC}_r}` (units :math:`\mathrm{rads}^{-1}`) and probability to scatter a photon during a :math:`\pi`-pulse :math:`P_\mathrm{sc}`
        """

        # Intensity/beam (W/m^2)
        Ip = 2.0 * Pp / (pi * wp**2)
        Ic = 2.0 * Pc / (pi * wc**2)
        # Electric field (V/m)
        Ep = np.sqrt(2.0 * Ip / (epsilon_0 * C_c))
        Ec = np.sqrt(2.0 * Ic / (epsilon_0 * C_c))

        # Excited State Properties

        # Reduced Matrix Element (au)
        ng = self.groundStateN
        lg = 0
        jg = 0.5
        rme_j = self.getReducedMatrixElementJ(ng, lg, jg, ne, le, je)

        # Rescale to (Cm)
        rme_j *= C_e * physical_constants["Bohr radius"][0]

        # Hyperfine Coefficients (Hz)
        [A, B] = self.getHFSCoefficients(ne, le, je)
        # Linewidth (rad s-1)
        Gamma = 1.0 / self.getStateLifetime(ne, le, je)

        # Rydberg State Reduced Matrix Element (au)
        rme_jRyd = self.getReducedMatrixElementJ(ne, le, je, nr, lr, jr)
        # Rescale to (Cm)
        rme_jRyd *= C_e * physical_constants["Bohr radius"][0]

        # Initialise Output Variables
        OmegaR = np.zeros(np.shape(Delta))
        ACg = np.zeros(np.shape(Delta))
        ACr = np.zeros(np.shape(Delta))
        Pe = np.zeros(np.shape(Delta))

        # Loop over excited state energylevels
        for fe in range(int(abs(je - self.I)), 1 + int(je + self.I)):
            # Hyperfine energy shift (rad s-1)
            Ehfs = 2.0 * np.pi * self.getHFSEnergyShift(je, fe, A, B)
            # range(max(-fe,min(mf1,mf0)-1),1+min(fe,max(mf1,mf0)+1)):
            for mfe in range(-fe, fe + 1):

                # Probe Rabi Frequency (rad s-1)
                OmP = Ep * rme_j / hbar * self.getSphericalDipoleMatrixElement(fg, mfg, fe, mfe, qp) \
                    * self._reducedMatrixElementFJ(jg, fg, je, fe)
                # Coupling Rabi Frequency (rad s-1)
                OmC = Ec * rme_jRyd / hbar * \
                    self.getSphericalMatrixElementHFStoFS(
                        je, fe, mfe, jr, mjr, qc)

                # AC Stark shift on ground state (rad s-1)
                ACg += (OmP**2) / (4 * (Delta - Ehfs))
                # AC Stark shift on Rydberg state (rad s-1)
                ACr += (OmC**2) / (4 * (Delta - Ehfs))

                # Two-Photon Rabi Frequency (rad s-1)
                OmegaR += OmP * OmC / (2 * (Delta - Ehfs))

                # Excitated state population Pe
                Pe += 0.5 * (OmP**2 + OmC**2) / (2 * (Delta - Ehfs)**2)

        # Pi-rotation time (s)
        tau_pi = pi / abs(OmegaR)
        # Spontaneous Emission Probability
        Psc = Gamma * tau_pi * Pe

        return OmegaR, ACg, ACr, Psc

    def _spinMatrices(self, j):
        # SPINMATRICES Generates spin-matrices for spin S
        #   [Sx,Sy,Sz]=SPINMATRICES(S) returns the Sx,Sy,Sz spin
        #   matrices calculated using raising and lowering operators
        mj = -np.arange(-j + 1, j + 1)
        jm = np.sqrt(j * (j + 1) - mj * (mj + 1))
        Jplus = np.matrix(np.diag(jm, 1))  # Raising Operator
        Jminus = np.matrix(np.diag(jm, -1))  # Lowering Operator
        Jx = (Jplus + Jminus) / 2.0
        Jy = (-Jplus + Jminus) * 1j / 2.0
        Jz = (Jplus * Jminus - Jminus * Jplus) / 2.0
        # J2=Jx**2+Jy**2+Jz**2
        return Jx, Jy, Jz

    def breitRabi(self, n, l, j, B):
        r"""
             Returns exact Zeeman energies math:`E_z` for states :math:`\vert F,m_f\rangle` in the :math:`\ell,j` manifold via exact diagonalisation of the Zeeman interaction :math:`\mathcal{H}_z` and the hyperfine interaction :math:`\mathcal{H}_\mathrm{hfs}` given by equations

                    :math:`\mathcal{H}_Z=\frac{\mu_B}{\hbar}(g_J J_z+g_I I_z)B_z`

                and

                    :math:`\mathcal{H}_\mathrm{hfs}=A_\mathrm{hfs}I\cdot J + B_\mathrm{hfs}\frac{3(I\cdot J)^2+3/2 I\cdot J -I^2J^2}{2I(2I-1)2J(2J-1)}`.

            Args:
                n,l,j: principal,orbital, total orbital quantum numbers
                B: Magnetic Field (units T)

            Returns:
                float: State energy :math:`E_z` in SI units (Hz), state f, state mf
        """

        Ahfs, Bhfs = self.getHFSCoefficients(n, l, j)

        # Bohr Magneton
        uB = physical_constants["Bohr magneton in Hz/T"][0]

        # Define Spin Matrices
        N = int((2 * j + 1) * (2 * self.I + 1))
        [jx, jy, jz] = self._spinMatrices(j)
        ji = np.eye(int(2.0 * j + 1.0))
        [ix, iy, iz] = self._spinMatrices(self.I)
        ii = np.eye(int(2.0 * self.I + 1.0))

        # Calculate Tensor Products
        Jx = np.kron(jx, ii)
        Jy = np.kron(jy, ii)
        Jz = np.kron(jz, ii)
        Ix = np.kron(ji, ix)
        Iy = np.kron(ji, iy)
        Iz = np.kron(ji, iz)
        J2 = Jx**2 + Jy**2 + Jz**2
        I2 = Ix**2 + Iy**2 + Iz**2
        IJ = Ix * Jx + Iy * Jy + Iz * Jz
        # F Basis
        Fx = Jx + Ix
        Fy = Jy + Iy
        Fz = Jz + Iz
        F2 = Fx**2 + Fy**2 + Fz**2

        # Hyperfine Interaction
        Hhfs = Ahfs * IJ
        if(Bhfs != 0):
            Hhfs += Bhfs * (3 * IJ * IJ + 3 / 2 * IJ - I2 * J2) / \
                (2 * self.I * (2 * self.I - 1) * 2 * j * (2 * j - 1))

        # Zeeman Interaction
        Hz = uB * (self.getLandegjExact(l, j) * Jz + self.gI * Iz)

        # Initialise Output
        en = np.zeros([B.size, N])

        ctr = -1
        for b in B:
            ctr = ctr + 1
            eVal, eVec = eigh(Hhfs + b * Hz)
            en[ctr, :] = eVal

        # Determine States
        eVal, eVec = eigh(Hhfs + 1e-4 * Hz)
        eVec = np.matrix(eVec)
        f = np.zeros(N)
        mf = np.zeros(N)
        for ctr in range(N):
            f2 = eVec[:, ctr].conj().T * F2 * eVec[:, ctr]
            f[ctr] = np.round(
                1 / 2 * (-1 + np.sqrt(1 + 4 * np.real(f2[0, 0]))))
            m = eVec[:, ctr].conj().T * Fz * eVec[:, ctr]
            mf[ctr] = np.round(np.real(m[0, 0]))

        return en, f, mf

### JDP EDITS FINISH ###


def NumerovBack(innerLimit, outerLimit, kfun, step, init1, init2):
    """
        Full Python implementation of Numerov integration

        Calculates solution function :math:`rad(r)` with descrete step in
        :math:`r` size of `step`, integrating from `outerLimit` towards the
        `innerLimit` (from outside, inwards) equation
        :math:`\\frac{\\mathrm{d}^2 rad(r)}{\\mathrm{d} r^2} = \
            kfun(r)\\cdot rad(r)`.




        Args:
            innerLimit (float): inner limit of integration
            outerLimit (flaot): outer limit of integration
            kfun (function(double)): pointer to function used in equation (see
                longer explanation above)
            step: descrete step size for integration
            init1 (float): initial value, `rad`(`outerLimit`+`step`)
            init2 (float): initial value,
                `rad`(`outerLimit`+:math:`2\\cdot` `step`)

        Returns:
            numpy array of float , numpy array of float, int : :math:`r` (a.u),
            :math:`rad(r)`;

        Note:
            Returned function is not normalized!

        Note:
            If :obj:`AlkaliAtom.cpp_numerov` swich is set to True (default),
            much faster C implementation of the algorithm will be used instead.
            That is recommended option. See documentation installation
            instructions for more details.

    """

    br = int((sqrt(outerLimit) - sqrt(innerLimit)) / step)
    # integrated wavefunction R(r)*r^{3/4}
    sol = np.zeros(br, dtype=np.dtype('d'))
    # radial coordinate for integration \sqrt(r)
    rad = np.zeros(br, dtype=np.dtype('d'))

    br = br - 1
    x = sqrt(innerLimit) + step * (br - 1)
    sol[br] = (2. * (1. - 5.0 / 12.0 * step**2 * kfun(x)) * init1
               - (1. + 1. / 12.0 * step**2 * kfun(x + step)) * init2) /\
        (1 + 1 / 12.0 * step**2 * kfun(x - step))
    rad[br] = x

    x = x - step
    br = br - 1

    sol[br] = (2. * (1. - 5.0 / 12.0 * step**2 * kfun(x)) * sol[br + 1]
               - (1. + 1. / 12.0 * step**2 * kfun(x + step)) * init1) /\
        (1 + 1 / 12.0 * step**2 * kfun(x - step))
    rad[br] = x

    # check if the function starts diverging  before the innerLimit
    # -> in that case break integration earlier

    maxValue = 0.

    checkPoint = 0
    fromLastMax = 0

    while br > checkPoint:
        br = br - 1
        x = x - step
        sol[br] = (2. * (1. - 5.0 / 12.0 * step**2 * kfun(x)) * sol[br + 1]
                   - (1. + 1. / 12.0 * step**2 * kfun(x + step)) * sol[br + 2]
                   ) /\
            (1. + 1. / 12.0 * step**2 * kfun(x - step))
        rad[br] = x
        if abs(sol[br] * sqrt(x)) > maxValue:
            maxValue = abs(sol[br] * sqrt(x))
        else:
            fromLastMax += 1
            if fromLastMax > 50:
                checkPoint = br
    # now proceed with caution - checking if the divergence starts
    # - if it does, cut earlier

    divergencePoint = 0

    while (br > 0)and(divergencePoint == 0):
        br = br - 1
        x = x - step
        sol[br] = (2. * (1. - 5.0 / 12.0 * step**2 * kfun(x)) * sol[br + 1]
                   - (1. + 1. / 12.0 * step**2 * kfun(x + step)) * sol[br + 2]
                   ) /\
            (1. + 1. / 12.0 * step**2 * kfun(x - step))
        rad[br] = x
        if (divergencePoint == 0)and(abs(sol[br] * sqrt(x)) > maxValue):
            divergencePoint = br
            while (abs(sol[divergencePoint]) > abs(sol[divergencePoint + 1])) \
                    and (divergencePoint < checkPoint):
                divergencePoint += 1
            if divergencePoint > checkPoint:
                print("Numerov error")
                exit()

    br = divergencePoint
    while (br > 0):
        rad[br] = rad[br + 1] - step
        sol[br] = 0
        br -= 1

    # convert R(r)*r^{3/4} to  R(r)*r
    sol = np.multiply(sol, np.sqrt(rad))
    # convert \sqrt(r) to r
    rad = np.multiply(rad, rad)

    return rad, sol


def _atomLightAtomCoupling(n, l, j, nn, ll, jj, n1, l1, j1, n2, l2, j2,
                           atom1, atom2=None, s=0.5, s2=None):
    """
        Calculates radial part of atom-light coupling

        This function might seem redundant, since similar function exist for
        each of the atoms. Function that is not connected to specific
        atomic species is provided in order to provides route to implement
        inter-species coupling.
    """
    if atom2 is None:
        # if not explicitly inter-species, assume it's the same species
        atom2 = atom1
    if s2 is None:
        s2 = s

    # determine coupling
    dl = abs(l - l1)
    dj = abs(j - j1)
    c1 = 0
    if dl == 1 and (dj < 1.1):
        c1 = 1  # dipole couplings1
    elif (dl == 0 or dl == 2 or dl == 1) and(dj < 2.1):
        c1 = 2  # quadrupole coupling
    else:
        return False
    dl = abs(ll - l2)
    dj = abs(jj - j2)
    c2 = 0
    if dl == 1 and (dj < 1.1):
        c2 = 1  # dipole coupling
    elif (dl == 0 or dl == 2 or dl == 1) and(dj < 2.1):
        c2 = 2  # quadrupole coupling
    else:
        return False

    radial1 = atom1.getRadialCoupling(n, l, j, n1, l1, j1, s=s)
    radial2 = atom2.getRadialCoupling(nn, ll, jj, n2, l2, j2, s=s2)

    # TO-DO: check exponent of the Boht radius (from where it comes?!)

    coupling = C_e**2 / (4.0 * pi * epsilon_0) * radial1 * radial2 *\
        (physical_constants["Bohr radius"][0])**(c1 + c2)
    return coupling


# ================== Saving and loading calculations (START) ==================

def saveCalculation(calculation, fileName):
    """
    Saves calculation for future use.

    Saves :obj:`calculations_atom_pairstate.PairStateInteractions` and
    :obj:`calculations_atom_single.StarkMap`
    calculations in compact binary format in file named `filename`. It uses
    cPickle serialization library in Python, and also zips the final file.

    Calculation can be retrieved and used with :obj:`loadSavedCalculation`

    Args:
        calculation: class instance of calculations (instance of
            :obj:`calculations_atom_pairstate.PairStateInteractions`
            or :obj:`calculations_atom_single.StarkMap`)
            to be saved.
        fileName: name of the file where calculation will be saved

    Example:
        Let's suppose that we did the part of the
        :obj:`calculation_atom_pairstate.PairStateInteractions`
        calculation that involves generation of the interaction
        matrix. After that we can save the full calculation in a single file::

            calc = PairStateInteractions(Rubidium(),
                    60,0,0.5,
                    60,0,0.5,
                    0.5,0.5)
            calc.defineBasis(0,0, 5,5, 25.e9)
            calc.diagonalise(np.linspace(0.5,10.0,200),150)
            saveCalculation(calc, "mySavedCalculation.pkl")

        Then, at a later time, and even on the another machine, we can load
        that file and continue with calculation. We can for example explore
        the calculated level diagram::

            calc = loadSavedCalculation("mySavedCalculation.pkl")
            calc.plotLevelDiagram()
            calc.showPlot()
            rvdw = calc.getVdwFromLevelDiagram(0.5,14,
                minStateContribution=0.5,
                showPlot = True)

        Or, we can do additional matrix diagonalization, in some new range,
        then and find C6 by fitting the obtained level diagram::

            calc = loadSavedCalculation("mySavedCalculation.pkl")
            calc.diagonalise(np.linspace(3,6.0,200),20)
            calc.getC6fromLevelDiagram(3,6.0,showPlot=True)

        Note that for all loading of saved calculations we've been using
        function :obj:`loadSavedCalculation` .


    Note:
        This doesn't save results of :obj:`plotLevelDiagram` for the
        corresponding calculations. Call the plot function before calling
        :obj:`showPlot` function for the corresponding calculation.

    """

    try:
        ax = calculation.ax
        fig = calculation.fig
        calculation.ax = 0
        calculation.fig = 0

        # close database connections
        atomNumber = 0
        if hasattr(calculation, 'atom'):
            atomNumber = 1
            atomDatabaseConn1 = calculation.atom.conn
            calculation.atom.conn = False
        elif hasattr(calculation, 'atom1'):
            atomNumber = 2
            atomDatabaseConn1 = calculation.atom1.conn
            calculation.atom1.conn = False
            atomDatabaseConn2 = calculation.atom2.conn
            calculation.atom2.conn = False

        output = gzip.GzipFile(fileName, 'wb')
        pickle.dump(calculation, output, pickle.HIGHEST_PROTOCOL)
        output.close()

        calculation.ax = ax
        calculation.fig = fig
        if atomNumber == 1:
            calculation.atom.conn = atomDatabaseConn1
        elif atomNumber == 2:
            calculation.atom1.conn = atomDatabaseConn1
            calculation.atom2.conn = atomDatabaseConn2
    except Exception as ex:
        print(ex)
        print("ERROR: saving of the calculation failed.")
        print(sys.exc_info())
        return 1
    return 0


def loadSavedCalculation(fileName):
    """
    Loads previously saved calculation.

    Loads :obj:`calculations_atom_pairstate.PairStateInteractions` and
    :obj:`calculations_atom_single.StarkMap`
    calculation instance from file named `filename` where it was previously
    saved with :obj:`saveCalculation` .

    Example:
        See example for :obj:`saveCalculation`.

    Args:
        fileName: name of the file where calculation will be saved

    Returns:
        saved calculation
    """

    calculation = False
    try:
        calcInput = gzip.GzipFile(fileName, 'rb')
        calculation = pickle.load(calcInput)
    except Exception as ex:
        print(ex)
        print("ERROR: loading of the calculation from '%s' failed" % fileName)
        print(sys.exc_info())
        return False
    print("Loading of " + calculation.__class__.__name__ + " from '"
          + fileName + "' successful.")

    # establish conneciton to the database
    if hasattr(calculation, 'atom'):
        calculation.atom._databaseInit()
    elif hasattr(calculation, 'atom'):
        calculation.atom1._databaseInit()
        calculation.atom2._databaseInit()

    return calculation

# =================== Saving and loading calculations (END) ===================

# =================== State generation and printing (START) ===================


def singleAtomState(j, m):
    a = np.zeros((int(round(2.0 * j + 1.0, 0)), 1), dtype=np.complex128)
    a[int(round(j + m, 0))] = 1
    return a


def compositeState(s1, s2):
    return np.kron(s1, s2).reshape((s1.shape[0] * s2.shape[0], 1))


def printState(n, l, j, s=None):
    """
        Prints state spectroscopic label for numeric :math:`n`,
        :math:`l`, :math:`s` label of the state

        Args:
            n (int): principal quantum number
            l (int): orbital angular momentum
            j (float): total angular momentum
            s (float): (optional) total spin momentum
    """
    print(printStateString(n, l, j, s=s))


def printStateString(n, l, j, s=None):
    """
        Returns state spectroscopic label for numeric :math:`n`,
        :math:`l`, :math:`j` label of the state.

        Optionally users can define :math:`s`, prompting printing :math:`2S+1`
        index too (commonly used for Alkaline Earth atoms, while it is usually
        omitted for Alkali atoms)

        Args:
            n (int): principal quantum number
            l (int): orbital angular momentum
            j (float): total angular momentum
            s (float): (optional) total spin momentum

        Returns:
            string: label for the state in standard spectroscopic notation
    """
    if s is None:
        return str(n) + " " + printStateLetter(l) + (" %.0d/2" % (j * 2))
    else:
        if abs(floor(j) - j) < 0.1:
            subscript = " %.0d" % (j)
        else:
            subscript = " %.0d/2" % (j * 2)
        return str(n) + (" %d" % (round(2 * s + 1))) + \
            printStateLetter(l) + subscript


def printStateStringLatex(n, l, j, s=None):
    """
        Returns latex code for spectroscopic label for numeric :math:`n`,
        :math:`l`, :math:`j` label of the state.

        Args:
            n (int): principal quantum number
            l (int): orbital angular momentum
            j (float): total angular momentum
            s (float): (optional) total spin momentum

        Returns:
            string: label for the state in standard spectroscopic notation
    """
    if s is None:
        return str(n) + printStateLetter(l) + ("_{%.0d/2}" % (j * 2))
    else:
        if abs(floor(j) - j) < 0.1:
            subscript = "_{%.0d}" % (j)
        else:
            subscript = "_{%.0d/2}" % (j * 2)
        return str(n) + (" ^{%d}" % (round(2 * s + 1))) + \
            printStateLetter(l) + subscript


def printStateLetter(l):
    let = ''
    if l == 0:
        let = "S"
    elif l == 1:
        let = "P"
    elif l == 2:
        let = "D"
    elif l == 3:
        let = "F"
    elif l == 4:
        let = "G"
    elif l == 5:
        let = "H"
    elif l == 6:
        let = "I"
    elif l == 7:
        let = "K"
    elif l == 8:
        let = "L"
    elif l == 9:
        let = "M"
    elif l == 10:
        let = "N"
    else:
        let = " l=%d" % l
    return let


def formatNumberSI(datum,precision=4):
# format datum with SI abbreviation to specified precision (# digits)

    exponent = np.floor(np.log10(np.abs(datum)))
    expInt   = np.floor(exponent/3).astype('int')
    expRange = (expInt * 3).astype('double')

    digitsLeftOfDecimal  = exponent - expRange + 1
    digitsRightOfDecimal = np.max((precision - digitsLeftOfDecimal,0))

    newDatum = datum * 10**(-expRange);

    sisym = ('y','z','a','f','p','n','\mu','m','','k','M','G','T','P','E','Z','Y')
    if np.abs(expRange) <= 24:
        sym = " " + sisym[expInt + 8]
    else:
        sym = " x 10^{%d}"%expRange

    if digitsLeftOfDecimal == precision: # if the last significant figure is in the
                                         # ones place, add the decimal to indicate
                                         # it as such
        sym = "." + sym


    # Formally, if digitsLeftOfDecimal > precision, newDatum should be rounded off
    # to requested precision, but since we are showing no more than 3 digits left
    # of the decimal, it's probably better not to round off

    fmtString = "%%%d.%df%s"%(digitsLeftOfDecimal,digitsRightOfDecimal,sym);

    return fmtString%(newDatum)

# =================== State generation and printing (END) ===================

# =================== E FIELD Coupling (START) ===================


class _EFieldCoupling:
    dataFolder = DPATH

    def __init__(self, theta=0., phi=0.):
        self.theta = theta
        self.phi = phi

        # STARK memoization
        self.conn = sqlite3.connect(os.path.join(self.dataFolder,
                                                 "precalculated_stark.db"))

        # ANGULAR PARTS
        c = self.conn.cursor()
        c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table'
                            AND name='eFieldCoupling_angular';''')
        if (c.fetchone()[0] == 0):
            # create table
            c.execute('''CREATE TABLE IF NOT EXISTS eFieldCoupling_angular
             (l1 TINYINT UNSIGNED, j1_x2 TINYINT UNSIGNED,
             j1_mj1 TINYINT UNSIGNED,
              l2 TINYINT UNSIGNED, j2_x2 TINYINT UNSIGNED,
              j2_mj2 TINYINT UNSIGNED, s_x2 TINYINT UNSIGNED,
             sumPart DOUBLE,
             PRIMARY KEY (l1,j1_x2,j1_mj1,l2,j2_x2,j2_mj2, s_x2)
            ) ''')
            self.conn.commit()

        # COUPLINGS IN ROTATED BASIS (depend on theta, phi)
        self.wgd = WignerDmatrix(self.theta, self.phi)

        c.execute('''DROP TABLE IF EXISTS eFieldCoupling''')
        c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='eFieldCoupling';''')
        if (c.fetchone()[0] == 0):
            # create table
            c.execute('''CREATE TABLE IF NOT EXISTS eFieldCoupling
             (l1 TINYINT UNSIGNED, j1_x2 TINYINT UNSIGNED,
             j1_mj1 TINYINT UNSIGNED,
              l2 TINYINT UNSIGNED, j2_x2 TINYINT UNSIGNED,
              j2_mj2 TINYINT UNSIGNED, s_x2 TINYINT_UNSIGNED,
             coupling DOUBLE,
             PRIMARY KEY (l1,j1_x2,j1_mj1,l2,j2_x2,j2_mj2, s_x2)
            ) ''')
            self.conn.commit()

    def getAngular(self, l1, j1, mj1, l2, j2, mj2, s=0.5):
        c = self.conn.cursor()
        c.execute('''SELECT sumPart FROM eFieldCoupling_angular WHERE
         l1= ? AND j1_x2 = ? AND j1_mj1 = ? AND
         l2 = ? AND j2_x2 = ? AND j2_mj2 = ? AND s_x2 = ?
         ''', (l1, 2 * j1, j1 + mj1, l2, j2 * 2, j2 + mj2, s * 2))
        answer = c.fetchone()
        if (answer):
            return answer[0]

        # calulates sum (See PRA 20:2251 (1979), eq.(10))
        sumPart = 0.

        for ml in np.linspace(mj1 - s, mj1 + s, round(2 * s + 1)):
            if (abs(ml) - 0.1 < l1)and(abs(ml) - 0.1 < l2):
                angularPart = 0.
                if (abs(l1 - l2 - 1) < 0.1):
                    angularPart = ((l1**2 - ml**2)
                                   / ((2. * l1 + 1.) * (2. * l1 - 1.)))**0.5
                elif(abs(l1 - l2 + 1) < 0.1):
                    angularPart = ((l2**2 - ml**2)
                                   / ((2. * l2 + 1.) * (2. * l2 - 1.)))**0.5
                sumPart += CG(l1, ml, s, mj1 - ml, j1, mj1) \
                    * CG(l2, ml, s, mj1 - ml, j2, mj2) \
                    * angularPart

        c.execute(''' INSERT INTO eFieldCoupling_angular
                            VALUES (?,?,?, ?,?,?, ?, ?)''',
                  [l1, 2 * j1, j1 + mj1, l2, j2 * 2, j2 + mj2,
                   s * 2, sumPart])
        self.conn.commit()

        return sumPart

    def getCouplingDivEDivDME(self, l1, j1, mj1, l2, j2, mj2, s=0.5):
        # returns angular coupling without radial part and electric field

        # if calculated before, retrieve from memory
        c = self.conn.cursor()
        c.execute('''SELECT coupling FROM eFieldCoupling WHERE
         l1= ? AND j1_x2 = ? AND j1_mj1 = ? AND
         l2 = ? AND j2_x2 = ? AND j2_mj2 = ? AND s_x2 = ?
         ''', (l1, 2 * j1, j1 + mj1, l2, j2 * 2, j2 + mj2, s * 2))
        answer = c.fetchone()
        if (answer):
            return answer[0]

        # if it is not calculated before, calculate now

        coupling = 0.

        # rotate individual states
        statePart1 = singleAtomState(j1, mj1)
        dMatrix = self.wgd.get(j1)
        statePart1 = np.conj(dMatrix.dot(statePart1))

        statePart2 = singleAtomState(j2, mj2)
        dMatrix = self.wgd.get(j2)
        statePart2 = dMatrix.dot(statePart2)

        # find first common index and start summation
        start = min(j1, j2)

        for mj in np.linspace(-start, start, floor(2 * start + 1)):
            coupling += (self.getAngular(l1, j1, mj, l2, j2, mj)
                         * (statePart1[j1 + mj] * statePart2[j2 + mj])[0].real)

        # save in memory for later use

        c.execute(''' INSERT INTO eFieldCoupling
                            VALUES (?,?,?, ?,?,?, ?, ?)''',
                  [l1, 2 * j1, j1 + mj1, l2, j2 * 2, j2 + mj2, s * 2,
                   coupling])
        self.conn.commit()

        # return result

        return coupling

    def _closeDatabase(self):
        self.conn.commit()
        self.conn.close()
        self.conn = False

# =================== E FIELD Coupling (END) ===================

# we copy the data files to the user home at first run. This avoids
# permission trouble.


setup_data_folder()
