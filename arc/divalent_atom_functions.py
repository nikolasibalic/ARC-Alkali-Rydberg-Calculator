# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
from .alkali_atom_functions import *

import os
import numpy as np
from math import sqrt

import sqlite3
sqlite3.register_adapter(np.float64, float)
sqlite3.register_adapter(np.float32, float)
sqlite3.register_adapter(np.int64, int)
sqlite3.register_adapter(np.int32, int)


class DivalentAtom(AlkaliAtom):
    """
        Implements general calculations for Alkaline Earths, and other divalent
        atoms.

        This class inherits :obj:`arc.alkali_atom_functions.AlkaliAtom` .
        Most of the methods can be directly used from there, and the source
        for them is provided in the base class. Few methods that are
        implemented differently for Alkaline Earths are defined here.

        Args:
            preferQuantumDefects (bool):
                Use quantum defects for energy level calculations. If False,
                uses NIST ASD values where available. If True, uses quantum
                defects for energy calculations for principal quantum numbers
                within the range specified in :obj:`defectFittingRange` which
                is specified for each element and series separately.
                For principal quantum numbers below this value, NIST ASD values
                are used if existing, since quantum defects. Default is True.
            cpp_numerov (bool):
                This switch for Alkaline Earths at the moment doesn't have
                any effect since wavefunction calculation function is not
                implemented (d.m.e. and quadrupole matrix elements are
                calculated directly semiclassically)
    """

    modelPotential_coef = dict()
    """
        Model potential parameters fitted from experimental observations for
        different l (electron angular momentum)
    """

    quantumDefect = [[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
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

    #: file with .csv data, each row is
    #: `[n, l, s, j, energy, source, absolute uncertanty]`
    levelDataFromNIST = ""

    #: Not used with DivalentAtom, see :obj:`defectFittingRange` instead.
    minQuantumDefectN = None

    #: Used for AlkalineEarths to define minimum and maximum principal quantum
    #: number for which quantum defects are valid. Ranges are stored under
    #: keys defined as state terms ({'stateLabel':[minN, maxN]}, e.g. '1S0').
    #: Dictionary returns array
    #: stating minimal and maximal principal quantum number for which quantum
    #: defects were fitted. For example::
    #:      limits = self.defectFittingRange['1S0']
    #:      print("Minimal n = %d" % limits[0])
    #:      print("Maximal n = %d" % limits[1]) 1
    defectFittingRange = {}

    #: flag that is turned to True if the energy levels of this atom were
    #: calculated by extrapolating with quantum defects values outside the
    #: quantum defect fitting range.
    energyLevelsExtrapolated = False

    def __init__(self, preferQuantumDefects=True, cpp_numerov=True):
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
                      + os.path.join(self.dataFolder,
                                     self.dipoleMatrixElementFile))
                print(e)
        # save to SQLite database
        try:
            c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='dipoleME';''')
            if (c.fetchone()[0] == 0):
                # create table
                c.execute('''CREATE TABLE IF NOT EXISTS dipoleME
                 (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED,
                 j1 TINYINT UNSIGNED,
                 n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED,
                 j2 TINYINT UNSIGNED, s TINYINT UNSIGNED,
                 dme DOUBLE,
                 PRIMARY KEY (n1,l1,j1,n2,l2,j2,s)
                ) ''')
                if (len(data) > 0):
                    c.executemany(
                        'INSERT INTO dipoleME VALUES (?,?,?,?,?,?,?,?)',
                        data)
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
                 j1 TINYINT UNSIGNED,
                 n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED,
                 j2 TINYINT UNSIGNED, s TINYINT UNSIGNED,
                 qme DOUBLE,
                 PRIMARY KEY (n1,l1,j1,n2,l2,j2,s)
                ) ''')
                if (len(data) > 0):
                    c.executemany(
                        'INSERT INTO quadrupoleME VALUES (?,?,?,?,?,?,?,?)',
                        data)
                self.conn.commit()
        except sqlite3.Error as e:
            print("Error while loading precalculated values into the database")
            print(e)
            exit()

        if (self.levelDataFromNIST == ""):
            print("NIST level data file not specified."
                  " Only quantum defects will be used.")
        else:
            levels = self._parseLevelsFromNIST(
                os.path.join(self.dataFolder, self.levelDataFromNIST)
            )
            br = 0
            while br < len(levels):
                self._addEnergy(*levels[br])
                br = br + 1
            try:
                self.conn.commit()
            except sqlite3.Error as e:
                print("Error while loading precalculated values"
                      "into the database")
                print(e)
                exit()

        self._readLiteratureValues()

    def _parseLevelsFromNIST(self, fileData):
        data = np.loadtxt(fileData, delimiter=",",
                          usecols=(0, 1, 3, 2, 4))
        return data

    def _addEnergy(self, n, l, j, s, energy):
        """
            Adds energy level relative to

            NOTE:
            Requres changes to be commited to the sql database afterwards!

            Args:
                n: principal quantum number
                l: orbital angular momentum quantum number
                j: total angular momentum quantum number
                s: spin quantum number
                energy: energy in cm^-1 relative to the ground state
        """
        c = self.conn.cursor()
        c.execute(
            'INSERT INTO energyLevel VALUES (?,?,?,?,?)',
            (int(n), int(l), int(j), int(s),
             energy * 1.e2
             * physical_constants["inverse meter-electron volt relationship"][0]
             - self.ionisationEnergy)
        )
        self.NISTdataLevels = max(self.NISTdataLevels, int(n))
        # saves energy in eV

    def _databaseInit(self):
        self.conn = sqlite3.connect(os.path.join(self.dataFolder,
                                                 self.precalculatedDB))
        c = self.conn.cursor()

        # create space for storing NIST/literature energy levels
        c.execute('''SELECT COUNT(*) FROM sqlite_master
                        WHERE type='table' AND name='energyLevel';''')
        if (c.fetchone()[0] != 0):
            c.execute('''DROP TABLE energyLevel''')
        # create fresh table
        c.execute('''CREATE TABLE IF NOT EXISTS energyLevel
             (n TINYINT UNSIGNED, l TINYINT UNSIGNED, j TINYINT UNSIGNED,
             s TINYINT UNSIGNED,
             energy DOUBLE,
             PRIMARY KEY (n, l, j, s)
            ) ''')

        self.conn.commit()

    def getEnergy(self, n, l, j, s=None):
        if s is None:
            raise ValueError("Spin state for DivalentAtom has to be "
                             "explicitly defined as a keyword argument "
                             "s=0 or s=1")
        if l >= n:
            raise ValueError(
                "Requested energy for state l=%d >= n=%d !" % (l, n))

        stateLabel = "%d%s%d" % (int(2*s+1), printStateLetter(l), j)
        minQuantumDefectN = 100000
        maxQuantumDefectN = 0

        if stateLabel in self.defectFittingRange.keys():
            minQuantumDefectN = self.defectFittingRange[stateLabel][0]
            maxQuantumDefectN = self.defectFittingRange[stateLabel][1]

        # use NIST data ?
        if (not self.preferQuantumDefects or n < minQuantumDefectN):
            savedEnergy = self._getSavedEnergy(n, l, j, s=s)
            if (abs(savedEnergy) > 1e-8):
                return savedEnergy
            else:
                if (n < minQuantumDefectN or n > maxQuantumDefectN):
                    self.energyLevelsExtrapolated = True

        # else, use quantum defects
        defect = self.getQuantumDefect(n, l, j, s=s)
        return -self.scaledRydbergConstant / ((n - defect)**2)

    def _getSavedEnergy(self, n, l, j, s=0):
        c = self.conn.cursor()
        c.execute('''SELECT energy FROM energyLevel WHERE
            n= ? AND l = ? AND j = ? AND
            s = ? ''', (n, l, j, s))
        energy = c.fetchone()
        if (energy):
            return energy[0]
        else:
            return 0      # there is no saved energy level measurement

    def getRadialMatrixElement(self,
                               n1, l1, j1,
                               n2, l2, j2,
                               s=None,
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
                s (float): is required argument, total spin angular momentum of
                    state. Specify `s=0` for singlet state or `s=1` for
                    triplet state.
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
        if s is None:
            raise ValueError("You must specify total angular momentum s "
                             "explicitly.")

        dl = abs(l1 - l2)
        dj = abs(j2 - j2)
        if not(dl == 1 and (dj < 1.1)):
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

        c = self.conn.cursor()

        if useLiterature:
            # is there literature value for this DME? If there is,
            # use the best one (smalles error)
            c.execute('''SELECT dme FROM literatureDME WHERE
             n1= ? AND l1 = ? AND j1 = ? AND
             n2 = ? AND l2 = ? AND j2 = ? AND s = ?
             ORDER BY errorEstimate ASC''',
             (n1, l1, j1, n2, l2, j2, s))

            answer = c.fetchone()
            if (answer):
                # we did found literature value
                return answer[0]
        # was this calculated before? If it was, retrieve from memory

        c.execute(
            '''SELECT dme FROM dipoleME WHERE
            n1= ? AND l1 = ? AND j1 = ? AND
            n2 = ? AND l2 = ? AND j2 = ? AND s = ?''',
            (n1, l1, j1, n2, l2, j2, s)
            )
        dme = c.fetchone()
        if (dme):
            return dme[0]

        dipoleElement = self._getRadialDipoleSemiClassical(
            n1, l1, j1, n2, l2, j2, s=s
            )

        c.execute(
            ''' INSERT INTO dipoleME VALUES (?,?,?, ?,?,?, ?, ?)''',
            [n1, l1, j1,
             n2, l2, j2, s,
             dipoleElement])
        self.conn.commit()

        return dipoleElement

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
             (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED, j1 TINYINT UNSIGNED,
             n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED, j2 TINYINT UNSIGNED,
             s TINYINT UNSIGNED,
             dme DOUBLE,
             typeOfSource TINYINT,
             errorEstimate DOUBLE,
             comment TINYTEXT,
             ref TINYTEXT,
             refdoi TINYTEXT
            );''')
            c.execute('''CREATE INDEX compositeIndex
            ON literatureDME (n1,l1,j1,n2,l2,j2,s); ''')
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
                    j1 = int(row[2])
                    s1 = int(row[3])

                    n2 = int(row[4])
                    l2 = int(row[5])
                    j2 = int(row[6])
                    s2 = int(row[7])
                    if (s1 != s2):
                        raise ValueError("Error reading litearture: database "
                                         "cannot accept spin changing "
                                         "transitions")
                    s = s1
                    if (
                        self.getEnergy(n1, l1, j1, s=s)
                            > self.getEnergy(n2, l2, j2, s=s)
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

                    # To-DO : see in what notation are Strontium literature elements saved
                    print("To-do (_readLiteratureValues): see in what notation are Sr literature saved (angular part)")
                    dme = float(row[8]) / (
                        (-1)**(int(l1 + s + j2 + 1.))
                        * sqrt((2. * j1 + 1.) * (2. * j2 + 1.))
                        * Wigner6j(j1, 1., j2, l2, s, l1)
                        * (-1)**l1 * sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0))
                        * Wigner3j(l1, 1, l2, 0, 0, 0))

                    comment = row[9]
                    typeOfSource = int(row[10])  # 0 = experiment; 1 = theory
                    errorEstimate = float(row[11])
                    ref = row[12]
                    refdoi = row[13]

                    literatureDME.append([n1, l1, j1,
                                          n2, l2, j2, s, dme,
                                          typeOfSource, errorEstimate,
                                          comment, ref,
                                          refdoi])
                i += 1
            fn.close()

            try:
                if i > 1:
                    c.executemany('''INSERT INTO literatureDME
                                        VALUES (?,?,?, ?,?,?,?
                                                ?,?,?,?,?,?)''',
                                       literatureDME)
                    self.conn.commit()

            except sqlite3.Error as e:
                print("Error while loading precalculated values "
                      "into the database")
                print(e)
                print(literatureDME)
                exit()

        except IOError as e:
            print("Error reading literature values File "
                  + self.literatureDMEfilename)
            print(e)

    def getLiteratureDME(self,
                         n1, l1, j1,
                         n2, l2, j2, s=0):
        """
            Returns literature information on requested transition.

            Args:
                n1,l1,j1: one of the states we are coupling
                n2,l2,j2: the other state to which we are coupling
                s: (optional) spin of the state. Default s=0.

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

                    * typeOfSource=1 if the value is theoretical\
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
                * s
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

        # is there literature value for this DME? If there is,
        # use the best one (wit the smallest error)

        c = self.conn.cursor()
        c.execute('''SELECT dme, typeOfSource,
                     errorEstimate ,
                     comment ,
                     ref,
                     refdoi FROM literatureDME WHERE
                     n1= ? AND l1 = ? AND j1 = ? AND
                     n2 = ? AND l2 = ? AND j2 = ? AND s = ?
                     ORDER BY errorEstimate ASC''',
                       (n1, l1, j1, n2, l2, j2, s))
        answer = c.fetchone()
        if (answer):
            # we did found literature value
            return True, answer[0], [answer[1], answer[2], answer[3],
                                     answer[4], answer[5]]

        # if we are here, we were unsucessfull in literature search
        # for this value
        return False, 0, []

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
                s (float): optional. Spin of the state. Default 0.5 is for
                    Alkali

            Returns:
                float: quadrupole matrix element (:math:`a_0^2 e`).
        """

        dl = abs(l1 - l2)
        dj = abs(j1 - j2)
        if not ((dl == 0 or dl == 2 or dl == 1)and (dj < 2.1)):
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

        # was this calculated before? If yes, retrieve from memory.
        c = self.conn.cursor()
        c.execute('''SELECT qme FROM quadrupoleME WHERE
            n1= ? AND l1 = ? AND j1 = ? AND
            n2 = ? AND l2 = ? AND j2 = ? AND s= ?''',
            (n1, l1, j1, n2, l2, j2, s))
        qme = c.fetchone()
        if (qme):
            return qme[0]

        # if it wasn't, calculate now

        quadrupoleElement = self._getRadialQuadrupoleSemiClassical(
            n1, l1, j1, n2, l2, j2, s=s
        )

        c.execute(''' INSERT INTO quadrupoleME VALUES (?,?,?, ?,?,?,?, ?)''',
                       [n1, l1, j1, n2, l2, j2, s, quadrupoleElement])
        self.conn.commit()

        return quadrupoleElement

    def radialWavefunction(self, l, s, j, stateEnergy,
                           innerLimit, outerLimit, step):
        """
            Not implemented for Alkaline earths
        """
        raise NotImplementedError("radialWavefunction calculation for alkaline"
                                  " earths has not been implemented yet.")

    def effectiveCharge(self, l, r):
        """
            Not implemented for Alkaline earths
        """
        raise NotImplementedError("effectiveCharge calculation for alkaline"
                                  " earths has not been implemented yet.")

    def corePotential(self, l, r):
        """
            Not implemented for Alkaline earths
        """
        raise NotImplementedError("corePotential calculation for alkaline"
                                  " earths has not been implemented yet.")

    def potential(self, l, s, j, r):
        """
            Not implemented for Alkaline earths
        """
        raise NotImplementedError("potential calculation for alkaline"
                                  " earths has not been implemented yet.")

    def getStateLifetime(self, n, l, j, temperature=0, includeLevelsUpTo=0,
                         s=0):
        print("WARNING:  For AlkalineEarths, lifetimes are observed to be "
              "significantly modified by inter-electron correlations that are "
              "not included in this code (see Vaillant et al., J. Phys B 47 "
              "155001 (2015) for examples).  Use with caution.")
        # after waring user, call method from the parent class
        # (parent of DivalentAtom is AlkaliAtom)
        return super(DivalentAtom, self).getStateLifetime(
            n, l, j,
            temperature=temperature,
            includeLevelsUpTo=includeLevelsUpTo,
            s=s)
