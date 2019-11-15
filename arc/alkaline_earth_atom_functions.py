# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
from .alkali_atom_functions import *


class AlkalineEarthAtom(AlkaliAtom):

    modelPotential_coef = dict()
    """
        Model potential parameters fitted from experimental observations for
        different l (electron angular momentum)
    """

    quantumDefect = [[[0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0]],
                     [[0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0]],
                     [[0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0]],
                     [[0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0],\
                      [0.0,0.0,0.0,0.0,0.0,0.0]]]
    """ Contains list of modified Rydberg-Ritz coefficients for calculating
        quantum defects for
        [[ :math:`^1S_{0},^1P_{1},^1D_{2},^1F_{3}`],
        [ :math:`^3S_{0},^3P_{0},^3D_{1},^3F_{2}`],
        [ :math:`^3S_{0},^3P_{1},^3D_{2},^3F_{3}`],
        [ :math:`^3S_{1},^3P_{2},^3D_{3},^3F_{4}`]]."""



    def __init__(self,preferQuantumDefects=True,cpp_numerov=True):
        self.cpp_numerov = cpp_numerov
        self.preferQuantumDefects = preferQuantumDefects

        self._databaseInit()

        if self.cpp_numerov:
            from .arc_c_extensions import NumerovWavefunction
            self.NumerovWavefunction = NumerovWavefunction

        # load dipole matrix elements previously calculated
        data=[]
        if (self.dipoleMatrixElementFile != ""):
            if (preferQuantumDefects == False):
                self.dipoleMatrixElementFile  = "NIST_"+self.dipoleMatrixElementFile

            try:
                data = np.load(os.path.join(self.dataFolder,\
                                            self.dipoleMatrixElementFile),\
                               encoding = 'latin1', allow_pickle=True)
            except IOError as e:
                print("Error reading dipoleMatrixElement File "+\
                    os.path.join(self.dataFolder,self.dipoleMatrixElementFile))
                print(e)
        # save to SQLite database
        try:
            self.c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='dipoleME';''')
            if (self.c.fetchone()[0] == 0):
                # create table
                self.c.execute('''CREATE TABLE IF NOT EXISTS dipoleME
                 (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED, j1 TINYINT UNSIGNED,
                 s1 TINYINT UNSIGNED,
                 n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED, j2 TINYINT UNSIGNED,
                 s2 TINYINT UNSIGNED,
                 dme DOUBLE,
                 PRIMARY KEY (n1,l1,j1,s1,n2,l2,j2,s2)
                ) ''')
                if (len(data)>0):
                    self.c.executemany('INSERT INTO dipoleME VALUES (?,?,?,?,?,?,?,?,?)',
                                        data)
                self.conn.commit()
        except sqlite3.Error as e:
            print("Error while loading precalculated values into the database")
            print(e)
            exit()

        # load quadrupole matrix elements previously calculated
        data=[]
        if (self.quadrupoleMatrixElementFile != ""):
            if (preferQuantumDefects == False):
                self.quadrupoleMatrixElementFile  = "NIST_"+self.quadrupoleMatrixElementFile
            try:
                data = np.load(os.path.join(self.dataFolder,\
                                            self.quadrupoleMatrixElementFile),\
                               encoding = 'latin1', allow_pickle=True)

            except IOError as e:
                print("Error reading quadrupoleMatrixElementFile File "+\
                    os.path.join(self.dataFolder,self.quadrupoleMatrixElementFile))
                print(e)
        # save to SQLite database
        try:
            self.c.execute('''SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='quadrupoleME';''')
            if (self.c.fetchone()[0] == 0):
                # create table
                self.c.execute('''CREATE TABLE IF NOT EXISTS quadrupoleME
                 (n1 TINYINT UNSIGNED, l1 TINYINT UNSIGNED, j1 TINYINT UNSIGNED,
                 s1 TINYINT UNSIGNED,
                 n2 TINYINT UNSIGNED, l2 TINYINT UNSIGNED, j2 TINYINT UNSIGNED,
                 s2 TINYINT UNSIGNED,
                 qme DOUBLE,
                 PRIMARY KEY (n1,l1,j1,s1,n2,l2,j2,s2)
                ) ''')
                if (len(data)>0):
                    self.c.executemany('INSERT INTO quadrupoleME VALUES (?,?,?,?,?,?,?,?,?)',
                                        data)
                self.conn.commit()
        except sqlite3.Error as e:
            print("Error while loading precalculated values into the database")
            print(e)
            exit()

        if (self.levelDataFromNIST == ""):
            print("NIST level data file not specified. Only quantum defects will be used.")
        else:
            levels = self._parseLevelsFromNIST(os.path.join(self.dataFolder,\
                                               self.levelDataFromNIST))
            br = 0
            while br<len(levels):
                self._addEnergy(*levels[br])
                br = br+1
            try:
                self.conn.commit()
            except sqlite3.Error as e:
                print("Error while loading precalculated values into the database")
                print(e)
                print(n," ",l," ",j," ",s)
                exit()


    def _parseLevelsFromNIST(self,fileData):
        print(fileData)
        data = np.loadtxt(fileData, delimiter=",",
                          usecols=(0,1,3,2,4))
        return data


    def _addEnergy(self, n, l ,j, s, energy):
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
        self.c.execute('INSERT INTO energyLevel VALUES (?,?,?,?,?)',
            (int(n), int(l), int(j), int(s),
            energy * 1.e2
            * physical_constants["inverse meter-electron volt relationship"][0]
            - self.ionisationEnergy)
            )
        self.NISTdataLevels = max(self.NISTdataLevels, int(n))
        # saves energy in eV


    def _databaseInit(self):
        self.conn = sqlite3.connect(os.path.join(self.dataFolder,\
                                                 self.precalculatedDB))
        self.c = self.conn.cursor()

        # create space for storing NIST/literature energy levels
        self.c.execute('''SELECT COUNT(*) FROM sqlite_master
                        WHERE type='table' AND name='energyLevel';''')
        if (self.c.fetchone()[0] != 0):
            self.c.execute('''DROP TABLE energyLevel''')
        # create fresh table
        self.c.execute('''CREATE TABLE IF NOT EXISTS energyLevel
             (n TINYINT UNSIGNED, l TINYINT UNSIGNED, j TINYINT UNSIGNED,
             s TINYINT UNSIGNED,
             energy DOUBLE,
             PRIMARY KEY (n, l, j, s)
            ) ''')

        self.conn.commit()



    def _getSavedEnergy(self, n, l, j, s=0):
        self.c.execute('''SELECT energy FROM energyLevel WHERE
            n= ? AND l = ? AND j = ? AND
            s = ? ''',(n, l, j, s))
        energy = self.c.fetchone()
        if (energy):
            return energy[0]
        else:
            return 0      # there is no saved energy level measurement
