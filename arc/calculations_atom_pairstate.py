# -*- coding: utf-8 -*-
# ruff: noqa: E741

"""
Pair-state basis level diagram calculations

Calculates Rydberg spaghetti of level diagrams, as well as pertubative C6
and similar properties. It also allows calculation of Foster resonances
tuned by DC electric fields.

Example:
    Calculation of the Rydberg eigenstates in pair-state basis for Rubidium
    in the vicinity of the
    :math:`|60~S_{1/2}~m_j=1/2,~60~S_{1/2}~m_j=1/2\\rangle` state. Colour
    highlights coupling strength from state :math:`6~P_{1/2}~m_j=1/2` with
    :math:`\\pi` (:math:`q=0`) polarized light.
    eigenstates::

        from arc import *
        calc1 = PairStateInteractions(Rubidium(), 60, 0, 0.5, 60, 0, 0.5,0.5, 0.5)
        calc1.defineBasis( 0., 0., 4, 5,10e9)
        # optionally we can save now results of calculation for future use
        saveCalculation(calc1,"mycalculation.pkl")
        calculation1.diagonalise(linspace(1,10.0,30),250,progressOutput = True,drivingFromState=[6,1,0.5,0.5,0])
        calc1.plotLevelDiagram()
        calc1.ax.set_xlim(1,10)
        calc1.ax.set_ylim(-2,2)
        calc1.showPlot()

"""

from __future__ import division, print_function, absolute_import

from arc._database import sqlite3, UsedModulesARC
from arc.wigner import Wigner6j, CG, WignerDmatrix
from arc.alkali_atom_functions import (
    _atomLightAtomCoupling,
    singleAtomState,
    compositeState,
)
from scipy.constants import physical_constants, pi
import gzip
import sys
import os
import datetime
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from arc.calculations_atom_single import StarkMap
from arc.alkali_atom_functions import (
    printStateStringLatex,
    printStateString,
    printStateLetter,
)
from arc.divalent_atom_functions import DivalentAtom
from scipy.special import factorial
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
from scipy.constants import e as C_e
from scipy.constants import h as C_h
from scipy.constants import c as C_c

import numpy as np
from math import exp, sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl
from inspect import getmodule as inspectgetmodule
import urllib.request
import h5py

mpl.rcParams["xtick.minor.visible"] = True
mpl.rcParams["ytick.minor.visible"] = True
mpl.rcParams["xtick.major.size"] = 8
mpl.rcParams["ytick.major.size"] = 8
mpl.rcParams["xtick.minor.size"] = 4
mpl.rcParams["ytick.minor.size"] = 4
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["font.family"] = "serif"


# for matrices


if sys.version_info > (2,):
    xrange = range


DPATH = os.path.join(os.path.expanduser("~"), ".arc-data")

__all__ = ["PairStateInteractions", "StarkMapResonances"]


class PairStateInteractions:
    """
        Calculates Rydberg level diagram (spaghetti) for the given pair state

        Initializes Rydberg level spaghetti calculation for the given atom
        species (or for two atoms of different species) in the vicinity
        of the given pair state. For details of calculation see
        Ref. [1]_. For a quick start point example see
        `interactions example snippet`_.
        For inter-species calculations see
        `inter-species interaction calculation snippet`_.

        .. _`interactions example snippet`:
            ./Rydberg_atoms_a_primer.html#Short-range-interactions

        .. _`inter-species interaction calculation snippet`:
           ./ARC_3_0_introduction.html#Inter-species-pair-state-calculations

        Parameters:
            atom (:obj:`arc.alkali_atom_functions.AlkaliAtom` or :obj:`arc.divalent_atom_functions.DivalentAtom`):
                = {
                :obj:`arc.alkali_atom_data.Lithium6`,
                :obj:`arc.alkali_atom_data.Lithium7`,
                :obj:`arc.alkali_atom_data.Sodium`,
                :obj:`arc.alkali_atom_data.Potassium39`,
                :obj:`arc.alkali_atom_data.Potassium40`,
                :obj:`arc.alkali_atom_data.Potassium41`,
                :obj:`arc.alkali_atom_data.Rubidium85`,
                :obj:`arc.alkali_atom_data.Rubidium87`,
                :obj:`arc.alkali_atom_data.Caesium`,
                :obj:`arc.divalent_atom_data.Strontium88`,
                :obj:`arc.divalent_atom_data.Calcium40`
                :obj:`arc.divalent_atom_data.Ytterbium174` }
                Select the alkali metal for energy level
                diagram calculation
            n (int):
                principal quantum number for the *first* atom
            l (int):
                orbital angular momentum for the *first* atom
            j (float):
                total angular momentum for the *first* atom
            nn (int):
                principal quantum number for the *second* atom
            ll (int):
                orbital angular momentum for the *second* atom
            jj (float):
                total angular momentum for the *second* atom
            m1 (float):
                projection of the total angular momentum on z-axis
                for the *first* atom
            m2 (float):
                projection of the total angular momentum on z-axis
                for the *second* atom
            interactionsUpTo (int):
                Optional. If set to 1, includes only
                dipole-dipole interactions. If set to 2 includes interactions
                up to quadrupole-quadrupole. Default value is 1.
            s (float):
                optional, spin state of the first atom. Default value
                of 0.5 is correct for :obj:`arc.alkali_atom_functions.AlkaliAtom`
                but for :obj:`arc.divalent_atom_functions.DivalentAtom`
                it has to be explicitly set to 0 or 1 for
                singlet and triplet states respectively.
                **If `s2` is not specified, it is assumed that the second
                atom is in the same spin state.**
            s2 (float):
                optinal, spin state of the second atom. If not
                specified (left to default value None) it will assume spin
                state of the first atom.
            atom2 (:obj:`arc.alkali_atom_functions.AlkaliAtom` or :obj:`arc.divalent_atom_functions.DivalentAtom`):
                optional,
                specifies atomic species for the second atom, enabeling
                calculation of **inter-species pair-state interactions**.
                If not specified (left to default value None) it will assume
                spin state of the first atom.

        References:
            .. [1] T. G Walker, M. Saffman, PRA **77**, 032723 (2008)
                https://doi.org/10.1103/PhysRevA.77.032723

        Examples:
            **Advanced interfacing of pair-state is2=None, atom2=Nonenteractions calculations
            (PairStateInteractions class).** This
            is an advanced example intended for building up extensions to the
            existing code. If you want to directly access the pair-state
            interaction matrix, constructed by :obj:`defineBasis`,
            you can assemble it easily from diagonal part
            (stored in :obj:`matDiagonal` ) and off-diagonal matrices whose
            spatial dependence is :math:`R^{-3},R^{-4},R^{-5}` stored in that
            order in :obj:`matR`. Basis states are stored in :obj:`basisStates`
            array.

            >>> from arc import *
            >>> calc = PairStateInteractions(Rubidium(), 60,0,0.5, \
                60,0,0.5, 0.5,0.5,interactionsUpTo = 1)
            >>> # theta=0, phi = 0, range of pqn, range of l, deltaE = 25e9
            >>> calc.defineBasis(0 ,0 , 5, 5, 25e9, progressOutput=True)
            >>> # now calc stores interaction matrix and relevant basis
            >>> # we can access this directly and generate interaction matrix
            >>> # at distance rval :
            >>> rval = 4  # in mum
            >>> matrix = calc.matDiagonal
            >>> rX = (rval*1.e-6)**3
            >>> for matRX in self.matR:
            >>>     matrix = matrix + matRX/rX
            >>>     rX *= (rval*1.e-6)
            >>> # matrix variable now holds full interaction matrix for
            >>> # interacting atoms at distance rval calculated in
            >>> # pair-state basis states can be accessed as
            >>> basisStates = calc.basisStates
    """

    dataFolder = DPATH

    # =============================== Methods ===============================

    def __init__(
        self,
        atom,
        n,
        l,
        j,
        nn,
        ll,
        jj,
        m1,
        m2,
        interactionsUpTo=1,
        s=0.5,
        s2=None,
        atom2=None,
    ):
        # alkali atom type, principal quantum number, orbital angular momentum,
        #  total angular momentum projections of the angular momentum on z axis
        self.atom1 = atom  #: the first atom type (isotope)
        if atom2 is None:
            self.atom2 = atom  #: the second atom type (isotope)
        else:
            self.atom2 = atom2  #: thge second atom type (isotope)
        self.n = n  # : pair-state definition: principal quantum number of the first atom
        self.l = l  # : pair-state definition: orbital angular momentum of the first atom
        self.j = j  # : pair-state definition: total angular momentum of the first atom
        self.nn = nn  # : pair-state definition: principal quantum number of the second atom
        self.ll = ll  # : pair-state definition: orbital angular momentum of the second atom
        self.jj = jj  # : pair-state definition: total angular momentum of the second atom
        self.m1 = m1  # : pair-state definition: projection of the total ang. momentum for the *first* atom
        self.m2 = m2  # : pair-state definition: projection of the total angular momentum for the *second* atom
        self.interactionsUpTo = interactionsUpTo
        """ Specifies up to which approximation we include in pair-state interactions.
            By default value is 1, corresponding to pair-state interactions up to
            dipole-dipole coupling. Value of 2 is also supported, corresponding
            to pair-state interactions up to quadrupole-quadrupole coupling.
        """

        if issubclass(type(atom), DivalentAtom) and not (s == 0 or s == 1):
            raise ValueError(
                "total angular spin s has to be defined explicitly "
                "for calculations, and value has to be 0 or 1 "
                "for singlet and tripplet states respectively."
            )
        self.s1 = s  #: total spin angular momentum, optional (default 0.5)

        if s2 is None:
            self.s2 = s
        else:
            self.s2 = s2

        # check that values of spin states are valid for entered atomic species

        if issubclass(type(self.atom1), DivalentAtom):
            if abs(self.s1) > 0.1 and abs(self.s1 - 1) > 0.1:
                raise ValueError(
                    "atom1 is DivalentAtom and its spin has to be "
                    "s=0 or s=1 (for singlet and triplet states "
                    "respectively)"
                )
        elif abs(self.s1 - 0.5) > 0.1:
            raise ValueError("atom1 is AlkaliAtom and its spin has to be s=0.5")
        if issubclass(type(self.atom2), DivalentAtom):
            if abs(self.s2) > 0.1 and abs(self.s2 - 1) > 0.1:
                raise ValueError(
                    "atom2 is DivalentAtom and its spin has to be "
                    "s=0 or s=1 (for singlet and triplet states "
                    "respectively)"
                )
        elif abs(self.s2 - 0.5) > 0.1:
            # we have divalent atom
            raise ValueError("atom2 is AlkaliAtom and its spin has to be s=0.5")
        if abs((self.s1 - self.m1) % 1) > 0.1:
            raise ValueError(
                "atom1 with spin s = %.1d cannot have m1 = %.1d"
                % (self.s1, self.m1)
            )
        if abs((self.s2 - self.m2) % 1) > 0.1:
            raise ValueError(
                "atom2 with spin s = %.1d cannot have m2 = %.1d"
                % (self.s2, self.m2)
            )

        # ====================== J basis (not resolving mj) ===================

        self.coupling = []
        """
            List of matrices defineing coupling strengths between the states in
            J basis (not resolving :math:`m_j` ). Basis is given by
            :obj:`PairStateInteractions.channel`. Used as intermediary for full
            interaction matrix calculation by
            :obj:`PairStateInteractions.defineBasis`.
        """
        self.channel = []
        """
            states relevant for calculation, defined in J basis (not resolving
            :math:`m_j`. Used as intermediary for full interaction matrix
            calculation by :obj:`PairStateInteractions.defineBasis`.
        """

        # ======================= Full basis (resolving mj) ===================

        self.basisStates = []
        """
            List of pair-states for calculation. In the form
            [[n1,l1,j1,mj1,n2,l2,j2,mj2], ...].
            Each state is an array [n1,l1,j1,mj1,n2,l2,j2,mj2] corresponding to
            :math:`|n_1,l_1,j_1,m_{j1},n_2,l_2,j_2,m_{j2}\\rangle` state.
            Calculated by :obj:`PairStateInteractions.defineBasis`.
        """
        self.matrixElement = []
        """
            `matrixElement[i]` gives index of state in
            :obj:`PairStateInteractions.channel` basis
            (that doesn't resolve :math:`m_j` states), for the given index `i`
            of the state in :obj:`PairStateInteractions.basisStates`
            ( :math:`m_j` resolving) basis.
        """

        # variuos parts of interaction matrix in pair-state basis
        self.matDiagonal = []
        """
            Part of interaction matrix in pair-state basis that doesn't depend
            on inter-atomic distance. E.g. diagonal elements of the interaction
            matrix, that describe energies of the pair states in unperturbed
            basis, will be stored here. Basis states are stored in
            :obj:`PairStateInteractions.basisStates`. Calculated by
            :obj:`PairStateInteractions.defineBasis`.
        """
        self.matR = []
        """
            Stores interaction matrices in pair-state basis
            that scale as :math:`1/R^3`, :math:`1/R^4` and :math:`1/R^5`
            with distance in  :obj:`matR[0]`, :obj:`matR[1]` and :obj:`matR[2]`
            respectively. These matrices correspond to dipole-dipole
            ( :math:`C_3`), dipole-quadrupole ( :math:`C_4`) and
            quadrupole-quadrupole ( :math:`C_5`) interactions
            coefficients. Basis states are stored in
            :obj:`PairStateInteractions.basisStates`.
            Calculated by :obj:`PairStateInteractions.defineBasis`.
        """
        self.originalPairStateIndex = 0
        """
            index of the original n,l,j,m1,nn,ll,jj,m2 pair-state in the
            :obj:`PairStateInteractions.basisStates` basis.
        """

        self.matE = []
        self.matB_1 = []
        self.matB_2 = []

        # ===================== Eigen states and plotting =====================

        # finding perturbed energy levels
        self.r = []  # detuning scale
        self.y = []  # energy levels
        self.highlight = []

        # pointers towards figure
        self.fig = 0
        self.ax = 0

        # for normalization of the maximum coupling later
        self.maxCoupling = 0.0

        # n,l,j,mj, drive polarization q
        self.drivingFromState = [0, 0, 0, 0, 0]

        # sam = saved angular matrix metadata
        self.angularMatrixFile = "angularMatrix.npy"
        self.angularMatrixFile_meta = "angularMatrix_meta.npy"
        # self.sam = []
        self.savedAngularMatrix_matrix = []

        # intialize precalculated values for factorial term
        # in __getAngularMatrix_M
        def fcoef(l1, l2, m):
            return (
                factorial(l1 + l2)
                / (
                    factorial(l1 + m)
                    * factorial(l1 - m)
                    * factorial(l2 + m)
                    * factorial(l2 - m)
                )
                ** 0.5
            )

        x = self.interactionsUpTo
        self.fcp = np.zeros((x + 1, x + 1, 2 * x + 1))
        for c1 in range(1, x + 1):
            for c2 in range(1, x + 1):
                for p in range(-min(c1, c2), min(c1, c2) + 1):
                    self.fcp[c1, c2, p + x] = fcoef(c1, c2, p)

        self.conn = False

    def __getAngularMatrix_M(self, l, j, ll, jj, l1, j1, l2, j2):
        # did we already calculated this matrix?
        c = self.conn.cursor()
        c.execute(
            """SELECT ind FROM pair_angularMatrix WHERE
             l1 = ? AND j1_x2 = ? AND
             l2 = ? AND j2_x2 = ? AND
             l3 = ? AND j3_x2 = ? AND
             l4 = ? AND j4_x2 = ?
             """,
            (l, j * 2, ll, jj * 2, l1, j1 * 2, l2, j2 * 2),
        )

        index = c.fetchone()
        if index:
            return self.savedAngularMatrix_matrix[index[0]]

        # determine coupling
        dl = abs(l - l1)
        dj = abs(j - j1)
        c1 = 0
        if dl == 1 and (dj < 1.1):
            c1 = 1  # dipole coupling
        elif dl == 0 or dl == 2 or dl == 1:
            c1 = 2  # quadrupole coupling
        else:
            raise ValueError("error in __getAngularMatrix_M")

        dl = abs(ll - l2)
        dj = abs(jj - j2)
        c2 = 0
        if dl == 1 and (dj < 1.1):
            c2 = 1  # dipole coupling
        elif dl == 0 or dl == 2 or dl == 1:
            c2 = 2  # quadrupole coupling
        else:
            raise ValueError("error in __getAngularMatrix_M")

        am = np.zeros(
            (
                round((2 * j1 + 1) * (2 * j2 + 1)),
                round((2 * j + 1) * (2 * jj + 1)),
            ),
            dtype=np.float64,
        )

        if (c1 > self.interactionsUpTo) or (c2 > self.interactionsUpTo):
            return am

        j1range = np.linspace(-j1, j1, round(2 * j1) + 1)
        j2range = np.linspace(-j2, j2, round(2 * j2) + 1)
        jrange = np.linspace(-j, j, round(2 * j) + 1)
        jjrange = np.linspace(-jj, jj, round(2 * jj) + 1)

        for m1 in j1range:
            for m2 in j2range:
                # we have chosen the first index
                index1 = round(
                    m1 * (2.0 * j2 + 1.0) + m2 + (j1 * (2.0 * j2 + 1.0) + j2)
                )
                for m in jrange:
                    for mm in jjrange:
                        # we have chosen the second index
                        index2 = round(
                            m * (2.0 * jj + 1.0)
                            + mm
                            + (j * (2.0 * jj + 1.0) + jj)
                        )

                        # angular matrix element from Sa??mannshausen, Heiner,
                        # Merkt, Fr??d??ric, Deiglmayr, Johannes
                        # PRA 92: 032505 (2015)
                        elem = (
                            (-1.0) ** (j + jj + self.s1 + self.s2 + l1 + l2)
                            * CG(l, 0, c1, 0, l1, 0)
                            * CG(ll, 0, c2, 0, l2, 0)
                        )
                        elem = (
                            elem
                            * sqrt((2.0 * l + 1.0) * (2.0 * ll + 1.0))
                            * sqrt((2.0 * j + 1.0) * (2.0 * jj + 1.0))
                        )
                        elem = (
                            elem
                            * Wigner6j(l, self.s1, j, j1, c1, l1)
                            * Wigner6j(ll, self.s2, jj, j2, c2, l2)
                        )

                        sumPol = 0.0  # sum over polarisations
                        limit = min(c1, c2)
                        for p in xrange(-limit, limit + 1):
                            sumPol = sumPol + self.fcp[
                                c1, c2, p + self.interactionsUpTo
                            ] * CG(j, m, c1, p, j1, m1) * CG(
                                jj, mm, c2, -p, j2, m2
                            )
                        am[index1, index2] = elem * sumPol

        index = len(self.savedAngularMatrix_matrix)

        c.execute(
            """ INSERT INTO pair_angularMatrix
                            VALUES (?,?, ?,?, ?,?, ?,?, ?)""",
            (l, j * 2, ll, jj * 2, l1, j1 * 2, l2, j2 * 2, index),
        )
        self.conn.commit()

        self.savedAngularMatrix_matrix.append(am)
        self.savedAngularMatrixChanged = True

        return am

    def __updateAngularMatrixElementsFile(self):
        if not (self.savedAngularMatrixChanged):
            return

        try:
            c = self.conn.cursor()
            c.execute("""SELECT * FROM pair_angularMatrix """)
            data = []
            for v in c.fetchall():
                data.append(v)

            data = np.array(data, dtype=np.float32)

            data[:, 1] /= 2.0  # 2 r j1 -> j1
            data[:, 3] /= 2.0  # 2 r j2 -> j2
            data[:, 5] /= 2.0  # 2 r j3 -> j3
            data[:, 7] /= 2.0  # 2 r j4 -> j4

            fileHandle = gzip.GzipFile(
                os.path.join(self.dataFolder, self.angularMatrixFile_meta), "wb"
            )
            np.save(fileHandle, data)
            fileHandle.close()
        except IOError:
            print(
                "Error while updating angularMatrix \
                data meta (description) File "
                + self.angularMatrixFile_meta
            )

        try:
            fileHandle = gzip.GzipFile(
                os.path.join(self.dataFolder, self.angularMatrixFile), "wb"
            )
            np.save(
                fileHandle,
                np.array(self.savedAngularMatrix_matrix, dtype=object),
            )
            fileHandle.close()
        except (IOError, ValueError) as e:
            print(
                "Error while updating angularMatrix \
                    data File "
                + self.angularMatrixFile
            )
            print(e)

    def __loadAngularMatrixElementsFile(self):
        try:
            fileHandle = gzip.GzipFile(
                os.path.join(self.dataFolder, self.angularMatrixFile_meta), "rb"
            )
            data = np.load(fileHandle, encoding="latin1", allow_pickle=True)
            fileHandle.close()
        except Exception as ex:
            print(ex)
            print("Note: No saved angular matrix metadata files to be loaded.")
            print(sys.exc_info())
            return

        data[:, 1] *= 2  # j1 -> 2 r j1
        data[:, 3] *= 2  # j2 -> 2 r j2
        data[:, 5] *= 2  # j3 -> 2 r j3
        data[:, 7] *= 2  # j4 -> 2 r j4

        data = np.array(np.rint(data), dtype=int)

        try:
            c = self.conn.cursor()
            c.executemany(
                """INSERT INTO pair_angularMatrix
                (l1, j1_x2 ,
                 l2 , j2_x2 ,
                 l3, j3_x2,
                 l4 , j4_x2 ,
                 ind)
                      VALUES (?,?,?,?,?,?,?,?,?)""",
                data,
            )

            self.conn.commit()

        except sqlite3.Error as e:
            print("Error while loading precalculated values into the database!")
            print(e)
            exit()
        if len(data) == 0:
            print("error")
            return

        try:
            fileHandle = gzip.GzipFile(
                os.path.join(self.dataFolder, self.angularMatrixFile), "rb"
            )
            self.savedAngularMatrix_matrix = np.load(
                fileHandle, encoding="latin1", allow_pickle=True
            ).tolist()
            fileHandle.close()
        except Exception as ex:
            print(ex)
            print("Note: No saved angular matrix files to be loaded.")
            print(sys.exc_info())

    def __isCoupled(self, n, l, j, nn, ll, jj, n1, l1, j1, n2, l2, j2, limit):
        if (
            (
                abs(
                    self.__getEnergyDefect(
                        n, l, j, nn, ll, jj, n1, l1, j1, n2, l2, j2
                    )
                )
                / C_h
                < limit
            )
            and not (
                n == n1
                and nn == n2
                and l == l1
                and ll == l2
                and j == j1
                and jj == j2
            )
            and not (
                (
                    abs(l1 - l) != 1
                    and (
                        (
                            abs(j - 0.5) < 0.1 and abs(j1 - 0.5) < 0.1
                        )  # j = 1/2 and j'=1/2 forbidden
                        or (
                            abs(j) < 0.1 and abs(j1 - 1) < 0.1
                        )  # j = 0 and j'=1 forbidden
                        or (
                            abs(j - 1) < 0.1 and abs(j1) < 0.1
                        )  # j = 1 and j'=0 forbidden
                    )
                )
                or (
                    abs(l2 - ll) != 1
                    and (
                        (
                            abs(jj - 0.5) < 0.1 and abs(j2 - 0.5) < 0.1
                        )  # j = 1/2 and j'=1/2 forbidden
                        or (
                            abs(jj) < 0.1 and abs(j2 - 1) < 0.1
                        )  # j = 0 and j'=1 forbidden
                        or (
                            abs(jj - 1) < 0.1 and abs(j2) < 0.1
                        )  # j = 1 and j'=0 forbidden
                    )
                )
            )
            and not (abs(j) < 0.1 and abs(j1) < 0.1)  # j = 0 and j'=0 forbiden
            and not (abs(jj) < 0.1 and abs(j2) < 0.1)
            and not (
                abs(l) < 0.1 and abs(l1) < 0.1
            )  # l = 0 and l' = 0 is forbiden
            and not (abs(ll) < 0.1 and abs(l2) < 0.1)
        ):
            # determine coupling
            dl = abs(l - l1)
            dj = abs(j - j1)
            c1 = 0
            if dl == 1 and (dj < 1.1):
                c1 = 1  # dipole coupling
            elif (
                (dl == 0 or dl == 2 or dl == 1)
                and (dj < 2.1)
                and (2 <= self.interactionsUpTo)
            ):
                c1 = 2  # quadrupole coupling
            else:
                return False
            dl = abs(ll - l2)
            dj = abs(jj - j2)
            c2 = 0
            if dl == 1 and (dj < 1.1):
                c2 = 1  # dipole coupling
            elif (
                (dl == 0 or dl == 2 or dl == 1)
                and (dj < 2.1)
                and (2 <= self.interactionsUpTo)
            ):
                c2 = 2  # quadrupole coupling
            else:
                return False
            return c1 + c2
        else:
            return False

    def __getEnergyDefect(self, n, l, j, nn, ll, jj, n1, l1, j1, n2, l2, j2):
        """
        Energy defect between |n,l,j>x|nn,ll,jj> state and |n1,l1,j1>x|n1,l1,j1>
        state of atom1 and atom2 in respective spins states s1 and s2

        Takes spin vales s1 and s2 as the one defined when defining calculation.

        Parameters:
            n (int): principal quantum number
            l (int): orbital angular momenutum
            j (float): total angular momentum
            nn (int): principal quantum number
            ll (int): orbital angular momenutum
            jj (float): total angular momentum
            n1 (int): principal quantum number
            l1 (int): orbital angular momentum
            j1 (float): total angular momentum
            n2 (int): principal quantum number
            l2 (int): orbital angular momentum
            j2 (float): total angular momentum

        Returns:
            float:  energy defect (SI units: J)
        """
        return C_e * (
            self.atom1.getEnergy(n1, l1, j1, s=self.s1)
            + self.atom2.getEnergy(n2, l2, j2, s=self.s2)
            - self.atom1.getEnergy(n, l, j, s=self.s1)
            - self.atom2.getEnergy(nn, ll, jj, s=self.s2)
        )

    def __makeRawMatrix2(
        self,
        n,
        l,
        j,
        nn,
        ll,
        jj,
        k,
        lrange,
        limit,
        limitBasisToMj,
        progressOutput=False,
        debugOutput=False,
    ):
        # limit = limit in Hz on energy defect
        # k defines range of n' = [n-k, n+k]
        dimension = 0

        # which states/channels contribute significantly in the second order perturbation?
        states = []

        # original pairstate index
        opi = 0

        # this numbers are conserved if we use only dipole-dipole interactions
        Lmod2 = (l + ll) % 2

        l1start = max(l - self.interactionsUpTo, 0)
        l2start = max(ll - self.interactionsUpTo, 0)

        if debugOutput:
            print("\n ======= Relevant states =======\n")

        for n1 in xrange(max(n - k, 1), n + k + 1):
            for n2 in xrange(max(nn - k, 1), nn + k + 1):
                l1max = max(l + self.interactionsUpTo, lrange) + 1
                l1max = min(l1max, n1 - 1)
                for l1 in xrange(l1start, l1max):
                    l2max = max(ll + self.interactionsUpTo, lrange) + 1
                    l2max = min(l2max, n2 - 1)
                    for l2 in xrange(l2start, l2max):
                        j1 = l1 - self.s1
                        while j1 < -0.1:
                            j1 += 2 * self.s1
                        while j1 <= l1 + self.s1 + 0.1:
                            j2 = l2 - self.s2
                            while j2 < -0.1:
                                j2 += 2 * self.s2

                            while j2 <= l2 + self.s2 + 0.1:
                                ed = (
                                    self.__getEnergyDefect(
                                        n,
                                        l,
                                        j,
                                        nn,
                                        ll,
                                        jj,
                                        n1,
                                        l1,
                                        j1,
                                        n2,
                                        l2,
                                        j2,
                                    )
                                    / C_h
                                )
                                if (
                                    abs(ed) < limit
                                    and (
                                        not (self.interactionsUpTo == 1)
                                        or (Lmod2 == ((l1 + l2) % 2))
                                    )
                                    and (
                                        (not limitBasisToMj)
                                        or (j1 + j2 + 0.1 > self.m1 + self.m2)
                                    )
                                    and (
                                        n1 >= self.atom1.groundStateN
                                        or [n1, l1, j1]
                                        in self.atom1.extraLevels
                                    )
                                    and (
                                        n2 >= self.atom2.groundStateN
                                        or [n2, l2, j2]
                                        in self.atom2.extraLevels
                                    )
                                ):
                                    if debugOutput:
                                        pairState = (
                                            "|"
                                            + printStateString(
                                                n1, l1, j1, s=self.s1
                                            )
                                            + ","
                                            + printStateString(
                                                n2, l2, j2, s=self.s2
                                            )
                                            + ">"
                                        )
                                        print(
                                            pairState
                                            + (
                                                "\t EnergyDefect = %.3f GHz"
                                                % (ed * 1.0e-9)
                                            )
                                        )

                                    states.append([n1, l1, j1, n2, l2, j2])

                                    if (
                                        n == n1
                                        and nn == n2
                                        and l == l1
                                        and ll == l2
                                        and j == j1
                                        and jj == j2
                                    ):
                                        opi = dimension

                                    dimension = dimension + 1
                                j2 = j2 + 1.0
                            j1 = j1 + 1.0

        if debugOutput:
            print("\tMatrix dimension\t=\t", dimension)

        # mat_value, mat_row, mat_column for each sparce matrix describing
        # dipole-dipole, dipole-quadrupole (and quad-dipole) and quadrupole-quadrupole
        couplingMatConstructor = [
            [[], [], []] for i in xrange(2 * self.interactionsUpTo - 1)
        ]

        # original pair-state (i.e. target pair state) Zeeman Shift
        opZeemanShift = (
            (
                self.atom1.getZeemanEnergyShift(
                    self.l, self.j, self.m1, self.Bz, s=self.s1
                )
                + self.atom2.getZeemanEnergyShift(
                    self.ll, self.jj, self.m2, self.Bz, s=self.s2
                )
            )
            / C_h
            * 1.0e-9
        )  # in GHz

        if debugOutput:
            print("\n ======= Coupling strengths (radial part only) =======\n")

        maxCoupling = "quadrupole-quadrupole"
        if self.interactionsUpTo == 1:
            maxCoupling = "dipole-dipole"
        if debugOutput:
            print(
                "Calculating coupling (up to ",
                maxCoupling,
                ") between the pair states",
            )

        for i in xrange(dimension):
            ed = (
                self.__getEnergyDefect(
                    states[opi][0],
                    states[opi][1],
                    states[opi][2],
                    states[opi][3],
                    states[opi][4],
                    states[opi][5],
                    states[i][0],
                    states[i][1],
                    states[i][2],
                    states[i][3],
                    states[i][4],
                    states[i][5],
                )
                / C_h
                * 1.0e-9
                - opZeemanShift
            )

            pairState1 = (
                "|"
                + printStateString(
                    states[i][0], states[i][1], states[i][2], s=self.s1
                )
                + ","
                + printStateString(
                    states[i][3], states[i][4], states[i][5], s=self.s2
                )
                + ">"
            )

            states[i].append(ed)  # energy defect of given state

            for j in xrange(i + 1, dimension):
                coupled = self.__isCoupled(
                    states[i][0],
                    states[i][1],
                    states[i][2],
                    states[i][3],
                    states[i][4],
                    states[i][5],
                    states[j][0],
                    states[j][1],
                    states[j][2],
                    states[j][3],
                    states[j][4],
                    states[j][5],
                    limit,
                )

                if states[i][0] == 24 and states[j][0] == 18:
                    print("\n")
                    print(states[i])
                    print(states[j])
                    print(coupled)

                if coupled and (
                    abs(states[i][0] - states[j][0]) <= k
                    and abs(states[i][3] - states[j][3]) <= k
                ):
                    if debugOutput:
                        pairState2 = (
                            "|"
                            + printStateString(
                                states[j][0],
                                states[j][1],
                                states[j][2],
                                s=self.s1,
                            )
                            + ","
                            + printStateString(
                                states[j][3],
                                states[j][4],
                                states[j][5],
                                s=self.s2,
                            )
                            + ">"
                        )
                        print(pairState1 + " <---> " + pairState2)

                    couplingStregth = (
                        _atomLightAtomCoupling(
                            states[i][0],
                            states[i][1],
                            states[i][2],
                            states[i][3],
                            states[i][4],
                            states[i][5],
                            states[j][0],
                            states[j][1],
                            states[j][2],
                            states[j][3],
                            states[j][4],
                            states[j][5],
                            self.atom1,
                            atom2=self.atom2,
                            s=self.s1,
                            s2=self.s2,
                        )
                        / C_h
                        * 1.0e-9
                    )

                    couplingMatConstructor[coupled - 2][0].append(
                        couplingStregth
                    )
                    couplingMatConstructor[coupled - 2][1].append(i)
                    couplingMatConstructor[coupled - 2][2].append(j)

                    exponent = coupled + 1
                    if debugOutput:
                        print(
                            (
                                "\tcoupling (C_%d/R^%d) = %.5f"
                                % (
                                    exponent,
                                    exponent,
                                    couplingStregth * (1e6) ** (exponent),
                                )
                            ),
                            "/R^",
                            exponent,
                            " GHz  (mu m)^",
                            exponent,
                            "\n",
                        )

        # coupling = [1,1] dipole-dipole, [2,1]  quadrupole dipole, [2,2] quadrupole quadrupole

        couplingMatArray = [
            csr_matrix(
                (
                    couplingMatConstructor[i][0],
                    (
                        couplingMatConstructor[i][1],
                        couplingMatConstructor[i][2],
                    ),
                ),
                shape=(dimension, dimension),
            )
            for i in xrange(len(couplingMatConstructor))
        ]
        return states, couplingMatArray

    def __initializeDatabaseForMemoization(self):
        # memoization of angular parts
        self.conn = sqlite3.connect(
            os.path.join(self.dataFolder, "precalculated_pair.db")
        )
        c = self.conn.cursor()

        # ANGULAR PARTS
        c.execute("""DROP TABLE IF EXISTS pair_angularMatrix""")
        c.execute(
            """SELECT COUNT(*) FROM sqlite_master
                            WHERE type='table' AND name='pair_angularMatrix';"""
        )
        if c.fetchone()[0] == 0:
            # create table
            try:
                c.execute(
                    """CREATE TABLE IF NOT EXISTS pair_angularMatrix
                 (l1 TINYINT UNSIGNED, j1_x2 TINYINT UNSIGNED,
                 l2 TINYINT UNSIGNED, j2_x2 TINYINT UNSIGNED,
                 l3 TINYINT UNSIGNED, j3_x2 TINYINT UNSIGNED,
                 l4 TINYINT UNSIGNED, j4_x2 TINYINT UNSIGNED,
                 ind INTEGER,
                 PRIMARY KEY (l1,j1_x2, l2,j2_x2, l3,j3_x2, l4,j4_x2)
                ) """
                )
            except sqlite3.Error as e:
                print(e)
            self.conn.commit()
        self.__loadAngularMatrixElementsFile()
        self.savedAngularMatrixChanged = False

    def __closeDatabaseForMemoization(self):
        self.conn.commit()
        self.conn.close()
        self.conn = False

    def getLeRoyRadius(self):
        """
        Returns Le Roy radius for initial pair-state.

        Le Roy radius [#leroy]_ is defined as
        :math:`2(\\langle r_1^2 \\rangle^{1/2} + \\langle r_2^2 \\rangle^{1/2})`,
        where :math:`r_1` and :math:`r_2` are electron coordinates for the
        first and the second atom in the initial pair-state.
        Below this radius, calculations are not valid since electron
        wavefunctions start to overlap.

        Returns:
            float: LeRoy radius measured in :math:`\\mu m`

        References:
            .. [#leroy] R.J. Le Roy, Can. J. Phys. **52**, 246 (1974)
                http://www.nrcresearchpress.com/doi/abs/10.1139/p74-035
        """
        step = 0.001
        r1, psi1_r1 = self.atom1.radialWavefunction(
            self.l,
            0.5,
            self.j,
            self.atom1.getEnergy(self.n, self.l, self.j, s=self.s1)
            / 27.211_386_245_981,
            self.atom1.alphaC ** (1 / 3.0),
            2.0 * self.n * (self.n + 15.0),
            step,
        )

        sqrt_r1_on2 = trapezoid(
            np.multiply(np.multiply(psi1_r1, psi1_r1), np.multiply(r1, r1)),
            x=r1,
        )

        r2, psi2_r2 = self.atom2.radialWavefunction(
            self.ll,
            0.5,
            self.jj,
            self.atom2.getEnergy(self.nn, self.ll, self.jj, s=self.s2)
            / 27.211_386_245_981,
            self.atom2.alphaC ** (1 / 3.0),
            2.0 * self.nn * (self.nn + 15.0),
            step,
        )

        sqrt_r2_on2 = trapezoid(
            np.multiply(np.multiply(psi2_r2, psi2_r2), np.multiply(r2, r2)),
            x=r2,
        )

        return (
            2.0
            * (sqrt(sqrt_r1_on2) + sqrt(sqrt_r2_on2))
            * (physical_constants["Bohr radius"][0] * 1.0e6)
        )

    def getC6perturbatively(
        self, theta, phi, nRange, energyDelta, degeneratePerturbation=False
    ):
        r"""
        Calculates :math:`C_6` from second order perturbation theory.

        Calculates
        :math:`C_6=\sum_{\rm r',r''}|\langle {\rm r',r''}|V|\
        {\rm r1,r2}\rangle|^2/\Delta_{\rm r',r''}`, where
        :math:`\Delta_{\rm r',r''}\equiv E({\rm r',r''})-E({\rm r1, r2})`
        When second order perturbation couples to multiple energy degenerate
        states, users shold use **degenerate perturbation calculations** by
        setting `degeneratePerturbation=True` .

        This calculation is faster then full diagonalization, but it is valid
        only far from the so called spaghetti region that occurs when atoms
        are close to each other. In that region multiple levels are strongly
        coupled, and one needs to use full diagonalization. In region where
        perturbative calculation is correct, energy level shift can be
        obtained as :math:`V(R)=-C_6/R^6`

        See `perturbative C6 calculations example snippet`_ and for
        degenerate perturbation calculation see
        `degenerate pertubation C6 calculation example snippet`_

        .. _`perturbative C6 calculations example snippet`:
            ./Rydberg_atoms_a_primer.html#Dispersion-Coefficients

       .. _`degenerate pertubation C6 calculation example snippet`:
           ./ARC_3_0_introduction.html#Pertubative-C6-calculation-in-the-manifold-of-degenerate-states

        Parameters:
            theta (float):
                orientation of inter-atomic axis with respect
                to quantization axis (:math:`z`) in Euler coordinates
                (measured in units of radian)
            phi (float):
                orientation of inter-atomic axis with respect
                to quantization axis (:math:`z`) in Euler coordinates
                (measured in units of radian)
            nRange (int):
                how much below and above the given principal quantum number
                of the pair state we should be looking
            energyDelta (float):
                what is maximum energy difference ( :math:`\Delta E/h` in Hz)
                between the original pair state and the other pair states that we are including in
                calculation
            degeneratePerturbation (bool):
                optional, default False. Should one
                use degenerate perturbation theory. This should be used whenever
                angle between quantisation and interatomic axis is non-zero,
                as well as when one considers non-stretched states.

        Returns:
            float: if **degeneratePerturbation=False**, returns
            :math:`C_6` measured in :math:`\text{GHz }\mu\text{m}^6`;
            if **degeneratePerturbation=True**, returns array of
            :math:`C_6` measured in :math:`\text{GHz }\mu\text{m}^6`
            AND array of corresponding eigenvectors in
            :math:`\{m_{j_1}=-j_1, \ldots, m_{j_1} = +j1\}\bigotimes \
            \{ m_{j_2}=-j_2, \ldots, m_{j_2} = +j2\}`
            basis


        Example:
            If we want to quickly calculate :math:`C_6` for two Rubidium
            atoms in state :math:`62 D_{3/2} m_j=3/2`, positioned in space
            along the shared quantization axis::

                from arc import *
                calculation = PairStateInteractions(Rubidium(), 62, 2, 1.5, 62, 2, 1.5, 1.5, 1.5)
                c6 = calculation.getC6perturbatively(0,0, 5, 25e9)
                print "C_6 = %.0f GHz (mu m)^6" % c6

            Which returns::

                C_6 = 767 GHz (mu m)^6

            Quick calculation of angular anisotropy of for Rubidium
            :math:`D_{2/5},m_j=5/2` states::

                # Rb 60 D_{2/5}, mj=2.5 , 60 D_{2/5}, mj=2.5 pair state
                calculation1 = PairStateInteractions(Rubidium(), 60, 2, 2.5, 60, 2, 2.5, 2.5, 2.5)
                # list of atom orientations
                thetaList = np.linspace(0,pi,30)
                # do calculation of C6 pertubatively for all atom orientations
                c6 = []
                for theta in thetaList:
                    value = calculation1.getC6perturbatively(theta,0,5,25e9)
                    c6.append(value)
                    print ("theta = %.2f * pi \tC6 = %.2f GHz  mum^6" % (theta/pi,value))
                # plot results
                plot(thetaList/pi,c6,"b-")
                title("Rb, pairstate  60 $D_{5/2},m_j = 5/2$, 60 $D_{5/2},m_j = 5/2$")
                xlabel(r"$\Theta /\pi$")
                ylabel(r"$C_6$ (GHz $\mu$m${}^6$")
                show()

        """
        self.__initializeDatabaseForMemoization()

        # ========= START OF THE MAIN CODE ===========

        # wigner D matrix allows calculations with arbitrary orientation of
        # the two atoms
        wgd = WignerDmatrix(theta, phi)

        # any conservation?
        # this numbers are conserved if we use only dipole-dipole interactions
        Lmod2 = (self.l + self.ll) % 2

        # find nearby states

        lmin1 = self.l - 1
        if lmin1 < -0.1:
            lmin1 = 1
        lmin2 = self.ll - 1
        if lmin2 < -0.1:
            lmin2 = 1

        interactionMatrix = np.zeros(
            (
                round((2 * self.j + 1) * (2 * self.jj + 1)),
                round((2 * self.j + 1) * (2 * self.jj + 1)),
            ),
            dtype=np.complex128,
        )

        for n1 in xrange(max(self.n - nRange, 1), self.n + nRange + 1):
            for n2 in xrange(max(self.nn - nRange, 1), self.nn + nRange + 1):
                lmax1 = min(self.l + 2, n1)
                for l1 in xrange(lmin1, lmax1, 2):
                    lmax2 = min(self.ll + 2, n2)
                    for l2 in xrange(lmin2, lmax2, 2):
                        if (l1 + l2) % 2 == Lmod2:
                            j1 = l1 - self.s1
                            while j1 < -0.1:
                                j1 += 2 * self.s1
                            while j1 <= l1 + self.s1 + 0.1:
                                j2 = l2 - self.s2
                                while j2 < -0.1:
                                    j2 += 2 * self.s2

                                while j2 <= l2 + self.s2 + 0.1:
                                    coupled = self.__isCoupled(
                                        self.n,
                                        self.l,
                                        self.j,
                                        self.nn,
                                        self.ll,
                                        self.jj,
                                        n1,
                                        l1,
                                        j1,
                                        n2,
                                        l2,
                                        j2,
                                        energyDelta,
                                    )
                                    if (
                                        coupled
                                        and (
                                            not (self.interactionsUpTo == 1)
                                            or (Lmod2 == ((l1 + l2) % 2))
                                        )
                                        and (
                                            n1 >= self.atom1.groundStateN
                                            or [n1, l1, j1]
                                            in self.atom1.extraLevels
                                        )
                                        and (
                                            n2 >= self.atom2.groundStateN
                                            or [n2, l2, j2]
                                            in self.atom2.extraLevels
                                        )
                                    ):
                                        energyDefect = (
                                            self.__getEnergyDefect(
                                                self.n,
                                                self.l,
                                                self.j,
                                                self.nn,
                                                self.ll,
                                                self.jj,
                                                n1,
                                                l1,
                                                j1,
                                                n2,
                                                l2,
                                                j2,
                                            )
                                            / C_h
                                        )
                                        energyDefect = (
                                            energyDefect * 1.0e-9
                                        )  # GHz
                                        if abs(energyDefect) < 1e-10:
                                            raise ValueError(
                                                "The requested pair-state "
                                                "is dipole coupled resonatly "
                                                "(energy defect = 0)"
                                                "to other pair-states"
                                                "Aborting pertubative "
                                                "calculation."
                                                "(This usually happens for "
                                                "high-L states for which "
                                                "identical quantum defects give "
                                                "raise to degeneracies, making "
                                                "total L ultimately not "
                                                "conserved quantum number) "
                                            )

                                        # calculate radial part
                                        couplingStregth = (
                                            _atomLightAtomCoupling(
                                                self.n,
                                                self.l,
                                                self.j,
                                                self.nn,
                                                self.ll,
                                                self.jj,
                                                n1,
                                                l1,
                                                j1,
                                                n2,
                                                l2,
                                                j2,
                                                self.atom1,
                                                atom2=self.atom2,
                                                s=self.s1,
                                                s2=self.s2,
                                            )
                                            * (1.0e-9 * (1.0e6) ** 3 / C_h)
                                        )  # GHz / mum^3

                                        d = self.__getAngularMatrix_M(
                                            self.l,
                                            self.j,
                                            self.ll,
                                            self.jj,
                                            l1,
                                            j1,
                                            l2,
                                            j2,
                                        )

                                        interactionMatrix += (
                                            d.conj().T.dot(d)
                                            * abs(couplingStregth) ** 2
                                            / energyDefect
                                        )

                                    j2 = j2 + 1.0
                                j1 = j1 + 1.0

        rotationMatrix = np.kron(
            wgd.get(self.j).toarray(), wgd.get(self.jj).toarray()
        )

        interactionMatrix = rotationMatrix.dot(
            interactionMatrix.dot(rotationMatrix.conj().T)
        )
        # ========= END OF THE MAIN CODE ===========
        self.__closeDatabaseForMemoization()

        value, vectors = np.linalg.eigh(interactionMatrix)
        vectors = vectors.T
        stateCom = compositeState(
            singleAtomState(self.j, self.m1), singleAtomState(self.jj, self.m2)
        ).T

        if not degeneratePerturbation:
            for i, v in enumerate(vectors):
                if abs(np.vdot(v, stateCom)) > 1 - 1e-9:
                    return value[i]
            #    else:
            #        print(np.vdot(v, stateCom))
            # if initial state is not eigen state print warning and return
            # results for eigenstates, and eigenstate composition
            """
            print("WARNING: Requested state is not eigenstate when dipole-dipole "
                  "interactions and/or relative position of atoms are "
                  "taken into account.\n"
                  "We will use degenerate pertubative theory to correctly "
                  "calculate C6.\n"
                  "Method will return values AND eigenvectors in basis \n"
                  "{mj1 = -j1, ... , mj1 = +j1} x {mj2 = -j2, ... , m2 = +j2}, "
                  "where x denotes Kronecker product\n"
                  "To not see this warning request explicitly "
                  "degeneratePerturbation=True in call of this method.\n")
            """
            # print(stateCom.conj().dot(interactionMatrix.dot(stateCom.T)))
            # print(stateCom.conj().dot(interactionMatrix.dot(stateCom.T)).shape)
            return np.real(
                stateCom.conj().dot(interactionMatrix.dot(stateCom.T))[0][0]
            )
        return np.real(value), vectors

    def _getd(self, l, j, ll, jj, l1, j1, l2, j2):
        r"""
        Gets the mj-resolved matrix for the transition weights.
        Note that this function is slow due to database initialisation.
        Only use if necessary.

        Args:
            l, j (floats) -   atom 1, initial state orbital and total angular momentum
            ll, jj (floats) - atom 2, initial state orbital and total angular momentum
            l1, j1 (floats) - atom 1, final state orbital and total angular momentum
            l2, j2 (floats) - atom 2, final state orbital and total angular momentum

        Output:
            d (ndarray) - mj-resolved matrix for transition weights
        """
        self.__initializeDatabaseForMemoization()
        d = self.__getAngularMatrix_M(l, j, ll, jj, l1, j1, l2, j2)
        self.__closeDatabaseForMemoization()
        return d

    def _getAngularBasisRotationMatrix(self, j1, j2):
        r"""
        Returns basis rotation matrix when interchanging state of atom1 and atom2, e.g. in hopping process.

        basis2 = rot * basis1 => rot = basis2 * (basis1)^{-1} = basis2 * (1)^{-1} = basis2

        Args:
            j1 (float) - total orbital angular momentum of atom 1, half integer or integer >= 0
            j2 (float) - total orbital angular momentum of atom 2, half integer or integer >= 0

        Output:
            basis2 (ndarray) -  (2*j1+1)*(2*j2+1) \times (2*j1+1)*(2*j2+1) fine-structure (mj)
                                basis rotation matrix
        """
        basis2 = np.zeros(
            (
                int((2 * j1 + 1) * (2 * j2 + 1)),
                int((2 * j1 + 1) * (2 * j2 + 1)),
            ),
            dtype=np.complex128,
        )
        for i, mj2 in enumerate(np.arange(-j2, j2 + 1, 1)):
            for j, mj1 in enumerate(np.arange(-j1, j1 + 1, 1)):
                basis2[int(i * (2 * j1 + 1) + j), :] = compositeState(
                    singleAtomState(j1, mj1), singleAtomState(j2, mj2)
                ).T
        return basis2

    def _ljCoupledCheck(self, l, j, l1, j1, s):
        r"""
        Checks if a pair of angular momentum quantum numbers (l,j) and (l1, j1) is
        coupled via dipole or up to quadrupole transitions
        (for `self.interactionsUpTo=1` and `self.interactionsUpTo=1` respectively).

        Args:
            l, j (floats) -   angular momentum quantum number of initial state
            l1, j1 (floats) - angular momentum quantum number of final state
            s (float) -       spin angular momentum of the atom

        Output:
            boolean - True or False
        """
        if (
            not (l == l1 and j == j1)
            and not (
                abs(j) < 0.1 and abs(j1) < 0.1
            )  # j = 0 and j'=0 always forbiden
            and not (
                abs(l) < 0.1 and abs(l1) < 0.1
            )  # l = 0 and l' = 0 always forbiden
            and (
                abs(round(l - j)) < s + 0.1
            )  # check for total angular momentum
            and (abs(round(l1 - j1)) < s + 0.1)
        ):
            # check if the pair is dipole coupled
            # if so, then all is okay and no further selection rules need be enforced
            if (abs(round(l - l1)) == 1) and (abs(round(j - j1)) in [0, 1]):
                return True
            # if not dipole coupled, check if quadrupole coupling was allowed
            # and enforce additional quadrupole selection rules
            elif (
                (self.interactionsUpTo == 2)
                and (abs(round(l - l1)) in [0, 2])
                and (abs(round(j - j1)) in [0, 1, 2])
                and not (
                    abs(j) < 0.1 and abs(round(j1 - 1)) < 0.1
                )  # j=0 to j1=1 forbidden
                and not (
                    abs(round(j - 1)) < 0.1 and abs(j1) < 0.1
                )  # j=1 to j1=0 forbidden
                and not (
                    abs(round(j - 0.5)) < 0.1 and abs(round(j1 - 0.5)) < 0.1
                )  # j=1/2 to j1=1/2 forbidden
                and not (
                    abs(l) < 0.1 and abs(round(l1 - 1)) < 0.1
                )  # l=0 to l1=1 forbidden
                and not (
                    abs(round(l - 1)) < 0.1 and abs(l1) < 0.1
                )  # l=1 to l1=0 forbidden
            ):
                return True
            else:
                return False
        else:
            return False

    def _findAllCoupledAngularMomentumStates(
        self, l, j, s1, ll, jj, s2, stateHopping=False
    ):
        r"""
        Finds all second-order coupled angular momentum states for an initial pair
        configuration (l,j; ll,jj) --> (l1,j1; l2,j2) --> (l',j'; ll',jj').
        If hopping == False, (l',j'; ll',jj') = (l,j; ll,jj)
        Elif hopping == True, (l',j'; ll',jj') = (ll,jj; l,j)

        Args:
            l, j, s1 (floats) -     angular momentum quantum numbers of first atom
            ll, jj, s2 (floats) -   angular momentum quantum numbers of second atom
            hopping (boolean) -     determines whether the final angMomentum configuration
                                    is equal to the initial one or if the states are
                                    interchanged between the atoms ('state hopping')

        Output:
            coupledStates (list) -  list of tuples containing the coupled angular
                                    momentum configurations in the form
                                    (l,j, ll,jj, l1,j1, l2,j2, l',j', ll',jj')
        """
        coupledStates = []

        # iterate through potential states for atom 1
        for l1 in range(
            max(0, l - self.interactionsUpTo), l + self.interactionsUpTo + 1
        ):
            for j1 in np.arange(
                max(s1, j - self.interactionsUpTo),
                j + self.interactionsUpTo + 1,
                1,
            ):
                # check if angular momentum coupling is valid
                if self._ljCoupledCheck(l, j, l1, j1, s1):
                    # iterate through potential states for atom 2
                    for l2 in range(
                        max(0, ll - self.interactionsUpTo),
                        ll + self.interactionsUpTo + 1,
                    ):
                        for j2 in np.arange(
                            max(s2, jj - self.interactionsUpTo),
                            jj + self.interactionsUpTo + 1,
                            1,
                        ):
                            # check if angular momentum coupling is valid
                            if self._ljCoupledCheck(ll, jj, l2, j2, s2):
                                j1, j2 = float(j1), float(j2)
                                # atoms each return into their initial state, respectively
                                if not stateHopping:
                                    # append state to list, finalState=initialState is certainly coupled
                                    coupledStates.append(
                                        (
                                            l,
                                            j,
                                            ll,
                                            jj,
                                            l1,
                                            j1,
                                            l2,
                                            j2,
                                            l,
                                            j,
                                            ll,
                                            jj,
                                        )
                                    )
                                # atoms return to swappd states, check if these do couple
                                elif (
                                    stateHopping
                                    # and (abs(j1-jj) <= 1) and (abs(j2-j) <= 1)
                                    and self._ljCoupledCheck(l1, j1, ll, jj, s1)
                                    and self._ljCoupledCheck(l2, j2, l, j, s2)
                                ):
                                    # append state to list, final state is swapped w.r.t. initial state
                                    coupledStates.append(
                                        (
                                            l,
                                            j,
                                            ll,
                                            jj,
                                            l1,
                                            j1,
                                            l2,
                                            j2,
                                            ll,
                                            jj,
                                            l,
                                            j,
                                        )
                                    )
        return coupledStates

    def _getC6contributions_lj(self, nRange, energyDelta, stateHopping=False):
        r"""
        Returns the interaction strengths for the different
        (l,j; ll,jj) --> (l1,j1; l2,j2) --> (l',j'; ll',jj') configurations.

        Args:
            nRange (int) -  how much below and above the given principal quantum number
                            of the pair state we should be looking
            energyDelta (float) -   what is maximum energy difference ( :math:`\Delta E/h` in Hz)
                                    between the original pair state and the other pair states that
                                    we are including in the calculation
            stateHopping (bool) -   whether or not the final state is interchanged ('hopped')
                                    w.r.t. the initial state
        Output:
            ljInteractions (list) - list containing entries of the form
                                    [(l,j, ll,jj, l1,j1, l2,j2, l',j' ll',jj'), V_{lj}]
                                    with V_{lj} the interaction strength for the given
                                    configuration in GHz(um)^6
        """

        ljInteractions = []

        # find all (l,j; ll,jj) --> (l1,j1; l2,j2) --> (l',j'; ll',jj') pairs
        coupledStates = self._findAllCoupledAngularMomentumStates(
            self.l,
            self.j,
            self.s1,
            self.ll,
            self.jj,
            self.s2,
            stateHopping=stateHopping,
        )

        for lj in coupledStates:
            V_lj = 0
            # unpack angular momentum info
            [l1, j1, ll1, jj1, l2, j2, ll2, jj2, l3, j3, ll3, jj3] = list(lj)
            # iterate through n1 states
            for n2 in range(max(self.n - nRange, 1), self.n + nRange + 1):
                # iterate through n2 states
                for nn2 in range(
                    max(self.nn - nRange, 1), self.nn + nRange + 1
                ):
                    # to check if nVals are okay
                    nCheck = (
                        n2 >= self.atom1.groundStateN
                        or [n2, l2, j2] in self.atom1.extraLevels
                    ) and (
                        nn2 >= self.atom2.groundStateN
                        or [nn2, ll2, jj2] in self.atom2.extraLevels
                    )
                    if stateHopping:
                        nCheck = (
                            nCheck
                            and (
                                self.nn >= self.atom1.groundStateN
                                or [self.nn, l3, j3] in self.atom1.extraLevels
                            )
                            and (
                                self.n >= self.atom2.groundStateN
                                or [self.n, ll3, jj3] in self.atom2.extraLevels
                            )
                        )

                    # calculate energy defect
                    energyDefect = (
                        self.__getEnergyDefect(
                            self.n,
                            l1,
                            j1,
                            self.nn,
                            ll1,
                            jj1,
                            n2,
                            l2,
                            j2,
                            nn2,
                            ll2,
                            jj2,
                        )
                        / C_h
                    )
                    energyDefect = energyDefect * 1e-9  # GHz
                    if abs(energyDefect) < 1e-10:
                        print(n2, l2, j2, nn2, ll2, jj2, stateHopping, "error")
                        raise ValueError(
                            "The requested pair-state "
                            "is dipole coupled resonatly "
                            "(energy defect = 0) "
                            "to other pair-states. "
                            "Aborting pertubative "
                            "calculation. "
                            "(This usually happens for "
                            "high-L states for which "
                            "identical quantum defects give "
                            "raise to degeneracies, making "
                            "total L ultimately not "
                            "conserved quantum number) "
                        )

                    # proceed only if energy defect is within limit and nCheck was positive
                    if (abs(energyDefect) < energyDelta * 10**-9) and nCheck:
                        # calculate radial overlaps
                        couplingStrength1 = _atomLightAtomCoupling(
                            self.n,
                            l1,
                            j1,
                            self.nn,
                            ll1,
                            jj1,
                            n2,
                            l2,
                            j2,
                            nn2,
                            ll2,
                            jj2,
                            self.atom1,
                            atom2=self.atom2,
                            s=self.s1,
                            s2=self.s2,
                        ) * (1.0e-9 * (1.0e6) ** 3 / C_h)  # GHz / mum^3
                        if not stateHopping:
                            couplingStrength2 = couplingStrength1
                        else:
                            couplingStrength2 = _atomLightAtomCoupling(
                                n2,
                                l2,
                                j2,
                                nn2,
                                ll2,
                                jj2,
                                self.nn,
                                l3,
                                j3,
                                self.n,
                                ll3,
                                jj3,
                                self.atom1,
                                atom2=self.atom2,
                                s=self.s1,
                                s2=self.s2,
                            ) * (1.0e-9 * (1.0e6) ** 3 / C_h)  # GHz / mum^3

                        V_lj += (
                            abs(couplingStrength1 * couplingStrength2)
                            / energyDefect
                        )  # GHz um^6
            ljInteractions.append([*list(lj), float(V_lj)])
        return ljInteractions

    def _getPerturbativeC6Matrix_lj(self, ljInteractions):
        r"""
        Construct full Imat from lj, V_{lj} information.

        Args:
            ljInteractions (list) - list contains entries of the form
                                    [l,j, ll,jj, l1,j1, l2,j2, l',j', ll',jj', V_{lj}]
                                    Only those angular channels contained in the list are
                                    included in the resulting Imat. So make sure you pass
                                    a complete list to this function.

        Output:
            Imat (ndarray) -    interaction matrix with mj-basis resolution as reconstructed
                                from ljInteraction list
        """
        # open database
        self.__initializeDatabaseForMemoization()
        Imat = 0
        # iterate through channels
        for vals in ljInteractions:
            d1 = self.__getAngularMatrix_M(*vals[0:8])
            d2 = self.__getAngularMatrix_M(*vals[4:12])
            # no need to take the hermitian conjugate of d2 here as this code uses the right order
            # of l,j, l1,j1 as opposed to the original code
            Imat += vals[-1] * d2.dot(d1)
        # close database
        self.__closeDatabaseForMemoization()
        return np.array(Imat)

    def _getInteractionMatrix_lj(self, nRange, energyDelta):
        r"""
        Small helper function to get the interaction matrix from the d_lj method, used for debugging
        and double-checking. Can also be used to call the interaction matrix via this method - but only
        for the case theta = phi = 0. Also, it automatically builds the full Imat from all four blocks
        if atomState1 != atomState2.

        Args:
            nRange (int) -  how much below and above the given principal quantum number
                            of the pair state we should be looking
            energyDelta (float) -   what is maximum energy difference ( :math:`\Delta E/h` in Hz)
                                    between the original pair state and the other pair states that
                                    we are including in the calculation

        Output:
            Imat (ndarray) - full interaction matrix calculated via d_lj-method for theta=phi=0.
        """
        interaction11 = self._getC6contributions_lj(
            nRange, energyDelta, stateHopping=False
        )
        Imat11 = self._getPerturbativeC6Matrix_lj(interaction11)
        if (
            self.n == self.nn
            and self.l == self.ll
            and self.j == self.jj
            and self.s1 == self.s2
        ):
            return Imat11
        else:
            # get rotation matrix for basis changes from [atomState2,atomState1] --> [atomState1,atomState2]
            basisRotationMatrix12 = self._getAngularBasisRotationMatrix(
                self.jj, self.j
            )
            # get second Imat block
            interaction21 = self._getC6contributions_lj(
                nRange, energyDelta, stateHopping=True
            )
            Imat21 = self._getPerturbativeC6Matrix_lj(interaction21)
            # put all together
            Imat = np.block(
                [
                    [
                        Imat11,
                        basisRotationMatrix12.dot(
                            Imat21.dot(basisRotationMatrix12)
                        ),
                    ],
                    [
                        Imat21,
                        basisRotationMatrix12.T.dot(
                            Imat11.dot(basisRotationMatrix12)
                        ),
                    ],
                ]
            )
            return Imat

    def getC6perturbativelyAngularChannel(
        self,
        theta,
        phi,
        nRange,
        energyDelta,
        degeneratePerturbation=False,
        returnInteractionMatrix=False,
    ):
        r"""
        Calculates :math:`C_6` from second order perturbation theory.

        Calculates
        :math:`C_6=\sum_{\rm r',r''}|\langle {\rm r',r''}|V|\
        {\rm r1,r2}\rangle|^2/\Delta_{\rm r',r''}`, where
        :math:`\Delta_{\rm r',r''}\equiv E({\rm r',r''})-E({\rm r1, r2})`
        When second order perturbation couples to multiple energy degenerate
        states, users shold use **degenerate perturbation calculations** by
        setting `degeneratePerturbation=True` .

        This calculation is faster then full diagonalization, but it is valid
        only far from the so called spaghetti region that occurs when atoms
        are close to each other. In that region multiple levels are strongly
        coupled, and one needs to use full diagonalization. In region where
        perturbative calculation is correct, energy level shift can be
        obtained as :math:`V(R)=-C_6/R^6`

        Args:
            theta (float) - azimuthal angular orientation of atomic pair
                            state in rad
            phi (float) -   polar angular orientation of atomic pair state
                            in rad
            nRange (int) -  how much below and above the given principal quantum
                            number of the pair state we should be looking
            energyDelta (float) -   what is maximum energy difference
                                    ( :math:`\Delta E/h` in Hz)
                                    between the original pair state and the other
                                    pair states that we are including in
                                    calculation
            degeneratePerturbation (bool) - optional, default False. Should one
                                            use degenerate perturbation theory. 
                                            This should be used whenever
                                            angle between quantisation and
                                            interatomic axis is non-zero,
                                            as well as when one considers
                                            non-stretched states.
            returnInteractionMatrix (bool) -    optional, default False.
                                                Option to return the interaction
                                                matrix V(r)*R^6 in [GHz]
        Output:
            C6 (float) -    C6 value in [GHz] for the [n1,l1,j1,mj1; n2,l2,j2,mj2]
                            state specified in the PairStateInteraction class
                            initialisation
            if degeneratePerturbation == False:
                C6hop (float) - C6 value in [GHz] for the
                                [n1,l1,j1,mj1; n2,l2,j2,mj2] -> 
                                [n2,l2,j2,mj2; n1,l1,j1,mj1] state hopping contribution
            elif degeneratePerturbation == True:
                C6 (ndarray) -  array of eigenvalues of the full
                                interaction matrix in [GHz]
                eigenvectors (ndarray) -    corresponding list of eigenvectors
                                            :math:`\{m_{j_1}=-j_1, \ldots,
                                            m_{j_1} = +j1\}\bigotimes \
                                            \{ m_{j_2}=-j_2, \ldots, m_{j_2} = +j2\}`
                                            basis
            if returnInteractionMatrix == True:
                Imat_rot (ndarray) -    interaction matrix, fine-structure basis resolved
                                        for atomState1 == atomState2:
                                        [n1,l1,j1, mj1; n2,l2,j2,mj2] with 
                                        :math:`\{m_{j_1}=-j_1, \ldots,
                                        m_{j_1} = +j1\}\bigotimes \
                                        \{ m_{j_2}=-j_2, \ldots, m_{j_2} = +j2\}`
                                        for aomState != atomState2:
                                        first basis above, then basis with the 
                                        atomStates interchanged.
        """
        UsedModulesARC.pairstate_angular_channels = True

        atomState1 = [self.n, self.l, self.j, self.s1]
        atomState2 = [self.nn, self.ll, self.jj, self.s2]

        if degeneratePerturbation:
            degenerateStates = [
                [self.n, self.l, self.j, self.nn, self.ll, self.jj],
                [self.nn, self.ll, self.jj, self.n, self.l, self.j],
            ]

        # calculate interaction matrix wthout any basis changes (top left, Imat11)
        interaction11 = self._getC6contributions_lj(
            nRange, energyDelta, stateHopping=False
        )
        Imat11 = self._getPerturbativeC6Matrix_lj(interaction11)

        if atomState1 != atomState2:
            # if pair states are not identical, also calculate bottom left Imat (Imat21)
            interaction21 = self._getC6contributions_lj(
                nRange, energyDelta, stateHopping=True
            )
            Imat21 = self._getPerturbativeC6Matrix_lj(interaction21)
            # get rotation matrix for basis changes from [atomState2,atomState1] --> [atomState1,atomState2]
            basisRotationMatrix12 = self._getAngularBasisRotationMatrix(
                self.jj, self.j
            )

        # wigner D matrix allows calculations with arbitrary orientation of
        # the two atoms
        wgd = WignerDmatrix(theta, phi)
        angRotationMatrix = np.kron(
            wgd.get(atomState1[2]).toarray(), wgd.get(atomState2[2]).toarray()
        )
        if atomState1 != atomState2:
            angRotationMatrix2 = np.kron(
                wgd.get(atomState2[2]).toarray(),
                wgd.get(atomState1[2]).toarray(),
            )

        # rotate Imat's into correct basis for angles theta, phi (angle1, angle2)
        Imat11_rot = angRotationMatrix.dot(
            Imat11.dot(angRotationMatrix.conj().T)
        )
        if atomState1 != atomState2:
            Imat21_rot = angRotationMatrix2.dot(
                Imat21.dot(angRotationMatrix.conj().T)
            )

        # if degeneratePerturbation == False, calculate C6 value for a given mj1,mj2 state
        if not degeneratePerturbation:
            # calculate C6 value for non-hopped case
            compositeState1 = compositeState(
                singleAtomState(self.j, self.m1),
                singleAtomState(self.jj, self.m2),
            ).T
            C6 = np.real(
                compositeState1.dot(Imat11_rot.dot(compositeState1.T))[0][0]
            )

            # if atom states are different, also calculate C6 contribution from hopping
            if atomState1 != atomState2:
                compositeState2 = compositeState(
                    singleAtomState(self.jj, self.m2),
                    singleAtomState(self.j, self.m1),
                ).T
                C6hop = np.real(
                    compositeState2.dot(Imat21_rot.dot(compositeState1.T))[0][0]
                )
            # if atom states are the same, then C6hop = C6
            else:
                C6hop = C6

        # if degeneratePerturbation == True, construct full interaction matrix from above two matrices
        if degeneratePerturbation or returnInteractionMatrix:
            # construct full Imat
            if atomState1 == atomState2:
                Imat_rot = Imat11_rot
            elif atomState1 != atomState2:
                # compose resulting interaction marix
                Imat_rot = np.block(
                    [
                        [
                            Imat11_rot,
                            basisRotationMatrix12.dot(
                                Imat21_rot.dot(basisRotationMatrix12)
                            ),
                        ],
                        [
                            Imat21_rot,
                            basisRotationMatrix12.T.dot(
                                Imat11_rot.dot(basisRotationMatrix12)
                            ),
                        ],
                    ]
                )

            if degeneratePerturbation:
                # calculate eigenvalues, eigenstates etc
                eigenvalues, eigenvectors = np.linalg.eigh(Imat_rot)
                eigenvectors = eigenvectors.T

        # return function output
        if (not degeneratePerturbation) and (not returnInteractionMatrix):
            return C6, C6hop
        elif (not degeneratePerturbation) and (returnInteractionMatrix):
            return C6, C6hop, Imat_rot
        elif (degeneratePerturbation) and (not returnInteractionMatrix):
            return eigenvalues, eigenvectors, degenerateStates
        elif degeneratePerturbation and returnInteractionMatrix:
            return eigenvalues, eigenvectors, degenerateStates, Imat_rot

    def getC6perturbatively_anglePairs(
        self,
        anglePairs,
        nRange,
        energyDelta,
        degeneratePerturbation=False,
        returnInteractionMatrix=False,
    ):
        r"""
        Calculates :math:`C_6` from second order perturbation theory.


        Calculates
        :math:`C_6=\sum_{\rm r',r''}|\langle {\rm r',r''}|V|\
        {\rm r1,r2}\rangle|^2/\Delta_{\rm r',r''}`, where
        :math:`\Delta_{\rm r',r''}\equiv E({\rm r',r''})-E({\rm r1, r2})`
        When second order perturbation couples to multiple energy degenerate
        states, users shold use **degenerate perturbation calculations** by
        setting `degeneratePerturbation=True` .

        This calculation is faster then full diagonalization, but it is valid
        only far from the so called spaghetti region that occurs when atoms
        are close to each other. In that region multiple levels are strongly
        coupled, and one needs to use full diagonalization. In region where
        perturbative calculation is correct, energy level shift can be
        obtained as :math:`V(R)=-C_6/R^6`

        Args:
            anglePairs (list/array) -   contains lists/arrays of pairs of (theta, phi),
                                        i.e. azimuthal and polar orientations of
                                        atomic pair state in rad

            nRange (int) -  how much below and above the given principal quantum
                            number of the pair state we should be looking

            energyDelta (float) -   what is maximum energy difference
                                    ( :math:`\Delta E/h` in Hz)
                                    between the original pair state and the other
                                    pair states that we are including in
                                    calculation

            degeneratePerturbation (bool) - optional, default False. Should one
                                            use degenerate perturbation theory. 
                                            This should be used whenever
                                            angle between quantisation and
                                            interatomic axis is non-zero,
                                            as well as when one considers
                                            non-stretched states.

            returnInteractionMatrix (bool) -    optional, default False.
                                                Option to return the interaction
                                                matrix V(r)*R^6 in [GHz]

        Output:
            C6 (list) - list of arrays of C6 values in [GHz] for the [n1,l1,j1,mj1; n2,l2,j2,mj2]
                        state specified in the PairStateInteraction class initialisation

            if degeneratePerturbation == False:
                C6hop (list) -  list of arrays containing C6 value in [GHz] for the
                                [n1,l1,j1,mj1; n2,l2,j2,mj2] -> 
                                [n2,l2,j2,mj2; n1,l1,j1,mj1] state hopping contribution

            elif degeneratePerturbation == True:

                C6 (list) - list of arrays containing eigenvalues of the full
                            interaction matrix in [GHz]

                eigenvectors (list) -   list of arrays containing the corresponding list of eigenvectors
                                        :math:`\{m_{j_1}=-j_1, \ldots, m_{j_1} = +j1\}\bigotimes \
                                        \{ m_{j_2}=-j_2, \ldots, m_{j_2} = +j2\}` basis


            if returnInteractionMatrix == True:

                Imat_rot (list) -   list of arrays containing interaction matrices, fine-structure
                                    basis resolved for atomState1 == atomState2:
                                    [n1,l1,j1, mj1; n2,l2,j2,mj2] with 
                                    :math:`\{m_{j_1}=-j_1, \ldots, m_{j_1} = +j1\}\bigotimes \
                                    \{ m_{j_2}=-j_2, \ldots, m_{j_2} = +j2\}`
                                    for aomState != atomState2:
                                    first basis above, then basis with the atomStates interchanged                          
        """
        UsedModulesARC.pairstate_angular_channels = True

        atomState1 = [self.n, self.l, self.j, self.s1]
        atomState2 = [self.nn, self.ll, self.jj, self.s2]

        # save outputVals

        C6 = []

        if degeneratePerturbation:
            eigenvectors = []
            degenerateStates = [
                [self.n, self.l, self.j, self.nn, self.ll, self.jj],
                [self.nn, self.ll, self.jj, self.n, self.l, self.j],
            ]

        else:  # degeneratePerturbation == False
            C6hop = []

        if returnInteractionMatrix:
            interactionMatrices = []

        # calculate interaction matrix wthout any basis changes (top left, Imat11)

        interaction11 = self._getC6contributions_lj(
            nRange, energyDelta, stateHopping=False
        )

        Imat11 = self._getPerturbativeC6Matrix_lj(interaction11)

        if atomState1 != atomState2:
            # if pair states are not identical, also calculate bottom left Imat (Imat21)
            interaction21 = self._getC6contributions_lj(
                nRange, energyDelta, stateHopping=True
            )
            Imat21 = self._getPerturbativeC6Matrix_lj(interaction21)

            # get rotation matrix for basis changes from [atomState2,atomState1] --> [atomState1,atomState2]
            basisRotationMatrix12 = self._getAngularBasisRotationMatrix(
                self.jj, self.j
            )

        for i in range(len(anglePairs)):
            angle1, angle2 = anglePairs[i][0], anglePairs[i][1]

            # wigner D matrix allows calculations with arbitrary orientation of
            # the two atoms
            wgd = WignerDmatrix(angle1, angle2)
            angRotationMatrix = np.kron(
                wgd.get(atomState1[2]).toarray(),
                wgd.get(atomState2[2]).toarray(),
            )

            if atomState1 != atomState2:
                angRotationMatrix2 = np.kron(
                    wgd.get(atomState2[2]).toarray(),
                    wgd.get(atomState1[2]).toarray(),
                )

            # rotate Imat's into correct basis for angles theta, phi (angle1, angle2)
            Imat11_rot = angRotationMatrix.dot(
                Imat11.dot(angRotationMatrix.conj().T)
            )

            if atomState1 != atomState2:
                Imat21_rot = angRotationMatrix2.dot(
                    Imat21.dot(angRotationMatrix.conj().T)
                )

            # if degeneratePerturbation == False, calculate C6 value for a given mj1,mj2 state
            if not degeneratePerturbation:
                # calculate C6 value for non-hopped case
                compositeState1 = compositeState(
                    singleAtomState(self.j, self.m1),
                    singleAtomState(self.jj, self.m2),
                ).T
                C6val = np.real(
                    compositeState1.dot(Imat11_rot.dot(compositeState1.T))
                )[0][0]
                C6.append(C6val)

                # if atom states are different, also calculate C6 contribution from hopping
                if atomState1 != atomState2:
                    compositeState2 = compositeState(
                        singleAtomState(self.jj, self.m2),
                        singleAtomState(self.j, self.m1),
                    ).T
                    C6hop_val = np.real(
                        compositeState2.dot(Imat21_rot.dot(compositeState1.T))
                    )[0][0]
                    C6hop.append(C6hop_val)

            # if degeneratePerturbation == True, construct full interaction matrix from above two matrices
            if degeneratePerturbation or returnInteractionMatrix:
                # construct full Imat
                if atomState1 == atomState2:
                    Imat_rot = Imat11_rot
                elif atomState1 != atomState2:
                    # compose resulting interaction marix
                    Imat_rot = np.block(
                        [
                            [
                                Imat11_rot,
                                basisRotationMatrix12.dot(
                                    Imat21_rot.dot(basisRotationMatrix12)
                                ),
                            ],
                            [
                                Imat21_rot,
                                basisRotationMatrix12.T.dot(
                                    Imat11_rot.dot(basisRotationMatrix12)
                                ),
                            ],
                        ]
                    )

                if degeneratePerturbation:
                    # calculate eigenvalues, eigenstates etc
                    eigenvalues, vectors = np.linalg.eigh(Imat_rot)
                    vectors = vectors.T

                    # append output value arrays by current loop result
                    C6.append(np.array(eigenvalues))
                    eigenvectors.append(np.array(vectors))

                if returnInteractionMatrix:
                    interactionMatrices.append(np.array(Imat_rot))

        # return function output
        if (not degeneratePerturbation) and (not returnInteractionMatrix):
            return C6, C6hop

        elif (not degeneratePerturbation) and (returnInteractionMatrix):
            return C6, C6hop, interactionMatrices

        elif (degeneratePerturbation) and (not returnInteractionMatrix):
            return C6, eigenvectors, degenerateStates

        elif degeneratePerturbation and returnInteractionMatrix:
            return C6, eigenvectors, interactionMatrices, degenerateStates

    def _calcLJcontribution_allParamsFree(
        self,
        pathway,
        atom1,
        atom2,
        nRange,
        energyDelta,
        stateHopping,
        interactionsUpTo=1,
    ):
        r"""
        Returns the interaction strengths for the different
        (l,j; ll,jj) --> (l1,j1; l2,j2) --> (l',j'; ll',jj') configurations.

        Args:
            pathway (list) -    list containing the lj coupling pathway
                                [l,j, ll,jj, l1,j1, l2,j2, l',j' ll',jj']
            atom1 (list) -  infos on init state of atom 1 [n,s, atomType (ARC, e.g. Rubidium())]
            atom2 (list) -  infos on init state of atom 2 [n,s, atomType]
            nRange (int) -  how much below and above the given principal quantum number
                            of the pair state we should be looking
            energyDelta (float) -   what is maximum energy difference ( :math:`\Delta E/h` in Hz)
                                    between the original pair state and the other pair states that
                                    we are including in the calculation
            stateHopping (bool) -   whether or not the final state is interchanged ('hopped')
                                    w.r.t. the initial state
        Output:
            ljInteractions (list) - list containing entries of the form
                                    [(l,j, ll,jj, l1,j1, l2,j2, l',j' ll',jj'), V_{lj}]
                                    with V_{lj} the interaction strength for the given
                                    configuration in GHz(um)^6
        """

        V_lj = 0
        # unpack angular momentum info
        [l1, j1, ll1, jj1, l2, j2, ll2, jj2, l3, j3, ll3, jj3] = list(pathway)
        # iterate through n1 states
        for n2 in range(max(atom1[0] - nRange, 1), atom1[0] + nRange + 1):
            # iterate through n2 states
            for nn2 in range(max(atom2[0] - nRange, 1), atom2[0] + nRange + 1):
                # to check if nVals are okay
                nCheck = (
                    n2 >= atom1[2].groundStateN
                    or [n2, l2, j2] in atom1[2].extraLevels
                ) and (
                    nn2 >= atom2[2].groundStateN
                    or [nn2, ll2, jj2] in atom2[2].extraLevels
                )
                if stateHopping:
                    nCheck = (
                        nCheck
                        and (
                            atom2[0] >= atom1[2].groundStateN
                            or [atom2[0], ll3, jj3] in atom1[2].extraLevels
                        )
                        and (
                            atom1[0] >= atom2[2].groundStateN
                            or [atom1[0], l3, j3] in atom2[2].extraLevels
                        )
                    )

                # calculate energy defect
                energyDefect = (
                    self.__getEnergyDefect(
                        atom1[0],
                        l1,
                        j1,
                        atom2[0],
                        ll1,
                        jj1,
                        n2,
                        l2,
                        j2,
                        nn2,
                        ll2,
                        jj2,
                    )
                    / C_h
                )
                energyDefect = energyDefect * 1e-9  # GHz
                if abs(energyDefect) < 1e-10:
                    print(n2, l2, j2, nn2, ll2, jj2, stateHopping, "error")
                    raise ValueError(
                        "The requested pair-state "
                        "is dipole coupled resonatly "
                        "(energy defect = 0) "
                        "to other pair-states. "
                        "Aborting pertubative "
                        "calculation. "
                        "(This usually happens for "
                        "high-L states for which "
                        "identical quantum defects give "
                        "raise to degeneracies, making "
                        "total L ultimately not "
                        "conserved quantum number) "
                    )

                # proceed only if energy defect is within limit and nCheck was positive
                if (abs(energyDefect) < energyDelta * 10**-9) and nCheck:
                    # calculate radial overlaps
                    couplingStrength1 = _atomLightAtomCoupling(
                        atom1[0],
                        l1,
                        j1,
                        atom2[0],
                        ll1,
                        jj1,
                        n2,
                        l2,
                        j2,
                        nn2,
                        ll2,
                        jj2,
                        atom1[2],
                        atom2=atom2[2],
                        s=atom1[1],
                        s2=atom2[1],
                    ) * (1.0e-9 * (1.0e6) ** 3 / C_h)  # GHz / mum^3
                    if not stateHopping:
                        couplingStrength2 = couplingStrength1
                    else:
                        couplingStrength2 = _atomLightAtomCoupling(
                            n2,
                            l2,
                            j2,
                            nn2,
                            ll2,
                            jj2,
                            atom2[0],
                            l3,
                            j3,
                            atom1[0],
                            ll3,
                            jj3,
                            atom2[2],
                            atom2=atom1[2],
                            s=atom2[1],
                            s2=atom1[1],
                        ) * (1.0e-9 * (1.0e6) ** 3 / C_h)  # GHz / mum^3

                    V_lj += (
                        abs(couplingStrength1 * couplingStrength2)
                        / energyDefect
                    )  # GHz um^6
        return V_lj

    def __isAngularChannelDataCache(self):
        """
        Checks if the angular channel precalc data file exists or not.
        If not, creates the directory (if that doesn't exist yet) and creates an empty hdf5 file.

        Output:
            arcpath (str)     - local ARC directory path
            datapath (str)    - local ARC angular channel datacache directory path
            fileExist (bool)  - does datafile exist or not?
        """
        # data filename
        filename = "angularChannel_precalcData.hdf5"

        # get path to local ARC dir & data cache
        try:
            arcpath = inspectgetmodule(PairStateInteractions).__file__.strip(
                "calculations_atom_pairstate.py"
            )
            cachepath = (
                arcpath + "data" + arcpath[-1] + "C6_angularChannels_cache"
            )
            if not os.path.exists(cachepath):
                # create cache directory if it doesn't exist yet
                os.makedirs(cachepath)

            # check if file exists in local ARC angular channel data cache: local arc/data/C6_angularChannels_cache
            fileExist = os.path.isfile(cachepath + arcpath[-1] + filename)
        except Exception as e:
            raise ValueError(
                "Local ARC directory cannot be determined. " + str(e)
            )

        return arcpath, cachepath, fileExist

    def __loadAngularChannelPrecalcDataFromZenodo(self):
        """
        Loads precalculated datafile from Zenodo and stores locally.
        """
        # get directories
        arcpath, cachepath, _ = self.__isAngularChannelDataCache()

        # try to load data from Zenodo
        try:
            # load from Zenodo
            urllib.request.urlretrieve(
                "https://zenodo.org/record/15006915/files/angularChannel_precalcData.hdf5",
                cachepath + arcpath[-1] + "angularChannel_precalcData.hdf5",
            )
            print("Loaded data from Zenodo.")
        except Exception as _:
            # create empty file
            f = h5py.File(
                cachepath + arcpath[-1] + "angularChannel_precalcData.hdf5", "w"
            )
            # close
            f.close()

    def _checkLocalPrecalcForPairStateData(
        self, f, atom1Vals, atom2Vals, stateHopping
    ):
        """
        Checks whether or not a precalculated dataset for the specified pair state exists locally.

        Args:
            f (hdf5 object)  - data file handle
            atom1Vals (list) - [l1, j1, s1, atom1Type (ARC, e.g. Rubidium())]
            atom2Vals (list) - [l2, j2, s2, atom2Type (ARC, e.g. Rubidium())]
            stateHopping (bool) -   whether or not the final state is interchanged ('hopped')
                                    w.r.t. the initial state

        Output:
            status (bool) - True or False
            filekey (str) - key to dataset
        """
        # set status flag to False, filekey to nan
        status = False
        filekey = np.nan

        # unpack atom data
        [l1, j1, s1, atom1] = atom1Vals
        [l2, j2, s2, atom2] = atom2Vals

        # get atom pair from atom1Vals, atom2Vals and atom species interchanged
        atompairkey = atom1.elementName[:2] + atom2.elementName[:2]

        ## get folder
        datakeys = f[atompairkey].keys()

        # iterate through datasets in file until (hopefully) found the matching one
        for key2 in datakeys:
            # check that pair state details are the same
            if (
                not status
                and f[atompairkey + "/" + key2].attrs["atom 1"][:2]
                == atom1.elementName[:2]
                and f[atompairkey + "/" + key2].attrs["atom 2"][:2]
                == atom2.elementName[:2]
                and f[atompairkey + "/" + key2].attrs["j1"] == j1
                and f[atompairkey + "/" + key2].attrs["j2"] == j2
                and f[atompairkey + "/" + key2].attrs["l1"] == l1
                and f[atompairkey + "/" + key2].attrs["l2"] == l2
                and f[atompairkey + "/" + key2].attrs["stateHopping"]
                == str(stateHopping)
            ):
                status = True
                filekey = atompairkey + "/" + key2
            elif (
                not status
                and f[atompairkey + "/" + key2].attrs["atom 1"][:2]
                == atom2.elementName[:2]
                and f[atompairkey + "/" + key2].attrs["atom 2"][:2]
                == atom1.elementName[:2]
                and f[atompairkey + "/" + key2].attrs["j1"] == j2
                and f[atompairkey + "/" + key2].attrs["j2"] == j1
                and f[atompairkey + "/" + key2].attrs["l1"] == l2
                and f[atompairkey + "/" + key2].attrs["l2"] == l1
                and f[atompairkey + "/" + key2].attrs["stateHopping"]
                == str(stateHopping)
            ):
                print(
                    "Data exists, but for the properties of atom 1 and atom 2 interchanged. Swap input order to access data."
                )

        return status, filekey

    def __getAngularChannelFilename(self, atom1, l, j, atom2, ll, jj, stateHop):
        """
        Returns the filename for the specified atom parameters.

        Naming scheme example:
            angularChannels_CsS0p5_RbD2p5_stateHopFalse.txt
        Means:
            atom 1: Cs, S-state (i.e. L=0), j=0.5
            atom 2: Rb, D-state (i.e. L=2), j=2.5

        Args:
            atom1    - (ARC atom type) e.g. Rubidium()
            atom2    - (ARC atom type) e.g. Rubidium()
            l        - (int) l-value of first atom
            j        - (float) j-value of first atom
            ll       - (int) l-value of second atom
            jj       - (float) j-value of second atom
            stateHop - (bool) state hopping true or false?

        Output:
            filename - (str) filename for given atom data
        """
        # define l and j dictionaries
        lDict = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G", 5: "H"}
        jDict = {
            0.5: "0p5",
            1.5: "1p5",
            2.5: "2p5",
            3.5: "3p5",
            4.5: "4p5",
            5.5: "5p5",
            6.5: "6p5",
        }

        # get name block for atoms 1 and 2
        atom1block = atom1.elementName[:2] + lDict[l] + jDict[j]
        atom2block = atom2.elementName[:2] + lDict[ll] + jDict[jj]

        # get filename
        filename = (
            "angularChannels_"
            + atom1block
            + "_"
            + atom2block
            + "_stateHop"
            + str(stateHop)
            + ".txt"
        )

        return filename

    def __progressBar(self, n, ntot):
        """
        Prints a progress bar to std output.

        Args:
            n (int)    - current n
            ntot (int) - max n
        """
        # write progress bar output
        sys.stdout.write("\r")
        sys.stdout.write(
            "[%-20s] %d%%"
            % ("=" * int(np.floor(n / ntot * 20)), n / ntot * 100)
        )
        sys.stdout.flush()

    def loadAngularChannelData(
        self, atom1Vals, atom2Vals, stateHopping=False, loadFromZenodo=False
    ):
        """
        Checks if the precalculated angular channel values exist in the local cache.

        If yes, returns the data.

        If not, checks on Zenodo () if bulk data exists,
            - if yes loads datafile and saves it locally, then returns the data.
            - if not, output prompt to user explainaing how to calculate the requested dataset with the :obj:`saveAngularChannelData` function.

        Args:
            atom1Vals (list) - [l1, j1, s1, atom1Type (ARC, e.g. Rubidium())]
            atom2Vals (list) - [l2, j2, s2, atom2Type (ARC, e.g. Rubidium())]
            stateHopping (bool) -   whether or not the final state is interchanged ('hopped')
                                    w.r.t. the initial state
            loadFromZenodo (bool) - whether or not to load data from Zenodo
                                        If True then this will overwrite locally existing data for the specified pair state
                                        with data from Zenodo if exists
                                        If False but no local data exists, then checks Zenodo if data exists there and loads if True

        Output:
            metadata (dict)       - file calculation infos (atom infos, calculation settings)
            coupledChannels (list) - list of coupled channels, in same order as in data
            data (ndarray)         - list of angular channel values of the form: [n1, n2, C_lj1, C_lj2, ...] with C_lj the channel values
        """
        UsedModulesARC.pairstate_angular_channels = True

        # unpack atom data
        [l1, j1, s1, atom1] = atom1Vals
        [l2, j2, s2, atom2] = atom2Vals

        # initialise data containers
        metadata, coupledChannels, data = np.nan, np.nan, np.nan
        dataExists = False

        ## DOES FILE EXIST LOCALLY?

        # get ARC and precalc data directories
        arcpath, cachepath, fileExist = self.__isAngularChannelDataCache()

        # if file doesn't exist or loadFromZenodo=True: load file from Zenodo
        if not fileExist or loadFromZenodo:
            self.__loadAngularChannelPrecalcDataFromZenodo()

        ## DOES REQUESTED ATOM DATA EXIST LOCALLY?

        # check if we can find relevant data
        # open file
        f = h5py.File(
            cachepath + arcpath[-1] + "angularChannel_precalcData.hdf5", "r"
        )
        # get groups (i.e. equivalent to folders)
        atompairkeys = f.keys()

        # check if atom1-atom2 combination data exists
        if atom1.elementName[:2] + atom2.elementName[:2] in atompairkeys:
            # check if the specific pair state data exists
            dataExists, datakey = self._checkLocalPrecalcForPairStateData(
                f, atom1Vals, atom2Vals, stateHopping
            )
            if dataExists:
                print("Data exists.")
        elif atom2.elementName[:2] + atom1.elementName[:2] in atompairkeys:
            # check if the specific pair state data exists
            dataExists, datakey = self._checkLocalPrecalcForPairStateData(
                f, atom2Vals, atom1Vals, stateHopping
            )
            if dataExists:
                print(
                    "Precalculated data exists, but for atomic species interchanged.\nLoaded data for species in order: atom 1: "
                    + atom2.elementName[:2]
                    + ", atom 2: "
                    + atom1.elementName[:2]
                )

        # if the data exists: load all relevant data
        if dataExists:
            # get file attributes, i.e. calculation metadata
            fileattrs = [x for x in f[datakey].attrs.keys()][1:]
            metadata = dict([(x, f[datakey].attrs[x]) for x in fileattrs])
            # get angular channels
            coupledChannels = [
                [x] for x in (f[datakey].dims[1].label)[14:-2].split("], [")
            ]
            for i, channel in enumerate(coupledChannels):
                coupledChannels[i] = [float(x) for x in channel[0].split(",")]
            # load data
            data = f[
                datakey
            ][
                :
            ]  # 'empty' slicing is required to get ndarray that remains accessible after hdf5 file closure
        # else: prompt user to calculate data with calcAngularChannelData function
        else:
            print(
                "The requested data does not exist in the local cache. Please calculate the data via the calcAngularChannelData function."
            )

        # close file again
        f.close()

        return metadata, coupledChannels, data

    def calculateAngularChannelData(
        self,
        atom1Vals,
        atom2Vals,
        nValueRange,
        nRange,
        energyDelta,
        stateHopping=False,
        overwriteLocalData=False,
    ):
        """
        Saves the angular channel values C_{lj} for the atom1 and atom2 pair-interaction. Data is stored in local cache.
        With this data, the full interaction matrix can be reconstructed by e.g. passing the angular channel values to the function :obj:`_getPerturbativeC6Matrix_lj`,
        and the angular channel values can be loaded with the function :obj:`loadAngularChannelData`.
        For more information on how to implement this, check the example Jupyter notebook on the angular channel code.

        Args:
            atom1Vals (list) - [l1, j1, s1, atom1Type (ARC, e.g. Rubidium())]
            atom2Vals (list) - [l2, j2, s2, atom2Type (ARC, e.g. Rubidium())]
            nValueRange (list) - [nMin, nMax]
            nRange (int) -  how much below and above the given principal quantum number
                            of the pair state we should be looking
            energyDelta (float) -   what is maximum energy difference ( :math:`Delta E/h` in Hz)
                                    between the original pair state and the other pair states that
                                    we are including in the calculation
            stateHopping (bool) -   whether or not the final state is interchanged ('hopped')
                                    w.r.t. the initial state
            overwriteLocalData (bool) - allow to overwrite the local dataset if it exists already?


        Output:
            status (bool) - status flag, True if calculation exited successfully
        """
        UsedModulesARC.pairstate_angular_channels = True

        # unpack atom data
        [l1, j1, s1, atom1] = atom1Vals
        [l2, j2, s2, atom2] = atom2Vals

        # get ARC and precalc data directories
        arcpath, cachepath, fileExist = self.__isAngularChannelDataCache()

        # data exists boolean
        dataExists = False

        ## DOES REQUESTED ATOM DATA EXIST LOCALLY?

        if not overwriteLocalData:
            # check if we can find relevant data
            # open file
            f = h5py.File(
                cachepath + arcpath[-1] + "angularChannel_precalcData.hdf5", "r"
            )
            # get groups (i.e. equivalent to folders)
            atompairkeys = f.keys()

            # check if atom1-atom2 combination data exists
            if atom1.elementName[:2] + atom2.elementName[:2] in atompairkeys:
                # check if the specific pair state data exists
                dataExists, datakey = self._checkLocalPrecalcForPairStateData(
                    f, atom1Vals, atom2Vals, stateHopping
                )
            elif atom2.elementName[:2] + atom1.elementName[:2] in atompairkeys:
                # check if the specific pair state data exists
                dataExists, datakey = self._checkLocalPrecalcForPairStateData(
                    f, atom2Vals, atom1Vals, stateHopping
                )
                if dataExists:
                    print(
                        "Precalculated data exists, but for atomic species interchanged.\nTry to load data for atom 1: "
                        + atom2.elementName[:2]
                        + ", atom 2: "
                        + atom1.elementName[:2]
                    )

            # close file again
            f.close()

            # if data exists, then print warning and exit function
            if dataExists:
                raise Warning(
                    "The data exists locally and overwriteLocal was set to False. No new data was computed."
                    "\nYou can enforce recalculation of the data by setting the function parameter overwriteLocalData=True. /"
                    "\nYou can call the existing data via the function loadAngularChannelData."
                )

        ## DATA DOES NOT EXIST YET OR OVERWRITE LOCAL=TRUE ##

        # calculate data
        # get coupling pathways
        if not stateHopping:
            coupledStates = self._findAllCoupledAngularMomentumStates(
                l1, j1, s1, l2, j2, s2, stateHopping=False
            )
        else:  # stateHopping == True
            coupledStates = self._findAllCoupledAngularMomentumStates(
                l1, j1, s1, l2, j2, s2, stateHopping=True
            )
        if coupledStates == []:
            raise ValueError(
                "No interaction pathways found for the specified conditions."
            )

        # total number of values to be calculated
        ntot = int((nValueRange[1] - nValueRange[0] + 1) ** 2)
        ncurr = 1
        # get angular channel values
        ljValues = np.zeros(
            (
                (int(nValueRange[1] - nValueRange[0]) + 1) ** 2,
                2 + len(coupledStates),
            )
        )
        # iterate through n1Vals
        for n1 in range(nValueRange[0], nValueRange[1] + 1):
            i = int(n1 - nValueRange[0])
            # iterate through n2Vals
            for n2 in range(nValueRange[0], nValueRange[1] + 1):
                j = int(n2 - nValueRange[0])
                # calculate channel values
                vals = []
                for pathway in coupledStates:
                    V_lj = self._calcLJcontribution_allParamsFree(
                        pathway,
                        [n1, s1, atom1],
                        [n2, s2, atom2],
                        nRange,
                        energyDelta,
                        stateHopping,
                        interactionsUpTo=self.interactionsUpTo,
                    )
                    vals.append(V_lj)
                ljValues[
                    i * int(nValueRange[1] - nValueRange[0] + 1) + j, :
                ] = [n1, n2, *vals]
                # print progress bar
                self.__progressBar(ncurr, ntot)
                ncurr += 1

        ## write new data to file

        # open file
        f = h5py.File(
            cachepath + arcpath[-1] + "angularChannel_precalcData.hdf5", "r+"
        )

        # check if current atom combination is already key
        atomSpeciesKey = atom1.elementName[:2] + atom2.elementName[:2]
        if atomSpeciesKey in f.keys():
            datafolder = f[atomSpeciesKey]
        # if not: create group (folder) and add attributes
        else:
            datafolder = f.create_group(atomSpeciesKey)
            # add attributes to group element
            datafolder.attrs["atom 1"] = atom1.elementName[:2]
            datafolder.attrs["atom 2"] = atom2.elementName[:2]

        # create dataset
        filename = self.__getAngularChannelFilename(
            atom1, l1, j1, atom2, l2, j2, stateHopping
        )
        dset = datafolder.create_dataset(
            filename, data=ljValues, dtype="f", compression="gzip"
        )

        # add 'axis' label to dataset columns, i.e. coupled state info
        dset.dims[1].label = str(
            ["n1", "n2", *[list(x) for x in coupledStates]]
        )

        # add atom pair state attributes to dataset
        dset.attrs["atom 1"] = atom1.elementName
        dset.attrs["l1"] = l1
        dset.attrs["j1"] = j1
        dset.attrs["atom 2"] = atom2.elementName
        dset.attrs["l2"] = l2
        dset.attrs["j2"] = j2
        dset.attrs["stateHopping"] = str(stateHopping)
        # add calculation attributes to dataset
        dset.attrs["nRange"] = str(nRange)
        dset.attrs["energyDelta"] = str(energyDelta)
        dset.attrs["nValueRange"] = str(nValueRange)
        dset.attrs["interactionsUpTo"] = str(self.interactionsUpTo)

        # close file
        f.close()

        print("\nData was successfully saved to local cache.")

    def defineBasis(
        self,
        theta,
        phi,
        nRange,
        lrange,
        energyDelta,
        Bz=0,
        progressOutput=False,
        debugOutput=False,
    ):
        r"""
        Finds relevant states in the vicinity of the given pair-state

        Finds relevant pair-state basis and calculates interaction matrix.
        Pair-state basis is saved in :obj:`basisStates`.
        Interaction matrix is saved in parts depending on the scaling with
        distance. Diagonal elements :obj:`matDiagonal`, correponding to
        relative energy defects of the pair-states, don't change with
        interatomic separation. Off diagonal elements can depend
        on distance as :math:`R^{-3}, R^{-4}` or :math:`R^{-5}`,
        corresponding to dipole-dipole (:math:`C_3` ), dipole-qudrupole
        (:math:`C_4` ) and quadrupole-quadrupole coupling (:math:`C_5` )
        respectively. These parts of the matrix are stored in
        :obj:`PairStateInteractions.matR`
        in that order. I.e. :obj:`matR[0]`
        stores dipole-dipole coupling
        (:math:`\propto R^{-3}`),
        :obj:`matR[1]` stores dipole-quadrupole
        couplings etc.

        Parameters:
            theta (float):  relative orientation of the two atoms
                (see figure on top of the page), range 0 to :math:`\pi`
            phi (float): relative orientation of the two atoms (see figure
                on top of the page), range 0 to :math:`2\pi`
            nRange (int): how much below and above the given principal
                quantum number of the pair state we should be looking?
            lrange (int): what is the maximum angular orbital momentum
                state that we are including in calculation
            energyDelta (float): what is maximum energy difference (
                :math:`\Delta E/h` in Hz)
                between the original pair state and the other pair states
                that we are including in calculation
            Bz (float): optional, magnetic field directed along z-axis in
                units of Tesla. Calculation will be correct only for weak
                magnetic fields, where paramagnetic term is much stronger
                then diamagnetic term. Diamagnetic term is neglected.
            progressOutput (bool): optional, False by default. If true,
                prints information about the progress of the calculation.
            debugOutput (bool): optional, False by default. If true,
                similarly to progressOutput=True, this will print
                information about the progress of calculations, but with
                more verbose output.

        See also:
            :obj:`arc.alkali_atom_functions.saveCalculation` and
            :obj:`arc.alkali_atom_functions.loadSavedCalculation` for
            information on saving intermediate results of calculation for
            later use.
        """

        self.__initializeDatabaseForMemoization()

        # save call parameters
        self.theta = theta
        self.phi = phi
        self.nRange = nRange
        self.lrange = lrange
        self.energyDelta = energyDelta
        self.Bz = Bz

        self.basisStates = []

        # wignerDmatrix
        wgd = WignerDmatrix(theta, phi)

        limitBasisToMj = False
        if theta < 0.001:
            limitBasisToMj = True  # Mj will be conserved in calculations

        originalMj = self.m1 + self.m2

        self.channel, self.coupling = self.__makeRawMatrix2(
            self.n,
            self.l,
            self.j,
            self.nn,
            self.ll,
            self.jj,
            nRange,
            lrange,
            energyDelta,
            limitBasisToMj,
            progressOutput=progressOutput,
            debugOutput=debugOutput,
        )

        self.atom1.updateDipoleMatrixElementsFile()
        self.atom2.updateDipoleMatrixElementsFile()

        # generate all the states (with mj principal quantum number)

        # opi = original pairstate index
        opi = 0

        # NEW FOR SPACE MATRIX
        self.index = np.zeros(len(self.channel) + 1, dtype=np.int16)

        for i in xrange(len(self.channel)):
            self.index[i] = len(self.basisStates)

            stateCoupled = self.channel[i]

            for m1c in np.linspace(
                stateCoupled[2],
                -stateCoupled[2],
                round(1 + 2 * stateCoupled[2]),
            ):
                for m2c in np.linspace(
                    stateCoupled[5],
                    -stateCoupled[5],
                    round(1 + 2 * stateCoupled[5]),
                ):
                    if (not limitBasisToMj) or (
                        abs(originalMj - m1c - m2c) < 0.1
                    ):
                        self.basisStates.append(
                            [
                                stateCoupled[0],
                                stateCoupled[1],
                                stateCoupled[2],
                                m1c,
                                stateCoupled[3],
                                stateCoupled[4],
                                stateCoupled[5],
                                m2c,
                            ]
                        )
                        self.matrixElement.append(i)

                        if (
                            abs(stateCoupled[0] - self.n) < 0.1
                            and abs(stateCoupled[1] - self.l) < 0.1
                            and abs(stateCoupled[2] - self.j) < 0.1
                            and abs(m1c - self.m1) < 0.1
                            and abs(stateCoupled[3] - self.nn) < 0.1
                            and abs(stateCoupled[4] - self.ll) < 0.1
                            and abs(stateCoupled[5] - self.jj) < 0.1
                            and abs(m2c - self.m2) < 0.1
                        ):
                            opi = len(self.basisStates) - 1
            if self.index[i] == len(self.basisStates):
                print(stateCoupled)
        self.index[-1] = len(self.basisStates)

        if progressOutput or debugOutput:
            print("\nCalculating Hamiltonian matrix...\n")

        dimension = len(self.basisStates)
        if progressOutput or debugOutput:
            print("\n\tmatrix (dimension ", dimension, ")\n")

        # INITIALIZING MATICES
        # all (sparce) matrices will be saved in csr format
        # value, row, column
        matDiagonalConstructor = [[], [], []]

        matRConstructor = [
            [[], [], []] for i in xrange(self.interactionsUpTo * 2 - 1)
        ]

        matRIndex = 0
        for c in self.coupling:
            progress = 0.0
            for ii in xrange(len(self.channel)):
                if progressOutput:
                    dim = len(self.channel)
                    progress += (dim - ii) * 2 - 1
                    sys.stdout.write(
                        "\rMatrix R%d %.1f %% (state %d of %d)"
                        % (
                            matRIndex + 3,
                            float(progress) / float(dim**2) * 100.0,
                            ii + 1,
                            len(self.channel),
                        )
                    )
                    sys.stdout.flush()

                ed = self.channel[ii][6]

                # solves problems with exactly degenerate basisStates
                degeneracyOffset = 0.00000001

                i = self.index[ii]
                dMatrix1 = wgd.get(self.basisStates[i][2])
                dMatrix2 = wgd.get(self.basisStates[i][6])

                for i in xrange(self.index[ii], self.index[ii + 1]):
                    statePart1 = singleAtomState(
                        self.basisStates[i][2], self.basisStates[i][3]
                    )
                    statePart2 = singleAtomState(
                        self.basisStates[i][6], self.basisStates[i][7]
                    )
                    # rotate individual states

                    statePart1 = dMatrix1.T.conjugate().dot(statePart1)
                    statePart2 = dMatrix2.T.conjugate().dot(statePart2)

                    stateCom = compositeState(statePart1, statePart2)

                    if matRIndex == 0:
                        zeemanShift = (
                            (
                                self.atom1.getZeemanEnergyShift(
                                    self.basisStates[i][1],
                                    self.basisStates[i][2],
                                    self.basisStates[i][3],
                                    self.Bz,
                                    s=self.s1,
                                )
                                + self.atom2.getZeemanEnergyShift(
                                    self.basisStates[i][5],
                                    self.basisStates[i][6],
                                    self.basisStates[i][7],
                                    self.Bz,
                                    s=self.s2,
                                )
                            )
                            / C_h
                            * 1.0e-9
                        )  # in GHz
                        matDiagonalConstructor[0].append(
                            ed + zeemanShift + degeneracyOffset
                        )
                        degeneracyOffset += 0.00000001
                        matDiagonalConstructor[1].append(i)
                        matDiagonalConstructor[2].append(i)

                    for dataIndex in xrange(c.indptr[ii], c.indptr[ii + 1]):
                        jj = c.indices[dataIndex]
                        radialPart = c.data[dataIndex]

                        j = self.index[jj]
                        dMatrix3 = wgd.get(self.basisStates[j][2])
                        dMatrix4 = wgd.get(self.basisStates[j][6])

                        if self.index[jj] != self.index[jj + 1]:
                            d = self.__getAngularMatrix_M(
                                self.basisStates[i][1],
                                self.basisStates[i][2],
                                self.basisStates[i][5],
                                self.basisStates[i][6],
                                self.basisStates[j][1],
                                self.basisStates[j][2],
                                self.basisStates[j][5],
                                self.basisStates[j][6],
                            )
                            secondPart = d.dot(stateCom)
                        else:
                            print(" - - - ", self.channel[jj])

                        for j in xrange(self.index[jj], self.index[jj + 1]):
                            statePart1 = singleAtomState(
                                self.basisStates[j][2], self.basisStates[j][3]
                            )
                            statePart2 = singleAtomState(
                                self.basisStates[j][6], self.basisStates[j][7]
                            )
                            # rotate individual states

                            statePart1 = dMatrix3.T.conjugate().dot(statePart1)
                            statePart2 = dMatrix4.T.conjugate().dot(statePart2)
                            # composite state of two atoms
                            stateCom2 = compositeState(statePart1, statePart2)

                            angularFactor = stateCom2.T.conjugate().dot(
                                secondPart
                            )
                            if abs(self.phi) < 1e-9:
                                angularFactor = angularFactor[0, 0].real
                            else:
                                angularFactor = angularFactor[0, 0]

                            if abs(angularFactor) > 1.0e-5:
                                matRConstructor[matRIndex][0].append(
                                    (radialPart * angularFactor).conj()
                                )
                                matRConstructor[matRIndex][1].append(i)
                                matRConstructor[matRIndex][2].append(j)

                                matRConstructor[matRIndex][0].append(
                                    radialPart * angularFactor
                                )
                                matRConstructor[matRIndex][1].append(j)
                                matRConstructor[matRIndex][2].append(i)
            matRIndex += 1
            if progressOutput or debugOutput:
                print("\n")

        self.matDiagonal = csr_matrix(
            (
                matDiagonalConstructor[0],
                (matDiagonalConstructor[1], matDiagonalConstructor[2]),
            ),
            shape=(dimension, dimension),
        )

        self.matR = [
            csr_matrix(
                (
                    matRConstructor[i][0],
                    (matRConstructor[i][1], matRConstructor[i][2]),
                ),
                shape=(dimension, dimension),
            )
            for i in xrange(self.interactionsUpTo * 2 - 1)
        ]

        self.originalPairStateIndex = opi

        self.__updateAngularMatrixElementsFile()
        self.__closeDatabaseForMemoization()

    def diagonalise(
        self,
        rangeR,
        noOfEigenvectors,
        drivingFromState=[0, 0, 0, 0, 0],
        eigenstateDetuning=0.0,
        sortEigenvectors=False,
        progressOutput=False,
        debugOutput=False,
    ):
        r"""
        Finds eigenstates in atom pair basis.

        ARPACK ( :obj:`scipy.sparse.linalg.eigsh`) calculation of the
        `noOfEigenvectors` eigenvectors closest to the original state. If
        `drivingFromState` is specified as `[n,l,j,mj,q]` coupling between
        the pair-states and the situation where one of the atoms in the
        pair state basis is in :math:`|n,l,j,m_j\rangle` state due to
        driving with a laser field that drives :math:`q` transition
        (+1,0,-1 for :math:`\sigma^-`, :math:`\pi` and :math:`\sigma^+`
        transitions respectively) is calculated and marked by the
        colourmaping these values on the obtained eigenvectors.

        Parameters:
            rangeR ( :obj:`array`): Array of values for distance between
                the atoms (in :math:`\mu` m) for which we want to calculate
                eigenstates.
            noOfEigenvectors (int): number of eigen vectors closest to the
                energy of the original (unperturbed) pair state. Has to be
                smaller then the total number of states.
            eigenstateDetuning (float, optional): Default is 0. This
                specifies detuning from the initial pair-state (in Hz)
                around which we want to find `noOfEigenvectors`
                eigenvectors. This is useful when looking only for couple
                of off-resonant features.
            drivingFromState ([int,int,float,float,int]): Optional. State
                of one of the atoms from the original pair-state basis
                from which we try to drive to the excited pair-basis
                manifold, **assuming that the first of the two atoms is
                already excited to the specified Rydberg state**.
                By default, program will calculate just
                contribution of the original pair-state in the eigenstates
                obtained by diagonalization, and will highlight it's
                admixure by colour mapping the obtained eigenstates plot.
                State is specified as :math:`[n,\ell,j,mj, d]`
                where :math:`d` is +1, 0 or
                -1 for driving :math:`\sigma^-` , :math:`\pi`
                and :math:`\sigma^+` transitions respectively.
            sortEigenvectors(bool): optional, False by default. Tries to
                sort eigenvectors so that given eigen vector index
                corresponds to adiabatically changing eigenstate, as
                detirmined by maximising overlap between old and new
                eigenvectors.
            progressOutput (bool): optional, False by default. If true,
                prints information about the progress of the calculation.
            debugOutput (bool): optional, False by default. If true,
                similarly to progressOutput=True, this will print
                information about the progress of calculations, but with
                more verbose output.
        """

        self.r = np.sort(rangeR)
        dimension = len(self.basisStates)

        self.noOfEigenvectors = noOfEigenvectors

        # energy of the state - to be calculated
        self.y = []
        # how much original state is contained in this eigenvector
        self.highlight = []
        # what are the dominant contributing states?
        self.composition = []

        if noOfEigenvectors >= dimension - 1:
            noOfEigenvectors = dimension - 1
            print(
                "Warning: Requested number of eigenvectors >=dimension-1\n \
                 ARPACK can only find up to dimension-1 eigenvectors, where\
                dimension is matrix dimension.\n"
            )
            if noOfEigenvectors < 1:
                return

        coupling = []
        self.maxCoupling = 0.0
        self.maxCoupledStateIndex = 0
        if drivingFromState[0] != 0:
            self.drivingFromState = drivingFromState
            if progressOutput:
                print("Finding coupling strengths")
            # get first what was the state we are calculating coupling with
            state1 = drivingFromState
            n1 = round(state1[0])
            l1 = round(state1[1])
            j1 = state1[2]
            m1 = state1[3]
            q = state1[4]

            for i in xrange(dimension):
                thisCoupling = 0.0

                if (
                    round(abs(self.basisStates[i][5] - l1)) == 1
                    and abs(
                        self.basisStates[i][0]
                        - self.basisStates[self.originalPairStateIndex][0]
                    )
                    < 0.1
                    and abs(
                        self.basisStates[i][1]
                        - self.basisStates[self.originalPairStateIndex][1]
                    )
                    < 0.1
                    and abs(
                        self.basisStates[i][2]
                        - self.basisStates[self.originalPairStateIndex][2]
                    )
                    < 0.1
                    and abs(
                        self.basisStates[i][3]
                        - self.basisStates[self.originalPairStateIndex][3]
                    )
                    < 0.1
                ):
                    state2 = self.basisStates[i]
                    n2 = round(state2[0 + 4])
                    l2 = round(state2[1 + 4])
                    j2 = state2[2 + 4]
                    m2 = state2[3 + 4]
                    if debugOutput:
                        print(
                            n1,
                            " ",
                            l1,
                            " ",
                            j1,
                            " ",
                            m1,
                            " ",
                            n2,
                            " ",
                            l2,
                            " ",
                            j2,
                            " ",
                            m2,
                            " q=",
                            q,
                        )
                        print(self.basisStates[i])
                    dme = self.atom2.getDipoleMatrixElement(
                        n1, l1, j1, m1, n2, l2, j2, m2, q, s=self.s2
                    )
                    thisCoupling += dme

                thisCoupling = abs(thisCoupling) ** 2
                if thisCoupling > self.maxCoupling:
                    self.maxCoupling = thisCoupling
                    self.maxCoupledStateIndex = i
                if (thisCoupling > 0.000001) and debugOutput:
                    print(
                        "original pairstate index = ",
                        self.originalPairStateIndex,
                    )
                    print("this pairstate index = ", i)
                    print("state itself ", self.basisStates[i])
                    print("coupling = ", thisCoupling)
                coupling.append(thisCoupling)

            print("Maximal coupling from a state")
            print("is to a state ", self.basisStates[self.maxCoupledStateIndex])
            print("is equal to %.3e a_0 e" % self.maxCoupling)

        if progressOutput:
            print("\n\nDiagonalizing interaction matrix...\n")

        rvalIndex = 0.0
        previousEigenvectors = []

        for rval in self.r:
            if progressOutput:
                sys.stdout.write(
                    "\r%d%%" % (rvalIndex / len(self.r - 1) * 100.0)
                )
                sys.stdout.flush()
            rvalIndex += 1.0

            # calculate interaction matrix

            m = self.matDiagonal

            rX = (rval * 1.0e-6) ** 3
            for matRX in self.matR:
                m = m + matRX / rX
                rX *= rval * 1.0e-6

            # uses ARPACK algorithm to find only noOfEigenvectors eigenvectors
            # sigma specifies center frequency (in GHz)
            ev, egvector = eigsh(
                m,
                noOfEigenvectors,
                sigma=eigenstateDetuning * 1.0e-9,
                which="LM",
                tol=1e-6,
            )

            if sortEigenvectors:
                # Find which eigenvectors overlap most with eigenvectors from
                # previous diagonalisatoin, in order to find "adiabatic"
                # continuation for the respective states

                if previousEigenvectors == []:
                    previousEigenvectors = np.copy(egvector)
                rowPicked = [False for i in range(len(ev))]
                columnPicked = [False for i in range(len(ev))]

                stateOverlap = np.zeros((len(ev), len(ev)))
                for i in range(len(ev)):
                    for j in range(len(ev)):
                        stateOverlap[i, j] = (
                            np.vdot(egvector[:, i], previousEigenvectors[:, j])
                            ** 2
                        )

                sortedOverlap = np.dstack(
                    np.unravel_index(
                        np.argsort(stateOverlap.ravel()), (len(ev), len(ev))
                    )
                )[0]

                sortedEigenvaluesOrder = np.zeros(len(ev), dtype=np.int32)
                j = len(ev) ** 2 - 1
                for i in range(len(ev)):
                    while (
                        rowPicked[sortedOverlap[j, 0]]
                        or columnPicked[sortedOverlap[j, 1]]
                    ):
                        j -= 1
                    rowPicked[sortedOverlap[j, 0]] = True
                    columnPicked[sortedOverlap[j, 1]] = True
                    sortedEigenvaluesOrder[sortedOverlap[j, 1]] = sortedOverlap[
                        j, 0
                    ]

                egvector = egvector[:, sortedEigenvaluesOrder]
                ev = ev[sortedEigenvaluesOrder]
                previousEigenvectors = np.copy(egvector)

            self.y.append(ev)

            if drivingFromState[0] < 0.1:
                # if we've defined from which state we are driving
                sh = []
                comp = []
                for i in xrange(len(ev)):
                    sh.append(
                        abs(egvector[self.originalPairStateIndex, i]) ** 2
                    )
                    comp.append(self._stateComposition(egvector[:, i]))
                self.highlight.append(sh)
                self.composition.append(comp)
            else:
                sh = []
                comp = []
                for i in xrange(len(ev)):
                    sumCoupledStates = 0.0
                    for j in xrange(dimension):
                        sumCoupledStates += (
                            abs(coupling[j] / self.maxCoupling)
                            * abs(egvector[j, i]) ** 2
                        )
                    comp.append(self._stateComposition(egvector[:, i]))
                    sh.append(sumCoupledStates)
                self.highlight.append(sh)
                self.composition.append(comp)

        # end of FOR loop over inter-atomic dinstaces

    def exportData(self, fileBase, exportFormat="csv"):
        """
        Exports PairStateInteractions calculation data.

        Only supported format (selected by default) is .csv in a
        human-readable form with a header that saves details of calculation.
        Function saves three files: 1) `filebase` _r.csv;
        2) `filebase` _energyLevels
        3) `filebase` _highlight

        For more details on the format, see header of the saved files.

        Parameters:
            filebase (string): filebase for the names of the saved files
                without format extension. Add as a prefix a directory path
                if necessary (e.g. saving outside the current working directory)
            exportFormat (string): optional. Format of the exported file. Currently
                only .csv is supported but this can be extended in the future.
        """
        fmt = "on %Y-%m-%d @ %H:%M:%S"
        ts = datetime.datetime.now().strftime(fmt)

        commonHeader = "Export from Alkali Rydberg Calculator (ARC) %s.\n" % ts
        commonHeader += (
            "\n *** Pair State interactions for %s %s m_j = %d/2 , %s %s m_j = %d/2 pair-state. ***\n\n"
            % (
                self.atom1.elementName,
                printStateString(self.n, self.l, self.j),
                round(2.0 * self.m1),
                self.atom2.elementName,
                printStateString(self.nn, self.ll, self.jj),
                round(2.0 * self.m2),
            )
        )
        if self.interactionsUpTo == 1:
            commonHeader += " - Pair-state interactions included up to dipole-dipole coupling.\n"
        elif self.interactionsUpTo == 2:
            commonHeader += " - Pair-state interactions included up to quadrupole-quadrupole coupling.\n"
        commonHeader += (
            " - Pair-state interactions calculated for manifold with spin angular momentum s1 = %.1d s2 = %.1d .\n"
            % (self.s1, self.s2)
        )

        if hasattr(self, "theta"):
            commonHeader += " - Atom orientation:\n"
            commonHeader += "      theta (polar angle) = %.5f x pi\n" % (
                self.theta / pi
            )
            commonHeader += "      phi (azimuthal angle) = %.5f x pi\n" % (
                self.phi / pi
            )
            commonHeader += " - Calculation basis includes:\n"
            commonHeader += (
                "      States with principal quantum number in range [%d-%d]x[%d-%d],\n"
                % (
                    self.n - self.nRange,
                    self.n + self.nRange,
                    self.nn - self.nRange,
                    self.nn + self.nRange,
                )
            )
            commonHeader += (
                "      AND whose orbital angular momentum (l) is in range [%d-%d] (i.e. %s-%s),\n"
                % (
                    0,
                    self.lrange,
                    printStateLetter(0),
                    printStateLetter(self.lrange),
                )
            )
            commonHeader += (
                "      AND whose pair-state energy difference is at most %.3f GHz\n"
                % (self.energyDelta / 1.0e9)
            )
            commonHeader += "      (energy difference is measured relative to original pair-state).\n"
        else:
            commonHeader += " ! Atom orientation and basis not yet set (this is set in defineBasis method).\n"

        if hasattr(self, "noOfEigenvectors"):
            commonHeader += (
                " - Finding %d eigenvectors closest to the given pair-state\n"
                % self.noOfEigenvectors
            )

            if self.drivingFromState[0] < 0.1:
                commonHeader += (
                    " - State highlighting based on the relative contribution \n"
                    + "   of the original pair-state in the eigenstates obtained by diagonalization.\n"
                )
            else:
                commonHeader += (
                    " - State highlighting based on the relative driving strength \n"
                    + "   to a given energy eigenstate (energy level) from state\n"
                    + "   %s m_j =%d/2 with polarization q=%d.\n"
                    % (
                        printStateString(*self.drivingFromState[0:3]),
                        round(2.0 * self.drivingFromState[3]),
                        self.drivingFromState[4],
                    )
                )

        else:
            commonHeader += " ! Energy levels not yet found (this is done by calling diagonalise method).\n"

        if exportFormat == "csv":
            print("Exporting StarkMap calculation results as .csv ...")

            commonHeader += " - Export consists of three (3) files:\n"
            commonHeader += "       1) %s,\n" % (
                fileBase + "_r." + exportFormat
            )
            commonHeader += "       2) %s,\n" % (
                fileBase + "_energyLevels." + exportFormat
            )
            commonHeader += "       3) %s.\n\n" % (
                fileBase + "_highlight." + exportFormat
            )

            filename = fileBase + "_r." + exportFormat
            np.savetxt(
                filename,
                self.r,
                fmt="%.18e",
                delimiter=", ",
                newline="\n",
                header=(
                    commonHeader
                    + " - - - Interatomic distance, r (\\mu m) - - -"
                ),
                comments="# ",
            )
            print("   Interatomic distances (\\mu m) saved in %s" % filename)

            filename = fileBase + "_energyLevels." + exportFormat
            headerDetails = " NOTE : Each row corresponds to eigenstates for a single specified interatomic distance"
            np.savetxt(
                filename,
                self.y,
                fmt="%.18e",
                delimiter=", ",
                newline="\n",
                header=(
                    commonHeader + " - - - Energy (GHz) - - -\n" + headerDetails
                ),
                comments="# ",
            )
            print(
                "   Lists of energies (in GHz relative to the original pair-state energy)"
                + (" saved in %s" % filename)
            )

            filename = fileBase + "_highlight." + exportFormat
            np.savetxt(
                filename,
                self.highlight,
                fmt="%.18e",
                delimiter=", ",
                newline="\n",
                header=(
                    commonHeader
                    + " - - - Highlight value (rel.units) - - -\n"
                    + headerDetails
                ),
                comments="# ",
            )
            print("   Highlight values saved in %s" % filename)

            print("... data export finished!")
        else:
            raise ValueError("Unsupported export format (.%s)." % format)

    def _stateComposition(self, stateVector):
        contribution = np.absolute(stateVector)
        order = np.argsort(contribution, kind="heapsort")
        index = -1
        totalContribution = 0
        value = "$"
        while (index > -5) and (totalContribution < 0.95):
            i = order[index]
            if index != -1 and (
                stateVector[i].real > 0 or abs(stateVector[i].imag) > 1e-9
            ):
                value += "+"
            if abs(self.phi) < 1e-9:
                value = (
                    value
                    + ("%.2f" % stateVector[i])
                    + self._addState(*self.basisStates[i])
                )
            else:
                value = (
                    value
                    + (
                        "(%.2f+i%.2f)"
                        % (stateVector[i].real, stateVector[i].imag)
                    )
                    + self._addState(*self.basisStates[i])
                )
            totalContribution += contribution[i] ** 2
            index -= 1

        if totalContribution < 0.999:
            value += "+\\ldots"
        return value + "$"

    def _addState(self, n1, l1, j1, mj1, n2, l2, j2, mj2):
        stateString = ""
        if abs(self.s1 - 0.5) < 0.1:
            # Alkali atom
            stateString += "|%s %d/2" % (
                printStateStringLatex(n1, l1, j1, s=self.s1),
                round(2 * mj1),
            )
        else:
            # divalent atoms
            stateString += "|%s %d" % (
                printStateStringLatex(n1, l1, j1, s=self.s1),
                round(mj1),
            )

        if abs(self.s2 - 0.5) < 0.1:
            # Alkali atom
            stateString += ",%s %d/2\\rangle" % (
                printStateStringLatex(n2, l2, j2, s=self.s2),
                round(2 * mj2),
            )
        else:
            # divalent atom
            stateString += ",%s %d\\rangle" % (
                printStateStringLatex(n2, l2, j2, s=self.s2),
                round(mj2),
            )
        return stateString

    def plotLevelDiagram(
        self, highlightColor="red", highlightScale="linear", units="GHz"
    ):
        """
        Plots pair state level diagram

        Call :obj:`showPlot` if you want to display a plot afterwards.

        Parameters:
            highlightColor (string): optional, specifies the colour used
                for state highlighting
            highlightScale (string): optional, specifies scaling of
                state highlighting. Default is 'linear'. Use 'log-2' or
                'log-3' for logarithmic scale going down to 1e-2 and 1e-3
                respectively. Logarithmic scale is useful for spotting
                weakly admixed states.
            units (:obj:`char`,optional): possible values {'*GHz*','cm','eV'};
                [case insensitive] if value is 'GHz' (default), diagram will
                be plotted as energy :math:`/h` in units of GHz; if the
                string contains 'cm' diagram will be plotted in energy units
                cm :math:`{}^{-1}`; if the value is 'eV', diagram
                will be plotted as energy in units eV.
        """

        rvb = LinearSegmentedColormap.from_list(
            "mymap", ["0.9", highlightColor]
        )

        if units.lower() == "ev":
            self.units = "eV"
            self.scaleFactor = 1e9 * C_h / C_e
            eLabel = ""
        elif units.lower() == "ghz":
            self.units = "GHz"
            self.scaleFactor = 1
            eLabel = "/h"
        elif "cm" in units.lower():
            self.units = "cm$^{-1}$"
            self.scaleFactor = 1e9 / (C_c * 100)
            eLabel = "/(h c)"

        if highlightScale == "linear":
            cNorm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)
        elif highlightScale == "log-2":
            cNorm = matplotlib.colors.LogNorm(vmin=1e-2, vmax=1)
        elif highlightScale == "log-3":
            cNorm = matplotlib.colors.LogNorm(vmin=1e-3, vmax=1)
        else:
            raise ValueError(
                "Only 'linear', 'log-2' and 'log-3' are valid "
                "inputs for highlightScale"
            )

        print(" Now we are plotting...")
        self.fig, self.ax = plt.subplots(1, 1, figsize=(11.5, 5.0))

        self.y = np.array(self.y)
        self.highlight = np.array(self.highlight)

        colorfulX = []
        colorfulY = []
        colorfulState = []

        for i in xrange(len(self.r)):
            for j in xrange(len(self.y[i])):
                colorfulX.append(self.r[i])
                colorfulY.append(self.y[i][j])
                colorfulState.append(self.highlight[i][j])

        colorfulState = np.array(colorfulState)
        sortOrder = colorfulState.argsort(kind="heapsort")
        colorfulX = np.array(colorfulX)
        colorfulY = np.array(colorfulY)

        colorfulX = colorfulX[sortOrder]
        colorfulY = colorfulY[sortOrder]
        colorfulState = colorfulState[sortOrder]

        self.ax.scatter(
            colorfulX,
            colorfulY * self.scaleFactor,
            s=10,
            c=colorfulState,
            linewidth=0,
            norm=cNorm,
            cmap=rvb,
            zorder=2,
            picker=5,
        )
        cax = self.fig.add_axes([0.91, 0.1, 0.02, 0.8])
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=rvb, norm=cNorm)

        if self.drivingFromState[0] == 0:
            # colouring is based on the contribution of the original pair state here
            label = ""
            if abs(self.s1 - 0.5) < 0.1:
                # Alkali atom
                label += r"$|\langle %s m_j=%d/2 " % (
                    printStateStringLatex(self.n, self.l, self.j),
                    round(2.0 * self.m1),
                )
            else:
                # divalent atom
                label += r"$|\langle %s m_j=%d " % (
                    printStateStringLatex(self.n, self.l, self.j, s=self.s1),
                    round(self.m1),
                )

            if abs(self.s2 - 0.5) < 0.1:
                # Alkali atom
                label += r", %s m_j=%d/2 | \mu \rangle |^2$" % (
                    printStateStringLatex(self.nn, self.ll, self.jj),
                    round(2.0 * self.m2),
                )
            else:
                # divalent atom
                label += r", %s m_j=%d | \mu \rangle |^2$" % (
                    printStateStringLatex(self.nn, self.ll, self.jj, s=self.s2),
                    round(self.m2, 0),
                )

            cb.set_label(label)
        else:
            # colouring is based on the coupling to different states
            cb.set_label(r"$(\Omega_\mu/\Omega)^2$")

        self.ax.set_xlabel(r"Interatomic distance, $R$ ($\mu$m)")
        self.ax.set_ylabel(
            r"Pair-state relative energy, $\Delta E %s$ (%s)"
            % (eLabel, self.units)
        )

    def savePlot(self, filename="PairStateInteractions.pdf"):
        """
        Saves plot made by :obj:`PairStateInteractions.plotLevelDiagram`

        Args:
            filename (:obj:`str`, optional): file location where the plot
                should be saved
        """
        if self.fig != 0:
            self.fig.savefig(filename, bbox_inches="tight")
        else:
            print("Error while saving a plot: nothing is plotted yet")
        return 0

    def showPlot(self, interactive=True):
        """
        Shows level diagram printed by
        :obj:`PairStateInteractions.plotLevelDiagram`

        By default, it will output interactive plot, which means that
        clicking on the state will show the composition of the clicked
        state in original basis (dominant elements)

        Args:
            interactive (bool): optional, by default it is True. If true,
                plotted graph will be interactive, i.e. users can click
                on the state to identify the state composition

        Note:
            interactive=True has effect if the graphs are explored in usual
            matplotlib pop-up windows. It doesn't have effect on inline
            plots in Jupyter (IPython) notebooks.


        """
        if interactive:
            self.ax.set_title("Click on state to see state composition")
            self.clickedPoint = 0
            self.fig.canvas.draw()
            self.fig.canvas.mpl_connect("pick_event", self._onPick)

        plt.show()
        return 0

    def _onPick(self, event):
        if isinstance(event.artist, matplotlib.collections.PathCollection):
            x = event.mouseevent.xdata
            y = event.mouseevent.ydata / self.scaleFactor

            i = np.searchsorted(self.r, x)
            if i == len(self.r):
                i -= 1
            if (i > 0) and (abs(self.r[i - 1] - x) < abs(self.r[i] - x)):
                i -= 1

            j = 0
            for jj in xrange(len(self.y[i])):
                if abs(self.y[i][jj] - y) < abs(self.y[i][j] - y):
                    j = jj

            # now choose the most higlighted state in this area
            distance = abs(self.y[i][j] - y) * 1.5
            for jj in xrange(len(self.y[i])):
                if abs(self.y[i][jj] - y) < distance and (
                    abs(self.highlight[i][jj]) > abs(self.highlight[i][j])
                ):
                    j = jj

            if self.clickedPoint != 0:
                self.clickedPoint.remove()

            (self.clickedPoint,) = self.ax.plot(
                [self.r[i]],
                [self.y[i][j] * self.scaleFactor],
                "bs",
                linewidth=0,
                zorder=3,
            )

            self.ax.set_title(
                "State = "
                + self.composition[i][j]
                + ("   Colourbar = %.2f" % self.highlight[i][j]),
                fontsize=11,
            )

            event.canvas.draw()

    def getC6fromLevelDiagram(
        self, rStart, rStop, showPlot=False, minStateContribution=0.0
    ):
        """
        Finds :math:`C_6` coefficient for original pair state.

        Function first finds for each distance in the range
        [ `rStart` , `rStop` ] the eigen state with highest contribution of
        the original state. One can set optional parameter
        `minStateContribution` to value in the range [0,1), so that function
        finds only states if they have contribution of the original state
        that is bigger then `minStateContribution`.

        Once original pair-state is found in the range of interatomic
        distances, from smallest `rStart` to the biggest `rStop`, function
        will try to perform fitting of the corresponding state energy
        :math:`E(R)` at distance :math:`R` to the function
        :math:`A+C_6/R^6` where :math:`A` is some offset.

        Args:
            rStart (float): smallest inter-atomic distance to be used for fitting
            rStop (float): maximum inter-atomic distance to be used for fitting
            showPlot (bool): If set to true, it will print the plot showing
                fitted energy level and the obtained best fit. Default is
                False
            minStateContribution (float): valid values are in the range [0,1).
                It specifies minimum amount of the original state in the given
                energy state necessary for the state to be considered for
                the adiabatic continuation of the original unperturbed
                pair state.

        Returns:
            float:
                :math:`C_6` measured in :math:`\\text{GHz }\\mu\\text{m}^6`
                on success; If unsuccessful returns False.

        Note:
            In order to use this functions, highlighting in
            :obj:`diagonalise` should be based on the original pair
            state contribution of the eigenvectors (that this,
            `drivingFromState` parameter should not be set, which
            corresponds to `drivingFromState` = [0,0,0,0,0]).
        """
        initialStateDetuning = []
        initialStateDetuningX = []

        fromRindex = -1
        toRindex = -1
        for br in xrange(len(self.r)):
            if (fromRindex == -1) and (self.r[br] >= rStart):
                fromRindex = br
            if self.r[br] > rStop:
                toRindex = br - 1
                break
        if (fromRindex != -1) and (toRindex == -1):
            toRindex = len(self.r) - 1
        if fromRindex == -1:
            print(
                "\nERROR: could not find data for energy levels for interatomic"
            )
            print("distances between %2.f and %.2f mu m.\n\n" % (rStart, rStop))
            return 0

        for br in xrange(fromRindex, toRindex + 1):
            index = -1
            maxPortion = minStateContribution
            for br2 in xrange(len(self.y[br])):
                if abs(self.highlight[br][br2]) > maxPortion:
                    index = br2
                    maxPortion = abs(self.highlight[br][br2])
            if index != -1:
                initialStateDetuning.append(abs(self.y[br][index]))
                initialStateDetuningX.append(self.r[br])

        initialStateDetuning = np.log(np.array(initialStateDetuning))
        initialStateDetuningX = np.array(initialStateDetuningX)

        def c6fit(r, c6, offset):
            return np.log(c6 / r**6 + offset)

        try:
            popt, pcov = curve_fit(
                c6fit, initialStateDetuningX, initialStateDetuning, [1, 0]
            )
        except Exception as ex:
            print(ex)
            print("ERROR: unable to find a fit for C6.")
            return False
        print("c6 = ", popt[0], " GHz /R^6 (mu m)^6")
        print("offset = ", popt[1])

        y_fit = []
        for val in initialStateDetuningX:
            y_fit.append(c6fit(val, popt[0], popt[1]))
        y_fit = np.array(y_fit)

        if showPlot:
            fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.0))
            ax.loglog(
                initialStateDetuningX,
                np.exp(initialStateDetuning),
                "b-",
                lw=2,
                zorder=1,
            )
            ax.loglog(
                initialStateDetuningX, np.exp(y_fit), "r--", lw=2, zorder=2
            )

            ax.legend(
                ("calculated energy level", "fitted model function"),
                loc=1,
                fontsize=10,
            )

            ax.set_xlim(np.min(self.r), np.max(self.r))
            ymin = np.min(initialStateDetuning)
            ymax = np.max(initialStateDetuning)
            ax.set_ylim(exp(ymin), exp(ymax))

            minorLocator = mpl.ticker.MultipleLocator(1)
            minorFormatter = mpl.ticker.FormatStrFormatter("%d")
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_minor_formatter(minorFormatter)
            ax.xaxis.set_major_formatter(plt.NullFormatter())
            ax.set_xlabel(r"Interatomic distance, $r$ ($\mu$m)")
            ax.set_ylabel(r"Pair-state energy, $|E|$ (GHz)")
            ax.set_title(r"$C_6$ fit")
            plt.show()

        self.fitX = initialStateDetuningX
        self.fitY = initialStateDetuning
        self.fittedCurveY = y_fit

        return popt[0]

    def getC3fromLevelDiagram(
        self,
        rStart,
        rStop,
        showPlot=False,
        minStateContribution=0.0,
        resonantBranch=+1,
    ):
        """
        Finds :math:`C_3` coefficient for original pair state.

        Function first finds for each distance in the range
        [`rStart` , `rStop`] the eigen state with highest contribution of
        the original state. One can set optional parameter
        `minStateContribution` to value in the range [0,1), so that function
        finds only states if they have contribution of the original state
        that is bigger then `minStateContribution`.

        Once original pair-state is found in the range of interatomic
        distances, from smallest `rStart` to the biggest `rStop`, function
        will try to perform fitting of the corresponding state energy
        :math:`E(R)` at distance :math:`R` to the function
        :math:`A+C_3/R^3` where :math:`A` is some offset.

        Args:
            rStart (float): smallest inter-atomic distance to be used for fitting
            rStop (float): maximum inter-atomic distance to be used for fitting
            showPlot (bool): If set to true, it will print the plot showing
                fitted energy level and the obtained best fit. Default is
                False
            minStateContribution (float): valid values are in the range [0,1).
                It specifies minimum amount of the original state in the given
                energy state necessary for the state to be considered for
                the adiabatic continuation of the original unperturbed
                pair state.
            resonantBranch (int): optional, default +1. For resonant
                interactions we have two branches with identical
                state contributions. In this case, we will select only
                positively detuned branch (for resonantBranch = +1)
                or negatively detuned branch (fore resonantBranch = -1)
                depending on the value of resonantBranch optional parameter
        Returns:
            float:
                :math:`C_3` measured in :math:`\\text{GHz }\\mu\\text{m}^6`
                on success; If unsuccessful returns False.

        Note:
            In order to use this functions, highlighting in
            :obj:`diagonalise` should be based on the original pair
            state contribution of the eigenvectors (that this,
            `drivingFromState` parameter should not be set, which
            corresponds to `drivingFromState` = [0,0,0,0,0]).
        """

        selectBranch = False
        if abs(self.l - self.ll) == 1:
            selectBranch = True
            resonantBranch = float(resonantBranch)

        initialStateDetuning = []
        initialStateDetuningX = []

        fromRindex = -1
        toRindex = -1
        for br in xrange(len(self.r)):
            if (fromRindex == -1) and (self.r[br] >= rStart):
                fromRindex = br
            if self.r[br] > rStop:
                toRindex = br - 1
                break
        if (fromRindex != -1) and (toRindex == -1):
            toRindex = len(self.r) - 1
        if fromRindex == -1:
            print(
                "\nERROR: could not find data for energy levels for interatomic"
            )
            print("distances between %2.f and %.2f mu m.\n\n" % (rStart, rStop))
            return False

        discontinuityDetected = False
        for br in xrange(toRindex, fromRindex - 1, -1):
            index = -1
            maxPortion = minStateContribution
            for br2 in xrange(len(self.y[br])):
                if (abs(self.highlight[br][br2]) > maxPortion) and (
                    not selectBranch or (self.y[br][br2] * selectBranch > 0.0)
                ):
                    index = br2
                    maxPortion = abs(self.highlight[br][br2])

            if len(initialStateDetuningX) > 2:
                slope1 = (
                    initialStateDetuning[-1] - initialStateDetuning[-2]
                ) / (initialStateDetuningX[-1] - initialStateDetuningX[-2])
                slope2 = (abs(self.y[br][index]) - initialStateDetuning[-1]) / (
                    self.r[br] - initialStateDetuningX[-1]
                )
                if abs(slope2) > 3.0 * abs(slope1):
                    discontinuityDetected = True
            if (index != -1) and (not discontinuityDetected):
                initialStateDetuning.append(abs(self.y[br][index]))
                initialStateDetuningX.append(self.r[br])

        initialStateDetuning = np.log(np.array(initialStateDetuning))  # *1e9
        initialStateDetuningX = np.array(initialStateDetuningX)

        def c3fit(r, c3, offset):
            return np.log(c3 / r**3 + offset)

        try:
            popt, pcov = curve_fit(
                c3fit, initialStateDetuningX, initialStateDetuning, [1, 0]
            )
        except Exception as ex:
            print(ex)
            print("ERROR: unable to find a fit for C3.")
            return False
        print("c3 = ", popt[0], " GHz /R^3 (mu m)^3")
        print("offset = ", popt[1])

        y_fit = []

        for val in initialStateDetuningX:
            y_fit.append(c3fit(val, popt[0], popt[1]))
        y_fit = np.array(y_fit)

        if showPlot:
            fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.0))
            ax.loglog(
                initialStateDetuningX,
                np.exp(initialStateDetuning),
                "b-",
                lw=2,
                zorder=1,
            )
            ax.loglog(
                initialStateDetuningX, np.exp(y_fit), "r--", lw=2, zorder=2
            )

            ax.legend(
                ("calculated energy level", "fitted model function"),
                loc=1,
                fontsize=10,
            )
            ax.set_xlim(np.min(self.r), np.max(self.r))
            ymin = np.min(initialStateDetuning)
            ymax = np.max(initialStateDetuning)
            ax.set_ylim(exp(ymin), exp(ymax))

            minorLocator = mpl.ticker.MultipleLocator(1)
            minorFormatter = mpl.ticker.FormatStrFormatter("%d")
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_minor_formatter(minorFormatter)
            ax.xaxis.set_major_formatter(plt.NullFormatter())
            ax.set_xlabel(r"Interatomic distance, $r$ ($\mu$m)")
            ax.set_ylabel(r"Pair-state energy, $|E|$ (GHz)")

            locatorStep = 1.0
            while (locatorStep > (ymax - ymin)) and locatorStep > 1.0e-4:
                locatorStep /= 10.0

            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(locatorStep))
            ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.3f"))
            ax.yaxis.set_minor_locator(
                mpl.ticker.MultipleLocator(locatorStep / 10.0)
            )
            ax.yaxis.set_minor_formatter(plt.NullFormatter())
            # ax.yaxis.set_minor_formatter(mpl.ticker.FormatStrFormatter('%.3f'))

            ax.set_title(r"$C_3$ fit")

            plt.show()

        self.fitX = initialStateDetuningX
        self.fitY = initialStateDetuning
        self.fittedCurveY = y_fit

        return popt[0]

    def getVdwFromLevelDiagram(
        self, rStart, rStop, showPlot=False, minStateContribution=0.0
    ):
        """
        Finds :math:`r_{\\rm vdW}` coefficient for original pair state.

        Function first finds for each distance in the range [`rStart`,`rStop`]
        the eigen state with highest contribution of the original state.
        One can set optional parameter `minStateContribution` to value in
        the range [0,1), so that function finds only states if they have
        contribution of the original state that is bigger then
        `minStateContribution`.

        Once original pair-state is found in the range of interatomic
        distances, from smallest `rStart` to the biggest `rStop`, function
        will try to perform fitting of the corresponding state energy
        :math:`E(R)` at distance :math:`R` to the function
        :math:`A+B\\frac{1-\\sqrt{1+(r_{\\rm vdW}/r)^6}}{1-\\sqrt{1+r_{\\rm vdW}^6}}`
         where :math:`A` and :math:`B` are some offset.

        Args:
            rStart (float): smallest inter-atomic distance to be used for fitting
            rStop (float): maximum inter-atomic distance to be used for fitting
            showPlot (bool): If set to true, it will print the plot showing
                fitted energy level and the obtained best fit. Default is
                False
            minStateContribution (float): valid values are in the range [0,1).
                It specifies minimum amount of the original state in the given
                energy state necessary for the state to be considered for
                the adiabatic continuation of the original unperturbed
                pair state.
        Returns:
            float: :math:`r_{\\rm vdW}`  measured in :math:`\\mu\\text{m}`
            on success; If unsuccessful returns False.

        Note:
            In order to use this functions, highlighting in
            :obj:`diagonalise` should be based on the original pair
            state contribution of the eigenvectors (that this,
            `drivingFromState` parameter should not be set, which
            corresponds to `drivingFromState` = [0,0,0,0,0]).
        """

        initialStateDetuning = []
        initialStateDetuningX = []

        fromRindex = -1
        toRindex = -1
        for br in xrange(len(self.r)):
            if (fromRindex == -1) and (self.r[br] >= rStart):
                fromRindex = br
            if self.r[br] > rStop:
                toRindex = br - 1
                break
        if (fromRindex != -1) and (toRindex == -1):
            toRindex = len(self.r) - 1
        if fromRindex == -1:
            print(
                "\nERROR: could not find data for energy levels for interatomic"
            )
            print("distances between %2.f and %.2f mu m.\n\n" % (rStart, rStop))
            return False

        discontinuityDetected = False
        for br in xrange(toRindex, fromRindex - 1, -1):
            index = -1
            maxPortion = minStateContribution
            for br2 in xrange(len(self.y[br])):
                if abs(self.highlight[br][br2]) > maxPortion:
                    index = br2
                    maxPortion = abs(self.highlight[br][br2])
            if len(initialStateDetuningX) > 2:
                slope1 = (
                    initialStateDetuning[-1] - initialStateDetuning[-2]
                ) / (initialStateDetuningX[-1] - initialStateDetuningX[-2])
                slope2 = (abs(self.y[br][index]) - initialStateDetuning[-1]) / (
                    self.r[br] - initialStateDetuningX[-1]
                )
                if abs(slope2) > 3.0 * abs(slope1):
                    discontinuityDetected = True
            if (index != -1) and (not discontinuityDetected):
                initialStateDetuning.append(abs(self.y[br][index]))
                initialStateDetuningX.append(self.r[br])

        initialStateDetuning = np.log(abs(np.array(initialStateDetuning)))
        initialStateDetuningX = np.array(initialStateDetuningX)

        def vdwFit(r, offset, scale, vdw):
            return np.log(
                abs(
                    offset
                    + scale
                    * (1.0 - np.sqrt(1.0 + (vdw / r) ** 6))
                    / (1.0 - np.sqrt(1 + vdw**6))
                )
            )

        noOfPoints = len(initialStateDetuningX)
        print("Data points to fit = ", noOfPoints)

        try:
            popt, pcov = curve_fit(
                vdwFit,
                initialStateDetuningX,
                initialStateDetuning,
                [
                    0,
                    initialStateDetuning[noOfPoints // 2],
                    initialStateDetuningX[noOfPoints // 2],
                ],
            )
        except Exception as ex:
            print(ex)
            print("ERROR: unable to find a fit for van der Waals distance.")
            return False

        if (initialStateDetuningX[0] < popt[2]) or (
            popt[2] < initialStateDetuningX[-1]
        ):
            print("WARNING: vdw radius seems to be outside the fitting range!")
            print(
                "It's estimated to be around %.2f mu m from the current fit."
                % popt[2]
            )

        print("Rvdw =  ", popt[2], " mu m")
        print("offset = ", popt[0], "\n scale = ", popt[1])

        y_fit = []

        for val in initialStateDetuningX:
            y_fit.append(vdwFit(val, popt[0], popt[1], popt[2]))
        y_fit = np.array(y_fit)

        if showPlot:
            fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.0))
            ax.loglog(
                initialStateDetuningX,
                np.exp(initialStateDetuning),
                "b-",
                lw=2,
                zorder=1,
            )
            ax.loglog(
                initialStateDetuningX, np.exp(y_fit), "r--", lw=2, zorder=2
            )

            ax.set_xlim(np.min(self.r), np.max(self.r))
            ymin = np.min(initialStateDetuning)
            ymax = np.max(initialStateDetuning)
            ax.set_ylim(exp(ymin), exp(ymax))

            ax.axvline(x=popt[2], color="k")
            ax.text(
                popt[2],
                exp((ymin + ymax) / 2.0),
                r"$R_{vdw} = %.1f$ $\mu$m" % popt[2],
            )

            minorLocator = mpl.ticker.MultipleLocator(1)
            minorFormatter = mpl.ticker.FormatStrFormatter("%d")
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_minor_formatter(minorFormatter)
            ax.xaxis.set_major_formatter(plt.NullFormatter())
            ax.set_xlabel(r"Interatomic distance, $r$ ($\mu$m)")
            ax.set_ylabel(r"Pair-state energy, $|E|$ (GHz)")
            ax.legend(
                ("calculated energy level", "fitted model function"),
                loc=1,
                fontsize=10,
            )

            plt.show()

        self.fitX = initialStateDetuningX
        self.fitY = initialStateDetuning
        self.fittedCurveY = y_fit

        return popt[2]


class StarkMapResonances:
    """
    Calculates pair state Stark maps for finding resonances

    Tool for finding conditions for Foster resonances. For a given pair
    state, in a given range of the electric fields, looks for the pair-state
    that are close in energy and coupled via dipole-dipole interactions
    to the original pair-state.

    See `Stark resonances example snippet`_.

    .. _`Stark resonances example snippet`:
            ././Rydberg_atoms_a_primer.html#Tuning-the-interaction-strength-with-electric-fields

    Parameters:
        atom1 (:obj:`arc.alkali_atom_functions.AlkaliAtom` or :obj:`arc.divalent_atom_functions.DivalentAtom`):
            ={ :obj:`arc.alkali_atom_data.Lithium6`,
            :obj:`arc.alkali_atom_data.Lithium7`,
            :obj:`arc.alkali_atom_data.Sodium`,
            :obj:`arc.alkali_atom_data.Potassium39`,
            :obj:`arc.alkali_atom_data.Potassium40`,
            :obj:`arc.alkali_atom_data.Potassium41`,
            :obj:`arc.alkali_atom_data.Rubidium85`,
            :obj:`arc.alkali_atom_data.Rubidium87`,
            :obj:`arc.alkali_atom_data.Caesium`,
            :obj:`arc.divalent_atom_data.Strontium88`,
            :obj:`arc.divalent_atom_data.Calcium40`
            :obj:`arc.divalent_atom_data.Ytterbium174` }
             the first atom in the pair-state
        state1 ([int,int,float,float,(float)]): specification of the state
            of the first state as an array of values :math:`[n,l,j,m_j]`.
            For :obj:`arc.divalent_atom_functions.DivalentAtom` and other divalent atoms, 5th value
            should be added specifying total spin angular momentum `s`.
            Full definition of state then has format
            :math:`[n,l,j,m_j,s]`.
        atom2 (:obj:`arc.alkali_atom_functions.AlkaliAtom` or :obj:`arc.divalent_atom_functions.DivalentAtom`):
            ={ :obj:`arc.alkali_atom_data.Lithium6`,
            :obj:`arc.alkali_atom_data.Lithium7`,
            :obj:`arc.alkali_atom_data.Sodium`,
            :obj:`arc.alkali_atom_data.Potassium39`,
            :obj:`arc.alkali_atom_data.Potassium40`,
            :obj:`arc.alkali_atom_data.Potassium41`,
            :obj:`arc.alkali_atom_data.Rubidium85`,
            :obj:`arc.alkali_atom_data.Rubidium87`,
            :obj:`arc.alkali_atom_data.Caesium`,
            :obj:`arc.divalent_atom_data.Strontium88`,
            :obj:`arc.divalent_atom_data.Calcium40`
            :obj:`arc.divalent_atom_data.Ytterbium174` }
             the second atom in the pair-state
        state2 ([int,int,float,float,(float)]): specification of the state
            of the first state as an array of values :math:`[n,l,j,m_j]`,
            For :obj:`arc.divalent_atom_functions.DivalentAtom` and other divalent atoms, 5th value
            should be added specifying total spin angular momentum `s`.
            Full definition of state then has format
            :math:`[n,l,j,m_j,s]`.

    Note:
        In checking if certain state is dipole coupled to the original
        state, only the highest contributing state is checked for dipole
        coupling. This should be fine if one is interested in resonances
        in weak fields. For stronger fields, one might want to include
        effect of coupling to other contributing base states.



    """

    def __init__(self, atom1, state1, atom2, state2):
        self.atom1 = atom1
        if issubclass(type(self.atom1), DivalentAtom) and (
            len(state1) != 5 or (state1[4] != 0 and state1[4] != 1)
        ):
            raise ValueError(
                "For divalent atoms state specification has to "
                "include total spin angular momentum s as the last "
                "number in the state specification [n,l,j,m_j,s]."
            )
        self.state1 = state1
        # add exlicitly total spin of the state for Alkaline atoms
        if len(self.state1) == 4:
            self.state1.append(0.5)

        self.atom2 = atom2
        if issubclass(type(self.atom2), DivalentAtom) and (
            len(state1) != 5 or (state1[4] != 0 and state1[4] != 1)
        ):
            raise ValueError(
                "For divalent atoms state specification has to "
                "include total spin angular momentum s as the last "
                "numbre in the state specification [n,l,j,m_j,s]."
            )
        self.state2 = state2
        # add exlicitly total spin of the state for Alkaline atoms
        if len(self.state2) == 4:
            self.state2.append(0.5)

        self.pairStateEnergy = (
            (
                atom1.getEnergy(
                    self.state1[0],
                    self.state1[1],
                    self.state1[2],
                    s=self.state1[4],
                )
                + atom2.getEnergy(
                    self.state2[0],
                    self.state2[1],
                    self.state2[2],
                    s=self.state2[4],
                )
            )
            * C_e
            / C_h
            * 1e-9
        )

    def findResonances(
        self,
        nMin,
        nMax,
        maxL,
        eFieldList,
        energyRange=[-5.0e9, +5.0e9],
        Bz=0,
        progressOutput=False,
    ):
        r"""
        Finds near-resonant dipole-coupled pair-states

        For states in range of principal quantum numbers [`nMin`,`nMax`]
        and orbital angular momentum [0,`maxL`], for a range of electric fields
        given by `eFieldList` function will find near-resonant pair states.

        Only states that are in the range given by `energyRange` will be
        extracted from the pair-state Stark maps.

        Args:
            nMin (int): minimal principal quantum number of the state to be
                included in the StarkMap calculation
            nMax (int): maximal principal quantum number of the state to be
                included in the StarkMap calculation
            maxL (int): maximum value of orbital angular momentum for the states
                to be included in the calculation
            eFieldList ([float]): list of the electric fields (in V/m) for
                which to calculate level diagram (StarkMap)
            Bz (float): optional, magnetic field directed along z-axis in
                units of Tesla. Calculation will be correct only for weak
                magnetic fields, where paramagnetic term is much stronger
                then diamagnetic term. Diamagnetic term is neglected.
            energyRange ([float,float]): optinal argument. Minimal and maximal
                energy of that some dipole-coupled state should have in order
                to keep it in the plot (in units of Hz). By default it finds
                states that are :math:`\pm 5` GHz
            progressOutput (:obj:`bool`, optional): if True prints the
                progress of calculation; Set to false by default.
        """

        self.eFieldList = eFieldList
        self.Bz = Bz
        eMin = energyRange[0] * 1.0e-9  # in GHz
        eMax = energyRange[1] * 1.0e-9

        # find where is the original pair state

        sm1 = StarkMap(self.atom1)
        sm1.defineBasis(
            self.state1[0],
            self.state1[1],
            self.state1[2],
            self.state1[3],
            nMin,
            nMax,
            maxL,
            Bz=self.Bz,
            progressOutput=progressOutput,
            s=self.state1[4],
        )
        sm1.diagonalise(eFieldList, progressOutput=progressOutput)
        if (
            (self.atom2 is self.atom1)
            and (self.state1[0] == self.state2[0])
            and (self.state1[1] == self.state2[1])
            and (abs(self.state1[2] - self.state2[2]) < 0.1)
            and (abs(self.state1[3] - self.state2[3]) < 0.1)
            and (abs(self.state1[4] - self.state2[4]) < 0.1)
        ):
            sm2 = sm1
        else:
            sm2 = StarkMap(self.atom2)
            sm2.defineBasis(
                self.state2[0],
                self.state2[1],
                self.state2[2],
                self.state2[3],
                nMin,
                nMax,
                maxL,
                Bz=self.Bz,
                progressOutput=progressOutput,
                s=self.state2[4],
            )
            sm2.diagonalise(eFieldList, progressOutput=progressOutput)

        self.originalStateY = []
        self.originalStateContribution = []
        for i in xrange(len(sm1.eFieldList)):
            jmax1 = 0
            jmax2 = 0
            for j in xrange(len(sm1.highlight[i])):
                if sm1.highlight[i][j] > sm1.highlight[i][jmax1]:
                    jmax1 = j
            for j in xrange(len(sm2.highlight[i])):
                if sm2.highlight[i][j] > sm2.highlight[i][jmax2]:
                    jmax2 = j

            self.originalStateY.append(
                sm1.y[i][jmax1] + sm2.y[i][jmax2] - self.pairStateEnergy
            )
            self.originalStateContribution.append(
                (sm1.highlight[i][jmax1] + sm2.highlight[i][jmax2]) / 2.0
            )

        # M= mj1+mj2 is conserved with dipole-dipole interaction

        dmlist1 = [1, 0]
        if self.state1[3] != 0.5:
            dmlist1.append(-1)
        dmlist2 = [1, 0]
        if self.state2[3] != 0.5:
            dmlist2.append(-1)

        n1 = self.state1[0]
        l1 = self.state1[1] + 1
        j1 = self.state1[2] + 1
        mj1 = self.state1[3]

        n2 = self.state2[0]
        l2 = self.state2[1] + 1
        j2 = self.state2[2] + 1
        mj2 = self.state2[3]

        self.fig, self.ax = plt.subplots(1, 1, figsize=(9.0, 6))
        cm = LinearSegmentedColormap.from_list("mymap", ["0.9", "red", "black"])
        cNorm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)

        self.r = []
        self.y = []
        self.composition = []

        for dm1 in dmlist1:
            sm1.defineBasis(
                n1,
                l1,
                j1,
                mj1 + dm1,
                nMin,
                nMax,
                maxL,
                Bz=self.Bz,
                progressOutput=progressOutput,
                s=self.state1[4],
            )
            sm1.diagonalise(eFieldList, progressOutput=progressOutput)

            for dm2 in dmlist2:
                sm2.defineBasis(
                    n2,
                    l2,
                    j2,
                    mj2 + dm2,
                    nMin,
                    nMax,
                    maxL,
                    Bz=self.Bz,
                    progressOutput=progressOutput,
                    s=self.state2[4],
                )
                sm2.diagonalise(eFieldList, progressOutput=progressOutput)

                for i in xrange(len(sm1.eFieldList)):
                    yList = []
                    compositionList = []
                    if progressOutput:
                        sys.stdout.write("\rE=%.2f V/m " % sm1.eFieldList[i])
                        sys.stdout.flush()
                    for j in xrange(len(sm1.y[i])):
                        for jj in xrange(len(sm2.y[i])):
                            energy = (
                                sm1.y[i][j]
                                + sm2.y[i][jj]
                                - self.pairStateEnergy
                            )
                            statec1 = sm1.basisStates[
                                sm1.composition[i][j][0][1]
                            ]
                            statec2 = sm2.basisStates[
                                sm2.composition[i][jj][0][1]
                            ]
                            if (
                                (energy > eMin)
                                and (energy < eMax)
                                and (abs(statec1[1] - self.state1[1]) == 1)
                                and (abs(statec2[1] - self.state2[1]) == 1)
                            ):
                                # add this to PairStateMap
                                yList.append(energy)
                                compositionList.append(
                                    [
                                        sm1._stateComposition(
                                            sm1.composition[i][j]
                                        ),
                                        sm2._stateComposition(
                                            sm2.composition[i][jj]
                                        ),
                                    ]
                                )

                    if len(self.y) <= i:
                        self.y.append(yList)
                        self.composition.append(compositionList)
                    else:
                        self.y[i].extend(yList)
                        self.composition[i].extend(compositionList)

                if progressOutput:
                    print("\n")

        for i in xrange(len(sm1.eFieldList)):
            self.y[i] = np.array(self.y[i])
            self.composition[i] = np.array(self.composition[i])
            self.ax.scatter(
                [sm1.eFieldList[i] / 100.0] * len(self.y[i]),
                self.y[i],
                c="k",
                s=5,
                norm=cNorm,
                cmap=cm,
                lw=0,
                picker=5,
            )
        self.ax.plot(sm1.eFieldList / 100.0, self.originalStateY, "r-", lw=1)
        self.ax.set_ylim(eMin, eMax)
        self.ax.set_xlim(
            min(self.eFieldList) / 100.0, max(self.eFieldList) / 100.0
        )
        self.ax.set_xlabel("Electric field (V/cm)")
        self.ax.set_ylabel(r"Pair-state relative energy, $\Delta E/h$ (GHz)")

    def showPlot(self, interactive=True):
        """
        Plots initial state Stark map and its dipole-coupled resonances

        Args:

            interactive (optional, bool): if True (by default) points on plot
                will be clickable so that one can find the state labels
                and their composition (if they are heavily admixed).

        Note:
            Zero is given by the initial states of the atom given in
            initialisation of calculations, calculated **in absence of
            magnetic field B_z**. In other words, for non-zero magnetic
            field the inital states will have offset from zero even
            for zero electric field due to Zeeman shift.
        """

        if self.fig != 0:
            if interactive:
                self.ax.set_title("Click on state to see state composition")
                self.clickedPoint = 0
                self.fig.canvas.draw()
                self.fig.canvas.mpl_connect("pick_event", self._onPick)
            plt.show()
        else:
            print("Error while showing a plot: nothing is plotted yet")

    def _onPick(self, event):
        if isinstance(event.artist, matplotlib.collections.PathCollection):
            x = event.mouseevent.xdata * 100.0
            y = event.mouseevent.ydata

            i = np.searchsorted(self.eFieldList, x)
            if i == len(self.eFieldList):
                i -= 1
            if (i > 0) and (
                abs(self.eFieldList[i - 1] - x) < abs(self.eFieldList[i] - x)
            ):
                i -= 1

            j = 0
            for jj in xrange(len(self.y[i])):
                if abs(self.y[i][jj] - y) < abs(self.y[i][j] - y):
                    j = jj

            if self.clickedPoint != 0:
                self.clickedPoint.remove()

            (self.clickedPoint,) = self.ax.plot(
                [self.eFieldList[i] / 100.0],
                [self.y[i][j]],
                "bs",
                linewidth=0,
                zorder=3,
            )

            atom1 = self.atom1.elementName
            atom2 = self.atom2.elementName
            composition1 = str(self.composition[i][j][0])
            composition2 = str(self.composition[i][j][1])
            self.ax.set_title(
                ("[%s,%s]=[" % (atom1, atom2))
                + composition1
                + ","
                + composition2
                + "]",
                fontsize=10,
            )

            event.canvas.draw()

    def _onPick2(self, xdata, ydata):
        x = xdata * 100.0
        y = ydata

        i = np.searchsorted(self.eFieldList, x)
        if i == len(self.eFieldList):
            i -= 1
        if (i > 0) and (
            abs(self.eFieldList[i - 1] - x) < abs(self.eFieldList[i] - x)
        ):
            i -= 1

        j = 0
        for jj in xrange(len(self.y[i])):
            if abs(self.y[i][jj] - y) < abs(self.y[i][j] - y):
                j = jj

        if self.clickedPoint != 0:
            self.clickedPoint.remove()

        (self.clickedPoint,) = self.ax.plot(
            [self.eFieldList[i] / 100.0],
            [self.y[i][j]],
            "bs",
            linewidth=0,
            zorder=3,
        )

        atom1 = self.atom1.elementName
        atom2 = self.atom2.elementName
        composition1 = str(self.composition[i][j][0])
        composition2 = str(self.composition[i][j][1])
        self.ax.set_title(
            ("[%s,%s]=[" % (atom1, atom2))
            + composition1
            + ","
            + composition2
            + "]",
            fontsize=10,
        )
