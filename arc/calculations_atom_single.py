# -*- coding: utf-8 -*-

"""
    This module provides calculations of single-atom properties.

    Included calculations are Stark maps, level plot visualisations,
    lifetimes and radiative decays.

"""

from __future__ import print_function

from math import exp,log,sqrt
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import re
from .wigner import Wigner6j,Wigner3j,CG
from scipy.constants import physical_constants, pi , epsilon_0, hbar
from scipy.constants import k as C_k
from scipy.constants import c as C_c
from scipy.constants import h as C_h
from scipy.constants import e as C_e
from scipy.optimize import curve_fit

# for matrices
from numpy import zeros,savetxt, complex64,complex128
from numpy.linalg import eigvalsh,eig,eigh
from numpy.ma import conjugate
from numpy.lib.polynomial import real

from scipy.sparse import lil_matrix,csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.special.specfun import fcoef

import sys
if sys.version_info > (2,):
    xrange = range
from .alkali_atom_functions import printStateString, _EFieldCoupling, printStateLetter,printStateStringLatex

from matplotlib.colors import LinearSegmentedColormap
import matplotlib

import sqlite3
sqlite3.register_adapter(np.float64, float)
sqlite3.register_adapter(np.float32, float)
sqlite3.register_adapter(np.int64, int)
sqlite3.register_adapter(np.int32, int)
import datetime

class StarkMap:
    """
        Calculates Stark maps for single atom in a field

        This initializes calculation for the atom of a given type. For details
        of calculation see Zimmerman [1]_. For a quick working example
        see `Stark map example snippet`_.



        Args:
            atom (:obj:`AlkaliAtom`): ={ :obj:`alkali_atom_data.Lithium6`,
                :obj:`alkali_atom_data.Lithium7`,
                :obj:`alkali_atom_data.Sodium`,
                :obj:`alkali_atom_data.Potassium39`,
                :obj:`alkali_atom_data.Potassium40`,
                :obj:`alkali_atom_data.Potassium41`,
                :obj:`alkali_atom_data.Rubidium85`,
                :obj:`alkali_atom_data.Rubidium87`,
                :obj:`alkali_atom_data.Caesium` }
                Select the alkali metal for energy level
                diagram calculation



        Examples:
            State :math:`28~S_{1/2}~|m_j|=0.5` polarizability calculation

            >>> from arc import *
            >>> calc = StarkMap(Caesium())
            >>> calc.defineBasis(28, 0, 0.5, 0.5, 23, 32, 20)
            >>> calc.diagonalise(np.linspace(00.,6000,600))
            >>> print("%.5f MHz cm^2 / V^2 " % calc.getPolarizability())
            0.76705 MHz cm^2 / V^2

            Stark map calculation

            >>> from arc import *
            >>> calc = StarkMap(Caesium())
            >>> calc.defineBasis(28, 0, 0.5, 0.5, 23, 32, 20)
            >>> calc.diagonalise(np.linspace(00.,60000,600))
            >>> calc.plotLevelDiagram()
            >>> calc.showPlot()
            << matplotlib plot will open containing a Stark map >>

        Examples:
            **Advanced interfacing of Stark map calculations (StarkMap class)**
            Here we show one easy way to obtain the Stark matrix (from diagonal
            :obj:`mat1` and off-diagonal part :obj:`mat2` ) and basis states
            (stored in :obj:`basisStates` ), if this middle-product of the
            calculation is needed for some code build on top of the existing
            ARC package.

            >>> from arc import *
            >>> calc = StarkMap(Caesium())
            >>> calc.defineBasis(28, 0, 0.5, 0.5, 23, 32, 20)
            >>> # Now we have matrix and basis states, that we can used in our own code
            >>> # Let's say we want Stark map at electric field of 0.2 V/m
            >>> eField = 0.2 # V/m
            >>> # We can easily extract Stark matrix
            >>> # as diagonal matrix (state detunings)
            >>> #  + off-diagonal matrix (propotional to electric field)
            >>> matrix = calc.mat1+calc.mat2*eField
            >>> # and the basis states as array [ [n,l,j,mj] , ...]
            >>> basisStates = calc.basisStates
            >>> # you can do your own calculation now...

        References:
            .. [1] M. L. Zimmerman et.al, PRA **20**:2251 (1979)
                https://doi.org/10.1103/PhysRevA.20.2251

        .. _`Stark map example snippet`:
            ./Rydberg_atoms_a_primer.html#Rydberg-Atom-Stark-Shifts
    """

    def __init__(self,atom):

        self.atom = atom

        self.basisStates = []
        """
            List of basis states for calculation in the form [ [n,l,j,mj], ...].
            Calculated by :obj:`defineBasis` .
        """
        self.mat1 = []
        """
            diagonal elements of Stark-matrix (detuning of states) calculated by
            :obj:`defineBasis` in the basis :obj:`basisStates`.
        """
        self.mat2 = []
        """
            off-diagonal elements of Stark-matrix divided by electric
            field value. To get off diagonal elemements multiply this matrix
            with electric field value. Full Stark matrix is obtained as
            `fullStarkMatrix` = :obj:`mat1` + :obj:`mat2` *`eField`. Calculated by
            :obj:`defineBasis` in the basis :obj:`basisStates`.
        """
        self.indexOfCoupledState = []
        """
            Index of coupled state (initial state passed to :obj:`defineBasis`)
            in :obj:`basisStates` list of basis states
        """


        # finding energy levels
        self.eFieldList = []
        """
        Saves electric field (in units of V/m) for which energy levels are calculated

        See also:
            :obj:`y`, :obj:`highlight`, :obj:`diagonalise`
        """
        self.y = []  # eigenValues
        """
        `y[i]` is an array of eigenValues corresponding to the energies of the
        atom states at the electric field `eFieldList[i]`. For example `y[i][j]` is
        energy of the `j` eigenvalue (energy of the state) measured in
        cm :math:`{}^{-1}` relative to the ionization threshold.

        See also:
            :obj:`eFieldList`, :obj:`highlight`, :obj:`diagonalise`
        """
        self.highlight = [] #contribution of initial state there (overlap |<original state | given state>|^2)
        """
        `highlight[i]` is an array of values measuring highlighted feature in the
        eigenstates at electric field intensity `eFieldList[i]`. E.g. `highlight[i][j]`
        measures highlighted feature of the state with energy `y[i][j]` at electric
        field `eFieldList[i]`. What will be highlighted feature is defined in the
        call of :obj:`diagonalise` (see that part of documentation for details).

        See also:
            :obj:`eFieldList`, :obj:`y`, :obj:`diagonalise`
        """

        # pointers towards figure
        self.fig = 0
        self.ax = 0

        # values used for fitting polarizability, and fit
        self.fitX = []
        self.fitY = []
        self.fittedCurveY = []

        self.drivingFromState = [0,0,0,0,0]
        self.maxCoupling = 0.

        # STARK memoization
        self.eFieldCouplingSaved = _EFieldCoupling()



    def _eFieldCouplingDivE(self,n1,l1,j1,mj1,n2,l2,j2,mj2):
        # eFied coupling devided with E (witout actuall multiplication to getE)
        # delta(mj1,mj2') delta(l1,l2+-1)
        if ( (abs(mj1-mj2)>0.1) or (abs(l1-l2) !=1) ):
            return 0

        # matrix element
        result = self.atom.getRadialMatrixElement(n1,l1,j1,n2,l2,j2)*\
                physical_constants["Bohr radius"][0]*C_e

        sumPart = self.eFieldCouplingSaved.getAngular(l1,j1,mj1,l2,j2,mj2)

        return result*sumPart

    def _eFieldCoupling(self,n1,l1,j1,mj1,n2,l2,j2,mj2,eField):
        return self._eFieldCouplingDivE(n1,l1,j1,mj1,n2,l2,j2,mj2)*eField


    def defineBasis(self,n,l,j,mj,nMin,nMax,maxL,progressOutput = False,\
                    debugOutput=False):
        """
            Initializes basis of states around state of interest

            Defines basis of states for further calculation. :math:`n,l,j,m_j`
            specify state whose neighbourhood and polarizability we want
            to explore. Other parameters specify basis of calculations.
            This method stores basis in :obj:`basisStates`, while corresponding
            interaction matrix is stored in two parts. First part is diagonal
            electric-field independent part stored in :obj:`mat1`, while the
            second part :obj:`mat2` corresponds to off-diagonal elements that are
            propotional to electric field. Overall interaction matrix for
            electric field `eField` can be then obtained as
            `fullStarkMatrix` = :obj:`mat1` + :obj:`mat2` *`eField`

            Args:
                n (int): principal quantum number of the state
                l (int): angular orbital momentum of the state
                j (flaot): total angular momentum of the state
                mj (float): projection of total angular momentum of the state
                nMin (int): *minimal* principal quantum number of the states to
                    be included in the basis for calculation
                nMax (int): *maximal* principal quantum number of the states to
                    be included in the basis for calculation
                maxL (int): *maximal* value of orbital angular momentum for the
                    states to be included in the basis for calculation
                progressOutput (:obj:`bool`, optional): if True prints the
                    progress of calculation; Set to false by default.
                debugOutput (:obj:`bool`, optional): if True prints additional
                    information usefull for debuging. Set to false by default.
        """
        global wignerPrecal
        wignerPrecal = True

        states = []

        # save calculation details START
        self.n = n; self.l =l; self.j=j
        self.mj = mj; self.nMin = nMin; self.nMax = nMax; self.maxL = maxL
        # save calculation details END


        for tn in xrange(nMin,nMax):

            for tl in xrange(min(maxL+1,tn)):
                if (abs(mj)-0.1<=float(tl)+0.5):
                    states.append([tn,tl,float(tl)+0.5,mj])

                if (tl>0) and  (abs(mj)-0.1<=float(tl)-0.5):
                    states.append([tn,tl,float(tl)-0.5,mj])

        dimension = len(states)
        if progressOutput:
            print("Found ",dimension," states.")
            if debugOutput:
                print(states)

        indexOfCoupledState = 0
        index = 0
        for s in states:
            if (s[0]==n) and (abs(s[1]-l)<0.1) and (abs(s[2]-j)<0.1) and\
                 (abs(s[3]-mj)<0.1):
                indexOfCoupledState = index
            index +=1
        if debugOutput:
            print("Index of initial state")
            print(indexOfCoupledState)
            print("Initial state = ")
            print(states[indexOfCoupledState])


        self.mat1 = np.zeros((dimension,dimension),dtype=np.double)
        self.mat2 = np.zeros((dimension,dimension),dtype=np.double)

        self.basisStates = states
        self.indexOfCoupledState = indexOfCoupledState

        if progressOutput:
            print("Generating matrix...")
        progress = 0.

        for ii in xrange(dimension):
            if progressOutput:
                progress += ((dimension-ii)*2-1)
                sys.stdout.write("\r%d%%" % (float(progress)/float(dimension**2)*100))
                sys.stdout.flush()

            # add diagonal element
            self.mat1[ii][ii] = self.atom.getEnergy(states[ii][0],\
                                               states[ii][1],states[ii][2])\
                            *C_e/C_h*1e-9
            # add off-diagonal element

            for jj in xrange(ii+1,dimension):
                coupling = self._eFieldCouplingDivE(states[ii][0]\
                                                    ,states[ii][1],\
                                                    states[ii][2],mj,\
                                                    states[jj][0],\
                                                    states[jj][1],\
                                                    states[jj][2],mj)*\
                            1.e-9/C_h
                self.mat2[jj][ii] = coupling
                self.mat2[ii][jj] = coupling

        if progressOutput:
            print("\n")
        if debugOutput:
            print(self.mat1+self.mat2)
            print(self.mat2[0])

        self.atom.updateDipoleMatrixElementsFile()
        return 0

    def diagonalise(self,eFieldList,drivingFromState = [0,0,0,0,0],
                        progressOutput=False,debugOutput=False):
        """
            Finds atom eigenstates in a given electric field

            Eigenstates are calculated for a list of given electric fields. To
            extract polarizability of the originaly stated state see
            :obj:`getPolarizability` method. Results are saved in
            :obj:`eFieldList`, :obj:`y` and :obj:`highlight`.

            Args:
                eFieldList (array): array of electric field strength (in V/m)
                    for which we want to know energy eigenstates

                progressOutput (:obj:`bool`, optional): if True prints the
                    progress of calculation; Set to false by default.
                debugOutput (:obj:`bool`, optional): if True prints additional
                    information usefull for debuging. Set to false by default.
        """

        # if we are driving from some state
        # ========= FIND LASER COUPLINGS (START) =======

        coupling = []
        dimension = len(self.basisStates)
        self.maxCoupling = 0.
        self.drivingFromState = drivingFromState
        if (self.drivingFromState[0] != 0):
            if progressOutput: print("Finding driving field coupling...")
            # get first what was the state we are calculating coupling with
            state1 = drivingFromState
            n1 = int(round(state1[0]))
            l1 = int(round(state1[1]))
            j1 = state1[2]
            m1 = state1[3]
            q = state1[4]


            for i in xrange(dimension):
                thisCoupling = 0.
                if progressOutput:
                    sys.stdout.write("\r%d%%" %  (i/float(dimension-1)*100.))
                    sys.stdout.flush()
                if (int(abs(self.basisStates[i][1]-l1))==1)and\
                    (int(abs(self.basisStates[i][2]-j1))<=1) and\
                    (int(abs(self.basisStates[i][3]-m1-q))==0):
                    state2 = self.basisStates[i]
                    n2 = int(state2[0])
                    l2 = int(state2[1])
                    j2 = state2[2]
                    m2 = state2[3]
                    if debugOutput:
                        print(n1," ",l1," ",j1," ",m1," < - ",q," - >",n2," ",\
                            l2," ",j2," ",m2,"\n")
                    dme = self.atom.getDipoleMatrixElement(n1, l1,j1,m1,\
                                                            n2,l2,j2,m2,\
                                                            q)
                    thisCoupling += dme
                thisCoupling = abs(thisCoupling)**2
                if thisCoupling > self.maxCoupling:
                    self.maxCoupling = thisCoupling
                if (thisCoupling >0.00000001) and debugOutput:
                    print("coupling = ",thisCoupling)
                coupling.append(thisCoupling)

            if progressOutput:
                print("\n")

            if self.maxCoupling<0.00000001:
                raise Exception("State that you specified in drivingFromState, for a "+\
                "given laser polarization, is uncoupled from the specified Stark "+\
                "manifold. If you just want to see the specified Stark manifold "+\
                "remove driveFromState optional argument from call of function "+\
                "diagonalise. Or specify state and driving that is coupled "+\
                "to a given manifold to see coupling strengths.")

        # ========= FIND LASER COUPLINGS (END) =======


        indexOfCoupledState = self.indexOfCoupledState
        self.eFieldList = eFieldList

        self.y = []
        self.highlight = []
        self.composition = []

        if progressOutput:
            print("Finding eigenvectors...")
        progress = 0.
        for eField in eFieldList:
            if progressOutput:
                progress += 1.
                sys.stdout.write("\r%d%%" % \
                                 (float(progress)/float(len(eFieldList))*100))
                sys.stdout.flush()

            m = self.mat1+self.mat2*eField

            ev,egvector = eigh(m)

            self.y.append(ev)
            if (drivingFromState[0]<0.1):
                sh = []
                comp = []
                for i in xrange(len(ev)):
                    sh.append(abs(egvector[indexOfCoupledState,i])**2)
                    comp.append(self._stateComposition2(egvector[:,i]))
                self.highlight.append(sh)
                self.composition.append(comp)
            else:
                sh = []
                comp = []
                for i in xrange(len(ev)):
                    sumCoupledStates = 0.
                    for j in xrange(dimension):
                        sumCoupledStates += abs(coupling[j]/self.maxCoupling)*\
                                                abs(egvector[j,i]**2)
                    comp.append(self._stateComposition2(egvector[:,i]))
                    sh.append(sumCoupledStates)
                self.highlight.append(sh)
                self.composition.append(comp)


        if progressOutput:
            print("\n")
        return

    def exportData(self,fileBase,exportFormat = "csv"):
        """
            Exports StarkMap calculation data.

            Only supported format (selected by default) is .csv in a
            human-readable form with a header that saves details of calculation.
            Function saves three files: 1) `filebase` _eField.csv;
            2) `filebase` _energyLevels
            3) `filebase` _highlight

            For more details on the format, see header of the saved files.

            Args:
                filebase (string): filebase for the names of the saved files
                    without format extension. Add as a prefix a directory path
                    if necessary (e.g. saving outside the current working directory)
                exportFormat (string): optional. Format of the exported file. Currently
                    only .csv is supported but this can be extended in the future.
        """

        fmt='on %Y-%m-%d @ %H:%M:%S'
        ts = datetime.datetime.now().strftime(fmt)

        commonHeader = "Export from Alkali Rydberg Calculator (ARC) %s.\n" % ts
        commonHeader += ("\n *** Stark Map for %s %s m_j = %d/2. ***\n\n" % (self.atom.elementName,
                        printStateString(self.n, self.l, self.j), int(round(2.*self.mj)) ) )
        commonHeader += (" - Included states - principal quantum number (n) range [%d-%d].\n" %\
                         (self.nMin, self.nMax))
        commonHeader += (" - Included states with orbital momentum (l) in range [%d,%d] (i.e. %s-%s).\n"%\
                         (0, self.maxL, printStateLetter(0), printStateLetter(self.maxL)))
        if self.drivingFromState[0]<0.1:
            commonHeader += " - State highlighting based on the relative contribution \n"+\
            "   of the original state in the eigenstates obtained by diagonalization."
        else:
            commonHeader += (" - State highlighting based on the relative driving strength \n"+\
            "   to a given energy eigenstate (energy level) from state\n"+\
            "   %s m_j =%d/2 with polarization q=%d.\n"%\
             ( printStateString(*self.drivingFromState[0:3]),\
             int(round(2.*self.drivingFromState[3])),
             self.drivingFromState[4]))


        if exportFormat=="csv":
            print("Exporting StarkMap calculation results as .csv ...")

            commonHeader += " - Export consists of three (3) files:\n"
            commonHeader += ("       1) %s,\n" % (fileBase+"_eField."+exportFormat))
            commonHeader += ("       2) %s,\n" % (fileBase+"_energyLevels."+exportFormat))
            commonHeader += ("       3) %s.\n\n" % (fileBase+"_highlight."+exportFormat))

            filename = fileBase+"_eField."+exportFormat
            np.savetxt(filename, \
                self.eFieldList, fmt='%.18e', delimiter=', ',\
                newline='\n', \
                header=(commonHeader + " - - - eField (V/m) - - -"),\
                comments='# ')
            print("   Electric field values (V/m) saved in %s" % filename)

            filename = fileBase+"_energyLevels."+exportFormat
            headerDetails = " NOTE : Each row corresponds to eigenstates for a single specified electric field"
            np.savetxt(filename, \
                self.y, fmt='%.18e', delimiter=', ',\
                newline='\n', \
                header=(commonHeader + ' - - - Energy (GHz) - - -\n' + headerDetails),\
                comments='# ')
            print("   Lists of energies (in GHz relative to ionisation) saved in %s" % filename)

            filename = fileBase+"_highlight."+exportFormat
            np.savetxt(filename, \
                self.highlight, fmt='%.18e', delimiter=', ',\
                newline='\n', \
                header=(commonHeader + ' - - - Highlight value (rel.units) - - -\n'+ headerDetails),\
                comments='# ')
            print("   Highlight values saved in %s" % filename)

            print("... data export finished!")
        else:
            raise ValueError("Unsupported export format (.%s)." % format)

    def plotLevelDiagram(self,units=1,highlighState=True,progressOutput=False,\
                        debugOutput=False,highlightColour='red'):
        """
            Makes a plot of a stark map of energy levels

            To save this plot, see :obj:`savePlot`. To print this plot see
            :obj:`showPlot`.

            Args:
                units (:obj:`int`,optional): possible values {1,2} ; if the
                    value is 1 (default) Stark diagram will be plotted in
                    energy units cm :math:`{}^{-1}`; if value is 2, Stark
                    diagram will be plotted as energy :math:`/h` in units of GHz
                highlightState (:obj:`bool`, optional): False by default. If
                    True, scatter plot colour map will map in red amount of
                    original state for the given eigenState
                progressOutput (:obj:`bool`, optional): if True prints the
                    progress of calculation; Set to False by default.
                debugOutput (:obj:`bool`, optional): if True prints additional
                    information usefull for debuging. Set to False by default.
        """
        rvb = LinearSegmentedColormap.from_list('mymap',\
                                               ['0.9', highlightColour,'black'])

        self.units = units

        if progressOutput:
            print("plotting...")

        originalState = self.basisStates[self.indexOfCoupledState]
        n = originalState[0]
        l = originalState[1]
        j = originalState[2]

        existingPlot = False
        if (self.fig == 0):
            self.fig, self.ax = plt.subplots(1,1,figsize=(11.,5))
        else:
            existingPlot = True


        eFieldList = []
        y =[]
        yState = []


        for br in xrange(len(self.y)):

            for i in xrange(len(self.y[br])):
                eFieldList.append(self.eFieldList[br])
                y.append(self.y[br][i])
                yState.append(self.highlight[br][i])

        yState = np.array(yState)
        sortOrder = yState.argsort(kind='heapsort')
        eFieldList = np.array(eFieldList)
        y = np.array(y)

        eFieldList = eFieldList[sortOrder]
        y = y[sortOrder]
        yState = yState[sortOrder]


        if (units==1):
            ## in cm^-1


            if not highlighState:
                self.ax.scatter(eFieldList/100.,y*0.03336,s=1,color="k",picker=5)
            else:
                cm = rvb
                cNorm  = matplotlib.colors.Normalize(vmin=0., vmax=1.)
                self.ax.scatter(eFieldList/100,y*0.03336,\
                                c=yState,s=5,norm=cNorm, cmap=cm,lw=0,picker=5)
                if not existingPlot:
                    cax = self.fig.add_axes([0.91, 0.1, 0.02, 0.8])
                    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cm, norm=cNorm)
                    if (self.drivingFromState[0]<0.1):
                        cb.set_label(r"$|\langle %s | \mu \rangle |^2$" % \
                                 printStateStringLatex(n,l,j))
                    else:
                        cb.set_label(r"$( \Omega_\mu | \Omega )^2$")


        else:
            ## in GHz

            if not highlighState:
                self.ax.scatter(eFieldList/100.,y,\
                                s=1,color="k",picker=5) # in GHz
            else:
                cm = rvb
                cNorm  = matplotlib.colors.Normalize(vmin=0., vmax=1.)
                self.ax.scatter(eFieldList/100.,y,c=yState,\
                                s=5,norm=cNorm, cmap=cm,lw=0,picker=5)
                if not existingPlot:
                    cax = self.fig.add_axes([0.91, 0.1, 0.02, 0.8])
                    cb = matplotlib.colorbar.ColorbarBase(cax, \
                                                          cmap=cm, norm=cNorm)
                    if (self.drivingFromState[0]<0.1):
                        cb.set_label(r"$|\langle %s | \mu \rangle |^2$" %\
                                  printStateStringLatex(n,l,j))
                    else:
                        cb.set_label(r"$(\Omega_\mu / \Omega )^2$")

        self.ax.set_xlabel("Electric field (V/cm)")


        if (units==1):
            ## in cm^{-1}
            uppery = self.atom.getEnergy(n,l,j)*C_e/C_h*1e-9*0.03336+10
            lowery = self.atom.getEnergy(n,l,j)*C_e/C_h*1e-9*0.03336-10
            self.ax.set_ylabel("State energy, $E/(h c)$ (cm$^{-1}$)")
        else:
            ## in GHz
            uppery = self.atom.getEnergy(n,l,j)*C_e/C_h*1e-9+5
            lowery = self.atom.getEnergy(n,l,j)*C_e/C_h*1e-9-5
            self.ax.set_ylabel(r"State energy, $E/h$ (GHz)")


        self.ax.set_ylim(lowery,uppery)
        ##
        self.ax.set_xlim(min(eFieldList)/100.,max(eFieldList)/100.)
        return 0

    def savePlot(self,filename="StarkMap.pdf"):
        """
            Saves plot made by :obj:`plotLevelDiagram`

            Args:
                filename (:obj:`str`, optional): file location where the plot
                    should be saved
        """
        if (self.fig != 0):
            self.fig.savefig(filename,bbox_inches='tight')
        else:
            print("Error while saving a plot: nothing is plotted yet")
        return 0

    def showPlot(self, interactive = True):
        """
            Shows plot made by :obj:`plotLevelDiagram`
        """
        if (self.fig != 0):
            if interactive:
                self.ax.set_title("Click on state to see state composition")
                self.clickedPoint = 0
                self.fig.canvas.draw()
                self.fig.canvas.mpl_connect('pick_event', self._onPick)
            plt.show()
            self.fig.clear()
            self.fig = 0
            self.ax = 0
        else:
            print("Error while showing a plot: nothing is plotted yet")
        return 0

    def _onPick(self,event):
        if isinstance(event.artist, matplotlib.collections.PathCollection):
            if (self.units==1):
                scaleFactor = 0.03336

            x = event.mouseevent.xdata*100.
            y = event.mouseevent.ydata/scaleFactor

            i = np.searchsorted(self.eFieldList,x)
            if ((i>0) and (abs(self.eFieldList[i-1]-x)<abs(self.eFieldList[i]-x))):
                i -=1

            j = 0
            for jj in xrange(len(self.y[i])):
                if (abs(self.y[i][jj]-y) < abs(self.y[i][j]-y)):
                    j = jj

            # now choose the most higlighted state in this area
            distance = abs(self.y[i][j]-y)*1.5
            for jj in xrange(len(self.y[i])):
                if (abs(self.y[i][jj]-y) < distance and \
                    (abs(self.highlight[i][jj])>abs(self.highlight[i][j]))):
                    j = jj

            if (self.clickedPoint!=0):
                self.clickedPoint.remove()

            self.clickedPoint, = self.ax.plot([self.eFieldList[i]/100.],\
                                               [self.y[i][j]*scaleFactor],"bs",\
                                                 linewidth=0,zorder=3)

            self.ax.set_title(("[%s] = " % self.atom.elementName)+\
                              self._stateComposition(self.composition[i][j])+\
                             ("   Colourbar value = %.2f"% self.highlight[i][j]),
                             fontsize=11)

            event.canvas.draw()

    def _stateComposition(self,stateVector):
        i = 0
        totalContribution = 0
        value = "$"
        while (i<len(stateVector)) and (totalContribution<0.95):
            if (i!=0 and stateVector[i][0]>0):
                value+= "+"
            value = value+ ("%.2f" % stateVector[i][0])+\
                    self._addState(*self.basisStates[stateVector[i][1]])
            totalContribution += abs(stateVector[i][0])**2
            i += 1

        if totalContribution<0.999:
            value+="+\\ldots"
        return value+"$"


    def _stateComposition2(self,stateVector,upTo=4):
        contribution = np.absolute(stateVector)
        order = np.argsort(contribution,kind='heapsort')
        index = -1
        totalContribution = 0
        mainStates = []  #[state Value, state index]
        while (index>-upTo) and (totalContribution<0.95):
            i = order[index]
            mainStates.append([stateVector[i],i])
            totalContribution += contribution[i]**2
            index -= 1
        return mainStates

    def _addState(self,n1,l1,j1,mj1):
        return "|%s m_j=%d/2\\rangle" %\
             (printStateStringLatex(n1, l1, j1),int(2*mj1))

    def getPolarizability(self, maxField=1.e10, showPlot = False,\
                           debugOutput = False, minStateContribution=0.0):
        """
            Returns the polarizability of the state (set during the
            initalization process)

            Args:
                maxField (:obj:`float`, optional): maximum field (in V/m) to be
                    used for fitting the polarizability. By default, max field
                    is very large, so it will use eigenvalues calculated in the
                    whole range.
                showPlot (:obj:`bool`, optional): shows plot of calculated
                    eigenValues of the given state (dots), and the fit (solid
                    line) for extracting polarizability
                debugOutput (:obj:`bool`, optional): if True prints additional
                    information usefull for debuging. Set to false by default.


            Returns:
                float: scalar polarizability in units of MHz cm :math:`^2` / V \
                :math:`^2`
        """
        if (self.drivingFromState[0]!=0):
            raise Exception("Program can only find Polarizability of the original "+\
            "state if you highlight original state. You can do so by NOT "+\
            "specifying drivingFromState in diagonalise function.")


        eFieldList = self.eFieldList
        yState = self.highlight
        y = self.y

        originalState = self.basisStates[self.indexOfCoupledState]
        n = originalState[0]
        l = originalState[1]
        j = originalState[2]
        energyOfOriginalState = self.atom.getEnergy(n,l,j)*C_e/C_h*1e-9 # in  GHz

        if debugOutput:
            print("finding original state for each electric field value")

        stopFitIndex = 0
        while stopFitIndex<len(eFieldList)-1 and \
            eFieldList[stopFitIndex]<maxField:
            stopFitIndex += 1

        xOriginalState = []
        yOriginalState = []

        for ii in xrange(stopFitIndex):

            maxPortion = 0.
            yval = 0.
            jj=0
            for jj in xrange(len(y[ii])):
                if yState[ii][jj]>maxPortion:
                    maxPortion = yState[ii][jj]
                    yval = y[ii][jj]
            # measure state energy relative to the original state
            if (minStateContribution<maxPortion):
                xOriginalState.append(eFieldList[ii])
                yOriginalState.append(yval-energyOfOriginalState)


        xOriginalState = np.array(xOriginalState)/100. # converts to V/cm
        yOriginalState = np.array(yOriginalState)   # in GHz


        ## in GHz
        uppery = 5.0
        lowery = -5.0


        if debugOutput:
            print("found ",len(xOriginalState))
        if showPlot:
            self.fig, self.ax = plt.subplots(1, 1,figsize=(6.5, 3))
            self.ax.scatter(xOriginalState,yOriginalState,s=2,color="k")

            self.ax.set_xlabel("E field (V/cm)")

            self.ax.set_ylim(lowery,uppery)
            self.ax.set_ylabel(r"Energy/$h$ (GHz)")
            self.ax.set_xlim(xOriginalState[0],\
                                xOriginalState[-1])


        def polarizabilityFit(eField,offset,alpha):
            return offset-0.5*alpha*eField**2

        try:
            popt,pcov = curve_fit(polarizabilityFit,\
                              xOriginalState,\
                              yOriginalState,\
                              [0,0])
        except:
            print("\nERROR: fitting energy levels for extracting polarizability\
                    of the state failed. Please check the range of electric \
                    fields where you are trying to fit polarizability and ensure\
                    that there is only one state with continuous energy change\
                    that has dominant contribution of the initial state.\n\n")
            return 0

        if debugOutput:
            print("Scalar polarizability = ",popt[1]*1.e3," MHz cm^2 / V^2 ")

        y_fit = []
        for val in xOriginalState:
            y_fit.append(polarizabilityFit(val,popt[0],popt[1]))
        y_fit = np.array(y_fit)

        if showPlot:
            self.ax.plot(xOriginalState,y_fit,"r--")
            self.ax.legend(("fitted model function","calculated energy level"),\
                      loc=1,fontsize=10)

            self.ax.set_ylim(min(yOriginalState),max(yOriginalState))

            plt.show()

        self.fitX = xOriginalState
        self.fitY = yOriginalState
        self.fittedCurveY = y_fit

        return popt[1]*1.e3 # returned value is in  MHz cm^2 / V^2



# ================= Level plots, decays, cascades etc =======================

class LevelPlot:
    """
        Single atom level plots and decays

        For an example see `Rydberg energy levels example snippet`_.

        .. _`Rydberg energy levels example snippet`:
            ./Rydberg_atoms_a_primer.html#Rydberg-Atom-Energy-Levels

        Args:
            atom (:obj:`AlkaliAtom`): ={ :obj:`alkali_atom_data.Lithium6`,
                :obj:`alkali_atom_data.Lithium7`,
                :obj:`alkali_atom_data.Sodium`,
                :obj:`alkali_atom_data.Potassium39`,
                :obj:`alkali_atom_data.Potassium40`,
                :obj:`alkali_atom_data.Potassium41`,
                :obj:`alkali_atom_data.Rubidium85`,
                :obj:`alkali_atom_data.Rubidium87`,
                :obj:`alkali_atom_data.Caesium` }
                Alkali atom type whose levels we
                want to examine
    """



    def __init__(self,atomType ):
        self.atom = atomType
        self.nFrom = 0
        self.nTo = 0
        self.lFrom = 0
        self.lTo = 0

        self.listX = [] # list of l
        self.listY = [] # list of energies
        self.levelLabel = []

        self.fig = 0
        self.ax = 0
        self.width =0.2
        self.state1=[0,0,0]
        self.state2 =[0,-1,0]
        self.transitionMatrix = []
        self.populations = []
        self.transitionMatrixWavelength3 = []

        # characterization of the graph
        self.spectraX = []
        self.spectraY = []
        self.spectraLine = []


    def makeLevels(self,nFrom,nTo,lFrom,lTo):
        """
            Constructs energy level diagram in a given range

            Args:
                nFrom (int): minimal principal quantum number of the
                    states we are interested in
                nTo (int): maximal principal quantum number of the
                    states we are interested in
                lFrom (int): minimal orbital angular momentum
                    of the states we are interested in
                lTo (int): maximal orbital angular momentum
                    of the states we are interested in
        """
        #save local copy of the space restrictions
        self.nFrom = nFrom
        self.nTo = nTo
        self.lFrom = lFrom
        self.lTo = lTo

        # find all the levels within this space restrictions
        nFrom = max(nFrom,self.atom.groundStateN)
        while nFrom<=nTo:
            l = lFrom
            while l<=min(lTo,4,nFrom-1):
                if (l>0.5):
                    self.listX.append(l)
                    self.listY.append(self.atom.getEnergy(nFrom,l,l-0.5))
                    self.levelLabel.append([nFrom, l, l-0.5])
                self.listX.append(l)
                self.listY.append(self.atom.getEnergy(nFrom,l,l+0.5))
                self.levelLabel.append([nFrom, l, l+0.5])
                l = l+1
            nFrom += 1
        # if user requested principal quantum nuber below the
        # ground state principal quantum number
        # add those L states that are higher in energy then the ground state
        for state in self.atom.extraLevels:
            if state[1]<=lTo and state[0]>=self.nFrom:
                self.listX.append(state[1])
                self.listY.append(self.atom.getEnergy(state[0],state[1],state[2]))
                self.levelLabel.append(state)

    def makeTransitionMatrix(self,environmentTemperature = 0.0,printDecays=True):
        self.transitionMatrix =[]

        for i in xrange(len(self.levelLabel)):
            state1 = self.levelLabel[i]
            transitionVector = []

            # decay of the stay
            decay = 0.0

            for state2 in self.levelLabel:
                dipoleAllowed = (abs(state1[1]-state2[1])==1)and\
                                (abs(state1[2]-state2[2])<=1.01)
                if (dipoleAllowed):
                    # decay to this state
                    rate = self.atom.getTransitionRate(state2[0],state2[1],state2[2],\
                                                    state1[0],state1[1],state1[2],\
                                                    temperature=environmentTemperature)


                    transitionVector.append(rate)

                    # decay from this state
                    rate = self.atom.getTransitionRate(state1[0],state1[1],state1[2],\
                                                    state2[0],state2[1],state2[2],\
                                                    temperature=environmentTemperature)

                    decay = decay-rate
                else:
                    transitionVector.append(0.0)

            transitionVector[i] = decay
            if printDecays:
                print("Decay time of ")
                printState(state1[0], state1[1], state1[2])
                if decay < -1e-20:
                    print("\t is\t",-1.e9/decay," ns")
            self.transitionMatrix.append(transitionVector)

        np.array(self.transitionMatrix)

        self.transitionMatrix = np.transpose(self.transitionMatrix)

    def drawSpectra(self):
        self.fig, self.ax = plt.subplots(1, 1,figsize=(16, 5))

        lineWavelength = []
        lineStrength = []
        lineName = []
        i = 0
        while i<len(self.levelLabel):
            j = 0
            while j<len(self.levelLabel):
                if (i!=j):
                    wavelength = self.atom.getTransitionWavelength(\
                                     self.levelLabel[i][0],\
                                     self.levelLabel[i][1],self.levelLabel[i][2],
                                    self.levelLabel[j][0],\
                                    self.levelLabel[j][1],self.levelLabel[j][2])

                    intensity = self.atom.getTransitionRate(self.levelLabel[i][0],\
                                     self.levelLabel[i][1],self.levelLabel[i][2],\
                                     self.levelLabel[j][0],\
                                     self.levelLabel[j][1],self.levelLabel[j][2])

                    lineWavelength.append(abs(wavelength)*1.e9)
                    lineStrength.append(abs(intensity))
                    lineName.append(printStateString(self.levelLabel[i][0],\
                                                         self.levelLabel[i][1],\
                                                         self.levelLabel[i][2])+\
                                        " -> "+
                                        printStateString(self.levelLabel[j][0],\
                                                         self.levelLabel[j][1],\
                                                         self.levelLabel[j][2]))

                j = j+1
            i = i+1

        self.spectraX = np.copy(lineWavelength)
        self.spectraY = np.copy(lineStrength)
        self.spectraLine = np.copy(lineName)

    def drawSpectraConvoluted(self,lowerWavelength, higherWavelength,points,gamma):
        wavelengths = linspace(lowerWavelength,higherWavelength,points)
        spectra = np.zeros(points)
        i = 0
        while i<len(wavelengths):
            value = 0
            j = 0
            while j<len(self.spectraX):
                value = value + self.spectraY[j]*gamma/\
                                ((self.spectraX[j]-wavelengths[i])**2+gamma**2)
                j = j+1
            spectra[i] = value
            i = i+1
        self.ax.plot(wavelengths,spectra,"g-")

    def showSpectra(self,saveInFile="",showTransitionPoints=True):
        if showTransitionPoints:
            self.ax.plot(self.spectraX,self.spectraY,"ro",picker=5)
        self.ax.set_xlabel("Wavelength (nm)")
        self.ax.set_ylabel("Intensity (arb.un)")
        self.fig.subplots_adjust(right=0.95,left=0.1)
        #self.ax.set_xlim(300,600)
        self.fig.canvas.mpl_connect('pick_event', self.onpick3)
        if (saveInFile != ""):
            self.fig.savefig(saveInFile)
        plt.show()


    def drawLevels(self):
        """
            Draws a level diagram plot
        """
        self.fig, self.ax = plt.subplots(1, 1,figsize=(9.0, 11.5))

        i = 0
        while i<len(self.listX):
            self.ax.plot([self.listX[i]-self.width,self.listX[i]+self.width], \
                         [self.listY[i],self.listY[i]],"b-",picker=4)
            if (i<len(self.populations) and (self.populations[i]>1e-3)):
                self.ax.plot([self.listX[i]],[self.listY[i]],"ro",alpha=self.populations[i])

            i = i+1


    def showPlot(self):
        """
            Shows a level diagram plot
        """
        self.ax.set_ylabel("Energy (eV)")
        self.ax.set_xlim(-0.5+self.lFrom,self.lTo+0.5)

        # X AXIS
        majorLocator   = MultipleLocator(1)

        self.ax.xaxis.set_major_locator(majorLocator)
        tickNames = [" "]
        for l in xrange(self.lFrom,self.lTo+1):
            tickNames.append(printStateLetter(l))
        tickNum = len(self.ax.get_xticklabels())

        self.fig.canvas.draw()
        self.ax.set_xticklabels(tickNames)
        self.fig.canvas.mpl_connect('pick_event', self.onpick2)
        plt.show()

    def findState(self,x,y):
        distance = 100000000.0
        state=[0,0,0]
        i = 0
        while i<len(self.listX):
            dx = self.listX[i]-x
            dy = self.listY[i]-y
            dist = sqrt(dx*dx+dy*dy)
            if (dist<distance):
                distance = dist
                state = self.levelLabel[i]
            i = i+1
        return state

    def findStateNo(self,state):
        # returns no of the given state in the basis
        i = 0
        while i<len(self.levelLabel):
            if (self.levelLabel[i][0] == state[0])and\
                (self.levelLabel[i][1] == state[1])and\
                (abs(self.levelLabel[i][2] - state[2])<0.01):
                return i
            i = i+1

        print("Error: requested state ")
        print(state)
        print("could not be found!")
        return -1

    def findLine(self,x,y):
        distance = 1.e19
        line=""
        i = 0
        while i<len(self.spectraLine):
            dx = self.spectraX[i]-x
            dy = self.spectraY[i]-y
            dist = sqrt(dx*dx+dy*dy)
            if (dist<distance):
                distance = dist
                line = self.spectraLine[i]
            i = i+1
        return line

    def onpick2(self,event):
        if isinstance(event.artist, matplotlib.lines.Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()

            state = self.findState((xdata[0]+xdata[0])/2., ydata[0])
            if (self.state1[0]==0 or (state[1]== self.state2[1])):
                self.state1 = state
                self.ax.set_title(printStateString(state[0],state[1],state[2])+" -> ")
                self.state2=[-1,-1,-1]
            else:
                title = ""
                if (state[1] != self.state1[1]) and (state[1]!= self.state2[1]):
                    title = printStateString(self.state1[0],\
                                             self.state1[1],\
                                             self.state1[2])+\
                            " -> "+\
                            printStateString(state[0],state[1],state[2])+" "
                    title = title+(" %.2f nm (%.3f GHz)" % \
                                   (self.atom.getTransitionWavelength(self.state1[0],\
                                                                      self.state1[1],\
                                                                      self.state1[2],\
                                                                      state[0],state[1],\
                                                                      state[2])*1e9,\
                                    self.atom.getTransitionFrequency(self.state1[0],\
                                                                     self.state1[1],\
                                                                     self.state1[2],\
                                                                     state[0],\
                                                                     state[1],\
                                                                     state[2])*1e-9))
                    self.ax.set_title(title)
                    self.state1=[0,0,0]

                self.state2[1] = state[1]
            event.canvas.draw()

    def onpick3(self,event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            print(ind[0])

            line = self.findLine(xdata[ind][0], ydata[ind][0])
            self.ax.set_title(line)
            event.canvas.draw()


def printState(n,l,j):
    """
        Prints state spectroscopic label for numeric :math:`n`,
        :math:`l`, :math:`s` label of the state

        Args:
            n (int): principal quantum number
            l (int): orbital angular momentum
            j (float): total angular momentum
    """

    print(n," ",printStateLetter(l),(" %.0d/2" % (j*2)))

def printStateString(n,l,j):
    """
        Returns state spectroscopic label for numeric :math:`n`,
        :math:`l`, :math:`s` label of the state

        Args:
            n (int): principal quantum number
            l (int): orbital angular momentum
            j (float): total angular momentum

        Returns:
            string: label for the state in standard spectroscopic notation
    """

    return str(n)+" "+printStateLetter(l)+(" %.0d/2" % (j*2))
