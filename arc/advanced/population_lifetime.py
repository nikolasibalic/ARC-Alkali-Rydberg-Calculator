# -*- coding: utf-8 -*-
from scipy.integrate import odeint
from lmfit import minimize, Parameters, report_fit
from ..alkali_atom_data import *
import matplotlib.pyplot as plt

"""
    **Contributors:**
    getPopulationLifetime - written by Alessandro Greco,
    Dipartimento di Fisica *E. Fermi*, Università di Pisa,
    Largo Bruno Pontecorvo 3, 56127 Pisa, Italy (alessandrogreco08 at gmail dot com),
    the simulations have been compared with experimental data [#greco2019]_
"""


def getPopulationLifetime(atom, n, l, j,
                          temperature=0, includeLevelsUpTo=0, period=1,
                          plotting=1, thresholdState=False, detailedOutput=False):
    r"""
    Calculates lifetime of atomic **population** taking into account
    redistribution of population to other states under spontaneous and
    black body induced transitions.

    It simulates the time evolution of a system in which all the states,
    from the fundamental one to the highest state which you want to include,
    are taken into account.
    The orbital angular momenta taken into account are only S,P,D,F.

    This function is based on getStateLifetime but it takes into account
    the re-population processess due to BBR-induced transitions.
    For this reason lifetimes of Rydberg states are slightly longer
    than those returned by getStateLifetime up to 5-10%.

    This function creates a .txt file, plots the time evolution of the
    population of the Rydberg states and yields the lifetime values by using
    the fitting method from Ref. [#fit]_ .

    **Contributed by:** Alessandro Greco (alessandrogreco08 at gmail dot com),
    Dipartimento di Fisica *E. Fermi*, Università di Pisa, Largo Bruno Pontecorvo 3, 56127 Pisa, Italy.
    The simulations have been compared with experimental data [#greco2019]_ .

    **Please cite as:** `original ARC paper`_ and paper introducing
    extension [#greco2019]_

    .. _`original ARC paper`:
        https://doi.org/10.1016/j.cpc.2017.06.015

    References:

    .. [#fit] https://people.duke.edu/~ccc14/sta-663/CalibratingODEs.html

    .. [#greco2019] M. Archimi, C. Simonelli, L. Di Virgilio, A. Greco,
        M. Ceccanti, E. Arimondo, D. Ciampini, I. I. Ryabtsev, I. I. Beterov,
        and O. Morsch, *Phys. Rev. A* **100**, 030501(R) (2019)
        https://doi.org/10.1103/PhysRevA.100.030501

    **Some definitions:**

    What are the **ensemble**, the **support**, the **ground**?
    According to https://arxiv.org/abs/1907.01254
    The sum of the populations of every state which is detected as Rydberg state
    (above the threshold state which must be set\) is called **ensemble**
    The sum of the populations of every state which is detected as Rydberg state,
    without the target state, is called **support**
    The sum of the populations of every state which cannot be detected as Rydberg state
    (under the threshold state which must be set) is called **ground**
    **gammaTargetSpont** is the rate which describes transitions
    from Target State towards all the levels under the threshold
    state, i.e. Ground State
    **gammaTargetBBR** is the rate which describes transitions towards all
    the levels above the threshold state, i.e. Support State
    **gammaSupporSpont** is the rate which describes transitions from the
    Support State towards all the levels under the threshold state,
    i.e. Ground State
    **gammaSupportBBR** is the rate which describes transitions from
    upport State towards all the levels above the threshold state,
    i.e. Target State)

    Args:
        n (int): principal quantum number of the state whose population lifetime
            we are calculating, it's called the *Target* state and its color is
            green in the plot.
        l (int): orbital angular momentum number of the state whose population lifetime
            we are calculating, it's called the *Target* state and its color is
            green in the plot.
        j (float): total angular momentum of the state whose population lifetime
            we are calculating, it's called the *Target* state and its color is
            green in the plot.
        temperature (float): Temperature at which the atom environment
            is, measured in K. If this parameter is non-zero, user has
            to specify transitions up to which state (due to black-body
            decay) should be included in calculation.
        includeLevelsUpTo (int): At non zero temperatures,
            this specifies maximum principal quantum number of the state
            to which black-body induced transitions will be included.
            Minimal value of the parameter in that case is =`n+1
        period: Specifies the period that you want to consider for
            the time evolution, in microseconds.
        plotting (int): optional. It is set to 1 by default. The options are
            (see also image at the bottom of documentation):
            **plotting=0** no plot;
            **plotting=1** plots the population of the target (n,l,j) state
            with its fit and it yields the value of the target lifetime
            in microseconds;
            **plotting=2** plots the whole system (Ensemble, Support, Target),
            no fit;
            **plotting=3** plots the whole system (Ensemble, Support, Target)
            and it fits the Ensemble and Target curves, it yields the values
            of the Ensemble lifetime and Target lifetime in microseconds;
            **plotting=4** it plots the whole system (Ensemble, Support, Target) +
            the Ground (which is the complementary of the ensemble).
            It considers the whole system like a three-level model (Ground
            State, Support State, Target State) and yields four transition
            rates.

        thresholdState (int): optional. It specifies the principal quantum
            number n of the lowest state (it's referred to S state!) which is
            detectable by your experimental apparatus, it directly modifies
            the *Ensemble* and the *Support* (whose colors are red and blue
            respectively in the plot). It is necessary to define a threshold
            state if plotting = 2, 3 or 4 has been selected. It is not necessary
            to define a threshold state if plotting = 0 or 1 has been selected.

        detailedOutput=True: optional. It writes a .txt file with the time
            evolution of all the states. It is set to false by default.
            (The first column is the time, the other are the population of all
            the states. The order is time, nS, nP0.5, nP1.5, nD1.5, nD2.5,
            nF2.5, nF3.5, and n is ordered from the lowest state to the highest one.
            For example: time, 4S, 5S ,6S ,ecc... includeLevelsUpToS, 4P0.5,
            5P0.5, 6P0.5, ecc... includeLevelsUpToP0.5, 4P1.5, 5P1.5, 6P1.5, ecc...)

    Returns:
        Plots and a .txt file.
        **plotting = 0,1** create a .txt file with two coloumns
        (time \t target population);
        **plotting = 2,3,4** create a .txt file with four coloumns
        (time \t ensemble population \t support population \t target population)


    Example:
        >>> from arc import *
        >>> from arc.advanced.population_lifetime import getPopulationLifetime
        >>> atom = Rubidium()
        >>> getPopulationLifetime(atom, 10, 1, 1.5, temperature =300,
                includeLevelsUpTo=15, detailedOutput=True, plotting=1)

    """

    if l > 3:
        print("Error: this function takes into account only S, P, D, F states.")
        return

    if plotting > 4:
        print("Error: plotting must be equal to 0, 1, 2, 3 or 4.")
        return

    if ((thresholdState == False) and  (plotting >1 )) or (thresholdState==True):
        print("Error: you need to specify the principal quantum number of the "
              "thresholdState if you use plotting=2, 3 or 4.")
        return

    if ((plotting == 0) or  (plotting ==1)):
        thresholdState = False

    import time
    start = time.time()

    # What state do you want to excite?
    STATE = n
    L = l
    J = j

    # Which states do you want to consider for the BBR width?
    if includeLevelsUpTo - STATE < 0:
        raise valueError("Error: includeLevelsUpTo must be >= n")
    WidthBBR = includeLevelsUpTo - STATE
    # What is the temperature?
    if temperature == 0:
        raise valueError("Error: if you don't want BBR-induced transition, use getStateLifetime")
    TEMP_BBR = temperature
    # What is the critical state for the ionization?
    if thresholdState - STATE >= 0:
        raise valueError("Error: thresholdState must be < n")
    CState = thresholdState
    # It creates the references for the ensemble population
    cutoffs = int(atom.getQuantumDefect(STATE, 0, 0.5) -
                  atom.getQuantumDefect(STATE,0,0.5))
    cutoffp05 = int(atom.getQuantumDefect(STATE, 0, 0.5 ) -
                    atom.getQuantumDefect(STATE,1,0.5))
    cutoffp15 = int(atom.getQuantumDefect(STATE, 0, 0.5) -
                    atom.getQuantumDefect(STATE,1,1.5))
    cutoffd15 = int(atom.getQuantumDefect(STATE, 0, 0.5) -
                    atom.getQuantumDefect(STATE,2,1.5))
    cutoffd25 = int(atom.getQuantumDefect(STATE, 0, 0.5) -
                    atom.getQuantumDefect(STATE,2,2.5))
    cutofff25 = int(atom.getQuantumDefect(STATE, 0, 0.5) -
                    atom.getQuantumDefect(STATE,3,2.5))
    cutofff35 = int(atom.getQuantumDefect(STATE, 0, 0.5) -
                    atom.getQuantumDefect(STATE,3,3.5))
    # Total time of the dynamics
    totaltime = period * 1e-6
    # Parts of gammamax that you take for time step
    partg = 2.0

    #########################################################

    # It takes into account of the extra levels
    extraL = atom.extraLevels[1][:]

    # It creates the references for the matrix
    riftot = (STATE + WidthBBR - extraL[0] + 1) * 7
    rifs = ((STATE + WidthBBR - extraL[0] + 1) * 0) - extraL[0]
    rifp05 = ((STATE + WidthBBR - extraL[0] + 1) * 1) - extraL[0]
    rifp15 = ((STATE + WidthBBR - extraL[0] + 1) * 2) - extraL[0]
    rifd15 = ((STATE + WidthBBR - extraL[0] + 1) * 3) - extraL[0]
    rifd25 = ((STATE + WidthBBR - extraL[0] + 1) * 4) - extraL[0]
    riff25 = ((STATE + WidthBBR - extraL[0] + 1) * 5) - extraL[0]
    riff35 = ((STATE + WidthBBR - extraL[0] + 1) * 6) - extraL[0]

    # It creates the matrix of the rates
    c = np.zeros(shape=(riftot, riftot))

    print("Creating the rates matrix:")

    for pqn in xrange(extraL[0], STATE + WidthBBR+1):
        for fpqn in xrange(extraL[0], STATE + WidthBBR+1):
            # rate from s
            c[pqn + rifs, fpqn + rifp05] = atom.getTransitionRate(
                pqn, 0, 0.5, fpqn, 1, 0.5, TEMP_BBR)  # rate s -> p0.5
            c[pqn + rifs, fpqn + rifp15] = atom.getTransitionRate(
                pqn, 0, 0.5, fpqn, 1, 1.5, TEMP_BBR)  # rate s -> p1.5
        # rate from p0.5
            c[pqn + rifp05, fpqn + rifs] = atom.getTransitionRate(
                pqn, 1, 0.5, fpqn, 0, 0.5, TEMP_BBR)  # rate p0.5 -> s
            c[pqn + rifp05, fpqn + rifd15] = atom.getTransitionRate(
                pqn, 1, 0.5, fpqn, 2, 1.5, TEMP_BBR)  # rate p0.5 -> d1.5
        # rate from p1.5
            c[pqn + rifp15, fpqn + rifs] = atom.getTransitionRate(
                pqn, 1, 1.5, fpqn, 0, 0.5, TEMP_BBR)  # rate p1.5 -> s
            c[pqn + rifp15, fpqn + rifd15] = atom.getTransitionRate(
                pqn, 1, 1.5, fpqn, 2, 1.5, TEMP_BBR)  # rate p1.5 -> d1.5
            c[pqn + rifp15, fpqn + rifd25] = atom.getTransitionRate(
                pqn, 1, 1.5, fpqn, 2, 2.5, TEMP_BBR)  # rate p1.5 -> d2.5
        # rate from d1.5
            c[pqn + rifd15, fpqn + rifp05] = atom.getTransitionRate(
                pqn, 2, 1.5, fpqn, 1, 0.5, TEMP_BBR)  # rate d1.5 -> p0.5
            c[pqn + rifd15, fpqn + rifp15] = atom.getTransitionRate(
                pqn, 2, 1.5, fpqn, 1, 1.5, TEMP_BBR)  # rate d1.5 -> p1.5
            c[pqn + rifd15, fpqn + riff25] = atom.getTransitionRate(
                pqn, 2, 1.5, fpqn, 3, 2.5, TEMP_BBR)  # rate d1.5 -> f2.5
        # rate from d2.5
            c[pqn + rifd25, fpqn + rifp15] = atom.getTransitionRate(
                pqn, 2, 2.5, fpqn, 1, 1.5, TEMP_BBR)  # rate d2.5 -> p1.5
            c[pqn + rifd25, fpqn + riff25] = atom.getTransitionRate(
                pqn, 2, 2.5, fpqn, 3, 2.5, TEMP_BBR)  # rate d2.5 -> f2.5
            c[pqn + rifd25, fpqn + riff35] = atom.getTransitionRate(
                pqn, 2, 2.5, fpqn, 3, 3.5, TEMP_BBR)  # rate d2.5 -> f3.5
        # rate from f2.5
            c[pqn + riff25, fpqn + rifd15] = atom.getTransitionRate(
                pqn, 3, 2.5, fpqn, 2, 1.5, TEMP_BBR)  # rate f2.5 -> d1.5
            c[pqn + riff25, fpqn + rifd25] = atom.getTransitionRate(
                pqn, 3, 2.5, fpqn, 2, 2.5, TEMP_BBR)  # rate f2.5 -> d2.5
        # rate from f3.5
            c[pqn + riff35, fpqn + rifd25] = atom.getTransitionRate(
                pqn, 3, 3.5, fpqn, 2, 2.5, TEMP_BBR)  # rate f3.5 -> d2.5
        print(pqn, end=' ')

    # It deletes all the gammas for states under the ground state which are not the extra levels

    if extraL[1] > 2:
        c[extraL[0] +rifd15, :] = 0
        c[:, extraL[0] +rifd15] = 0
        c[extraL[0] +rifd25, :] = 0
        c[:, extraL[0] +rifd25] = 0
        if extraL[1] > 3:
            c[extraL[0] +riff25, :] = 0
            c[:, extraL[0] +riff25] = 0
            c[extraL[0] +riff35, :] = 0
            c[:, extraL[0] +riff35] = 0

    c[extraL[0] +rifs, :] = 0
    c[:, extraL[0] +rifs] = 0
    c[extraL[0] +rifp05, :] = 0
    c[:, extraL[0] +rifp05] = 0
    c[extraL[0] +rifp15, :] = 0
    c[:, extraL[0] +rifp15] = 0

    # It finds the maximum rate in the matrix
    gammamax = c.max()  # is from the 5P1.5 towards the 5S0.5
    # It defines Dtmin
    Dtmin = round(1 /(partg *gammamax), 9)
    print('\n', Dtmin)

    #########################################################

    # It inizialites the population and the auxiliry population vectors
    pop = np.zeros(shape=(1, riftot))
    popaus = np.zeros(shape=(1, riftot))

    # It inizializes the reference for the population vector
    if L == 0:
        rifinitial = rifs
    if L == 1:
        if J == 0.5:
            rifinitial = rifp05
        if J == 1.5:
            rifinitial = rifp15
    if L == 2:
        if J == 1.5:
            rifinitial = rifd15
        if J == 2.5:
            rifinitial = rifd25
    if L == 3:
        if J == 2.5:
            rifinitial = riff25
        if J == 3.5:
            rifinitial = riff35

    pop[0, (rifinitial +STATE)] = 1

    #########################################################

    # It inizializes the time and the time step
    t = 0.0
    Dt = 0.0

    #########################################################
    # References for the name of the .txt file
    if L == 0:
        StrL = 'S'
    elif L == 1:
        StrL = 'P'
    elif L == 2:
        StrL = 'D'
    elif L == 3:
        StrL = 'F'

    if J == 0.5:
        StrJ = '05'
    elif J == 1.5:
        StrJ = '15'
    elif J == 2.5:
        StrJ = '25'
    elif J == 3.5:
        StrJ = '35'

    # It creates the file for the three curves
    with open("Lifetime" + str(STATE) +StrL+StrJ+".txt", 'w') as fi:
        fi.writelines("")

    if detailedOutput == True:
        # It creates the file for the all states
        with open("Lifetime" + str(STATE) +StrL+StrJ+"All.txt", 'w') as fiall:
            fiall.writelines("")
    #########################################################

    # It creates four lists to quickly write the results to the file
    ListTime = []
    if thresholdState != False:
        ListRed = []
        ListBlue = []
    ListGreen = []

    # The core of the program starts
    while t < (totaltime):
        if detailedOutput == True:
            ListStates = []
            ListStates.append(t * 1e+6)
        for a in range(0, riftot):
            popaus[0, a] = 0.0
            for b in range(0, riftot):
                popaus[0, a] += -c[a, b]*pop[0, a] + c[b,a]*pop[0,b]
            popaus[0, a] = popaus[0, a] *Dt
        pop += popaus
        if t == 0:
            Dt = Dtmin
        if detailedOutput == True:
            ListStates.extend(pop[0, :])
            with open("Lifetime" + str(STATE) +StrL+StrJ+"All.txt", 'a') as fall:
                fall.writelines("%.5f \t" % (ListStates[ind]) for ind in range(0, len(ListStates)))
                fall.writelines("\n")
        ListTime.append(t * 1e+6)
        if thresholdState != False:
            popall = 0.0
            for k in range(0, riftot):
                if ((CState + rifs-cutoffs <= k < rifp05+extraL[0])
                    or (CState+rifp05-cutoffp05 <= k < rifp15+extraL[0])
                    or (CState+rifp15-cutoffp15 <= k < rifd15+extraL[0])
                    or (CState+rifd15-cutoffd15 <= k < rifd25+extraL[0])
                    or (CState+rifd25-cutoffd25 <= k < riff25+extraL[0])
                    or (CState+riff25-cutofff25 <= k < riff35+extraL[0])
                    or (CState+riff35-cutofff35 <= k < riftot)):
                    # above the threshold state
                    popall += pop[0, k]
            ListRed.append(popall)
            ListBlue.append(popall - pop[0, (rifinitial + STATE)])
        ListGreen.append(pop[0, (rifinitial + STATE)])
        sys.stdout.write("\rProgress: %d%%" % ((t / totaltime) * 100))
        sys.stdout.flush()
        t = t + Dt

    if thresholdState == False:
        with open("Lifetime" + str(STATE) +StrL+StrJ+".txt", 'a') as f:
            f.writelines("%.4f \t %.5f \n" %
                (ListTime[index], ListGreen[index]) for index in range(0, len(ListTime)))
    else:
        with open("Lifetime" + str(STATE) +StrL+StrJ+".txt", 'a') as f:
            f.writelines("%.4f \t %.5f \t %.5f \t %.5f \n" %
                (ListTime[index], ListRed[index], ListBlue[index],
                 ListGreen[index]) for index in range(0,len(ListTime)))

    #########################################################

    if plotting == 1:

        def f(xs, t, ps):
            """Lotka-Volterra predator-prey model."""
            try:
                gammaTarget = ps['gammaTarget'].value
            except Exception:
                gammaTarget = ps

            x, y = xs
            return [- gammaTarget * x, - gammaTarget *y]

        def g(t, x0, ps):
            """
            Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
            """
            x = odeint(f, x0, t, args=(ps,))
            return x

        def residual(ps, ts, data):
            x0 = ps['x0'].value, ps['y0'].value
            model = g(ts, x0, ps)
            return (model - data).ravel()

        t = np.array(ListTime)
        x0 = np.array([0, 0])

        data = np.zeros(shape=(len(t), 2))

        data[:, 0] = np.array(ListGreen)
        data[:, 1] = np.array(ListGreen)

        # set parameters incluing bounds
        params = Parameters()
        params.add('x0', value=1, vary=False)
        params.add('y0', value=1, vary=False)
        params.add('gammaTarget', value=0.01, min=0, max=1)

        # fit model and find predicted values
        result = minimize(residual, params, args=(
            t, data), method='leastsq')
        final = data + result.residual.reshape(data.shape)

        LifetimeTarget = 1. / (result.params['gammaTarget'].value)

        # Grafico
        fig, axes = plt.subplots(1, 1, figsize=(10, 6))

        axes.plot(t, data[:, 0], 'g*',label=r"Target")
        axes.plot(t, final[:, 0], 'k-', linewidth=2, label=r"Fit Target")

        axes.set_ylim(0, max(ListGreen))
        axes.set_xlim(0, max(ListTime))
        axes.legend(loc=0, fontsize=12)
        axes.set_ylabel("Number of Rydberg atoms", fontsize=12)
        axes.set_xlabel("Time, $\mu s$", fontsize=12)
        axes.grid()
        plt.legend
        plt.show()

        # display fitted statistics
        print("\n")
        report_fit(result)
        print("\n")
        print("Lifetime Target: %.6f us" % (LifetimeTarget))

    if plotting == 2:

        # Make the plot of the three curves
        fig, axes = plt.subplots(1, 1, figsize=(10, 6))

        axes.plot(ListTime, ListRed, 'r.', label=r"Ensemble")
        axes.plot(ListTime, ListBlue, 'b.', label=r"Other")
        axes.plot(ListTime, ListGreen, 'g.', label=r"Target")

        axes.set_ylim(0, 1)
        axes.set_xlim(0, ListTime[-1])
        axes.legend(loc=0, fontsize=12)
        axes.set_ylabel("Number of Rydberg atoms", fontsize=12)
        axes.set_xlabel("Time [$\mu s$]", fontsize=12)
        axes.grid()
        plt.legend
        plt.show()

    if plotting == 3:

        def f(xs, t, ps):
            """Lotka-Volterra predator-prey model."""
            try:
                gammaEnsemble = ps['gammaEnsemble'].value
                gammaTarget = ps['gammaTarget'].value
            except Exception:
                gammaEnsemble, gammaTarget = ps

            x, y = xs
            return [-gammaEnsemble * x, -gammaTarget *y]

        def g(t, x0, ps):
            """
            Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
            """
            x = odeint(f, x0, t, args=(ps,))
            return x

        def residual(ps, ts, data):
            x0 = ps['x0'].value, ps['y0'].value
            model = g(ts, x0, ps)
            return (model - data).ravel()

        t = np.array(ListTime)
        x0 = np.array([0, 0])

        dataAll = np.zeros(shape=(len(t), 3))
        dataAll[:, 0] = np.array(ListRed)
        dataAll[:, 1] = np.array(ListBlue)
        dataAll[:, 2] = np.array(ListGreen)

        data = np.zeros(shape=(len(t), 2))

        data[:, 0] = dataAll[:, 0]
        data[:, 1] = dataAll[:, 2]

        # set parameters incluing bounds
        params = Parameters()
        params.add('x0', value=max(ListRed), vary=False)
        params.add('y0', value=max(ListGreen), vary=False)
        params.add('gammaEnsemble',  value=0.005, min=0., max=1.)
        params.add('gammaTarget', value=0.01, min=0., max=1.)

        # fit model and find predicted values
        result = minimize(residual, params, args=(
            t, data), method='leastsq')
        final = data + result.residual.reshape(data.shape)

        LifetimeEnsemble = 1. / (result.params['gammaEnsemble'].value)
        LifetimeTarget = 1. / (result.params['gammaTarget'].value)

        # Grafico
        fig, axes = plt.subplots(1, 1, figsize=(10, 6))

        axes.plot(t, dataAll[:, 0], 'r*',label=r"Ensemble")
        axes.plot(t, dataAll[:, 1], 'b*',label=r"Support")
        axes.plot(t, dataAll[:, 2], 'g*',label=r"Target")
        axes.plot(t, final[:, 0], 'k-', linewidth=2, label=r"Fit Ensemble")
        axes.plot(t, final[:, 1], 'k-', linewidth=2, label=r"Fit Target")

        axes.set_ylim(0, max(ListRed))
        axes.set_xlim(0, max(ListTime))
        axes.legend(loc=0, fontsize=12)
        axes.set_ylabel("Number of Rydberg atoms", fontsize=12)
        axes.set_xlabel("Time, $\mu s$", fontsize=12)
        axes.grid()
        plt.legend
        plt.show()

        # display fitted statistics
        print("\n")
        report_fit(result)
        print("\n")
        print("Lifetime Ensemble: %.6f us \nLifetime Target: %.6f us"
              % (LifetimeEnsemble, LifetimeTarget))

    if plotting == 4:
        def f(xs, t, ps):
            """Lotka-Volterra predator-prey model."""
            try:
                gammaTargetSpont = ps['gammaTargetSpont'].value
                gammaTargetBBR = ps['gammaTargetBBR'].value
                gammaSupportSpont = ps['gammaSupportSpont'].value
                gammaSupportBBR = ps['gammaSupportBBR'].value
            except Exception:
                gammaTargetSpont, gammaTargetBBR, gammaSupportSpont, gammaSupportBBR = ps

            x, y, z = xs
            return [+gammaTargetSpont * z + gammaSupportSpont *y,
                    -gammaSupportSpont*y - gammaSupportBBR*y + gammaTargetBBR*z,
                    -gammaTargetSpont*z - gammaTargetBBR*z + gammaSupportBBR*y]

        def g(t, x0, ps):
            """
            Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
            """
            x = odeint(f, x0, t, args=(ps,))
            return x

        def residual(ps, ts, data):
            x0 = ps['x0'].value, ps['y0'].value, ps['z0'].value
            model = g(ts, x0, ps)
            return (model - data).ravel()

        ListRedAus = np.zeros(shape=(len(ListRed)))

        for i in range(0, len(ListRed)):
            ListRedAus[i] = max(ListRed) - ListRed[i]

        t = np.array(ListTime)

        data = np.zeros(shape=(len(t), 3))

        data[:, 0] = np.array(ListRedAus)
        data[:, 1] = np.array(ListBlue)
        data[:, 2] = np.array(ListGreen)

        # set parameters incluing bounds
        params = Parameters()
        params.add('x0', value=0, vary=False)
        params.add('y0', value=0, vary=False)
        params.add('z0', value=max(ListGreen), vary=False)
        params.add('gammaTargetSpont', value=0.02, min=0., max=1.)
        params.add('gammaTargetBBR', value=0.02, min=0., max=1.)
        params.add('gammaSupportSpont', value=0.02, min=0., max=1.)
        params.add('gammaSupportBBR', value=0.001, min=0., max=1.)

        # fit model and find predicted values
        result = minimize(residual, params, args=(
            t, data), method='leastsq')
        final = data + result.residual.reshape(data.shape)

        # Grafico
        fig, axes = plt.subplots(1, 1, figsize=(10, 6))

        axes.plot(t, data[:, 0], 'm*',label=r"Ground")
        axes.plot(t, ListRed, 'r*', label=r"Ensemble")
        axes.plot(t, data[:, 1], 'b*',label=r"Support")
        axes.plot(t, data[:, 2], 'g*',label=r"Target")
        axes.plot(t, final[:, 0], 'k-', linewidth=2, label=r"Fit Ground")
        axes.plot(t, final[:, 1], 'k-', linewidth=2, label=r"Fit Support")
        axes.plot(t, final[:, 2], 'k-', linewidth=2, label=r"Fit Target")

        axes.set_ylim(0, max(ListRed))
        axes.set_xlim(0, max(ListTime))
        axes.legend(loc=0, fontsize=12)
        axes.set_ylabel("Number of Rydberg atoms", fontsize=12)
        axes.set_xlabel("Time, $\mu s$", fontsize=12)
        axes.grid()
        plt.legend()
        plt.show()

        print("\n")
        # display fitted statistics
        report_fit(result)

    # It returns the time elapsed
    print('\nIt took', time.time() - start, 'seconds.')

    return
