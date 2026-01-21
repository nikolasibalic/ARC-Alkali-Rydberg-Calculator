from arc import StarkMap, Rubidium
import numpy as np


def testStark_Rb():
    calc = StarkMap(Rubidium())
    n = 75
    minEfield = 0.003e2
    maxEfield = 0.2e2
    calc.defineBasis(n, 0, 0.5, 0.5, n - 5, n + 5, 20)
    calc.diagonalise(np.linspace(minEfield, maxEfield, 100))
    assert (
        abs(calc.getPolarizability(minStateContribution=0.9) - 8.532e02) < 0.1
    )
