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
        pass;
