# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
from .alkali_atom_functions import *


class AlkalineEarthAtom(AlkaliAtom):

    modelPotential_coef = dict()
    """
        Model potential parameters fitted from experimental observations for
        different l (electron angular momentum)
    """


    def __init__(self,preferQuantumDefects=True,cpp_numerov=True):
        pass;
