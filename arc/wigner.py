# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from scipy import floor, sqrt
from scipy.misc import factorial
from numpy import arange
import numpy as np
import os
from scipy.sparse import csr_matrix
import sys
if sys.version_info > (2,):
    xrange = range

wignerPrecal = True  # use precalculated values - tested only for the main algorithm calls
wignerPrecalJmax = 23
wignerPrecal3j = np.load(os.path.join(os.path.dirname(os.path.realpath(__file__)),\
                                      "data","precalculated3j.npy"), encoding = 'latin1')

wignerPrecal6j = np.load(os.path.join(os.path.dirname(os.path.realpath(__file__)),\
                         "data","precalculated6j.npy"), encoding = 'latin1')

def Wigner3j(j1,j2,j3,m1,m2,m3):

    # use precalculated values
    if wignerPrecal and ((j2<2.1) and abs(m2)<2.1 and (j1<wignerPrecalJmax)):
        # we shoud have precalculated value
        if ((abs(j1-j2)-0.1<j3) and (j3<j1+j2+0.1) and abs(m1+m2+m3)<0.1):
            # return precalculated value
            return wignerPrecal3j[int(round(2*j1)),int(round(2*(wignerPrecalJmax+m1))),\
                                  int(round(2.*j2)),int(round(m2+j2)),int(round(2-j3+j1,0))]
        else:
            # that value is 0
            return 0

    ##print "unknown %.1f %.1f %.1f %.1f %.1f %.1f " % (j1,j2,j3,m1,m2,m3)
#======================================================================
# Wigner3j.m by David Terr, Raytheon, 6-17-04
#
# Compute the Wigner 3j symbol using the Racah formula [1].
#
# Usage:
# from wigner import Wigner3j
# wigner = Wigner3j(j1,j2,j3,m1,m2,m3)
#
#  / j1 j2 j3 \
#  |          |
#  \ m1 m2 m3 /
#
# Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld:
# http://mathworld.wolfram.com/Wigner3j-Symbol.html
#======================================================================

    # Error checking
    if ( ( 2*j1 != floor(2*j1) ) | ( 2*j2 != floor(2*j2) ) | \
         ( 2*j3 != floor(2*j3) ) | ( 2*m1 != floor(2*m1) ) | \
         ( 2*m2 != floor(2*m2) ) | ( 2*m3 != floor(2*m3) ) ):
        print('All arguments must be integers or half-integers.')
        return -1

    # Additional check if the sum of the second row equals zero
    if ( m1+m2+m3 != 0 ):
        #print('3j-Symbol unphysical')
        return 0

    if ( j1 - m1 != floor ( j1 - m1 ) ):
        print('2*j1 and 2*m1 must have the same parity')
        return 0

    if ( j2 - m2 != floor ( j2 - m2 ) ):
        print('2*j2 and 2*m2 must have the same parity')
        return; 0

    if ( j3 - m3 != floor ( j3 - m3 ) ):
        print('2*j3 and 2*m3 must have the same parity')
        return 0

    if ( j3 > j1 + j2)  | ( j3 < abs(j1 - j2) ):
        print('j3 is out of bounds.')
        return 0

    if abs(m1) > j1:
        print('m1 is out of bounds.')
        return 0

    if abs(m2) > j2:
        print('m2 is out of bounds.')
        return 0

    if abs(m3) > j3:
        print('m3 is out of bounds.')
        return 0

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max( 0, max( t1, t2 ) )
    tmax = min( t3, min( t4, t5 ) )
    tvec = arange(tmin, tmax+1, 1)

    wigner = 0

    for t in tvec:
        wigner += (-1)**t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) *\
                              factorial(t3-t) * factorial(t4-t) * factorial(t5-t) )

    return wigner * (-1)**(j1-j2-m3) * sqrt( factorial(j1+j2-j3) *\
                                              factorial(j1-j2+j3) * \
                                              factorial(-j1+j2+j3) / \
                                              factorial(j1+j2+j3+1) *\
                                               factorial(j1+m1) * \
                                               factorial(j1-m1) * \
                                               factorial(j2+m2) * \
                                               factorial(j2-m2) *\
                                                factorial(j3+m3) *\
                                                 factorial(j3-m3) )


def Wigner6j(j1,j2,j3,J1,J2,J3):
    # if possible, use precalculated values
    global wignerPrecal
    if wignerPrecal and ((round(2*j2,0)==1) and (J2==1 or J2==2)and \
                         (j1<wignerPrecalJmax) and (J3<wignerPrecalJmax)):
        # we have precalculated value
        return wignerPrecal6j[j1,J3,int(round(0.5 +j1-j3)),\
                              int(round(0.5+J3-J1)),J2-1]
    ##print("not in database %1.f %1.f %1.f %1.f %1.f %1.f" % (j1,j2,j3,J1,J2,J3))

#======================================================================
# Calculating the Wigner6j-Symbols using the Racah-Formula
# Author: Ulrich Krohn
# Date: 13th November 2009
#
# Based upon Wigner3j.m from David Terr, Raytheon
# Reference: http://mathworld.wolfram.com/Wigner6j-Symbol.html
#
# Usage:
# from wigner import Wigner6j
# WignerReturn = Wigner6j(j1,j2,j3,J1,J2,J3)
#
#  / j1 j2 j3 \
# <            >
#  \ J1 J2 J3 /
#
#======================================================================

    # Check that the js and Js are only integer or half integer
    if ( ( 2*j1 != round(2*j1) ) | ( 2*j2 != round(2*j2) ) | ( 2*j2 != round(2*j2) ) | ( 2*J1 != round(2*J1) ) | ( 2*J2 != round(2*J2) ) | ( 2*J3 != round(2*J3) ) ):
        print('All arguments must be integers or half-integers.')
        return -1

# Check if the 4 triads ( (j1 j2 j3), (j1 J2 J3), (J1 j2 J3), (J1 J2 j3) ) satisfy the triangular inequalities
    if ( ( abs(j1-j2) > j3 ) | ( j1+j2 < j3 ) | ( abs(j1-J2) > J3 ) | ( j1+J2 < J3 ) | ( abs(J1-j2) > J3 ) | ( J1+j2 < J3 ) | ( abs(J1-J2) > j3 ) | ( J1+J2 < j3 ) ):
        print('6j-Symbol is not triangular!')
        return 0

    # Check if the sum of the elements of each traid is an integer
    if ( ( 2*(j1+j2+j3) != round(2*(j1+j2+j3)) ) | ( 2*(j1+J2+J3) != round(2*(j1+J2+J3)) ) | ( 2*(J1+j2+J3) != round(2*(J1+j2+J3)) ) | ( 2*(J1+J2+j3) != round(2*(J1+J2+j3)) ) ):
        print('6j-Symbol is not triangular!')
        return 0

    # Arguments for the factorials
    t1 = j1+j2+j3
    t2 = j1+J2+J3
    t3 = J1+j2+J3
    t4 = J1+J2+j3
    t5 = j1+j2+J1+J2
    t6 = j2+j3+J2+J3
    t7 = j1+j3+J1+J3

    # Finding summation borders
    tmin = max(0, max(t1, max(t2, max(t3,t4))))
    tmax = min(t5, min(t6,t7))
    tvec = arange(tmin,tmax+1,1)

    # Calculation the sum part of the 6j-Symbol
    WignerReturn = 0
    for t in tvec:
        WignerReturn += (-1)**t*factorial(t+1)/( factorial(t-t1)*factorial(t-t2)*factorial(t-t3)*factorial(t-t4)*factorial(t5-t)*factorial(t6-t)*factorial(t7-t) )

    # Calculation of the 6j-Symbol
    return WignerReturn*sqrt( TriaCoeff(j1,j2,j3)*TriaCoeff(j1,J2,J3)*TriaCoeff(J1,j2,J3)*TriaCoeff(J1,J2,j3) )


def TriaCoeff(a,b,c):
    # Calculating the triangle coefficient
    return factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/(factorial(a+b+c+1))

# copied from https://sites.google.com/site/theodoregoetz/notes/wignerdfunction
# Jojann Goetz

from scipy.special import jv, legendre, sph_harm, jacobi
from scipy.misc import  comb
from numpy import floor, sqrt, sin, cos, exp, power
from numpy import conj as conjugate
from math import pi

def wignerd(j,m,n=0,approx_lim=10):
    '''
        Wigner "small d" matrix. (Euler z-y-z convention)
        example:
            j = 2
            m = 1
            n = 0
            beta = linspace(0,pi,100)
            wd210 = wignerd(j,m,n)(beta)

        some conditions have to be met:
             j >= 0
            -j <= m <= j
            -j <= n <= j

        The approx_lim determines at what point
        bessel functions are used. Default is when:
            j > m+10
              and
            j > n+10

        for integer l and n=0, we can use the spherical harmonics. If in
        addition m=0, we can use the ordinary legendre polynomials.
    '''

    if (j < 0) or (abs(m) > j) or (abs(n) > j):
        raise ValueError("wignerd(j = {0}, m = {1}, n = {2}) value error.".format(j,m,n) \
            + " Valid range for parameters: j>=0, -j<=m,n<=j.")

    if (j > (m + approx_lim)) and (j > (n + approx_lim)):
        #print('bessel (approximation)')
        return lambda beta: jv(m-n, j*beta)

    if (floor(j) == j) and (n == 0):
        if m == 0:
            #print('legendre (exact)')
            return lambda beta: legendre(j)(cos(beta))
        elif False:
            #print('spherical harmonics (exact)')
            a = sqrt(4.*pi / (2.*j + 1.))
            return lambda beta: a * conjugate(sph_harm(m,j,beta,0.))

    jmn_terms = {
        j+n : (m-n,m-n),
        j-n : (n-m,0.),
        j+m : (n-m,0.),
        j-m : (m-n,m-n),
        }

    k = min(jmn_terms)
    a, lmb = jmn_terms[k]

    b = 2.*j - 2.*k - a

    if (a < 0) or (b < 0):
        raise ValueError("wignerd(j = {0}, m = {1}, n = {2}) value error.".format(j,m,n) \
            + " Encountered negative values in (a,b) = ({0},{1})".format(a,b))

    coeff = power(-1.,lmb) * sqrt(comb(2.*j-k,k+a)) * (1./sqrt(comb(k+b,b)))

    #print('jacobi (exact)')
    return lambda beta: coeff \
        * power(sin(0.5*beta),a) \
        * power(cos(0.5*beta),b) \
        * jacobi(k,a,b)(cos(beta))


def wignerD(j,m,n=0,approx_lim=10):
    '''
        Wigner D-function. (Euler z-y-z convention)

        This returns a function of 2 to 3 Euler angles:
            (alpha, beta, gamma)

        gamma defaults to zero and does not need to be
        specified.

        The approx_lim determines at what point
        bessel functions are used. Default is when:
            j > m+10
              and
            j > n+10

        usage:
            from numpy import linspace, meshgrid
            a = linspace(0, 2*pi, 100)
            b = linspace(0,   pi, 100)
            aa,bb = meshgrid(a,b)
            j,m,n = 1,1,1
            zz = wignerD(j,m,n)(aa,bb)
    '''

    return lambda alpha,beta,gamma=0: \
          exp(-1j*m*alpha) \
        * wignerd(j,m,n,approx_lim)(beta) \
        * exp(-1j*n*gamma)

def CG(j1,m1,j2,m2,j3,m3):
    """
        returns < j1,m1,j2,m2 | j1,j2, j3, m3 >
    """
    return Wigner3j(j1,j2,j3,m1,m2,-m3)*sqrt(2*j3+1)*(-1)**(j1-j2+m3)


class wignerDmatrix:

    def __init__(self,theta,phi,gamma=0.):
        self.matSaved = []
        self.matLoc = np.zeros(100,dtype=np.int8)
        self.theta = theta
        self.phi = phi
        self.gamma = gamma
        if (abs(self.theta)<1e-5 and abs(self.phi)<1e-5 and\
             abs(self.gamma)<1e-5):
            self.trivial = True
        else:
            self.trivial = False

    def get(self,j):
        if self.trivial:
            return np.eye(int(round(2.*j+1.,0)),int(round(2.*j+1.,0)), dtype=np.complex128)
        savedIndex = self.matLoc[int(round(2*j,0))]
        if savedIndex != 0:
            return self.matSaved[savedIndex-1]
            # bacause 0 marks no entry; but matrix numbers starts from zero,
            # saved Index array is actually offsetted by 1
        # else
        mat = np.zeros((int(round(2.*j+1.,0)),int(round(2.*j+1.,0))),dtype=np.complex128)
        jrange = np.linspace(-j,j,int(2*j)+1)
        maxIndex = int(2*j)+1

        for index1 in xrange(maxIndex):
            for index2 in xrange(maxIndex):
                mat[index1,index2] = wignerD(j,jrange[index1],jrange[index2])\
                                            (self.phi,self.theta,self.gamma)

        mat = csr_matrix(mat)
        self.matSaved.append(mat)
        self.matLoc[int(round(2*j,0))] = len(self.matSaved)
        return mat
