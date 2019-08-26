# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:24:34 2019

@author: DaniJ
"""
"""
    Trying the example of 4.1 of "Physics and Chemistry of Interfaces" (2003), Jürgen Butt, Karlheinz Graf, Michael Kappl.
"""

import numpy as np

def inv_Debye_length (epsilon, T, Caq,Z):
    '''
    epsilon is te relative permitivity
    T is the temperture in Kelvins
    Caq is the vector of aqueous concentrations in mol
    Z is the vector of ionic charge
    Note: Z and Caq should be coherent in order
    '''
    e = 1.60218e-19 # C that is the elementary charge units are Coulombs
    kb = 1.38066e-23    # Boltzmann constant units are J/K
    epsilon0 =  8.854187871e-12  # vacuum permitivity units (A s)/(V m)
    Navo = 6.022141e23 
    
   # C = Navo*Caq
    C = Caq
    CZ = np.multiply(C,np.multiply(Z,Z))
    sumCZ = np.sum(CZ)
    A = e*e*sumCZ
    B = epsilon0*epsilon*T*kb
    k = np.sqrt(A/B)
    return k

Caq = np.array([861e23, 15e23,620e23,6e23,30e23,6e23,163e23,3e23])
Z = np.array([1,2,-1,-2,1,2,-1,-2])
T = 273.15 + 36
epsilon = 74.5

k = inv_Debye_length(epsilon, T, Caq, Z)
#k = k*1e9
d = 1/k