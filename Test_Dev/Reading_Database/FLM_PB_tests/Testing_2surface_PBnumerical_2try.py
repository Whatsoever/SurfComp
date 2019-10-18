# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:07:18 2019

@author: DaniJ
"""

import PB_coup_fourlayer_3try as flmPB
import four_layer_model_2try_withFixSpeciesOption_Scaling_2surface as flm1
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

def Boltzman_factor_2_psi (x,temp):
    '''
        Transforms the equation from Xb = exp(-psi*F/RT) to psi = -ln(Xb)RT/F
        from Boltzman factor to electrostatic potential
        The units of "temp" (short for temperature) should be Kelvin
    '''
    R =  8.314472                                       # J/(K*mol)
    F = 96485.3328959                                   # C/mol                                    
    D = R*temp
    psi = - np.log(x)*(D/F)
    return psi   


idx_fix_species=[0]

#idx_fix_species=None
#T_H = np.linspace(-3,-11.2,42)
#T_H = 10**T_H

X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9,8.7e-7, 0.9, 0.9, 0.9])

# A in the columns we have X, and in the rows with have the concentration (vector C)
# Remark A_transpose is the U or component matrix.

A = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,0],[1,0,0,1,0,1,0,0,0,0,0,0,0], \
              [-1,0,0,1,0,-1,0,0,0,0,0,0,0], [1,1,0,1,0,1,0,-1,0,0,0,0,0], [-1,0,1,1,0,-1,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0],[1,0,0,0,1,0,0,0,0,1,0,0,0],\
              [-1,0,0,0,1,0,0,0,0,-1,0,0,0],[1,1,0,0,1,0,0,0,0,1,0,-1,0],[-1,0,1,0,1,0,0,0,0,-1,1,0,0]])

# Z, is  a vector with the ionic charge of the aqueous species. It should be somehow in order with the aqueous species. Looking at C, we see only 3 aqueous species H+, Cl-, and Na+
Z=np.array([1,-1,1])
# Log_k the log of reactions. Primary species have log_k = 0, I assume that it follows the same order than C

log_k = np.array([0, 0, 0, 0, 4, -8.8, 5.82, -7,0,4,-8.8,5.82,-7])


idx_Aq=np.array([0,1,2])

pos_psi_S1_vec = np.array([5,6,7,8])
pos_psi_S2_vec = np.array([9,10,11,12])


# Temperature
temp=273.15+25
#s             is the specific surface area
sS1=1         # m2/g
aS1=1         # g/l
sS2=1         # m2/g
aS2=1         # g/l
e = 78.45203739768931
CapacitancesS1=[1.05, 3.36, 0.27] 
CapacitancesS2=[1.05, 3.36, 0.27] 

#

d0=0        # initial position
df=400e-9   # final position
h = np.linspace(d0,df,100)
psi_pb_0 = np.zeros((2, h.size))
T=np.array([1e-3, 1e-2, 1.1e-2, 9.9635e-6, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9,8.7e-7, 0.9, 0.9, 0.9]) 
#(T, X_guess, distance, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, x, idx_fix_species = None, zel=1, tolerance = 1e-6, max_iterations = 100,scalingRC = True):
[X_guess,C1] = flm1.four_layer_two_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, idx_fix_species, zel=1, tolerance = 1e-12, max_iterations = 100,scalingRC = True)
X_guess[8] = Boltzman_factor_2_psi (X_guess[8],temp)
X_guess[12] = Boltzman_factor_2_psi (X_guess[12],temp)


#def PB_and_fourlayermodel (distance, psi_pb_0, log_K, X, A, T,pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, idx_aq, Z_aq_vec, temp, epsilon,tolerance = 1e-8, max_iterations=80, idx_fix_species = None,scalingRC = True)
#[X,C] = flmPB.PB_and_fourlayermodel(T, X_guess, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, d0,df, idx_fix_species, zel=1, tolerance_NR = 1e-8, max_iterations = 10,scalingRC = True, tol_PB=1e-4)
[X,C] = flmPB.PB_and_fourlayermodel(h, psi_pb_0, log_k, X_guess, A, T, pos_psi_S1_vec, pos_psi_S2_vec, CapacitancesS1, CapacitancesS2, aS1,aS2,sS1,sS2, idx_Aq, Z, temp, e, idx_fix_species=idx_fix_species)
