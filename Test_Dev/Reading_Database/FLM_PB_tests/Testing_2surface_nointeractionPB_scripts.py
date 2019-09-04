# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 16:32:57 2019

@author: DaniJ
"""

import four_layer_model_2try_withFixSpeciesOption_Scaling_2surface as flm1
import four_layer_model_LNX_withFixSpeciesOption_Scaling_2surface as flm2
import numpy as np
import scipy as sp

from matplotlib import pyplot as plt



idx_fix_species=[0]
#idx_fix_species=None
T_H = np.linspace(-3,-11.2,42)
T_H = 10**T_H


#
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
ln_k = log_k/np.log10(np.e)      # Changing the base from log_10 to ln (log_e)
# - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
'I must check if the first position of idx_Aq is 0 or 1. Since it is python it is  0'
idx_Aq=np.array([0,1,2])

pos_psi_S1_vec = np.array([5,6,7,8])
pos_psi_S2_vec = np.array([9,10,11,12])


# Temperature
temp=273.15+25
#s             is the specific surface area
sS1=1         # m2/l
aS1=1         # g/l
sS2=1         # m2/l
aS2=1         # g/l
e = 78.45203739768931
CapacitancesS1=[1.05, 3.36, 0.27] 
CapacitancesS2=[1.05, 3.36, 0.27] 

Array_X1 = []
Array_C1 = []
Array_X2 = []
Array_C2 = []

for i in range(0,len(T_H)):
    T=np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9,8.7e-7, 0.9, 0.9, 0.9]) 
    lnX_guess = np.log(X_guess)
    [X1,C1] = flm1.four_layer_two_surface_speciation(T, X_guess, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, idx_fix_species , zel=1, tolerance = 1e-12, max_iterations = 100, scalingRC = True)
    [X2,C2] = flm2.four_layer_two_surface_speciation(T, lnX_guess, A, Z, ln_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, idx_fix_species, tolerance = 1e-12, max_iterations = 100, scalingRC = True, debug_flm = None)
    lnX_guess=np.log(X2)
    X_guess = X1
    if i == 0:
        Array_X1 = X1
        Array_C1 = C1
        Array_X2 = X2
        Array_C2 = C2
    else:
        Array_X1 = np.vstack([Array_X1, X1])
        Array_C1 = np.vstack([Array_C1, C1])
        Array_X2 = np.vstack([Array_X2, X2])
        Array_C2 = np.vstack([Array_C2, C2])
np.save('X_arr_X',Array_X1)
np.save('C_arr_X',Array_C1)
np.save('X_arr_lnX',Array_X2)
np.save('C_arr_lnX',Array_C2)