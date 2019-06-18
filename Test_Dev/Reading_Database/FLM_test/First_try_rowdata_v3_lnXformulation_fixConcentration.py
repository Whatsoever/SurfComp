# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 09:45:17 2019
"""

import four_layer_model_LNX_withFixSpeciesOption as flm
import numpy as np
import scipy as sp

from matplotlib import pyplot as plt


def funky (T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, tolerance_B, idx_fix_species):
    try:
        #[X,C]=four_layer_one_surface_speciation ( T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, tolerance = 1e-6, max_iterations = 100, debug_flm = None):
        [X,C]=flm.four_layer_one_surface_speciation ( T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, idx_fix_species, tolerance = tolerance_B,max_iterations = 200)
        T_error = tolerance_B
        F.write(str(tolerance_B))
        F.write("\n")
    except:
        tolerance_B = tolerance_B*10
        [X,C, T_error]=funky (T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, tolerance_B, idx_fix_species)
        F.write(str(tolerance_B))
        F.write("\n")
    
    return X,C, T_error

T_H = np.linspace(-3,-11.2,42)
T_H = 10**T_H

idx_fix_species=[0]
#X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 2, 1.5, 2e-14])
X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 0.9, 0.8, 0.8, 0.7])
lnX_guess = np.log(X_guess)

# A in the columns we have X, and in the rows with have the concentration of vector C --> C=[ H^+ Cl- Na+ SOH SOH2+ SO- SOH_2Cl SONa ]
# Remark A_transpose is the U or component matrix.

A = np.array([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[1,0,0,1,1,0,0,0],[-1,0,0,1,-1,0,0,0], [1,1,0,1,1,0,-1,0], [-1,0,1,1,-1,1,0,0]])

# Z, is  a vector with the ionic charge of the aqueous species. It should be somehow in order with the aqueous species. Looking at C, we see only 3 aqueous species H+, Cl-, and Na+
Z=np.array([1,-1,1])

# Log_k the log of reactions. Primary species have log_k = 0, I assume that it follows the same order than C

log_k = np.array([0, 0, 0, 0, 4, -8.8, 5.82, -7])
ln_k = log_k/np.log10(np.e)      # Changing the base from log_10 to ln (log_e)


# - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
'I must check if the first position of idx_Aq is 0 or 1. Since it is python I guess it will be 0'
idx_Aq=np.array([0,1,2])

# pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma basically the same thing that idx_Aq, but only scalar. (Somehow specified). Problably position should agree with X, or T.
'Maybe I can change the name of the psialpha for psiC, psibeta for psiA, and psigamma for psid. Right now, not really the aim.'
pos_eb_0=4
pos_eb_c=5
pos_eb_a=6
pos_eb_d=7

# Temperature
temp=273.15+25
#s             is the specific surface area
s=1         # m2/l
a=1         # g/l
epsilon = 78.45203739768931
C_vector=[1.05, 3.36, 0.27] 

tolerance_vector=[]
Array_X = []
Array_C = []



#T=np.array([T_H[0], 1e-3, 1e-3, 9.9635e-6, 1, 1, 1, 1]) 
#[X,C]=flm.four_layer_one_surface_speciation ( T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, max_iterations = 500)
#[X,C]=flm.four_layer_one_surface_speciation ( T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, max_iterations = 500)
global F 
F= open("T_efile_lnX.txt","w") 
for i in range(0,len(T_H)):
    # T_transpose =[ T_H+ T_Cl- T_Na+ T_SOH T_σ_0 T_σ_c T_σ_A T_σ_d ] 
    T=np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 1, 1, 1, 1]) 
    X_guess = np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 0.9, 0.8, 0.9])
    print(i)
    tolerance_B=1e-8
    [X,C, T_e]= funky (T, lnX_guess, A, Z, ln_k, idx_Aq, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, tolerance_B,idx_fix_species)
    tolerance_vector.append(T_e)
    if i == 0:
        Array_X = X
        Array_C = C
    else:
        Array_X = np.vstack([Array_X, X])
        Array_C = np.vstack([Array_C, C])
F.close()
np.save('tol_vec_v3_lnX_fix',tolerance_vector)
np.save('X_arr_v3_lnX_fix',Array_X)
np.save('C_arr_v3_lnX_fix',Array_C)

# The plotting is done by Data_plotting file
