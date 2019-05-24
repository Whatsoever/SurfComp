# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:38:31 2019

@author: DaniJ
"""
import four_layer_model_2try as flm
import numpy as np
import scipy as sp
"""
    NOTE: This excersice is the same that Westall(1980) but with one layer more. Rudzinksi et al. have the values, year of the paper???
"""
# The info of this file will be in the daniel_jara_FLM_adapting .docx file. Although, I will try to document everything as good as possible.
# As usual, we try to follow Westall formulation

# The arguments of the four layer should follow a determine logical order, such as the one that can be found in the Westall(1980) paper.

# T_transpose =[ T_H+ T_Cl- T_Na+ T_SOH T_σ_0 T_σ_c T_σ_A T_σ_d ] 
T=np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 1, 1, 1, 1])        # The values of T related to the electrostatic potentials are set to 0, it should not matter which values is given.

# X_guess
# X = [ H+ Cl- Na+ SOH (e^(-Fψ_0/RT) (e^(-Fψ_c/RT) (e^(-Fψ_A/RT) (e^(-Fψ_d/RT) ) ]
#X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 2, 1.5, 2e-14])
X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9])

# A in the columns we have X, and in the rows with have the concentration of vector C --> C=[ H^+ Cl- Na+ SOH SOH2+ SO- SOH_2Cl SONa ]
# Remark A_transpose is the U or component matrix.

A = np.array([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[1,0,0,1,1,0,0,0],[-1,0,0,1,-1,0,0,0], [1,1,0,1,1,0,-1,0], [-1,0,1,1,-1,1,0,0]])

# Z, is  a vector with the ionic charge of the aqueous species. It should be somehow in order with the aqueous species. Looking at C, we see only 3 aqueous species H+, Cl-, and Na+
Z=np.array([1,-1,1])

# Log_k the log of reactions. Primary species have log_k = 0, I assume that it follows the same order than C

log_k = np.array([0, 0, 0, 0, 4, -8.8, 5.82, -7])

# - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
'I must check if the first position of idx_Aq is 0 or 1. Since it is python I guess it will be 0'
idx_Aq=np.array([0,1,2])

# pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma basically the same thing that idx_Aq, but only scalar. (Somehow specified). Problably position should agree with X, or T.
'Maybe I can change the name of the psialpha for psiC, psibeta for psiA, and psigamma for psid. Right now, not really the aim.'
pos_psi0=4
pos_psialpha=5
pos_psibeta=6
pos_psigamma=7

# Temperature
temp=273.15+25
#s             is the specific surface area
s=1         # m2/l
a=1         # g/l
e = 78.45203739768931
Capacitances=[1.05, 3.36, 0.27] 

[X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-7,max_iterations = 300)