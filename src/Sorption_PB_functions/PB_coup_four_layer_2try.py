# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:58:31 2019

@author: DaniJ
"""

import numpy as np
import scipy as sp
from bvp import solve_bvp
from four_layer_model_2try_withFixSpeciesOption_Scaling import four_layer_model_2try_withFixSpeciesOption_Scaling as flm

'''
    In this first try we will assume that the vector of unknowns is composed in the following order:
'''


def PB_and_fourlayermodel (T, X_guess, distance, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, idx_fix_species = None, zel=1, tolerance = 1e-6, max_iterations = 100,scalingRC = True):
    counter_iterations = 0
    abs_err = tolerance + 1
    if idx_fix_species != None:
        X_guess [idx_fix_species] = T [idx_fix_species]
    while abs_err>tolerance and counter_iterations < max_iterations:
        # Calculate Y
        [Y, T] = func_NR_FLM (X_guess, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, idx_fix_species)
        

def func_NR_FLM (X, A, log_k, temp, idx_Aq,  sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, idx_fix_species=None):
    """
        This function is supossed to be linked to the four_layer_two_surface_speciation function.
        It just gave the evaluated vector of Y, and T for the Newton-raphson procedure.
        The formulation of Westall (1980) is followed.
        FLM = four layer model
    """
    # Speciation - mass action law
    log_C = log_k + np.matmul(A,np.log10(X))
    # transf
    C = 10**(log_C)
    # Update T - "Electrostatic parameters"
    psi_S1_v = [Boltzman_factor_2_psi(X[pos_psi_S1_vec[0]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[1]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[3]], temp)]
    psi_S2_v = [Boltzman_factor_2_psi(X[pos_psi_S2_vec[0]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[1]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[3]], temp)] 
    C_aq = C[idx_Aq]
    I = Calculate_ionic_strength(Z, C_aq)