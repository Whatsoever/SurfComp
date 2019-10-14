# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:29:41 2019

@author: DaniJ

Note for my self: The aim of this script is to be extremely readable and well organized.
    From scratch?
"""
import numpy as np
from scipy import linalg
######################### Auxiliar functions #####################################################
def mass_action_law (log_K, X, A):
    '''
        Based on Westall (1980) paper    
    
        log_K   vector of log Ki
        X       vector of primary species 
        A       matrix of stoichiometric coefficients
        
        Note: the returned value is C. Ci = Xi should match the values. Namely if Ci is species Ca+2, and Ca+2 is a primary species, which states that Ca+2 is in X. Xi that belongs to Ca+2, should match  
    '''
    log_C = log_K+np.matmul(A,np.log(X))
    C = 10**(log_C)
    return C

def u_componentvector(A,C):
    '''
       - A --> stoichiometrix matrix [columns=X, rows = C_i] 
       - C --> vector of concentrations
    '''
    u = np.matmul(A.transpose(),C)
    return u

def surface_charge_edgelayer_flm(C,psi_L0,psi_L1):
    '''
       A generic way to calculate the surface charge for layers on the edges in the flm i.e. the O layer
       and the d layer. Using the flm theory.
       - C --> the capacitance (i.e. C1 or C3 in the flm)
       - psi_L0 --> the electrostatic potential in the reference layer (i.e. psi_O or psi_d in the flm model)
       - psi_L1 --> the electrostatic potential away from the reference layer (i.e. the psi_C or psi_A in the flm model)
       Note: The user must be sure of the units, in general electrostatic potential is in volts and 
             the capacitance is in farrads.
    '''
    sigma = C*(psi_L0-psi_L1)
    return sigma

def surface_charge_between_layer_flm(C_left, C_right, psi_mid, psi_left, psi_right):
    '''
       A generic way to calculate the surface charge for the inbetween layers in the flm i.e. the C layer
       and the A layer. Using the flm theory.
       - C_left --> The capacitance between the psi_mid and the psi_left (i.e. C1,C2 or C3)
       - C_right --> The capacitance between the psi_mid and the psi_right
       - psi_mid --> the electrostatic potential of the middle (i.e. the layer reference electrostatic potential. So, psi_C or psi_A in the flm model)
       - psi_left --> the electrostatic potential on the left (i.e. psi_0 or psi_C in the flm model)
       - psi_right --> the electrostatic potential on the right (i.e. psi_A or psi_d in the flm model)
       Note: The user must be sure of the units, in general electrostatic potential is in volts and 
             the capacitance is in farrads.
    '''
    sigma = C_left*(psi_mid-psi_left) + C_right*(psi_mid-psi_right)
    return sigma

def charge_2_mol (charge, s, a, F):
    '''
       The surface charge is multiplyed by specific surface area (or area), solid concentration (or grams) depending what is desired
       the units should be coherent and agree with the whole problem.
       - s is the solid concentration (or grams)
       - a is the specific surface area (or area)
       - F is the Faraday constant
    '''
    Tmol = (charge*s*a)/F
    return Tmol

def boltzman_2_psi(X, R, T, F):
    '''
        - X is the boltzman factor
        - R is the universal gas constant
        - T is the temperature 
        - F is the Faraday constant
        As usual every constant should be coherent
    '''
    partA = (-R*T)/F
    partB = np.log(X)
    psi= partA*partB
    return psi

def calculate_ionicstrength(Z,C):
    '''
        It is supossed to be numpy format vector 
        Z is the vector of charge 
    '''
    # Multiplication must be pointwise for the vector
    # multiply function of numpy. Multiplies pointwise according to the documentation and own experience.
    I = np.matmul(np.multiply(Z,Z),C)
    I = I/2
    return I

def diagonal_row(J):
    num_rows = J.shape[0]
    D = np.zeros((num_rows,num_rows))
    for i in range(0,num_rows):
        D[i,i]=np.sqrt(linalg.norm(J[i,:], np.inf))
    return D

def diagonal_col(J):
    num_cols = J.shape[1]
    D = np.zeros((num_cols,num_cols))
    for i in range(0,num_cols):
        D[i,i]=np.sqrt(linalg.norm(J[:,i], np.inf))
    return D
##############################################################################################
    



#################################### Main function ###########################################
def PB_and_fourlayermodel (log_K, X, A,tolerance = 1e-8, max_iterations=80, idx_fix_species = None):
    """
    -Implements the fours layer model for two surfaces that interact between them.
    
    Arguments:
        A               stoichiometric matrix. Define as Westall (1980)
        X               vector of primary species. Define as Westall (1980)
        log_K           vector of log Ki. Define as Westall (1980)
    """
    abs_err = tolerance+1
    counter_iterations = 0
    abs_err = tolerance + 1
    if idx_fix_species != None:
        X_guess [idx_fix_species] = T [idx_fix_species]
    while abs_err>tolerance and counter_iterations < max_iterations: