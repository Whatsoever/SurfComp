# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:29:41 2019

@author: DaniJ

Note for my self: The aim of this script is to be extremely readable and well organized.
    From scratch?
"""
import numpy as np
from scipy import linalg
import scipy as sp
from bvp import solve_bvp
from scipy import linalg
import matplotlib.pyplot as plt
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
def bc_PB(ya, yb, args):   
    y0 = args[4]
    return np.array([ya[0]-y0[0,0] , yb[0]-y0[0,-1]])

def fun_PB(x,y,args):
    me_ee0 = args[0]
    Z_aq = args[1]
    ni =  args[2] 
    e_kbt =  args[3]
    for i in range(len(ni)):
        summation = Z_aq[i]*ni[i]*np.exp(-Z_aq[i]*e_kbt*y[0])
    pb_righthandside = me_ee0*summation
    return num.vstack((y[1] ,  pb_righthandside))    

def integrate_PB_solution(psi_PB, epsilon_0, epsilon, C_aq, Z_aq):
    '''
    psi_PB is supose to be an object of bvp solver of scipy.
    '''
    kb = 1.38064852e-23                                 # J/K other units --> kb=8,6173303e-5  eV/K
    Na = 6.022140857e23                                 # 1/mol
    elec_charge = 1.60217662e-19 #electron charge in C
    ni = Na*1000*C_aq
    e_kbt = elec_charge/(kb*temp)
    me_ee0=-elec_charge/(epsilon_0*epsilon)
    psi_vec = psi_PB.y
    n_mesh = len(psi_PB.x)
    y = np.zeros(n_mesh)
    # create vector y before numerical integration
    for i in range(len(n_mesh)):
        for j in range(len(ni)):
            summation = Z_aq[j]*ni[j]*np.exp(-Z_aq[j]*e_kbt*psi_vec[0,i])
        y[i] = me_ee0*summation
    sigma_PB = np.trapz(y,psi_PB.x)
    sigma_PB_S1 = (epsilon_0*epsilon)*psi_vec[1,0]
    sigma_PB_S2 = (epsilon_0*epsilon)*psi_vec[1,-1]
    return sigma_PB, sigma_PB_S1, sigma_PB_S2
    
def totals_sigmas(T, pos_ele_S1, pos_ele_S2, sigma_FLM_S1, sigma_FLM_S2,F,s1,s2,a1,a2):
    for i in range(len(pos_ele_S1)-1):
        T[pos_ele_S1[i]] = charge_2_mol (sigma_FLM_S1[i], s1, a1, F)
        T[pos_ele_S2[i]] = charge_2_mol (sigma_FLM_S2[i], s2, a2, F)
        
def calculate_psi_PB(epsilon_0, epsilon, C_aq, Z_aq_vec, psi_pb_0, distance, temp):    
    kb = 1.38064852e-23                                 # J/K other units --> kb=8,6173303e-5  eV/K
    Na = 6.022140857e23                                 # 1/mol
    elec_charge = 1.60217662e-19 #electron charge in C
    ni = Na*1000*C_aq
    e_kbt = elec_charge/(kb*temp)
    me_ee0=-elec_charge/(epsilon_0*epsilon)
    args=[me_ee0, Z_aq, ni, e_kbt, psi_pb_0]
    psi_PB= solve_bvp(fun_PB, bc_PB, distance, psi_pb_0, tol = 1e-4, args = args)
    return psi_PB
    
def modification_initialvalue_PBsolver(X_psi_S1, X_psi_S2, sigma_FLM_S1, sigma_FLM_S2, epsilon_0, epsilon, psi_pb_0):
    psi_pb_0[0,0] = X_psi_S1
    psi_pb_0[0,-1] = X_psi_S2
    # The following part I am not sure, but anyway we thing it is not really affecting or its effect is supposed to be negligible.
    #Gauss part
    psi_pb_0[1,0] = sigma_FLM_S1/(epsilon_0*epsilon)
    psi_pb_0[0,-1] = -sigma_FLM_S2/(epsilon_0*epsilon)
    return psi_pb_0
    
def calculate_surface_charge_FLM(X_psi, C_vec):
    sigma_0 = surface_charge_edgelayer_flm(C_vec[0],X_psi[0],X_psi[1])
    sigma_C = surface_charge_between_layer_flm(C_vec[0], C_vec[1], X_psi[1], X_psi[0], X_psi[2])
    sigma_A = surface_charge_between_layer_flm(C_vec[1], C_vec[2], X_psi[2], X_psi[1], X_psi[3])
    sigma_D = - sigma_0 - sigma_C - sigma_A
    return [sigma_0, sigma_C, sigma_A,sigma_D]

def T_electrostatic_update (X, T, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1, a2, s1, s2, temp,R, F, epsilon_0, epsilon, C_aq, Z_aq_vec,psi_pb_0, distance):
    X_psi_S1 = X[pos_ele_S1]
    X_psi_S2 = X[pos_ele_S2]
    X_psi_S1[:-1] = boltzman_2_psi(X_psi_S1[:-1], R, temp, F)
    X_psi_S2[:-1] = boltzman_2_psi(X_psi_S2[:-1], R, temp, F)
    sigma_FLM_S1 = calculate_surface_charge_FLM(X_psi_S1, C_vec_S1)
    sigma_FLM_S2 = calculate_surface_charge_FLM(X_psi_S2, C_vec_S2)
    psi_pb_0 = modification_initialvalue_PBsolver(X_psi_S1[-1], X_psi_S2[-1], sigma_FLM_S1[-1], sigma_FLM_S2[-1], epsilon_0, epsilon, psi_pb_0)
    psi_PB = calculate_psi_PB(epsilon_0, epsilon, C_aq, Z_aq_vec,psi_pb_0, distance, temp)
    [sigma_PB, sigma_PB_S1, sigma_PB_S2] = integrate_PB_solution(psi_PB, epsilon_0, epsilon, C_aq, Z_aq_vec)
    # Assigning some Ts
    T = totals_sigmas(T, pos_ele_S1, pos_ele_S2, sigma_FLM_S1, sigma_FLM_S2,F,s1,s2,a1,a2)
    # The T for sigma must be different, since the part of psi_d related to u is equal to 0, here we should have already 0.
    T[pos_ele_S1[-1]]= sigma_FLM_S1[-1] - sigma_PB_S1
    T[pos_ele_S2[-1]]= sigma_FLM_S2[-1] - sigma_PB_S2
    return T

def calculate_jacobian_AC_part(log_K, X, A):
    C = mass_action_law (log_K, X, A)
    n = X.size
    Z = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, n):
            Z[i,j]= np.matmul(np.multiply(A[:,i], A[:,j]), (C/X[j]))
    return Z, C

def calculate_jacobian_FLM_part(J, X, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, temp, F,R):
    X_psi_S1 = X[pos_ele_S1]
    X_psi_S2 = X[pos_ele_S2]
    
    RT_F = -(R*temp)/F
    p = [pos_ele_S1, pos_ele_S2]
    sa_F = [(s1*a1)/F, (s2*a2)/F]
    Cap = [C_vec_S1, C_vec_S2]
    X_psi = [X_psi_S1, X_psi_S2]
    # plane 0
    for i in [0,1]:
        J[p[i][0], p[i][0]] = J[p[i][0], p[i][0]] - sa_F[i]*Cap[i][0]*RT_F*(1/X_psi[i][0])
        J[p[i][0], p[i][1]] = J[p[i][0], p[i][1]] + sa_F[i]*Cap[i][0]*RT_F*(1/X_psi[i][1])
    
    # plane C
    for i in [0,1]:
        J[p[i][1], p[i][0]] = J[p[i][1], p[i][0]] + sa_F[i]*Cap[i][0]*RT_F*(1/X_psi[i][0])
        J[p[i][1], p[i][1]] = J[p[i][1], p[i][1]] - sa_F[i]*(Cap[i][0]+Cap[i][1])*RT_F*(1/X_psi[i][1])
        J[p[i][1], p[i][2]] = J[p[i][1], p[i][2]] + sa_F[i]*Cap[i][1]*RT_F*(1/X_psi[i][2])
    
    # plane A
    for i in [0,1]:
        J[p[i][2], p[i][1]] = J[p[i][2], p[i][1]] + sa_F[i]*Cap[i][1]*RT_F*(1/X_psi[i][1])
        J[p[i][2], p[i][2]] = J[p[i][2], p[i][2]] - sa_F[i]*(Cap[i][1]+Cap[i][2])*RT_F*(1/X_psi[i][2])
        J[p[i][2], p[i][3]] = J[p[i][2], p[i][3]] + sa_F[i]*Cap[i][2]
    # plane D
    for i in [0,1]:
        J[p[i][3],p[i][2]] = J[p[i][3],p[i][2]] - Cap[i][2]*RT_F*(1/X_psi[i][2])
        J[p[i][3],p[i][3]] = J[p[i][3],p[i][3]] + Cap[i][2]
    return J

def calculate_term_PB_jac (J, psi_pb_0, psi_pb_1, psi_pb_2, delta_psi, epsilon_0, epsilon, pos_d_S1, pos_d_S2):
    J[pos_d_S1, pos_d_S1] = J[pos_d_S1, pos_d_S1] - (epsilon_0*epsilon)*((psi_pb_1[1,0]-psi_pb_0[1,0])/delta_psi)
    J[pos_d_S2, pos_d_S2] = J[pos_d_S2, pos_d_S2] - (epsilon_0*epsilon)*((psi_pb_2[1,0]-psi_pb_0[1,0])/delta_psi)
    return J

def calculate_jacobian_PB_part(J, psi_pb_0, distance, X, C_aq, pos_ele_S1, pos_ele_S2,C_vec_S1, C_vec_S2, temp, epsilon_0, epsilon, Z_aq_vec):
    delta_psi = 0.001
    X_psi_S1 = X[pos_ele_S1]
    X_psi_S2 = X[pos_ele_S2]
    X_psi_S1[:-1] = boltzman_2_psi(X_psi_S1[:-1], R, temp, F)
    X_psi_S2[:-1] = boltzman_2_psi(X_psi_S2[:-1], R, temp, F)
    sigma_FLM_S1 = calculate_surface_charge_FLM(X_psi_S1, C_vec_S1)
    sigma_FLM_S2 = calculate_surface_charge_FLM(X_psi_S2, C_vec_S2)
    psi_pb_00 = modification_initialvalue_PBsolver(X_psi_S1[-1], X_psi_S2[-1], sigma_FLM_S1[-1], sigma_FLM_S2[-1], epsilon_0, epsilon, psi_pb_0)
    psi_pb_01 = modification_initialvalue_PBsolver(X_psi_S1[-1] + delta_psi, X_psi_S2[-1], sigma_FLM_S1[-1], sigma_FLM_S2[-1], epsilon_0, epsilon, psi_pb_0)
    psi_pb_02 = modification_initialvalue_PBsolver(X_psi_S1[-1], X_psi_S2[-1]+delta_psi, sigma_FLM_S1[-1], sigma_FLM_S2[-1], epsilon_0, epsilon, psi_pb_0)
    
    psi_pb_0 =calculate_psi_PB(epsilon_0, epsilon, C_aq, Z_aq_vec, psi_pb_00, distance, temp)
    psi_pb_1 =calculate_psi_PB(epsilon_0, epsilon, C_aq, Z_aq_vec, psi_pb_01, distance, temp)
    psi_pb_2 =calculate_psi_PB(epsilon_0, epsilon, C_aq, Z_aq_vec, psi_pb_02, distance, temp)
    
    J = calculate_term_PB_jac (J,psi_pb_0, psi_pb_1, psi_pb_2, delta_psi, epsilon_0, epsilon, pos_ele_S1[-1], pos_ele_S2[-1])


def residual_function_calculation(log_K, X, A, T, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, temp, F, R, epsilon_0, epsilon, idx_aq, Z_aq_vec,psi_pb_0, distance, idx_fix_species):
    # Calculate species
    C = mass_action_law (log_K, X, A)
    C_aq = C[idx_aq]
    # Components
    u = u_componentvector(A,C)
    # The electrostatic part of T must be update
    T = T_electrostatic_update(X, T, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1, a2, s1, s2, temp, R, F, epsilon_0, epsilon, C_aq, Z_aq_vec,psi_pb_0, distance)
    # return residual function
    Y = u - T
    if idx_fix_species != None:
        Y[idx_fix_species]=0
    return Y,T

def calculate_jacobian_function(log_K, X, A, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, temp, F,R,epsilon_0, epsilon, idx_aq, Z_aq_vec, psi_pb_0, distance, idx_fix_species):
    [J,C] = calculate_jacobian_AC_part(log_K, X, A)
    J = calculate_jacobian_FLM_part(J, X, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, temp, F,R)
    C_aq = C[idx_aq]
    J = calculate_jacobian_PB_part(J, psi_pb_0, distance, X, C_aq, pos_ele_S1, pos_ele_S2,C_vec_S1, C_vec_S2, temp, epsilon_0, epsilon, Z_aq_vec)
    if idx_fix_species != None:
        for d in idx_fix_species:
            v=np.zeros(X.size)
            v[d]=1
            J[d,:] = v
    return J
#################################### Main function ###########################################
def PB_and_fourlayermodel (distance, psi_pb_0, log_K, X, A, T,pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, idx_aq, Z_aq_vec, temp, epsilon,tolerance = 1e-8, max_iterations=80, idx_fix_species = None,scalingRC = True):
    """
    -Implements the fours layer model for two surfaces that interact between them.
    
    Arguments:
        A               stoichiometric matrix. Define as Westall (1980)
        X               vector of primary species. Define as Westall (1980)
        log_K           vector of log Ki. Define as Westall (1980)
        T               vector of totals --> THE POSITION OF T must be the same THAT THE POSITON of X!!!!!!
        pos_ele_S1      vector of the indexes of the electrostatic psi_0, psi_C, psi_A, psi_D. They have the same position at X and at T
        pos_ele_S2      idem but for surface2
        a1              surface area of the solid per mass [m2/g]
        a2              idem but for surface2
        s1
        s2              idem but for surface2
        idx_aq          vector of indexes of aqueous species that are found in C (log_C = log_K + A*log(X))
        Z_aq_vec        is a vector of the valancies of the aqueous species. The order of the values on this vector is given by the order of aqueous species in C (log_C = log_K + A*log(X))
    """
    F = 96485.3328959                                   # C/mol
    R =  8.314472                                       # J/(K*mol)
    epsilon_0 = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    abs_err = tolerance+1
    counter_iterations = 0
    abs_err = tolerance + 1
    if idx_fix_species != None:
        X [idx_fix_species] = T [idx_fix_species]
    while abs_err>tolerance and counter_iterations < max_iterations:
        # First step is to calculate the Residual function
        [Y,T] = residual_function_calculation(log_K, X, A, T, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, temp, F,R,epsilon_0, epsilon, idx_aq, Z_aq_vec, psi_pb_0, distance, idx_fix_species)
        # Second step is to calculate the Jacaobian
        J = calculate_jacobian_function(log_K, X, A, pos_ele_S1, pos_ele_S2, C_vec_S1, C_vec_S2, a1,a2,s1,s2, temp, F,R,epsilon_0, epsilon, idx_aq, Z_aq_vec, psi_pb_0, distance, idx_fix_species)
        # solve
        if scalingRC == True:
            D1 = diagonal_row(J)
            D2 = diagonal_col(J)
        
            J_new = np.matmul(D1,np.matmul(J, D2))
            Y_new = np.matmul(D1, Y)
            delta_X_new = linalg.solve(J_new,-Y_new)
            delta_X = np.matmul(D2, delta_X_new)
        else:
            # Calculating the diff, Delta_X
            delta_X = linalg.solve(J,-Y)
        max_1 = 1
        max_2 =np.amax(-2*np.multiply(delta_X, 1/X))
        Max_f = np.amax([max_1, max_2])
        Del_mul = 1/Max_f
        X=X + Del_mul*delta_X
        
        if idx_fix_species != None:
            Y[idx_fix_species] =0
        abs_err = max(abs(Y))
         counter_iterations += 1
         
    if counter_iterations >= max_iterations or np.isnan(abs_err):
            raise ValueError('Max number of iterations exceed.')      
    C = mass_action_law (log_K, X, A)
    return X, C