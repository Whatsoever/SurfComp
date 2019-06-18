# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 21:38:42 2019

"""

import numpy as np
from scipy import linalg

# try to keep it in block

##################### basic functions ################################################

def mass_action_law (ln_X, ln_K, A):
    '''
       all inputs are numpy arrays!!!
       NO_activity!!!!
       ln [C_i] = log_K_i + Sum(aij*ln_Xj)
       ln_C = A*ln_X+ln_K
       parameters:
           - ln_X --> vector of primary variables
           - A --> stoichiometrix matrix [columns=X, rows = C_i] 
           - ln_K --> vector of equilibrium constant
    '''
    ln_C = ln_K+np.matmul(A,ln_X)
    return ln_C

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

def surface_charge_diffusive_monovalentelectrolyte (R, T, epsilon, epsilon_0, ionic_strength, F, psi_d):
    '''
       If previously the units were important, here the coherence between units is even more important
       sigma_d =〖-(8*1000*RTε_o εI)〗^(1/2) sinh((Fψ_d)/2RT)
    '''
    partA = np.sqrt(8*1000*R*T*epsilon*epsilon_0*ionic_strength)
    inner_B = (F*psi_d)/(2*R*T)
    partB = np.sinh(inner_B)
    sigma_d = -partA*partB
    return sigma_d
    
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
####################### functions of basic functions ###############################
'relative to residual function'
def calculate_T (X, C, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, epsilon_0, C_vector, R, T, F,Z):
    X_0 = X[pos_eb_0]
    X_C = X[pos_eb_c]
    X_A = X[pos_eb_a]
    X_D = X[pos_eb_d]
    'Now the psi'
    psi_0 = boltzman_2_psi(X_0, R, temp, F)
    psi_C = boltzman_2_psi(X_C, R, temp, F)
    psi_A = boltzman_2_psi(X_A, R, temp, F)
    psi_D = boltzman_2_psi(X_D, R, temp, F)
    'Now the charge'
    #previous to charge
    Caq = C[idx_Aq]
    ionic_strength = calculate_ionicstrength(Z, Caq)
    # charge
    sigma_0 = surface_charge_edgelayer_flm(C_vector[0],psi_0,psi_C)
    sigma_C = surface_charge_between_layer_flm(C_vector[0], C_vector[1], psi_C, psi_0, psi_A)
    sigma_A = surface_charge_between_layer_flm(C_vector[1], C_vector[2], psi_A, psi_C, psi_D)
    sigma_d_flm = surface_charge_edgelayer_flm(C_vector[2],psi_D,psi_A)
    sigma_d_pb = surface_charge_diffusive_monovalentelectrolyte (R, temp, epsilon, epsilon_0, ionic_strength, F, psi_D)
    'Change value in the electrostatic positions'
    T[pos_eb_0] = charge_2_mol(sigma_0, s, a, F)
    T[pos_eb_c] = charge_2_mol(sigma_C, s, a, F)
    T[pos_eb_a] = charge_2_mol(sigma_A, s, a, F)
    T[pos_eb_d] = charge_2_mol(sigma_d_pb, s, a, F) - charge_2_mol(sigma_d_flm, s, a, F)
    return T


def calculate_residual_function(T,ln_X, ln_K, A, idx_Aq, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, epsilon_0, C_vector, R, F,Z, idx_fix_species = None):
    ln_C = mass_action_law (ln_X, ln_K, A)
    C = np.exp(ln_C)
    u = u_componentvector(A,C)
    X = np.exp(ln_X)
    T = calculate_T (X, C, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, epsilon_0, C_vector, R, T, F, Z) 
    Y = u-T
    if idx_fix_species != None:
        Y[idx_fix_species]=0
    return Y,T

'relative to Jacobian'

def calculate_J_classicalPart(ln_X, ln_K, A):
    ln_C = mass_action_law (ln_X, ln_K, A)
    C = np.exp(ln_C)
    n = len(ln_X)
    Z = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, n):
            Z[i,j]= np.matmul(np.multiply(A[:,i], A[:,j]), C)
    return Z, C

def calculate_electrostatic_part (J, s, a, R, T, C_vector, Caq, Z, F, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, epsilon, epsilon_0, psi_d):
    # plane 0
    J[pos_eb_0,pos_eb_0] = J[pos_eb_0, pos_eb_0] + ((C_vector[0]*R*T*s*a)/(F*F))
    J[pos_eb_0, pos_eb_c] = J[pos_eb_0, pos_eb_c] - ((C_vector[0]*R*T*s*a)/(F*F))
    
    # plane C
    J[pos_eb_c, pos_eb_0] = J[pos_eb_c, pos_eb_0] - ((C_vector[0]*R*T*s*a)/(F*F))
    J[pos_eb_c, pos_eb_c] = J[pos_eb_c, pos_eb_c] + (((C_vector[0]+C_vector[1])*R*T*s*a)/(F*F))
    J[pos_eb_c, pos_eb_a] = J[pos_eb_c, pos_eb_a] - ((C_vector[1]*R*T*s*a)/(F*F))
    
    # plane A
    J[pos_eb_a, pos_eb_c] = J[pos_eb_a, pos_eb_c] - ((C_vector[1]*R*T*s*a)/(F*F))
    J[pos_eb_a, pos_eb_a] = J[pos_eb_a, pos_eb_a] + (((C_vector[1]+C_vector[2])*R*T*s*a)/(F*F))
    J[pos_eb_a, pos_eb_d] = J[pos_eb_a, pos_eb_d] - ((C_vector[2]*R*T*s*a)/(F*F))
    
    #plane D
    J[pos_eb_d, pos_eb_a] = J[pos_eb_d, pos_eb_a] - ((R*T*s*a*C_vector[2])/(F*F))
    J[pos_eb_d, pos_eb_d] = calculate_derivative_Td (C_vector[2], R, T, F, Caq, Z, epsilon, epsilon_0, psi_d,s,a)
    return J

def calculate_derivative_Td (C, R, T, F, Caq, Z, epsilon, epsilon_0, psi_d,s,a):
    ionic_strength = calculate_ionicstrength(Z, Caq)
    #
    DT_Dpsid = -np.sqrt(8*1000*R*T*epsilon*epsilon_0*ionic_strength)*np.cosh((F*psi_d)/(2*R*T))*(F/(2*R*T)) - C
    Dpsid_DlnXpsid = (-R*T)/F
    j_d = DT_Dpsid*Dpsid_DlnXpsid*((s*a)/F)
    return j_d

def calculate_jacobian_function(ln_X, ln_K, A, idx_Aq, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, epsilon_0, C_vector, R, F,Z,idx_fix_species = None):
    length_X=len(ln_X)
    #
    [J,C] = calculate_J_classicalPart(ln_X, ln_K, A)
    Caq = C[idx_Aq]
    X = np.exp(ln_X)
    psi_d = boltzman_2_psi(X[pos_eb_d], R, temp, F)
    J = calculate_electrostatic_part (J, s, a, R, temp, C_vector, Caq, Z, F, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, epsilon, epsilon_0, psi_d)
    # finally just return Z
    if idx_fix_species != None:
        for d in idx_fix_species:
            v=np.zeros(length_X)
            v[d]=1
            J[d,:] = v
    return J

###################### SOLVING ####################################################

def four_layer_one_surface_speciation ( T, lnX_guess, A, Z, ln_k, idx_Aq,pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, C_vector, idx_fix_species = None, tolerance = 1e-6, max_iterations = 100, debug_flm = None):
    '''
        - T --> The vector of Total values (The electrostatic values will be recalculated, so it does not matter what has been introduced)
        - lnX_guess --> The vector of primary vairables, it might be preconditioned in the future.
        - A --> stoichiometrix and component matrix (i.e. transpose). Number of rows = number species, Number of columns = number of primary variables
        - ln_k --> A vector of log(Konstant equilibrium). Primary species of aquoues and sorption have a log_k=0
        - idx_Aq --> An index vector with the different aqueous species position. It must coincide with the rows of "A".
        - Z --> The vector of charge of the different ion. The order is determine by the rows of "A" for aqueous species. That means that it is link to idx_Aq somehow.
        - pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d -->  This is basically the position of the boltzman factor for the different planes
        - s --> concentration of suspended solid. 
        - a --> is the specific surface area
        - epsilon --> relative permittivity
        - C_vector --> [C1, C2, C3] 
        - temp --> Temperature of the chemical system in Kelvins.
        - debug_flm --> the class is given, only if important information about a problem is desired.
    '''
    # Instantiation of parameters that are constant
    F = 96485.3328959                                   # C/mol
    R =  8.314472                                       # J/(K*mol)
    epsilon_0 = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    if idx_fix_species != None:
        lnX_guess [idx_fix_species] = np.log(T [idx_fix_species])
    ln_X = lnX_guess
    #X = np.exp(ln_X)
    # instantiation variables for loop
    counter_iterations = 0
    abs_err = tolerance + 1
    while abs_err>tolerance and counter_iterations < max_iterations:
        # Calculate Residual function
        [Y,T] = calculate_residual_function(T,ln_X, ln_k, A, idx_Aq, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, epsilon_0, C_vector, R, F,Z,idx_fix_species)
        # Calculate Jacobian Residual function
        J = calculate_jacobian_function(ln_X, ln_k, A, idx_Aq, pos_eb_0, pos_eb_c, pos_eb_a,  pos_eb_d, temp, s, a, epsilon, epsilon_0, C_vector, R, F,Z, idx_fix_species)
        #print(J)
        # Here the precondition techniques can be implemented
        # solve
        delta_ln_X = linalg.solve(J,-Y)
        #print(delta_ln_X)
        #update X
        #X = X*np.exp(delta_ln_X)
        ln_X = ln_X + delta_ln_X
        ln_C = mass_action_law (ln_X, ln_k, A)
        C = np.exp(ln_C)
        u = u_componentvector(A,C)
        
       # Vector_error = 
        # error
        d = u-T
        if idx_fix_species != None:
            d[idx_fix_species] =0
        abs_err = max(abs(d))     
        # Relaxation factor borrow from Craig M.Bethke to avoid negative values
        #max_1 = 1
        #max_2 =np.amax(-2*np.multiply(delta_ln_X, 1/ln_X))
        #Max_f = np.amax([max_1, max_2])
        #Del_mul = 1/Max_f
        #ln_X = Del_mul*delta_ln_X
        #ln_X = ln_X+delta_ln_X
        
        counter_iterations += 1
    if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.') 
    # things to do if goes well
    X = np.exp(ln_X)
    ln_C = mass_action_law (ln_X, ln_k, A)
    C = np.exp(ln_C)
    if debug_flm is not None:
        return X, C, debug_flm
    else:
        return X, C
    
    
############################## DEBUG CLASS ############################################################
    
class Debug_flm:
    def __init__(self):
        self.array_Jacobians = []
        self.array_residuals = []
        self.array_tolerance = []
        self.n_iterations = 0
    def append_Jacobian (self,a):
        self.array_Jacobians.append(a)
    def append_Residual (self,a):
        self.array_Jacobians.append(a)
    def append_tolerance (self,a):
        self.array_Jacobians.append(a)
    def inc_iteration(self):
        self.n_iterations += 1