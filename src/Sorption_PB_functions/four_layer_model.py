# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 15:28:36 2019

@author: 
"""


import numpy as np
import scipy as sp


"""
Implementation of the fourth layer model with just one surface. In a function manner,
allowing its utilization in other codes and this code.
"""


def four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp,  zel=1):
    """
    -The implementation of these algorithm is based on Westall (1980), but slightly modified in order to allow a 4th electrostatic layer.
        Arguments:
            - X_guess       A vector containing the initial guesses of the primary aqueous species, primary sorption species, electrostatic species
            - A             A matrix containing the stoichiometric values of the mass balance parameters
            - log_k         A vector of log(Konstant equilibrium)
            - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
            - Z             The vector of charge of the different ion. The order is determine by the rows of "A" for aqueous species. That means that it is link to idx_Aq somehow.
            - pos_boltz0, pos_boltzalpha, pos_boltzbeta, pos_boltzgamma  This is basically the position of the boltzman factor for the different planes (gamma == diffusive)
        Outputs: the outputs right now are:
                        - C the vector of species concentrations (aqueous and surface species). The order of the species will depend on the given matrix A, so it is user dependent.
                        - The vector X of primary unknowns. The value of the primary species of aqueous and surface species should be equivalent to the C vector. Here we can find the values
                          of the boltzman factors, which are related to psi values.
        Preconditions: 
                        1) The order of the rows of matrix "A" must agree with the order of the unknowns in the vector X_guess.
                           Namely, if the first row correspond to the species "H+", the first unknow in X_guess must be "H+"
                           This also implies that number of rows of A equals the length of the vector of unknows.
                        2) Since the order of the species is not known, the positions in the "X_guess" of  the electrostatic species
                           is needed to update vector "T", and also the Jacobian matrix.
                        3) It is also assumed that T has the same order than X_guess. Namely, if the first components is "H+" in "T", 
                           it should also be in "H+" in "X_guess".
                        4) log_k is the vector of the logarithm (equilibrium constant). For each species a logK is given, if the species
                           is a primary species, the value would be zero. The K must be coherent with matrix A.
        Note: There are extra arguments related to the fsolve method of scipy.optimize class. I have not add them, they can be add later.
        The vectors, and matrix are suppossed to be in a numpy 'format', due to the fact that we are using its libraries. It should be like that.
        
        The plane gamma is place at the "same lcation" that the diffusion plane.
    """
    # Speciation
    #C_guess = log_k + A*log(X_guess)

    # scipy.optimize.fsolve(func, x0, args=(), fprime=None, full_output=0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)[source]¶
    X = sp.optimize.fsolve(func_NR_FLM, X_guess, args = (A, log_k, temp, idx_Aq, s, a, e, Capacitances, T, Z, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma), fprime = Jacobian_NR_FLM)
    #Speciation
    C = log_k + A*np.log10(X)
    return X, C

def func_NR_FLM (X, A, log_k, temp, idx_Aq, s, a, e, Capacitances, T, Z, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma):
    """
        This function is supossed to be linked to the fourth_layer_one_surface_speciation function.
        It just gave the evaluated vector of Y, for the Newton-raphson procedure.
        The formulation of Westall (1980) is followed.
        FLM = four layer model
    """
    # Speciation - mass action law
    log_C = log_k + A*np.log10(X)
    # transf
    C = 10**(log_C)
    # Update T - "Electrostatic parameters"
    psi_v = [Boltzman_factor_2_psi(X[pos_psi0], temp), Boltzman_factor_2_psi(X[pos_psialpha], temp), Boltzman_factor_2_psi(X[pos_psibeta], temp), Boltzman_factor_2_psi(X[pos_psigamma], temp)]
    C_aq = C[idx_Aq]
    I = Calculate_ionic_strength(Z, C_aq)
    T = Update_T_FLM(T, s, e, I, temp, a, Capacitances, psi_v, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma)
    # Calculation of Y
    Y= np.matmul(A.transpose(),C)-T
    return Y

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

def Calculate_ionic_strength(Z,C):
    '''
        It is supossed to be numpy format vector 
        Z is the vector of charge 
    '''
    # Multiplication must be pointwise for the vector
    # multiply function of numpy. Multiplies pointwise according to the documentation and own experience.
    I = np.multiply(np.multiply(Z,Z),C)
    I = I/2
    return I

def Update_T_FLM(T, s,e, I, temp, a, C, psi_v,zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma):
    """
        This equation is linked to func_NR_FLM. It updates the values of T for the electrostatic parameters.
        C       is the vector of capacitances. The units are supossed to be F/m². Do not mistake it with the species vector.
        psi_v   is the vector of electrostatic potentials. The units are supposed to be volts
        s       is the specific surface area
        a       area
        temp    is the temperature in Kelvins
        e       is the relative permittivity
        I       is the ionic strenght mol/m³
        zel     is the background electrolyte value. It is needed for the PB analytical solution.
    """
    #  constant
    F = 96485.3328959                                   # C/mol
    R =  8.314472                                       # J/(K*mol)
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    # Update of T
    # T = (sa/F)*sigma
    # sigma = C*(psi-psi_0) <-- The sigma depends on the plane
    
    
    # NOTE: I JUST REALISED I DID SOMETHING WRONG THE WHOLE TIME. In the T vector, It must be checked in the other approaches (TLM)
    # The equations os the layers that are not in contact with the diffusion layer must be mol/kg or mol/l
    # while the third component is the total charge balance.
    
    sigma_0 = C[0]*(psi_v[0]-psi_v[1])
    sigma_alpha = -sigma_0 + C[1]*(psi_v[1]-psi_v[2])
    sigma_beta = -sigma_0-sigma_alpha+C[2]*(psi_v[2]-psi_v[3])
    sigma_gamma = -sigma_0 - sigma_alpha - sigma_beta
    
    sigma_d = np.sqrt(8*R*temp*eo*e*I)*np.sinh((zel*psi_v[3]*F)/(2*R*temp))
    
    # T
    T_0 = ((s*a)/F)*sigma_0;                    # units mol/L or mol/kg
    T_alpha = ((s*a)/F)*sigma_alpha;            # units mol/L or mol/kg
    T_beta = ((s*a)/F)*sigma_beta;              # units mol/L or mol/kg
    #!! Important!!
    T_gammad = -sigma_gamma+sigma_d             # This part should be equal to C[2]*(psi_beta-psi_dorgamma)+sigma_d
    
    # Now the values must be put in T
    T[pos_psi0] = T_0
    T[pos_psialpha] = T_alpha
    T[pos_psibeta] = T_beta
    T[pos_psigamma] = T_gammad
    
    return T
    
def Jacobian_NR_FLM (X, A, log_k, temp, idx_Aq, s, a, e, Capacitances, T, Z, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma):
    '''
        This function should give the Jacobian. Here The jacobian is calculated as Westall (1980), except the electrostatic terms that are slightly different.
        The reason is because there seems to be some typos in Westall paper.
    '''
    #  constant
    F = 96485.3328959                                   # C/mol [Faraday constant]
    R =  8.314472                                       # J/(K*mol) [universal constant gas]
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    # Speciation - mass action law
    log_C = log_k + A*np.log10(X)
    # transf
    C = 10**(log_C)
    C_aq = C[idx_Aq]
    I = Calculate_ionic_strength(Z, C_aq)
    # instantiate Jacobian
    length_X = X.size
    Z = np.zeros((length_X,length_X))
    # First part is the common of the Jacbian derivation
    for i in range(0, length_X):
            for j in range(0, length_X):
                Z[i,j]= np.matmul(np.multiply(A[:,i], A[:,j]), (C/X[j]))
    # Now the electrostatic part must be modified, one question hang on the air:
    # Should we check that the electrostatic part is as we expected?
    sa_F2 = (s*a)/(F*F)
    #### plane 0
    C1_sa_F2_RT = sa_F2*Capacitances[0]*R*temp
    # Assigning in Jacobian (plane 0)
    Z[pos_psi0,pos_psi0]=Z[pos_psi0,pos_psi0] + C1_sa_F2_RT/X[pos_psi0]
    Z[pos_psi0,pos_psialpha]=Z[pos_psi0, pos_psialpha] - C1_sa_F2_RT/X[pos_psialpha]
    
    #### plane alpha
    C1C2_sa_F2_RT = sa_F2*R*temp*(Capacitances[0]+Capacitances[1])
    C2_sa_F2_RT = sa_F2*Capacitances[1]*R*temp
    # Assigning in Jacobian (plane alpha)
    Z[pos_psialpha,pos_psi0]=Z[pos_psialpha, pos_psi0] - C1_sa_F2_RT/X[pos_psi0]
    Z[pos_psialpha,pos_psialpha]=Z[pos_psialpha, pos_psialpha] + C1C2_sa_F2_RT/X[pos_psialpha]
    Z[pos_psialpha,pos_psibeta]= Z[pos_psialpha,pos_psibeta] - C2_sa_F2_RT/X[pos_psibeta]
    
    #### plane beta
    C3C2_sa_F2_RT = sa_F2*R*temp*(Capacitances[1]+Capacitances[2])
    C3_sa_F2_RT = sa_F2*Capacitances[2]*R*temp
    # Assigning in Jacobian (plane beta)
    Z[pos_psibeta,pos_psialpha] = Z[pos_psibeta,pos_psialpha] - C2_sa_F2_RT/X[pos_psialpha]
    Z[pos_psibeta, pos_psibeta] = Z[pos_psibeta, pos_psibeta] + C3C2_sa_F2_RT/X[pos_psibeta]
    Z[pos_psibeta, pos_psigamma] = Z[pos_psibeta, pos_psigamma] - C3_sa_F2_RT/X[pos_psigamma]
    
    #### plane gamma [diffusive plane]
    gb_term = (R*temp*Capacitances[2])/F
    # PB solution of diffusive layer
    # F/2RT (8RTε_o εI)^(1/2) cosh((Fψ_d)/2RT)
    dif_term = ((zel*F)/(2*R*temp))*np.sqrt(8*R*temp*eo*e*I)*np.cosh((zel*Boltzman_factor_2_psi(X[pos_psigamma], temp)*F)/(2*R*temp))
    dif_term = dif_term+Capacitances[2]
    dif_term = dif_term*((-R*temp)/F)
    # Assigning in Jacobian (plane beta)
    Z[pos_psigamma,pos_psibeta] = Z[pos_psigamma, pos_psibeta] - gb_term/X[pos_psibeta]
    Z[pos_psigamma, pos_psigamma] = Z[pos_psigamma, pos_psigamma] +  dif_term/X[pos_psigamma] #Z[pos_bgamma, pos_bgamma] should be equal to  0
    # finally just return Z
    return Z

