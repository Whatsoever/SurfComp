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

def fourth_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_boltz0, pos_boltzalpha, pos_boltzbeta, pos_boltzgamma):
    """
    -The implementation of these algorithm is based on Westall (1980), but slightly modified in order to allow a 4th electrostatic layer.
        Arguments:
            - X_guess       A vector containing the initial guesses of the primary aqueous species, primary sorption species, electrostatic species
            - A             A matrix containing the stoichiometric values of the mass balance parameters
            - log_k         A vector of log(Konstant equilibrium)
            - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
            - Z             The vector of charge of the different ion. The order is determine by the rows of "A" for aqueous species. That means that it is link to idx_Aq somehow.
            - pos_boltz0, pos_boltzalpha, pos_boltzbeta, pos_boltzgamma  This is basically the position of the boltzman factor for the different planes (gamma == diffusive)
        Outpus:
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
    X = sp.optimize.fsolve(func_NR_FLM, X_guess, args = (), fprime = Jacobian_NR_FLM)
    #Speciation
    C = log_k + A*log(X)
    return X, C

def func_NR_FLM (X, ):
    """
        This function is supossed to be linked to the fourth_layer_one_surface_speciation function.
        It just gave the evaluated vector of Y, for the newton rapshon procedure.
        Westall (1980) look at Y
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
    T = Update_T_FLM(T, s,e, I, temp, a, Capacitances, psi_v, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma)
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

def Update_T_FLM(T, s,e, I, temp, a, C, psi_v, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma):
    """
        This equation is linked to func_NR_FLM. It updates the values of T for the electrostatic parameters.
        C       is the vector of capacitances. The units are supossed to be F/m²
        psi_v   is the vector of electrostatic potentials. The units are supposed to be volts
        s       is the specific surface area
        a       area
        temp    is the temperature in Kelvins
        e       is the relative permittivity
        I       is the ionic strenght mol/m³
    """
    # Faraday constant
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
    
def Jacobian_NR_FLM (X):
    '''
        This function should give the Jacobian. Here The jacobian is calculated as Westall (1980), except the electrostatic terms that are slightly different.
        The reason is because there seems to be some typos in Westall paper.
    '''
    # Speciation - mass action law
    log_C = log_k + A*np.log10(X)
    
    
    
    
        def Jacobian_Speciation_Westall1980_func (self, x, T_chem, pos_start_elec, pos_end_elec, S1, S2):
        log_c2 = np.matmul(linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(x)))      # Update secondary
        c2 =10**log_c2
        c_n = np.concatenate ((x, c2))
        return self.Jacobian_Speciation_Westall1980(c_n, pos_start_elec, pos_end_elec)
    
    
    def Jacobian_Speciation_Westall1980 (self, C, n_aq_plus_n_sorpt, n_primaryspecies):
        '''
            The jacobian matrix following an implementation based on the algorithm of  Westall (1980) 
            "Chemical equilibrium Including Adsorption on Charged Surfaces"
            Pages 37-to-39
            It is assumed that C is order first with the primary species and then with the secondary species such as C = [C1 C2]
        '''
        # The first part treats all terms as it was a normal speciation
        Z = np.zeros((n_primaryspecies, n_primaryspecies))
        for i in range(0, n_primaryspecies):
            for j in range(0, n_primaryspecies):
                Z[i,j]= np.matmul(np.multiply(self.U[i,:], self.U[j,:]), (C/C[j]))
    
        # According to the point 2 of Table III of Westall the term C*sa/F*RT/Funknwon must be added to the electrostatic part
        # I am supposing here that all the sorption phases are CCM