# -*- coding: utf-8 -*-
"""
Created on Fri May 24 15:17:13 2019

@author: DaniJ
"""

'''
It is like the four_layer_model_2try_withFixSpeciesOption_Scaling.py, but with two surfaces. There is not Poisson-Boltzman interaction between the two surfaces.
'''
import numpy as np
from scipy import linalg

def four_layer_two_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, idx_fix_species = None, zel=1, tolerance = 1e-6, max_iterations = 100,scalingRC = True):
    """
    -The implementation of these algorithm is based on Westall (1980), but slightly modified in order to allow a 4th electrostatic layer and 2 surface which its diffuse layers does not interact 
     Arguments:
            - T             A vector needed for creating the residual function for the Newthon-Raphson. The vector has the same size of X_guess and contains values like the total number of moles or mol/L of an aquoeus component
            - X_guess       A vector containing the initial guesses of the primary aqueous species, primary sorption species, electrostatic species
            - A             A matrix containing the stoichiometric values of the mass balance parameters
            - Z             The vector of charge of the different ion. The order is determine by the rows of "A" for aqueous species. That means that it is link to idx_Aq somehow.
            - log_k         A vector of log(Konstant equilibrium). Primary species of aquoues and sorption have a log_k=0
            - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
            - pos_psi_S1_vec Is a vector that contains the position of the boltzmann factor of each plane for surface 1 such as [pos_boltz0, pos_boltzalpha, pos_boltzbeta, pos_boltzgamma](gamma == diffusive of S1)
            - pos_psi_S2_vec It is like pos_psi_S2_vec but for the surface 2
            - temp           Temperature of the chemical system in Kelvins.
            - sS1             is the specific surface area for the surface 1
            - aS1             concentration of suspended solid for the surface 1
            - sS2             is the specific surface area for the surface 2
            - aS2             concentration of suspended solid for the surface 2
            - e             relative permittivity
            - CapacitancesS1  [C1, C2, C3] for surface 1
            - CapacitancesS2  [C1, C2, C3] for surface 2
            - scalingRC       If true a scaling stp will be done if false not scaling step is done (based on Marinoni et al. 2017) [Default = true]
            - idx_fix_species Index of the primary species that have a fixed value, it must coincide with X_guess.
    Outputs: the outputs right now are:
            - C the vector of species concentrations (aqueous and surface species, not electrostatic). The order of the species will depend on the given matrix A, so it is user dependent.
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
                        5) The vectors, and matrix are suppossed to be in a numpy 'format', due to the fact that we are using its libraries, it should be like that.
                        6) The plane gamma is place at the "same lcation" that the diffusion plane. SO, basically is the same.
    """
    counter_iterations = 0
    abs_err = tolerance + 1
    if idx_fix_species != None:
        X_guess [idx_fix_species] = T [idx_fix_species]
    while abs_err>tolerance and counter_iterations < max_iterations:
        # Calculate Y
        [Y, T] = func_NR_FLM (X_guess, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, idx_fix_species)
        # Calculate Z
        J = Jacobian_NR_FLM (X_guess, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, idx_fix_species)
        # Calculating the diff, Delta_X
        # Scaling technique is the RC technique from "Thermodynamic Equilibrium Solutions Through a Modified Newton Raphson Method"-Marianna Marinoni, Jer^ome Carrayrou, Yann Lucas, and Philippe Ackerer (2016)
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
        #print(delta_X))
        # Relaxation factor borrow from Craig M.Bethke to avoid negative values
        max_1 = 1
        max_2 =np.amax(-2*np.multiply(delta_X, 1/X_guess))
        Max_f = np.amax([max_1, max_2])
        Del_mul = 1/Max_f
        X_guess=X_guess + Del_mul*delta_X
        
        log_C = log_k + np.matmul(A,np.log10(X_guess))
        # transf
        C = 10**(log_C)
        u = np.matmul(A.transpose(),C)
        
       # Vector_error 
        d = u-T
        #print(d)
        if idx_fix_species != None:
            d[idx_fix_species] =0
        abs_err = max(abs(d))
        
        counter_iterations += 1
    if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.') 
    # Speciation - mass action law
    log_C = log_k + np.matmul(A,np.log10(X_guess))
    # transf
    C = 10**(log_C)
    return X_guess, C
    
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
    
    T = Update_T_FLM(T, sS1, sS2, e, I, temp, aS1, aS2, Z,CapacitancesS1, CapacitancesS2, psi_S1_v, psi_S2_v, zel, pos_psi_S1_vec, pos_psi_S2_vec, C_aq)
    # Calculation of Y
    Y= np.matmul(A.transpose(),C)-T
    # fix?
    if idx_fix_species != None:
        Y[idx_fix_species]=0
    return Y,T   
    
    
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
    I = np.matmul(np.multiply(Z,Z),C)
    I = I/2
    return I


def Update_T_FLM(T, sS1, sS2, e, I, temp, aS1, aS2, Z,CapacitancesS1, CapacitancesS2, psi_S1_v, psi_S2_v, zel, pos_psi_S1_vec, pos_psi_S2_vec, C_aq):
    """
        This equation is linked to func_NR_FLM. It updates the values of T for the electrostatic parameters.
         - All the arguments of the function have been stated in four_layer_two_surface_speciation function.
    """
    #  constant
    F = 96485.3328959                                   # C/mol
    R =  8.314472                                       # J/(K*mol)
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    #e  = 1.602176620898e-19                             # C
    kb = 1.38064852e-23                                 # J/K other units --> kb=8,6173303e-5  eV/K
    Na = 6.022140857e23                                 # 1/mol

    ########## S1  #####################
    sigma_S1_0 = CapacitancesS1[0]*(psi_S1_v[0]-psi_S1_v[1])
    sigma_S1_alpha = -sigma_S1_0 + CapacitancesS1[1]*(psi_S1_v[1]-psi_S1_v[2])
    sigma_S1_beta = -sigma_S1_0-sigma_S1_alpha+CapacitancesS1[2]*(psi_S1_v[2]-psi_S1_v[3])
    sigma_S1_gamma = -sigma_S1_0 - sigma_S1_alpha - sigma_S1_beta
    # Now the diffusive layer surface potential (sigma_d) is calculated. Using the formula given by Bethke in his book Geochemical Modeling Reactions
    
    sigma_S1_d = np.sqrt(8*1000*R*temp*eo*e*I)*np.sinh((zel*psi_S1_v[3]*F)/(2*R*temp))
    
     # T
    T_S1_0 = ((sS1*aS1)/F)*sigma_S1_0;                    # units mol/L or mol/kg
    T_S1_alpha = ((sS1*aS1)/F)*sigma_S1_alpha;            # units mol/L or mol/kg
    T_S1_beta = ((sS1*aS1)/F)*sigma_S1_beta;              # units mol/L or mol/kg
    #!! Important!!
    #T_gammad = ((s*a)/F)*(-sigma_gamma+sigma_d)             # This part should be equal to C[2]*(psi_beta-psi_dorgamma)+sigma_d
    T_S1_gammad = ((sS1*aS1)/F)*(sigma_S1_gamma+sigma_S1_d) 
    
    ########## S2  #####################
    sigma_S2_0 = CapacitancesS2[0]*(psi_S2_v[0]-psi_S2_v[1])
    sigma_S2_alpha = -sigma_S2_0 + CapacitancesS2[1]*(psi_S2_v[1]-psi_S2_v[2])
    sigma_S2_beta = -sigma_S2_0-sigma_S2_alpha+CapacitancesS2[2]*(psi_S2_v[2]-psi_S2_v[3])
    sigma_S2_gamma = -sigma_S2_0 - sigma_S2_alpha - sigma_S2_beta
    # Now the diffusive layer surface potential (sigma_d) is calculated. Using the formula given by Bethke in his book Geochemical Modeling Reactions
    
    sigma_S2_d = np.sqrt(8*1000*R*temp*eo*e*I)*np.sinh((zel*psi_S2_v[3]*F)/(2*R*temp))
    
     # T
    T_S2_0 = ((sS2*aS2)/F)*sigma_S2_0;                    # units mol/L or mol/kg
    T_S2_alpha = ((sS2*aS2)/F)*sigma_S2_alpha;            # units mol/L or mol/kg
    T_S2_beta = ((sS2*aS2)/F)*sigma_S2_beta;              # units mol/L or mol/kg
    #!! Important!!
    #T_gammad = ((s*a)/F)*(-sigma_gamma+sigma_d)             # This part should be equal to C[2]*(psi_beta-psi_dorgamma)+sigma_d
    T_S2_gammad = ((sS2*aS2)/F)*(sigma_S2_gamma+sigma_S2_d) 
    
    # Now the values must be put in T
    T[pos_psi_S1_vec[0]] = T_S1_0
    T[pos_psi_S1_vec[1]] = T_S1_alpha
    T[pos_psi_S1_vec[2]] = T_S1_beta
    T[pos_psi_S1_vec[3]] = T_S1_gammad
    
    T[pos_psi_S2_vec[0]] = T_S2_0
    T[pos_psi_S2_vec[1]] = T_S2_alpha
    T[pos_psi_S2_vec[2]] = T_S2_beta
    T[pos_psi_S2_vec[3]] = T_S2_gammad
    
    return T

def Jacobian_NR_FLM (X, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, idx_fix_species=None):
    '''
        This function should give the Jacobian. Here The jacobian is calculated as Westall (1980), except the electrostatic terms that are slightly different.
        The reason is because there seems to be some typos in Westall paper.
        Also, if idx_fix_species is given then the rows of the unknown will be 1 for the unknown and 0 for the other points.
    '''
    #  constant
    F = 96485.3328959                                   # C/mol [Faraday constant]
    R =  8.314472                                       # J/(K*mol) [universal constant gas]
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    # Speciation - mass action law
    #log_C = log_k + A*np.log10(X)
    log_C = log_k + np.matmul(A,np.log10(X))
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
    ############S1#######################
    sa_F2S1 = (sS1*aS1)/(F*F)
    C1_sa_F2_RTS1 = sa_F2S1*CapacitancesS1[0]*R*temp
    # Assigning in Jacobian (plane 0)
    Z[pos_psi_S1_vec[0],pos_psi_S1_vec[0]]=Z[pos_psi_S1_vec[0],pos_psi_S1_vec[0]] + C1_sa_F2_RTS1/X[pos_psi_S1_vec[0]]
    Z[pos_psi_S1_vec[0],pos_psi_S1_vec[1]]=Z[pos_psi_S1_vec[0],pos_psi_S1_vec[1]] - C1_sa_F2_RTS1/X[pos_psi_S1_vec[1]]
    #### plane alpha
    C1C2_sa_F2_RTS1 = sa_F2S1*R*temp*(CapacitancesS1[0]+CapacitancesS1[1])
    C2_sa_F2_RTS1 = sa_F2S1*CapacitancesS1[1]*R*temp
    # Assigning in Jacobian (plane alpha)
    Z[pos_psi_S1_vec[1],pos_psi_S1_vec[0]]=Z[pos_psi_S1_vec[1],pos_psi_S1_vec[0]] - C1_sa_F2_RTS1/X[pos_psi_S1_vec[0]]
    Z[pos_psi_S1_vec[1],pos_psi_S1_vec[1]]=Z[pos_psi_S1_vec[1],pos_psi_S1_vec[1]] + C1C2_sa_F2_RTS1/X[pos_psi_S1_vec[1]]
    Z[pos_psi_S1_vec[1],pos_psi_S1_vec[2]]= Z[pos_psi_S1_vec[1],pos_psi_S1_vec[2]] - C2_sa_F2_RTS1/X[pos_psi_S1_vec[2]]
    #### plane beta
    C3C2_sa_F2_RTS1 = sa_F2S1*R*temp*(CapacitancesS1[1]+CapacitancesS1[2])
    C3_sa_F2_RTS1 = sa_F2S1*CapacitancesS1[2]*R*temp
    # Assigning in Jacobian (plane beta)
    Z[pos_psi_S1_vec[2],pos_psi_S1_vec[1]] = Z[pos_psi_S1_vec[2],pos_psi_S1_vec[1]] - C2_sa_F2_RTS1/X[pos_psi_S1_vec[1]]
    Z[pos_psi_S1_vec[2], pos_psi_S1_vec[2]] = Z[pos_psi_S1_vec[2],pos_psi_S1_vec[2]] + C3C2_sa_F2_RTS1/X[pos_psi_S1_vec[2]]
    Z[pos_psi_S1_vec[2], pos_psi_S1_vec[3]] = Z[pos_psi_S1_vec[2],pos_psi_S1_vec[3]] - C3_sa_F2_RTS1/X[pos_psi_S1_vec[3]]
    #### plane gamma [diffusive plane]
    Z[pos_psi_S1_vec[3],pos_psi_S1_vec[2]] = Z[pos_psi_S1_vec[3],pos_psi_S1_vec[2]] - C3_sa_F2_RTS1/X[pos_psi_S1_vec[2]]  
    # d_d plane
    psi_d = Boltzman_factor_2_psi(X[pos_psi_S1_vec[3]], temp)
    DY_Dpsid = -np.sqrt(8*1000*R*temp*e*eo*I)*np.cosh((zel*F*psi_d)/(2*R*temp))*((zel*F)/(2*R*temp)) - CapacitancesS1[2]
    Dpsid_DpsidB = (-R*temp)/(F*X[pos_psi_S1_vec[3]])
    Z[pos_psi_S1_vec[3], pos_psi_S1_vec[3]] = Z[pos_psi_S1_vec[3], pos_psi_S1_vec[3]] + (DY_Dpsid*Dpsid_DpsidB*((sS1*aS1)/F))

#(Problably S1 and S2 can be enclosed in a for loop, reducing lines of code. If I have time and will, I will look at it.)
    ############S1#######################
    sa_F2S2 = (sS2*aS2)/(F*F)
    C1_sa_F2_RTS2 = sa_F2S2*CapacitancesS2[0]*R*temp
    # Assigning in Jacobian (plane 0)
    Z[pos_psi_S2_vec[0],pos_psi_S2_vec[0]]=Z[pos_psi_S2_vec[0],pos_psi_S2_vec[0]] + C1_sa_F2_RTS2/X[pos_psi_S2_vec[0]]
    Z[pos_psi_S2_vec[0],pos_psi_S2_vec[1]]=Z[pos_psi_S2_vec[0],pos_psi_S2_vec[1]] - C1_sa_F2_RTS2/X[pos_psi_S2_vec[1]]
    #### plane alpha
    C1C2_sa_F2_RTS2 = sa_F2S2*R*temp*(CapacitancesS2[0]+CapacitancesS2[1])
    C2_sa_F2_RTS2 = sa_F2S2*CapacitancesS2[1]*R*temp
    # Assigning in Jacobian (plane alpha)
    Z[pos_psi_S2_vec[1],pos_psi_S2_vec[0]]=Z[pos_psi_S2_vec[1],pos_psi_S2_vec[0]] - C1_sa_F2_RTS2/X[pos_psi_S2_vec[0]]
    Z[pos_psi_S2_vec[1],pos_psi_S2_vec[1]]=Z[pos_psi_S2_vec[1],pos_psi_S2_vec[1]] + C1C2_sa_F2_RTS2/X[pos_psi_S2_vec[1]]
    Z[pos_psi_S2_vec[1],pos_psi_S2_vec[2]]= Z[pos_psi_S2_vec[1],pos_psi_S2_vec[2]] - C2_sa_F2_RTS2/X[pos_psi_S2_vec[2]]
    #### plane beta
    C3C2_sa_F2_RTS2 = sa_F2S2*R*temp*(CapacitancesS2[1]+CapacitancesS2[2])
    C3_sa_F2_RTS2 = sa_F2S2*CapacitancesS2[2]*R*temp
    # Assigning in Jacobian (plane beta)
    Z[pos_psi_S2_vec[2],pos_psi_S2_vec[1]] = Z[pos_psi_S2_vec[2],pos_psi_S2_vec[1]] - C2_sa_F2_RTS2/X[pos_psi_S2_vec[1]]
    Z[pos_psi_S2_vec[2], pos_psi_S2_vec[2]] = Z[pos_psi_S2_vec[2],pos_psi_S2_vec[2]] + C3C2_sa_F2_RTS2/X[pos_psi_S2_vec[2]]
    Z[pos_psi_S2_vec[2], pos_psi_S2_vec[3]] = Z[pos_psi_S2_vec[2],pos_psi_S2_vec[3]] - C3_sa_F2_RTS2/X[pos_psi_S2_vec[3]]
    #### plane gamma [diffusive plane]
    Z[pos_psi_S2_vec[3],pos_psi_S2_vec[2]] = Z[pos_psi_S2_vec[3],pos_psi_S2_vec[2]] - C3_sa_F2_RTS2/X[pos_psi_S2_vec[2]]  
    # d_d plane
    psi_dS2 = Boltzman_factor_2_psi(X[pos_psi_S2_vec[3]], temp)
    DY_Dpsid = -np.sqrt(8*1000*R*temp*e*eo*I)*np.cosh((zel*F*psi_dS2)/(2*R*temp))*((zel*F)/(2*R*temp)) - CapacitancesS2[2]
    Dpsid_DpsidB = (-R*temp)/(F*X[pos_psi_S2_vec[3]])
    Z[pos_psi_S2_vec[3], pos_psi_S2_vec[3]] = Z[pos_psi_S2_vec[3], pos_psi_S2_vec[3]] + (DY_Dpsid*Dpsid_DpsidB*((sS2*aS2)/F))

    # finally just return Z
    if idx_fix_species != None:
        for d in idx_fix_species:
            v=np.zeros(length_X)
            v[d]=1
            Z[d,:] = v
    return Z
 

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