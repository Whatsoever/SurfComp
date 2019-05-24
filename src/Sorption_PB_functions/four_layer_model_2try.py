# -*- coding: utf-8 -*-
"""
Created on Fri May 24 15:17:13 2019

@author: DaniJ
"""

'''
Firt version of the four layer model seem to not work because the appeareance of negative values of the poisson boltzman factor. The next implementation might have the same problems, but the later one,
was using a modified Maxwell method. Here, we use a rough Newton-Rapshon and probably will be adepted in the future.
'''
import numpy as np
from scipy import linalg

def four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-6, max_iterations = 100):
    """
    -The implementation of these algorithm is based on Westall (1980), but slightly modified in order to allow a 4th electrostatic layer.
        Arguments:
            - T             A vector needed for creating the residual function for the Newthon-Raphson. The vector has the same size of X_guess and contains values like the total number of moles or mol/L of an aquoeus component
            - X_guess       A vector containing the initial guesses of the primary aqueous species, primary sorption species, electrostatic species
            - A             A matrix containing the stoichiometric values of the mass balance parameters
            - log_k         A vector of log(Konstant equilibrium). Primary species of aquoues and sorption have a log_k=0
            - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
            - Z             The vector of charge of the different ion. The order is determine by the rows of "A" for aqueous species. That means that it is link to idx_Aq somehow.
            - pos_boltz0, pos_boltzalpha, pos_boltzbeta, pos_boltzgamma  This is basically the position of the boltzman factor for the different planes (gamma == diffusive)
            - s             is the specific surface area
            - a             concentration of suspended solid
            - e             relative permittivity
            - Capacitances  [C1, C2, C3] 
            - temp           Temperature of the chemical system in Kelvins.
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
        Note: There are extra arguments related to the fsolve method of scipy.optimize class. I have not add them, they can be add later.
        The vectors, and matrix are suppossed to be in a numpy 'format', due to the fact that we are using its libraries, it should be like that.
        
        The plane gamma is place at the "same lcation" that the diffusion plane. SO, basically is the same.
    """
    # Speciation
    #C_guess = log_k + A*log(X_guess)

    # instantiation variables for loop
    counter_iterations = 0
    err = tolerance + 1
    while err>tolerance and counter_iterations < max_iterations:
        # Calculate Y
        Y = func_NR_FLM (X_guess, A, log_k, temp, idx_Aq, s, a, e, Capacitances, T, Z, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma)
        # Calculate Z
        J = Jacobian_NR_FLM (X_guess, A, log_k, temp, idx_Aq, s, a, e, Capacitances, T, Z, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma)
        # Calculating the diff, Delta_X
        delta_X = linalg.solve(J,-Y)
        # The error will be equal to the maximum increment
        err = max(abs(delta_X))
        # Relaxation factor borrow from Craig M.Bethke to avoid negative values
        max_1 = 1
        max_2 =np.amax(-2*np.multiply(delta_X, 1/X_guess))
        Max_f = np.amax([max_1, max_2])
        Del_mul = 1/Max_f
        X_guess=X_guess + Del_mul*delta_X
        counter_iterations += 1
    if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.') 
    # Speciation - mass action law
    log_C = log_k + np.matmul(A,np.log10(X_guess))
    # transf
    C = 10**(log_C)
    return X_guess, C
    

def func_NR_FLM (X, A, log_k, temp, idx_Aq, s, a, e, Capacitances, T, Z, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma):
    """
        This function is supossed to be linked to the fourth_layer_one_surface_speciation function.
        It just gave the evaluated vector of Y, for the Newton-raphson procedure.
        The formulation of Westall (1980) is followed.
        FLM = four layer model
    """
    # Speciation - mass action law
    log_C = log_k + np.matmul(A,np.log10(X))
    # transf
    C = 10**(log_C)
    # Update T - "Electrostatic parameters"
    psi_v = [Boltzman_factor_2_psi(X[pos_psi0], temp), Boltzman_factor_2_psi(X[pos_psialpha], temp), Boltzman_factor_2_psi(X[pos_psibeta], temp), Boltzman_factor_2_psi(X[pos_psigamma], temp)]
    C_aq = C[idx_Aq]
    I = Calculate_ionic_strength(Z, C_aq)
    T = Update_T_FLM(T, s, e, I, temp, a, Z,Capacitances, psi_v, zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma, C_aq)
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
    I = np.matmul(np.multiply(Z,Z),C)
    I = I/2
    return I

def Update_T_FLM(T, s,e, I, temp, a, Z,C, psi_v,zel, pos_psi0, pos_psialpha, pos_psibeta, pos_psigamma, C_aq):
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
    #e  = 1.602176620898e-19                             # C
    kb = 1.38064852e-23                                 # J/K other units --> kb=8,6173303e-5  eV/K
    Na = 6.022140857e23                                 # 1/mol
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
    
    # Now the diffusive layer surface potential (sigma_d) is calculated. Due to electroneutrality conditions (Gaus Boundary electrostatic) sigma_gamma = sigma_d
    # The equation for calculatin the surface potential is equation (46) of Chapter 3 "Diffuse double layer equations for use in surface complexation models: Approximations and limits" Hiroyuki Ohshima
    # The book is called Surface Complexation Models of Adsorption and author is Johannes Lützenkirchen
    
    sigma_d = np.sqrt(8*1000*R*temp*eo*e*I)*np.sinh((zel*psi_v[3]*F)/(2*R*temp))
    
    # Calculating Debye-Huckel paramerter (https://de.wikipedia.org/wiki/Debye-H%C3%BCckel-Theorie)
    #k = F*np.sqrt((2*I*1000)/(e*eo*R*temp))
    #yd= (psi_v[3]*F)/(R*temp)
    #if yd <0:
      #  sign=-1
    #else:
     #   sign=1
    #ni = C_aq*Na*1000
    #partA = np.sum(np.multiply(ni,(np.exp(-yd*Z)-1)))
    #partB = np.sum(np.multiply(ni,np.multiply(Z,Z)))
    #print(partA)
    #print(partB)
    #sigma_d = -sign*k*np.sqrt(partA/partB)
    
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
    dif_term = ((zel*F)/(2*R*temp))*np.sqrt(8*R*1000*temp*eo*e*I)*np.cosh((zel*Boltzman_factor_2_psi(X[pos_psigamma], temp)*F)/(2*R*temp))
    dif_term = dif_term+Capacitances[2]
    dif_term = dif_term*((-R*temp)/F)
    # Assigning in Jacobian (plane beta)
    Z[pos_psigamma,pos_psibeta] = Z[pos_psigamma, pos_psibeta] - gb_term/X[pos_psibeta]
    Z[pos_psigamma, pos_psigamma] = Z[pos_psigamma, pos_psigamma] +  dif_term/X[pos_psigamma] #Z[pos_bgamma, pos_bgamma] should be equal to  0
    # finally just return Z
    return Z