# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:58:31 2019

@author: DaniJ
"""

import numpy as np
import scipy as sp
from bvp import solve_bvp
from scipy import linalg
#from four_layer_model_2try_withFixSpeciesOption_Scaling import four_layer_model_2try_withFixSpeciesOption_Scaling as flm
from matplotlib import pyplot as plt
'''
    In this first try we will assume that the vector of unknowns is composed in the following order:
'''


def PB_and_fourlayermodel (T, X_guess, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, d0,df, idx_fix_species = None, zel=1, tolerance_NR = 1e-6, max_iterations = 100,scalingRC = True, tol_PB=1e-6):
    counter_iterations = 0
    abs_err = tolerance_NR + 1
    if idx_fix_species != None:
        X_guess [idx_fix_species] = T [idx_fix_species]
        tempv = (np.linspace(d0,df,100));
    y_0 = np.zeros((2,tempv.shape[0]))
    bvp_class = [y_0]
    while abs_err>tolerance_NR and counter_iterations < max_iterations:
        # Calculate Y
        [Y, T, bvp_class] = func_NR_FLM (X_guess, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, d0,df, bvp_class, idx_fix_species,tol_PB)
        # Calculate Z
        J = Jacobian_NR_FLM (X_guess, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, d0,df, bvp_class, idx_fix_species,tol_PB)
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
        #print(X_guess)
        Xmod=X_guess.copy()
        for i in range(len(X_guess)):
            if X_guess[i]<=0:
                Xmod[i]=1
        log_C = log_k + np.matmul(A,np.log10(Xmod))
        
        # transf
        C = 10**(log_C)
        u = np.matmul(A.transpose(),C)
        #print(C)
       # Vector_error 
        d = u-T
        print(d)
        if idx_fix_species != None:
            d[idx_fix_species] =0
        abs_err = max(abs(d))
        
        counter_iterations += 1
    if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.') 
            #return X_guess, C
    # Speciation - mass action law
    Xmod=X_guess.copy()
    for i in range(len(X_guess)):
        if X_guess[i]<=0:
            Xmod[i]=1
    log_C = log_k + np.matmul(A,np.log10(Xmod))
    # transf
    C = 10**(log_C)
    return X_guess, C

def func_NR_FLM (X, A, log_k, temp, idx_Aq,  sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Z, zel, pos_psi_S1_vec, pos_psi_S2_vec, d0,df,  bvp_solver, idx_fix_species=None, tol_PB=1e-6):
    """
        This function is supossed to be linked to the four_layer_two_surface_speciation function.
        It just gave the evaluated vector of Y, and T for the Newton-raphson procedure.
        The formulation of Westall (1980) is followed.
        FLM = four layer model
    """
    # Speciation - mass action law
    Xmod=X.copy()
    for i in range(len(X)):
        if X[i]<=0:
            Xmod[i]=1
    log_C = log_k + np.matmul(A,np.log10(Xmod))
    # transf
    C = 10**(log_C)
    # Update T - "Electrostatic parameters"
    'Notice that the last term of the psi_S1/2_v is not transformed from the boltzmann factor to the electrostatic potential!!!!!!!!!!!!'
    psi_S1_v = [Boltzman_factor_2_psi(X[pos_psi_S1_vec[0]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[1]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp), X[pos_psi_S1_vec[3]]]
    psi_S2_v = [Boltzman_factor_2_psi(X[pos_psi_S2_vec[0]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[1]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp), X[pos_psi_S2_vec[3]]] 
    C_aq = C[idx_Aq]
    #
    [T, y0] = Update_T_FLM(T, sS1, sS2, e, temp, aS1, aS2, Z,CapacitancesS1, CapacitancesS2, psi_S1_v, psi_S2_v, zel, pos_psi_S1_vec, pos_psi_S2_vec, C_aq, d0,df, bvp_solver, tol_PB)
    
    # Calculation of Y
    Y= np.matmul(A.transpose(),C)-T
    return Y, T, y0


def Update_T_FLM(T, sS1, sS2, e, temp, aS1, aS2, Z,CapacitancesS1, CapacitancesS2, psi_S1_v, psi_S2_v, zel, pos_psi_S1_vec, pos_psi_S2_vec, C_aq, d0,df, bvp_solver, tol_PB):
    
    #  constant
    F = 96485.3328959                                   # C/mol
    R =  8.314472                                       # J/(K*mol)
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    #e  = 1.602176620898e-19                             # C
    kb = 1.38064852e-23                                 # J/K other units --> kb=8,6173303e-5  eV/K
    Na = 6.022140857e23                                 # 1/mol
    elec_charge = 1.60217662e-19 #electron charge in C
    ########## S1  #####################
    sigma_S1_0 = CapacitancesS1[0]*(psi_S1_v[0]-psi_S1_v[1])
    sigma_S1_alpha = -sigma_S1_0 + CapacitancesS1[1]*(psi_S1_v[1]-psi_S1_v[2])
    sigma_S1_beta = -sigma_S1_0-sigma_S1_alpha+CapacitancesS1[2]*(psi_S1_v[2]-psi_S1_v[3])
    sigma_S1_gamma = -sigma_S1_0 - sigma_S1_alpha - sigma_S1_beta
    ########## S2  #####################
    sigma_S2_0 = CapacitancesS2[0]*(psi_S2_v[0]-psi_S2_v[1])
    sigma_S2_alpha = -sigma_S2_0 + CapacitancesS2[1]*(psi_S2_v[1]-psi_S2_v[2])
    sigma_S2_beta = -sigma_S2_0-sigma_S2_alpha+CapacitancesS2[2]*(psi_S2_v[2]-psi_S2_v[3])
    sigma_S2_gamma = -sigma_S2_0 - sigma_S2_alpha - sigma_S2_beta
     ########## T  S1  #####################
    T_S1_0 = ((sS1*aS1)/F)*sigma_S1_0;                    # units mol/L or mol/kg
    T_S1_alpha = ((sS1*aS1)/F)*sigma_S1_alpha;            # units mol/L or mol/kg
    T_S1_beta = ((sS1*aS1)/F)*sigma_S1_beta;              # units mol/L or mol/kg
    ########## T  S2  #####################
    T_S2_0 = ((sS2*aS2)/F)*sigma_S2_0;                    # units mol/L or mol/kg
    T_S2_alpha = ((sS2*aS2)/F)*sigma_S2_alpha;            # units mol/L or mol/kg
    T_S2_beta = ((sS2*aS2)/F)*sigma_S2_beta;              # units mol/L or mol/kg
    
    ################## PB part starts heres ########################################################################
    ew = eo*e
    Q = Z*elec_charge                                         # Q is the charve of the aqueous elements times the electron charge 
    C = C_aq
    A =Na* 1000/ew                # a prefactor = Avogadro * 1000 /ew
    #A = Na/ew
    kbt = 1.38064852e-23 *temp # kb (J/K) * T in K
    #y0 = np.zeros((2, x.size))
    if type(bvp_solver) == list:
        y_0 = bvp_solver[0].copy()
        x = np.linspace(d0,df,y_0.shape[1])
        #deltaing this
        #y1= bvp_solver[0].copy()
        #y2= bvp_solver[0].copy()
        #y3= bvp_solver[0].copy()
    else :
        y_0 = bvp_solver.y
        x = bvp_solver.x
    'I think that  y0[1,0] and y0[1,-1] are not necessary to solve the problem, I would say that its values do not have implications. Although I am not 100% sure.'
    y_0[1,0] = sigma_S1_gamma/ew           # The negative value that I am given here is extremely arbitrary I am not sure why. IT MUST BE DISCUSSED
#    dpsi_d = -(sig_0 + sig_b + sig_d)/ew   # electric field at diffuse layer, x>d
    y_0[1,-1] = -sigma_S2_gamma/ew             
    y_0[0,0] = psi_S1_v[3]
    y_0[0,-1]= psi_S2_v[3]
    #y0[1,-1]= sigma_S2_gamma/ew
    args=[Q,C,A,kbt,y_0]
    
    result = solve_bvp(fun_PB, bc_PB, x, y_0,  args = args, tol=tol_PB)
    # Looking further
    #y1[0,0] = psi_S1_v[0]
    #y1[0,-1] = psi_S2_v[0]
    #y1[1,0] = sigma_S1_0/ew
    #y1[1,-1] = -sigma_S2_0/ew
    #args=[Q,C,A,kbt,y1]
    
    #result1 = solve_bvp(fun_PB, bc_PB, x, y1,  args = args, tol=tol_PB)
    #y2[0,0] = psi_S1_v[1]
    #y2[0,-1] = psi_S2_v[1]
    #y2[1,0] = sigma_S1_alpha/ew
    #y2[1,-1] = -sigma_S2_alpha/ew
    #args=[Q,C,A,kbt,y2]
    
    #result2 = solve_bvp(fun_PB, bc_PB, x, y2,  args = args, tol=tol_PB)
    #y3[0,0] = psi_S1_v[2]
    #y3[0,-1] = psi_S2_v[2]
    #y3[1,0] = sigma_S1_beta/ew
    #y3[1,-1] = -sigma_S2_beta/ew
    #args=[Q,C,A,kbt,y3]
    
    #result3 = solve_bvp(fun_PB, bc_PB, x, y3,  args = args, tol=tol_PB)
    
    #si1=-result1.y[1][0]*ew
    #si2=-result2.y[1][0]*ew
    #si3=-result3.y[1][0]*ew
    
    plt.figure(3)
    plt.plot(result.x, result.y[0])
    #
    #assert 5==3
    sigma_S1_d=-result.y[1][0]*ew
    sigma_S2_d=result.y[1][-1]*ew
    #sigma_S1_d=-result.y[1][0]*ew + CapacitancesS1[2]*(psi_S1_v[2]-psi_S1_v[3])
    #sigma_S2_d=result.y[1][-1]*ew + CapacitancesS2[2]*(psi_S2_v[2]-psi_S2_v[3])
    #
    T_S1_gammad = sigma_S1_gamma+sigma_S1_d
    T_S2_gammad = sigma_S2_gamma+sigma_S2_d
    #
    #print([sigma_S1_gamma, sigma_S2_gamma, sigma_S1_d, sigma_S2_d])
    #print([psi_S1_v[3], psi_S2_v[3]])
     # Now the values must be put in T
    T[pos_psi_S1_vec[0]] = T_S1_0
    T[pos_psi_S1_vec[1]] = T_S1_alpha
    T[pos_psi_S1_vec[2]] = T_S1_beta
    T[pos_psi_S1_vec[3]] = T_S1_gammad
    
    T[pos_psi_S2_vec[0]] = T_S2_0
    T[pos_psi_S2_vec[1]] = T_S2_alpha
    T[pos_psi_S2_vec[2]] = T_S2_beta
    T[pos_psi_S2_vec[3]] = T_S2_gammad
    
    return T, result

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



def Jacobian_NR_FLM (X, A, log_k, temp, idx_Aq, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, T, Zi, zel, pos_psi_S1_vec, pos_psi_S2_vec, d0,df,  bvp_class, idx_fix_species=None,tol_PB=1e-6):
    '''
        This function should give the Jacobian. Here The jacobian is calculated as Westall (1980), except the electrostatic terms that are slightly different.
        The reason is because there seems to be some typos in Westall paper.
        Also, if idx_fix_species is given then the rows of the unknown will be 1 for the unknown and 0 for the other points.
    '''
    #  constant
    F = 96485.3328959                                   # C/mol [Faraday constant]
    R =  8.314472                                       # J/(K*mol) [universal constant gas]
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    elec_charge = 1.60217662e-19 #electron charge in C
    # Speciation - mass action law
    #log_C = log_k + A*np.log10(X)
    Xmod=X.copy()
    for i in range(len(X)):
        if X[i]<=0:
            Xmod[i]=1
    log_C = log_k + np.matmul(A,np.log10(Xmod))
    # transf
    C = 10**(log_C)
    C_aq = C[idx_Aq]
    #I = Calculate_ionic_strength(Z, C_aq)
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
    #Z[pos_psi_S1_vec[3],pos_psi_S1_vec[2]] = Z[pos_psi_S1_vec[3],pos_psi_S1_vec[2]] - C3_sa_F2_RTS1/X[pos_psi_S1_vec[2]]  
    # d_d plane
    #psi_d = Boltzman_factor_2_psi(X[pos_psi_S1_vec[3]], temp)
    #DY_Dpsid = -np.sqrt(8*1000*R*temp*e*eo*I)*np.cosh((zel*F*psi_d)/(2*R*temp))*((zel*F)/(2*R*temp)) - CapacitancesS1[2]
    #Dpsid_DpsidB = (-R*temp)/(F*X[pos_psi_S1_vec[3]])
    #Z[pos_psi_S1_vec[3], pos_psi_S1_vec[3]] = Z[pos_psi_S1_vec[3], pos_psi_S1_vec[3]] + (DY_Dpsid*Dpsid_DpsidB*((sS1*aS1)/F))

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
    #Z[pos_psi_S2_vec[3],pos_psi_S2_vec[2]] = Z[pos_psi_S2_vec[3],pos_psi_S2_vec[2]] - C3_sa_F2_RTS2/X[pos_psi_S2_vec[2]]  
    # d_d plane
    #psi_dS2 = Boltzman_factor_2_psi(X[pos_psi_S2_vec[3]], temp)
    #DY_Dpsid = -np.sqrt(8*1000*R*temp*e*eo*I)*np.cosh((zel*F*psi_dS2)/(2*R*temp))*((zel*F)/(2*R*temp)) - CapacitancesS2[2]
    #Dpsid_DpsidB = (-R*temp)/(F*X[pos_psi_S2_vec[3]])
    #Z[pos_psi_S2_vec[3], pos_psi_S2_vec[3]] = Z[pos_psi_S2_vec[3], pos_psi_S2_vec[3]] + (DY_Dpsid*Dpsid_DpsidB*((sS2*aS2)/F))
    
    
    #### plane gamma [diffusive plane]
    dpsiA_dXpsiA = (-R*temp)/(F*X[pos_psi_S2_vec[2]])
    Z[pos_psi_S1_vec[3],pos_psi_S1_vec[2]] =Z[pos_psi_S1_vec[3],pos_psi_S1_vec[2]] - CapacitancesS1[2]*dpsiA_dXpsiA
    Z[pos_psi_S2_vec[3],pos_psi_S2_vec[2]] =Z[pos_psi_S2_vec[3],pos_psi_S2_vec[2]] - CapacitancesS2[2]*dpsiA_dXpsiA
    #
    ew = eo*e
    Q = Zi*elec_charge                                         # Q is the charve of the aqueous elements times the electron charge 
    Cb = C_aq
    A =6.02214e23 * 1000/ew                # a prefactor = Avogadro * 1000 /ew
    kbt = 1.38064852e-23 *temp # kb (J/K) * T in K
    delta_psi = 0.0001
    x = bvp_class.x
    y0 = bvp_class.y
    #y0 = np.zeros((2, x.size))
    y_0 = y0.copy()
    y_0[0,0] = X[pos_psi_S1_vec[3]]
    'I think that  y0[1,0] and y0[1,-1] are not necessary to solve the problem, I would say that its values do not have implications. Although I am not 100% sure.'
    y_0[1,0] = (CapacitancesS1[2]*(X[pos_psi_S1_vec[3]]-Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp)))/ew           # The negative value that I am given here is extremely arbitrary I am not sure why. IT MUST BE DISCUSSED
    # y0[1,0] = dpsi_d                dpsi_d = -(sig_0 + sig_b + sig_d)/ew   # electric field at diffuse layer, x>d    
    y_0[0,-1]= X[pos_psi_S2_vec[3]]
    y_0[1,-1]= -(CapacitancesS2[2]*(X[pos_psi_S2_vec[3]]-Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp)))/ew
    
    y1 =y0.copy()
    y1[0,0] = X[pos_psi_S1_vec[3]]+delta_psi
    y1[1,0] = (CapacitancesS1[2]*((X[pos_psi_S1_vec[3]]+delta_psi)-Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp)))/ew
    y1[0,-1]= X[pos_psi_S2_vec[3]]
    y1[1,-1]= -(CapacitancesS2[2]*(X[pos_psi_S2_vec[3]]-Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp)))/ew
    
    y2 = y0.copy()
    y2[0,0] = X[pos_psi_S1_vec[3]]
    y2[1,0] = (CapacitancesS1[2]*((X[pos_psi_S1_vec[3]])-Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp)))/ew
    y2[0,-1]= X[pos_psi_S2_vec[3]]+delta_psi
    y2[1,-1]= -(CapacitancesS2[2]*((X[pos_psi_S2_vec[3]]+delta_psi)-Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp)))/ew
    
    args0=[Q,Cb,A,kbt,y_0]
    args1=[Q,Cb,A,kbt,y1]
    args2=[Q,Cb,A,kbt,y2]
    
    y1_p =y0.copy()
    y1_p[0,0] = X[pos_psi_S1_vec[3]]-delta_psi
    y1_p[1,0] = (CapacitancesS1[2]*((X[pos_psi_S1_vec[3]]-delta_psi)-Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp)))/ew
    y1_p[0,-1]= X[pos_psi_S2_vec[3]]
    y1_p[1,-1]= -(CapacitancesS2[2]*(X[pos_psi_S2_vec[3]]-Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp)))/ew
    
    y2_p = y0.copy()
    y2_p[0,0] = X[pos_psi_S1_vec[3]]
    y2_p[1,0] = (CapacitancesS1[2]*((X[pos_psi_S1_vec[3]])-Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp)))/ew
    y2_p[0,-1]= X[pos_psi_S2_vec[3]]-delta_psi
    y2_p[1,-1]= -(CapacitancesS2[2]*((X[pos_psi_S2_vec[3]]-delta_psi)-Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp)))/ew
    args3=[Q,Cb,A,kbt,y1_p]
    args4=[Q,Cb,A,kbt,y2_p]
    #PB solving  
    result0 = solve_bvp(fun_PB, bc_PB, x, y_0, tol = tol_PB, args = args0)
    result1 = solve_bvp(fun_PB, bc_PB, x, y1, tol = tol_PB, args = args1)
    result2 = solve_bvp(fun_PB, bc_PB, x, y2, tol = tol_PB, args = args2)
    result3 = solve_bvp(fun_PB, bc_PB, x, y1_p, tol = tol_PB, args = args3)
    result4 = solve_bvp(fun_PB, bc_PB, x, y2_p, tol = tol_PB, args = args4)
    
    
    d_sigma_d_psi_S1= -(ew/delta_psi)*(result1.y[1,0]-result0.y[1, 0])
    d_sigma_d_psi_S2= (ew/delta_psi)*(result2.y[1,-1]-result0.y[1, -1])
    d_sigma_d_psi_S1= -(ew/2*delta_psi)*(result1.y[1,0]-result3.y[1, 0])
    d_sigma_d_psi_S2= (ew/2*delta_psi)*(result2.y[1,-1]-result4.y[1, -1])
    
    Z[pos_psi_S1_vec[3], pos_psi_S1_vec[3]] = CapacitancesS1[2] + d_sigma_d_psi_S1
    Z[pos_psi_S2_vec[3], pos_psi_S2_vec[3]] = CapacitancesS2[2] + d_sigma_d_psi_S2
    

    # finally just return Z
    if idx_fix_species != None:
        for d in idx_fix_species:
            v=np.zeros(length_X)
            v[d]=1
            Z[d,:] = v
    return Z
 
def fun_PB(x, y, args):
    Q = args[0]
    C = args[1]
    A = args[2]
    kbt = args[3]
    arg1 = np.zeros((x.size))
    for i in range(len(Q)):
        arg1 += Q[i]*C[i]*np.exp(-Q[i]*y[0]/kbt)
    arg1 = -A*arg1
    return np.vstack((y[1] , arg1))

def bc_PB(ya, yb, args):
    y0 = args[4]
    return np.array([ya[0]-y0[0,0] , yb[0]-y0[0,-1]])


    
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