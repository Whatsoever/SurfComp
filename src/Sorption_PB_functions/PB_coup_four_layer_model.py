# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:08:29 2019

@author: 
"""

from four_layer_model import four_layer_model
import numpy as np
import scipy as sp
from bvp import solve_bvp

def PB_and_fourlayermodel (X_guess, log_k, A, T, sS1, sS2,e, temp, aS1, aS2, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq):
    """
        The following function tries to implement a specific coupling between the local chemistry given at two colloidal surface that interact between them through the Possion-Boltzman relationship.
        For the local chemistry a four layer model is used, based on the notation and methodology presented by Westall (1980). The fact that we use Westall (1980) formulation is important in order to know
        some assumptions that tend to be slightly different between codes.
    """
    # The initial guess should be done using the four layer model alone, but right now it is quite open to different possibilities
    # Hence, we assume that the initial guess is given.
    
    
    
    # scipy.optimize.fsolve(func, x0, args=(), fprime=None, full_output=0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)[source]¶
    X = sp.optimize.fsolve(func_NR_PB_FLM, X_guess, args = ( A, log_k, T, sS1, sS2,e, temp, aS1, aS2, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq), fprime = Jacobian_NR_PB_FLM)
    #Speciation
    log_C = log_k + A*np.log10(X)
    C = 10**(log_C)
    return X, C

def func_NR_PB_FLM(X, A, log_k, T, sS1, sS2,e, temp, aS1, aS2, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq):
    """
        This function is supossed to be linked to the PB_and_fourlayermodel function.
        It just gave the evaluated vector of Y, for the Newton-Raphson procedure.
        The formulation of Westall (1980) is followed. 
        FLM = four layer model. There are two colloids which means 4 electrostatic unknowns.
    """
     # Speciation - mass action law
    log_C = log_k + A*np.log10(X)
    # transf
    C = 10**(log_C)
    # Update T - "Electrostatic parameters"
    psi_vS1 = [Boltzman_factor_2_psi(X[pos_psi0S1], temp), Boltzman_factor_2_psi(X[pos_psialphaS1], temp), Boltzman_factor_2_psi(X[pos_psibetaS1], temp), X[pos_psigammaS1]]  # [e^((-F〖ψ_0〗_s1)/Rtemp) ],[e^((-F〖ψ_α〗_s1)/Rtemp) ],[e^((-F〖ψ_β〗_s1)/Rtemp) ],〖ψ_d〗_s1
    psi_vS2 = [Boltzman_factor_2_psi(X[pos_psi0S2], temp), Boltzman_factor_2_psi(X[pos_psialphaS2], temp), Boltzman_factor_2_psi(X[pos_psibetaS2], temp), X[pos_psigammaS2]]
    #I = Calculate_ionic_strength(Z, C_aq)
    T = Update_T_PB_FLM(T, sS1, sS2,e, temp, aS1, aS2, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq)
    # Calculation of Y
    Y= np.matmul(A.transpose(),C)-T
    return Y

def Update_T_PB_FLM(T, sS1, sS2,e, temp, aS1, aS2, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq):
    """
        This equation is linked to func_NR_PB_FLM. It updates the values of T for the electrostatic parameters.
        CapS1       is the vector of capacitances for surface 1. The units are supossed to be F/m². 
        CapS2       is the vector of capacitances for surface 2. The units are supossed to be F/m². 
        psi_vS1   is the vector of electrostatic potentials for the surface1. The units are supposed to be volts
        psi_vS2   is the vector of electrostatic potentials for the surface2. The units are supposed to be volts
        s       is the specific surface area
        a       concentration of suspended solid
        temp    is the temperature in Kelvins
        e       is the relative permittivity
        I       is the ionic strenght mol/m³
    """
    #  constant
    F = 96485.3328959                                   # C/mol
    R =  8.314472                                       # J/(K*mol)
    elec_charge = 1.60217662e-19 #electron charge in C
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    # Update of T
    # T = (sa/F)*sigma
    # sigma = C*(psi-psi_0) <-- The sigma depends on the plane
    
    
    # NOTE: I JUST REALISED I DID SOMETHING WRONG THE WHOLE TIME. In the T vector, It must be checked in the other approaches (TLM)
    # The equations os the layers that are not in contact with the diffusion layer must be mol/kg or mol/l
    # while the third component is the total charge balance.
    
    # Surface 1
    sigma_0_S1 = CapS1[0]*(psi_vS1[0]-psi_vS1[1])
    sigma_alpha_S1= -sigma_0_S1 + CapS1[1]*(psi_vS1[1]-psi_vS1[2])
    sigma_beta_S1 = -sigma_0_S1-sigma_alpha_S1+CapS1[2]*(psi_vS1[2]-psi_vS1[3])
    sigma_gamma_S1 = -sigma_0_S1 - sigma_alpha_S1 - sigma_beta_S1
    #Surface 2
    sigma_0_S2 = CapS2[0]*(psi_vS2[0]-psi_vS2[1])
    sigma_alpha_S2 = -sigma_0_S2 + CapS2[1]*(psi_vS2[1]-psi_vS2[2])
    sigma_beta_S2 = -sigma_0_S2-sigma_alpha_S2+CapS2[2]*(psi_vS2[2]-psi_vS2[3])
    sigma_gamma_S2 = -sigma_0_S2 - sigma_alpha_S2 - sigma_beta_S2
    
    # T
    T_0S1 = ((sS1*aS1)/F)*sigma_0_S1;                    # units mol/L or mol/kg
    T_alphaS1 = ((sS1*aS1)/F)*sigma_alpha_S1;            # units mol/L or mol/kg
    T_betaS1 = ((sS1*aS1)/F)*sigma_beta_S1;              # units mol/L or mol/kg
    
    T_0S2 = ((sS2*aS2)/F)*sigma_0_S2;                    # units mol/L or mol/kg
    T_alphaS2 = ((sS2*aS2)/F)*sigma_alpha_S2;            # units mol/L or mol/kg
    T_betaS2 = ((sS2*aS2)/F)*sigma_beta_S2;              # units mol/L or mol/kg
    #!! Important!! Here starts the PB part
    ew = eo*e
    #
    Q = Z*elec_charge                                         # Q is the charve of the aqueous elements times the electron charge 
    C = C_aq
    A =6.02214e23 * 1000/ew                # a prefactor = Avogadro * 1000 /ew
    kbt = 1.38064852e-23 *temp # kb (J/K) * T in K
    y0[0,0] = psi_vS1[3]
    'I think that  y0[1,0] and y0[1,-1] are not necessary to solve the problem, I would say that its values do not have implications. Although I am not 100% sure.'
    y0[1,0] = -sigma_gamma_S1/ew           # The negative value that I am given here is extremely arbitrary I am not sure why. IT MUST BE DISCUSSED
    # y0[1,0] = dpsi_d                dpsi_d = -(sig_0 + sig_b + sig_d)/ew   # electric field at diffuse layer, x>d
    y0[0,-1]= psi_vS2[3]
    y0[1,-1]= sigma_gamma_S2/ew
    args=[Q,C,A,kbt,y0]
    result = solve_bvp(fun_PB, bc_PB, x, y0, tol = 1e-8, args = args)
    #
    sigma_dS1=result.y[1][1]*ew
    sigma_dS2=result.y[1][-1]*ew
    #
    
    
    
    T_gammadS1 = -sigma_gamma_S1+sigma_dS1             # This part should be equal to C[2]*(psi_beta-psi_dorgamma)+sigma_d
    T_gammadS2 = -sigma_gamma_S2+sigma_dS2             # This part should be equal to C[2]*(psi_beta-psi_dorgamma)+sigma_d
    # Now the values must be put in T
    T[pos_psi0S1] = T_0S1
    T[pos_psialphaS1] = T_alphaS1
    T[pos_psibetaS1] = T_betaS1
    T[pos_psi0S2] = T_0S2
    T[pos_psialphaS2] = T_alphaS2
    T[pos_psibetaS2] = T_betaS2
    
    
    
    T[pos_psigammaS1] = T_gammadS1
    T[pos_psigammaS2] = T_gammadS2
    return T

def Jacobian_NR_PB_FLM(X, A, log_k, T, sS1, sS2,e, temp, aS1, aS2, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq):
    #  constant
    F = 96485.3328959                                   # C/mol [Faraday constant]
    R =  8.314472                                       # J/(K*mol) [universal constant gas]
    elec_charge = 1.60217662e-19 #electron charge in C
    eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum
    # Speciation - mass action law
    log_C = log_k + A*np.log10(X)
    # transf
    C = 10**(log_C)
    # First part is the common of the Jacbian derivation
    length_X = X.size
    for i in range(0, length_X):
            for j in range(0, length_X):
                Z[i,j]= np.matmul(np.multiply(A[:,i], A[:,j]), (C/X[j]))
    

    # Now the electrostatic part must be modified, one question hang on the air:
    # Should we check that the electrostatic part is as we expected?
    sa_F2_S1 = (sS1*aS1)/(F*F)
    sa_F2_S2 = (sS2*aS2)/(F*F)
    #### plane 0
    C1_sa_F2_RT_S1 = sa_F2_S1*CapS1[0]*R*temp
    C1_sa_F2_RT_S2 = sa_F2_S2*CapS2[0]*R*temp
    # Assigning in Jacobian (plane 0)
    #S1
    Z[pos_psi0S1,pos_psi0S1]=Z[pos_psi0S1,pos_psi0S1] + C1_sa_F2_RT_S1/X[pos_psi0S1]
    Z[pos_psi0S1,pos_psialphaS1]=Z[pos_psi0S1, pos_psialphaS1] - C1_sa_F2_RT_S1/X[pos_psialphaS1]    
    #S2        
    Z[pos_psi0S2,pos_psi0S2]=Z[pos_psi0S2,pos_psi0S2] + C1_sa_F2_RT_S2/X[pos_psi0S2]
    Z[pos_psi0S2,pos_psialphaS2]=Z[pos_psi0S2, pos_psialphaS2] - C1_sa_F2_RT_S2/X[pos_psialphaS2]    
                
    #### plane alpha
    C1C2_sa_F2_RT_S1 = sa_F2_S1*R*temp*(CapS1[0]+CapS1[1])
    C2_sa_F2_RT_S1 = sa_F2_S1*CapS1[1]*R*temp
    C1C2_sa_F2_RT_S2 = sa_F2_S2*R*temp*(CapS2[0]+CapS2[1])
    C2_sa_F2_RT_S2 = sa_F2_S2*CapS2[1]*R*temp
    # Assigning in Jacobian (plane alpha)
    #S1
    Z[pos_psialphaS1,pos_psi0S1]=Z[pos_psialphaS1, pos_psi0S1] - C1_sa_F2_RT_S1/X[pos_psi0S1]
    Z[pos_psialphaS1,pos_psialphaS1]=Z[pos_psialphaS1, pos_psialphaS1] + C1C2_sa_F2_RT_S1/X[pos_psialphaS1]
    Z[pos_psialphaS1,pos_psibetaS1]= Z[pos_psialphaS1,pos_psibetaS1] - C2_sa_F2_RT_S1/X[pos_psibetaS1] 
    #S2           
    Z[pos_psialphaS2,pos_psi0S2]=Z[pos_psialphaS2, pos_psi0S2] - C1_sa_F2_RT_S2/X[pos_psi0S2]
    Z[pos_psialphaS2,pos_psialphaS2]=Z[pos_psialphaS2, pos_psialphaS2] + C1C2_sa_F2_RT_S2/X[pos_psialphaS2]
    Z[pos_psialphaS2,pos_psibetaS2]= Z[pos_psialphaS2,pos_psibetaS2] - C2_sa_F2_RT_S2/X[pos_psibetaS2]    

    #### plane beta
    C3C2_sa_F2_RT_S1 = sa_F2_S1*R*temp*(CapS1[1]+CapS1[2])
    C3_sa_F2_RT_S1 = sa_F2_S1*CapS1[2]*R*temp
    C3C2_sa_F2_RT_S2 = sa_F2_S2*R*temp*(CapS2[1]+CapS2[2])
    C3_sa_F2_RT_S2 = sa_F2_S2*CapS2[2]*R*temp
    # Assigning in Jacobian (plane beta)
    #S1
    Z[pos_psibetaS1, pos_psialphaS1] = Z[pos_psibetaS1, pos_psialphaS1] - C2_sa_F2_RT_S1/X[pos_psialphaS1]
    Z[pos_psibetaS1, pos_psibetaS1] = Z[pos_psibetaS1, pos_psibetaS1] + C3C2_sa_F2_RT_S1/X[pos_psibetaS1]
    Z[pos_psibetaS1, pos_psigammaS1] = Z[pos_psibetaS1, pos_psigammaS1] - C3_sa_F2_RT_S1/X[pos_psigammaS1]
    #S2
    Z[pos_psibetaS2, pos_psialphaS2] = Z[pos_psibetaS2, pos_psialphaS2] - C2_sa_F2_RT_S2/X[pos_psialphaS2]
    Z[pos_psibetaS2, pos_psibetaS2] = Z[pos_psibetaS2, pos_psibetaS2] + C3C2_sa_F2_RT_S2/X[pos_psibetaS2]
    Z[pos_psibetaS2, pos_psigammaS2] = Z[pos_psibetaS2, pos_psigammaS2] - C3_sa_F2_RT_S2/X[pos_psigammaS2]               
                
    #### plane gamma [diffusive plane]      
    
    Z[pos_psigammaS1,pos_psibetaS1] = Z[pos_psigammaS1, pos_psibetaS1] - CapS1[2]
    Z[pos_psigammaS2,pos_psibetaS2] = Z[pos_psigammaS2, pos_psibetaS2] - CapS2[2]
    #
    ew = eo*e
    #
    Q = Z*elec_charge                                         # Q is the charve of the aqueous elements times the electron charge 
    Cb = C_aq
    A =6.02214e23 * 1000/ew                # a prefactor = Avogadro * 1000 /ew
    kbt = 1.38064852e-23 *temp # kb (J/K) * T in K
    delta_psi = 0.001
    
    
    y0 = np.zeros((2, x.size))
    y0[0,0] = X[pos_psigammaS1]
    'I think that  y0[1,0] and y0[1,-1] are not necessary to solve the problem, I would say that its values do not have implications. Although I am not 100% sure.'
    y0[1,0] = -(CapS1[2]*(X[pos_psigammaS1]-Boltzman_factor_2_psi(X[pos_psibetaS1], temp)))/ew           # The negative value that I am given here is extremely arbitrary I am not sure why. IT MUST BE DISCUSSED
    # y0[1,0] = dpsi_d                dpsi_d = -(sig_0 + sig_b + sig_d)/ew   # electric field at diffuse layer, x>d    
    y0[0,-1]= X[pos_psigammaS2]
    y0[1,-1]= (CapS2[2]*(X[pos_psigammaS2]-Boltzman_factor_2_psi(X[pos_psibetaS2], temp)))/ew
    
    y1 = np.zeros((2, x.size))
    y1[0,0] = X[pos_psigammaS1]*delta_psi
    y1[1,0] = -(CapS1[2]*((X[pos_psigammaS1]*delta_psi)-Boltzman_factor_2_psi(X[pos_psibetaS1], temp)))/ew
    y1[0,-1]= X[pos_psigammaS2]
    y1[1,-1]= (CapS2[2]*(X[pos_psigammaS2]-Boltzman_factor_2_psi(X[pos_psibetaS2], temp)))/ew
    
    y2 = np.zeros((2, x.size))
    y2[0,0] = X[pos_psigammaS1]
    y2[1,0] = -(CapS1[2]*((X[pos_psigammaS1]*delta_psi)-Boltzman_factor_2_psi(X[pos_psibetaS1], temp)))/ew
    y2[0,-1]= X[pos_psigammaS2]*delta_psi
    y2[1,-1]= (CapS2[2]*((X[pos_psigammaS2]*delta_psi)-Boltzman_factor_2_psi(X[pos_psibetaS2], temp)))/ew
    
    args0=[Q,Cb,A,kbt,y0]
    args1=[Q,Cb,A,kbt,y1]
    args2=[Q,Cb,A,kbt,y2]
    #PB solving  
    result0 = solve_bvp(fun_PB, bc_PB, x, y0, tol = 1e-8, args = args0)
    result1 = solve_bvp(fun_PB, bc_PB, x, y1, tol = 1e-8, args = args1)
    result2 = solve_bvp(fun_PB, bc_PB, x, y2, tol = 1e-8, args = args2)
    
    d_sigma_d_psi_S1= (ew/(delta_psi*X[pos_psigammaS1]))*(result1.y[1,0]-result0.y[1, 0])
    d_sigma_d_psi_S2= (ew/(delta_psi*X[pos_psigammaS2]))*(result2.y[1,-1]-result0.y[1, -1])
    
    Z[pos_psigammaS1, pos_psigammaS1] = CapS1[2] + d_sigma_d_psi_S1
    Z[pos_psigammaS2, pos_psigammaS2] = CapS2[2] + d_sigma_d_psi_S2
    
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