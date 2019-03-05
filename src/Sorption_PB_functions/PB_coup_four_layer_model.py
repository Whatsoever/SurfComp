# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:08:29 2019

@author: 
"""

from four_layer_model import four_layer_model
import numpy as np
import scipy as sp


def PB_and_fourlayermodel ():
    """
        The following function tries to implement a specific coupling between the local chemistry given at two colloidal surface that interact between them through the Possion-Boltzman relationship.
        For the local chemistry a four layer model is used, based on the notation and methodology presented by Westall (1980). The fact that we use Westall (1980) formulation is important in order to know
        some assumptions that tend to be slightly different between codes.
    """
    # The initial guess should be done using the four layer model alone, but right now it is quite open to different possibilities
    # Hence, we assume that the initial guess is given.
    
    
    
    # scipy.optimize.fsolve(func, x0, args=(), fprime=None, full_output=0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)[source]¶
    X = sp.optimize.fsolve(func_NR_PB_FLM, X_guess, args = (), fprime = Jacobian_NR_PB_FLM)
    #Speciation
    C = log_k + A*np.log10(X)
    return X, C

def func_NR_PB_FLM():
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
    C_aq = C[idx_Aq]
    #I = Calculate_ionic_strength(Z, C_aq)
    T = Update_T_PB_FLM(T, s,e, I, temp, a, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq)
    # Calculation of Y
    Y= np.matmul(A.transpose(),C)-T
    return Y

def Update_T_PB_FLM(T, s,e, I, temp, a, CapS1, CapS2, psi_vS1, psi_vS2, pos_psi0S1, pos_psialphaS1, pos_psibetaS1, pos_psigammaS1,pos_psi0S2, pos_psialphaS2, pos_psibetaS2, pos_psigammaS2, x, Z, C_aq):
    """
        This equation is linked to func_NR_PB_FLM. It updates the values of T for the electrostatic parameters.
        CapS1       is the vector of capacitances for surface 1. The units are supossed to be F/m². 
        CapS2       is the vector of capacitances for surface 2. The units are supossed to be F/m². 
        psi_vS1   is the vector of electrostatic potentials for the surface1. The units are supposed to be volts
        psi_vS2   is the vector of electrostatic potentials for the surface2. The units are supposed to be volts
        s       is the specific surface area
        a       area
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
    sigma_alpha_S1= -sigma_0 + CapS1[1]*(psi_vS1[1]-psi_vS1[2])
    sigma_beta_S1 = -sigma_0-sigma_alpha+CapS1[2]*(psi_vS1[2]-psi_vS1[3])
    sigma_gamma_S1 = -sigma_0 - sigma_alpha - sigma_beta
    #Surface 2
    sigma_0_S2 = CapS2[0]*(psi_vS2[0]-psi_vS2[1])
    sigma_alpha_S2 = -sigma_0 + CapS2[1]*(psi_vS2[1]-psi_vS2[2])
    sigma_beta_S2 = -sigma_0-sigma_alpha+CapS2[2]*(psi_vS2[2]-psi_vS2[3])
    sigma_gamma_S2 = -sigma_0 - sigma_alpha - sigma_beta
    
    # T
    T_0S1 = ((s*a)/F)*sigma_0_S1;                    # units mol/L or mol/kg
    T_alphaS1 = ((s*a)/F)*sigma_alpha_S1;            # units mol/L or mol/kg
    T_betaS1 = ((s*a)/F)*sigma_beta_S1;              # units mol/L or mol/kg
    
    T_0S2 = ((s*a)/F)*sigma_0_S2;                    # units mol/L or mol/kg
    T_alphaS2 = ((s*a)/F)*sigma_alpha_S2;            # units mol/L or mol/kg
    T_betaS2 = ((s*a)/F)*sigma_beta_S2;              # units mol/L or mol/kg
    #!! Important!! Here starts the PB part
    ew = eo*e
    #
    Q = Z*elec_charge                                         # Q is the charve of the aqueous elements times the electron charge 
    C = C_aq
    A =6.02214e23 * 1000/ew                # a prefactor = Avogadro * 1000 /ew
    kbt = 1.38064852e-23 *temp # kb (J/K) * T in K
    y0[0,0] = psi_vS1[3]
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

def fun_PB(x, y, args):
    Q = args[0]
    C = args[1]
    A = args[2]
    kbt = args[3]
    arg1 = num.zeros((x.size))
    for i in range(len(Q)):
        arg1 += Q[i]*C[i]*num.exp(-Q[i]*y[0]/kbt)
    arg1 = -A*arg1
    return num.vstack((y[1] , arg1))

def bc_PB(ya, yb, args):
    y0 = args[4]
    return num.array([ya[0]-y0[0,0] , yb[0]-y0[0,-1]])



def Jacobian_NR_PB_FLM():
    
    
    
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