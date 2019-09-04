# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:24:34 2019

@author: DaniJ
"""
"""

"""


import PB_coup_four_layer_2try as flmPB
import four_layer_model_2try_withFixSpeciesOption_Scaling_2surface as flm1
import numpy as np
import scipy as sp
from bvp import solve_bvp
import matplotlib.pyplot as plt

from matplotlib import pyplot as plt

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
def inv_Debye_length (epsilon, T, Caq,Z):
    '''
    epsilon is te relative permitivity
    T is the temperture in Kelvins
    Caq is the vector of aqueous concentrations in mol
    Z is the vector of ionic charge
    Note: Z and Caq should be coherent in order
    '''
    e = 1.60218e-19 # C that is the elementary charge units are Coulombs
    kb = 1.38066e-23    # Boltzmann constant units are J/K
    epsilon0 =  8.854187871e-12  # vacuum permitivity units (A s)/(V m)
    Navo = 6.022141e23 
    
    C = Navo*Caq
    CZ = np.multiply(C,np.multiply(Z,Z))
    sumCZ = np.sum(CZ)
    A = e*e*sumCZ
    B = epsilon0*epsilon*T*kb
    k = np.sqrt(A/B)
    return k
def create_y(psi_D_S2, T):
    e = 1.60218e-19 # C that is the elementary charge units are Coulombs
    kb = 1.38066e-23    # Boltzmann constant units are J/K
    A = e/(kb*T)
    y = A*psi_D_S2
    return y
def get_inival (kappa, boltzmann_psid, x):
    A = np.exp(boltzmann_psid/2)
    B = np.exp(-kappa*x)
    Numerator = A+1+(A-1)*B
    Denominator =  A+1-(A-1)*B
    y = 2*np.log(Numerator/Denominator)
    return y
def return_y(y, T):
    e = 1.60218e-19 # C that is the elementary charge units are Coulombs
    kb = 1.38066e-23    # Boltzmann constant units are J/K
    A = e/(kb*T)
    psi_D_S2 = y/A
    return psi_D_S2

idx_fix_species=[0]

#idx_fix_species=None
#T_H = np.linspace(-3,-11.2,42)
#T_H = 10**T_H

X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9,8.7e-7, 0.9, 0.9, 0.9])

# A in the columns we have X, and in the rows with have the concentration (vector C)
# Remark A_transpose is the U or component matrix.

A = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,0],[1,0,0,1,0,1,0,0,0,0,0,0,0], \
              [-1,0,0,1,0,-1,0,0,0,0,0,0,0], [1,1,0,1,0,1,0,-1,0,0,0,0,0], [-1,0,1,1,0,-1,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0],[1,0,0,0,1,0,0,0,0,1,0,0,0],\
              [-1,0,0,0,1,0,0,0,0,-1,0,0,0],[1,1,0,0,1,0,0,0,0,1,0,-1,0],[-1,0,1,0,1,0,0,0,0,-1,1,0,0]])

# Z, is  a vector with the ionic charge of the aqueous species. It should be somehow in order with the aqueous species. Looking at C, we see only 3 aqueous species H+, Cl-, and Na+
Z=np.array([1,-1,1])
# Log_k the log of reactions. Primary species have log_k = 0, I assume that it follows the same order than C

log_k = np.array([0, 0, 0, 0, 4, -8.8, 5.82, -7,0,4,-8.8,5.82,-7])


idx_Aq=np.array([0,1,2])

pos_psi_S1_vec = np.array([5,6,7,8])
pos_psi_S2_vec = np.array([9,10,11,12])


# Temperature
temp=273.15+25
#s             is the specific surface area
sS1=1         # m2/l
aS1=1         # g/l
sS2=1         # m2/l
aS2=1         # g/l
e = 78.45203739768931
CapacitancesS1=[1.05, 3.36, 0.27] 
CapacitancesS2=[1.05, 3.36, 0.27] 
#
x = np.linspace(0,800e-9,100)
T=np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9,8.7e-7, 0.9, 0.9, 0.9]) 
#(T, X_guess, distance, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, x, idx_fix_species = None, zel=1, tolerance = 1e-6, max_iterations = 100,scalingRC = True):
[X,C1] = flm1.four_layer_two_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq, pos_psi_S1_vec, pos_psi_S2_vec, temp, sS1, aS1, sS2, aS2, e, CapacitancesS1, CapacitancesS2, idx_fix_species, zel=1, tolerance = 1e-6, max_iterations = 100,scalingRC = True)
boltzmann_psid = X[8]
X[8] = Boltzman_factor_2_psi (X[8],temp)
X[12] = Boltzman_factor_2_psi (X[12],temp)
# How should look the boltzman profile?????

# Speciation - mass action law
log_C = log_k + np.matmul(A,np.log10(X))
# transf
C = 10**(log_C)

psi_S1_v = [Boltzman_factor_2_psi(X[pos_psi_S1_vec[0]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[1]], temp), Boltzman_factor_2_psi(X[pos_psi_S1_vec[2]], temp), X[pos_psi_S1_vec[3]]]
psi_S2_v = [Boltzman_factor_2_psi(X[pos_psi_S2_vec[0]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[1]], temp), Boltzman_factor_2_psi(X[pos_psi_S2_vec[2]], temp), X[pos_psi_S2_vec[3]]] 


sigma_S1_0 = CapacitancesS1[0]*(psi_S1_v[0]-psi_S1_v[1])
sigma_S1_alpha = -sigma_S1_0 + CapacitancesS1[1]*(psi_S1_v[1]-psi_S1_v[2])
sigma_S1_beta = -sigma_S1_0-sigma_S1_alpha+CapacitancesS1[2]*(psi_S1_v[2]-psi_S1_v[3])
sigma_S1_dTLM = -sigma_S1_0 - sigma_S1_alpha - sigma_S1_beta
    ########## S2  #####################
sigma_S2_0 = CapacitancesS2[0]*(psi_S2_v[0]-psi_S2_v[1])
sigma_S2_alpha = -sigma_S2_0 + CapacitancesS2[1]*(psi_S2_v[1]-psi_S2_v[2])
sigma_S2_beta = -sigma_S2_0-sigma_S2_alpha+CapacitancesS2[2]*(psi_S2_v[2]-psi_S2_v[3])
sigma_S2_dTLM = -sigma_S2_0 - sigma_S2_alpha - sigma_S2_beta










C_aq = C[idx_Aq]
eo = 8.854187871e-12                                # Farrads = F/m   - permittivity in vaccuum                         # C
kb = 1.38064852e-23                                 # J/K other units --> kb=8,6173303e-5  eV/K
Na = 6.022140857e23                                 # 1/mol
elec_charge = 1.60217662e-19 #electron charge in C

kappa = inv_Debye_length (e, temp, C_aq,Z)
y = get_inival (kappa, create_y(X[8], temp), x)

psi_y = return_y(y,temp)

#plt.figure()
#plt.plot(x, y)
#plt.figure()
#plt.plot(x, psi_y)

ew = eo*e
Q = Z*elec_charge                                         # Q is the charve of the aqueous elements times the electron charge 
A =Na* 1000/ew                # a prefactor = Avogadro * 1000 /ew
kbt = 1.38064852e-23 *temp # kb (J/K) * T in K






y0 = np.zeros((2, x.size))
y0[0,0] = X[8]
y0[0,-1] = X[12]
y0[1,0] = -sigma_S1_dTLM/ew 
y0[1,-1]= sigma_S2_dTLM/ew
args=[Q,C_aq,A,kbt,y0]
y02 = np.zeros((2, x.size))
y02[0,0] = X[8]
y02[0,-1] = X[12]
args2=[Q,C_aq,A,kbt,y02]
tol_PB = 1e-3
result = solve_bvp(fun_PB, bc_PB, x, y0,  args = args, tol=tol_PB,verbose = 2)
result2 = solve_bvp(fun_PB, bc_PB, x, y02,  args = args2, tol=tol_PB,verbose = 2)
plt.figure()
plt.plot(result.x, result.y[0])



y03 = np.zeros((2, x.size))
y03[0,0:50] = psi_y[0:50]
ff=psi_y[::-1]
y03[0,50:] = ff[50:]
args3=[Q,C_aq,A,kbt,y03]
result3 = solve_bvp(fun_PB, bc_PB, x, y03,  args = args3, tol=tol_PB,verbose = 2)
plt.figure()
plt.plot(result3.x, result3.y[0])


x = np.linspace(0,800e-9,62)
y04 = np.zeros((2, 62))
y04[0,:] = [0.110684801354927,0.0351366357009602,0.0120739252548629,0.00928365673741179,0.00895773089151665,0.00891967765053806,0.00891523479923835,0.00891471608058322,0.00891465551833957,0.00891464844748273,0.00891464762193510,0.00891464752554948,0.00891464751429612,0.00891464751298225,0.00891464751283962,0.00891464751281413,0.00891464751280957,0.00891464751280875,0.00891464751280862,0.00891464751280859,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280858,0.00891464751280864,0.00891464751280910,0.00891464751281089,0.00891464751281879,0.00891464751287998,0.00891464751330760,0.00891464751501557,0.00891464752256918,0.00891464755597560,0.00891464770371810,0.00891464803050745,0.00891464891667836,0.00891465732403850,0.00891471608058598,0.00891523479926005,0.00891967765067112,0.00895773089293939,0.00928365674581903,0.0120739252983962,0.0351366357457432,0.110684801354927];
y04[1,:]=[-28895446.8131893,-3435965.76963390,-390917.699171017,-45622.4022592364,-5326.54765625128,-621.893292647840,-72.6082486050005,-8.47727067211872,-0.989751432176056,-0.115556991884568,-0.0134916887535240,-0.00157520264094332,-0.000183910314805805,-2.14722885876059e-05,-3.83811078947183e-06,-6.86176765102564e-07,-1.22653670862259e-07,-2.20152607109814e-08,-4.81946345927546e-09,-1.70188443770758e-09,-7.62869182982919e-10,-4.62778452563986e-10,-3.15792381813486e-10,-3.65915733941120e-11,-1.28523079951723e-11,-4.52918864898380e-12,-1.63425355904211e-12,-3.68053906133982e-13,-4.25912401535252e-13,-2.01386401554055e-12,-1.08405344940722e-11,-3.02456950522095e-11,-8.58820358267262e-11,-8.50365591420824e-11,-5.04290170466442e-11,-1.24030311039289e-11,1.03999024837637e-11,4.46517768005715e-11,4.25250112692256e-11,3.04775488822333e-10,8.58518889865602e-10,7.59477842255297e-09,6.45728810260557e-08,2.85523544882341e-07,1.26321312181308e-06,8.82794196549841e-06,6.16966346041813e-05,0.000272858282927239,0.00120673676200361,0.00533688794167901,0.0236028043895202,0.0640048985344280,0.173565267680231,1.21299620817432,8.47727098353075,72.6082507908549,621.893311275133,5326.54780926240,45622.4030228492,390917.704252827,3435965.76365446,28895446.8110238];
args4=[Q,C_aq,A,kbt,y04]
result4 = solve_bvp(fun_PB, bc_PB, x, y04,  args = args4, tol=tol_PB, max_nodes=100000,verbose = 2)
plt.figure()
plt.plot(result4.x, result4.y[0])