# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:51:02 2019

@author: DaniJ
"""

"""
I would like first to draw the diffusive curve and later two know for the different potentials at which distance in z is 0.
Such thing is necessary for the first benchmark of the numerical PB approach.
"""
import numpy as np
import scipy as sp

from matplotlib import pyplot as plt

# Equation 6.20
def PBdistance (y,y0,k):
    """
        Equation 4.21 of the book of Jürgen butt, Karlheinz Graf and Michael Kappl is used.
        "Physics and Chemistry of Interfaces" (2003)
        y0 is the boltzman factor at x=0 (or where the diffusive layer starts)
        y is the boltzman factor at x
        k is the debye length
    """
    a1 = np.exp(y/2)
    a2 = np.exp(y0/2)
    r1 = np.log((a1-1)/(a1+1))
    r2 = np.log((a2-1)/(a2+1))
    D = r1-r2
    x = (-1/k)*D
    return x



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

def inv_Debye_length_quick(c):
    '''
    Equation 4.10
    '''
    k = (np.sqrt(c)/3.04)  # result in 1/nm
    return k

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

def vec_inv_Debye_length(epsilon, T, Caq,Z):
    k = np.zeros(Caq.shape[0])
    for i in range(0,Caq.shape[0]):
        k[i] = inv_Debye_length(epsilon, T, Caq[i,:],Z)
    return k

def create_y(psi_D_S2, T):
    e = 1.60218e-19 # C that is the elementary charge units are Coulombs
    kb = 1.38066e-23    # Boltzmann constant units are J/K
    A = e/(kb*T)
    y = A*psi_D_S2
    return y
def plot_distance(y0,k):
    '''
    y is the psi*e/kb*t = y
    It is defined in the section 4.2.3 of book of Jürgen butt, Karlheinz Graf and Michael Kappl is used.
        "Physics and Chemistry of Interfaces" (2003)
    '''
    y = np.linspace(y0, 0, num=20)
    f=len(y)
    dist=np.zeros(f)
    for j in range(0,f):
        dist[j] = PBdistance (y[j],y0,k)
    plt.figure()
    plt.plot(dist*1e9,y)

    
# What are the values of my potentials?????? Is it High or Low potential
Xarr1 = np.load('X_arr_X.npy')
Carr1 = np.load('C_arr_X.npy')

Xarr2 = np.load('X_arr_lnX.npy')
Carr2 = np.load('C_arr_lnX.npy')

F = 96485.3328959                                   # C/mol
R =  8.314472                                       # J/(K*mol)
T = 273.15 + 25
epsilon = 78.45203739768931

pos_v1 = [5,6,7,8]
pos_v2 = [9,10,11,12]
psi_0_S1 = boltzman_2_psi(Xarr1[:,pos_v1[0]], R, T, F)
psi_C_S1 = boltzman_2_psi(Xarr1[:,pos_v1[1]], R, T, F)
psi_A_S1 = boltzman_2_psi(Xarr1[:,pos_v1[2]], R, T, F)
psi_D_S1 = boltzman_2_psi(Xarr1[:,pos_v1[3]], R, T, F)

psi_0_S2 = boltzman_2_psi(Xarr1[:,pos_v1[0]], R, T, F)
psi_C_S2 = boltzman_2_psi(Xarr1[:,pos_v1[1]], R, T, F)
psi_A_S2 = boltzman_2_psi(Xarr1[:,pos_v1[2]], R, T, F)
psi_D_S2 = boltzman_2_psi(Xarr1[:,pos_v1[3]], R, T, F)


#R = psi_D_S1-psi_D_S2
plt.plot(psi_D_S2)

Caq = Xarr1[:,0:3]
Z = [1,-1,1]

kappa_vec = vec_inv_Debye_length(epsilon, T, Caq,Z)
y = create_y(psi_D_S2,T)

for i in range(0, len(y)):
    plot_distance(y[i], kappa_vec[i])