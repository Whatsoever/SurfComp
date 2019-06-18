# -*- coding: utf-8 -*-
"""
Created on Sun May 26 08:50:16 2019

@author: DaniJ
"""
import four_layer_model_2try as flm
import numpy as np
import scipy as sp

from matplotlib import pyplot as plt


def funky (T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, tolerance_B):
    try:
        [X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = tolerance_B,max_iterations = 200)
        T_error = tolerance_B
        F.write(str(tolerance_B))
        F.write("\n")
    except:
        tolerance_B = tolerance_B*10
        [X,C, T_error]=funky (T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, tolerance_B)
        F.write(str(tolerance_B))
        F.write("\n")
    
    return X,C, T_error





T_H = np.linspace(-3,-11.2,42)
T_H = 10**T_H


#X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 2, 1.5, 2e-14])
X_guess = np.array([1e-3, 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9])

# A in the columns we have X, and in the rows with have the concentration of vector C --> C=[ H^+ Cl- Na+ SOH SOH2+ SO- SOH_2Cl SONa ]
# Remark A_transpose is the U or component matrix.

A = np.array([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[1,0,0,1,1,0,0,0],[-1,0,0,1,-1,0,0,0], [1,1,0,1,1,0,-1,0], [-1,0,1,1,-1,1,0,0]])

# Z, is  a vector with the ionic charge of the aqueous species. It should be somehow in order with the aqueous species. Looking at C, we see only 3 aqueous species H+, Cl-, and Na+
Z=np.array([1,-1,1])

# Log_k the log of reactions. Primary species have log_k = 0, I assume that it follows the same order than C

log_k = np.array([0, 0, 0, 0, 4, -8.8, 5.82, -7])

# - idx_Aq        An index vector with the different aqueous species position. It must coincide with the rows of "A".
'I must check if the first position of idx_Aq is 0 or 1. Since it is python I guess it will be 0'
idx_Aq=np.array([0,1,2])

# pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma basically the same thing that idx_Aq, but only scalar. (Somehow specified). Problably position should agree with X, or T.
'Maybe I can change the name of the psialpha for psiC, psibeta for psiA, and psigamma for psid. Right now, not really the aim.'
pos_psi0=4
pos_psialpha=5
pos_psibeta=6
pos_psigamma=7

# Temperature
temp=273.15+25
#s             is the specific surface area
s=1         # m2/l
a=1         # g/l
e = 78.45203739768931
Capacitances=[1.05, 3.36, 0.27] 

n=0
ng=0
nn=0
#for i in range(0,len(T_H)):
#    # T_transpose =[ T_H+ T_Cl- T_Na+ T_SOH T_σ_0 T_σ_c T_σ_A T_σ_d ] 
#    T=np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 1, 1, 1, 1]) 
#    X_guess = np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9])
#    print(i)
#    n+=1
#    try:
#        ng+=1
#        [X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-6,max_iterations = 200)
#    except:
#        nn+=1
#        [X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-5,max_iterations = 200)

tolerance_vector=[]
Array_X = []
Array_C = []
global F 
F= open("T_efile.txt","w") 
for i in range(0,len(T_H)):
    # T_transpose =[ T_H+ T_Cl- T_Na+ T_SOH T_σ_0 T_σ_c T_σ_A T_σ_d ] 
    T=np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 1, 1, 1, 1]) 
    X_guess = np.array([T_H[i], 1e-3, 1e-3, 9.9635e-6, 8.7e-7, 0.9, 0.9, 0.9])
    print(i)
    n+=1
    tolerance_B=1e-8
    [X,C, T_e]= funky (T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, tolerance_B)
    tolerance_vector.append(T_e)
    if i == 0:
        Array_X = X
        Array_C = C
    else:
        Array_X = np.vstack([Array_X, X])
        Array_C = np.vstack([Array_C, C])
F.close()
np.save('tol_vec_v2',tolerance_vector)
np.save('X_arr_v2',Array_X)
np.save('C_arr_v2',Array_C)


#    try:
#        ng+=1
#        [X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-6,max_iterations = 200)
#    except:
#        nn+=1
#        [X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-5,max_iterations = 200)


#def funky (T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance_B):
#    try:
#        [X,C]=flm.four_layer_one_surface_speciation ( T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance = 1e-6,max_iterations = 200)
#    except:
#        tolerance_B = tolerance_B*10
#        [X,C, T_error]=funky (T, X_guess, A, Z, log_k, idx_Aq,pos_psi0, pos_psialpha, pos_psibeta,  pos_psigamma,temp, s, a, e, Capacitances, zel=1, tolerance_B)
    
#    return X,C, T_error
def plotElement(X, Y, Sep, xlab, ylab, Title):
    "It is somehow specific for this example"
    plt.figure()
    for i in range(0, len(X)):
        s,t = get_colorline_thiscase(Sep[i])
        plt.plot(X[i],Y[i],s,label=t)
        plt.legend(loc='best')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(Title)

def get_colorline_thiscase(value):
    if value-1e8<np.finfo(float).eps:
        col="g-"
        tol=str(1e8)
    elif value-1e7<np.finfo(float).eps:
        col="b-"
        tol=str(1e7)
    elif value-1e6<np.finfo(float).eps:
        col="k-"
        tol=str(1e6)
    elif value-1e5<np.finfo(float).eps:
        col="c-"
        tol=str(1e5)
    elif value-1e4<np.finfo(float).eps:
        col="r-"
        tol=str(1e4)
    return col, tol


#plotElement(np.log10(T_H), np.log10(Array_X[:,0]), tolerance_vector, "H component (mol/L)", "H+ (mol/L)", "H+ vs. Hcomp")








