# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:00:55 2019

@author: DaniJ
"""
import numpy as np
import scipy as sp

from matplotlib import pyplot as plt
######################################## functions ##############################################
def plotElement(X, Y, Sep, xlab, ylab, Title):
    "It is somehow specific for this example"
    plt.figure()
    for i in range(0, len(X)):
        s,t = get_colorline_thiscase(Sep[i])
        plt.plot(X[i],Y[i],s,label=t,markersize=12)
        #plt.legend(loc='best')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(Title)

def get_colorline_thiscase(value):
    if value-1e-8<np.finfo(float).eps:
        col="g."
        tol=str(1e-8)
    elif value-1e-7<np.finfo(float).eps:
        col="b."
        tol=str(1e-7)
    elif value-1e-6 <np.finfo(float).eps:
        col="k."
        tol=str(1e-6)
    elif value-1e-5<np.finfo(float).eps:
        col="c."
        tol=str(1e-5)
    elif value-1e-4<np.finfo(float).eps:
        col="r."
        tol=str(1e-4)
    elif value-1e-3<np.finfo(float).eps:
        col="y."
        tol=str(1e-3)
    elif value-1e-2<np.finfo(float).eps:
        col="m."
        tol=str(1e-2)
    elif value-1e-1<np.finfo(float).eps:
        col="w."
        tol=str(1e-1)
    return col, tol


################################# Param #############################################################
    
#tol_v = np.load('tol_vec_v3_lnX.npy')
#Xarr = np.load('X_arr_v3_lnX.npy')
#Carr = np.load('C_arr_v3_lnX.npy')
tol_v = np.load('tol_vec_v2.npy')
Xarr = np.load('X_arr_v2.npy')
Carr = np.load('C_arr_v2.npy')
T_H = np.linspace(-3,-11.2,42)
T_H = 10**T_H



################### calling functions ##########################################################

plotElement(-np.log10(T_H), -np.log10(Xarr[:,0]), tol_v, "-log[H] component (mol/L)", "-log[H+] (mol/L)", "H+ vs. Hcomp")
plotElement(-np.log10(T_H), np.log(Xarr[:,4]), tol_v, "-log[H] component (mol/L)", "Boltzmanf_psi0", "Bf_psi0 vs. Hcomp")
plotElement(-np.log10(T_H), np.log(Xarr[:,5]), tol_v, "-log[H] component (mol/L)", "Boltzmanf_psiC", "Bf_psiC vs. Hcomp")
plotElement(-np.log10(T_H), np.log(Xarr[:,6]), tol_v, "-log[H] component (mol/L)", "Boltzmanf_psiA", "Bf_psiA vs. Hcomp")
plotElement(-np.log10(T_H), np.log(Xarr[:,7]), tol_v, "-log[H] component (mol/L)", "Boltzmanf_psid", "Bf_psid vs. Hcomp")


#plt.figure()
#plt.plot(0,0,"g.",label="1e-8",markersize=12)
#plt.plot(1,1,"b.",label="1e-7",markersize=12)
#plt.plot(2,2,"k.",label="1e-6",markersize=12)
#plt.plot(3,3,"c.",label="1e-5",markersize=12)
#plt.plot(4,4,"r.",label="1e-4",markersize=12)
#plt.legend(loc='best')


#plt.figure()
#plt.plot(-np.log10(T_H), Carr[:,1],label="Cl+",markersize=12)
#plt.figure()
#plt.plot(-np.log10(T_H), Carr[:,6],label='SOH2Cl',markersize=12)


#plotElement(-np.log10(T_H), Carr[:,1], tol_v, "-log[H] component (mol/L)", "Cl- (mol/l)", "Bf_psi0 vs. Hcomp")
#plotElement(-np.log10(T_H), Carr[:,6], tol_v, "-log[H] component (mol/L)", "SOH2Cl", "Bf_psiC vs. Hcomp")