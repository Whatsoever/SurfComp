# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 10:30:04 2019

@author: DaniJ
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def plotElement(X, Y, Sep, xlab, ylab, Title):
    "It is somehow specific for this example"
    plt.figure()
    for i in range(0, len(X)):
        s="g."
        plt.plot(X[i],Y[i],s,markersize=12)
        #plt.legend(loc='best')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(Title)

def plot_2_Element(X, Y1, Y2,xlab=None, ylab=None, Title=None):
    plt.figure()
    plt.plot(X, Y2, "b.",markersize=12)
    plt.plot(X, Y1, "g.",markersize=12)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(Title)

T_H = np.linspace(-3,-11.2,42)
T_H = 10**T_H

Xarr1 = np.load('X_arr_X.npy')
Carr1 = np.load('C_arr_X.npy')

Xarr2 = np.load('X_arr_lnX.npy')
Carr2 = np.load('C_arr_lnX.npy')

#plotElement(-np.log10(T_H), -np.log10(Xarr2[:,0]), 'A', "-log[H] component (mol/L)", "-log[H+] (mol/L)", "H+ vs. Hcomp")
#plot_2_Element(-np.log10(T_H), -np.log10(Xarr1[:,2]), -np.log10(Xarr2[:,2]))

# Comparing between the two different formulation
#for i in range(0,12):
 #   plot_2_Element(-np.log10(T_H), -np.log10(Xarr1[:,i]), -np.log10(Xarr2[:,i]))
 
# Comparing between surface
#d1=[5,6,7,8]
#d2=[9,10,11,12]
#for i in range(0,len(d1)):
 #   plot_2_Element(-np.log10(T_H), -np.log10(Xarr1[:,d1[i]]), -np.log10(Xarr1[:,d2[i]])) 
  #  plot_2_Element(-np.log10(T_H), -np.log10(Xarr2[:,d1[i]]), -np.log10(Xarr2[:,d2[i]])) 