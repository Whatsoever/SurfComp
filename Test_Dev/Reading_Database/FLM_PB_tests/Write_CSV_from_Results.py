# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:57:32 2019

@author: DaniJ
"""

import numpy
Xarr2 = np.load('X_arr_lnX.npy')
Carr2 = np.load('C_arr_lnX.npy')
numpy.savetxt("X_values.csv", Xarr2, delimiter=",")

numpy.savetxt("C_values.csv", Carr2, delimiter=",")
