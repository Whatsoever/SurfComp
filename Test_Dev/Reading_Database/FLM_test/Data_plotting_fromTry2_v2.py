# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:17:55 2019

@author: DaniJ
"""

import numpy as np
import scipy as sp

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


tol_v = np.load('tol_vec.npy')
Xarr = np.load('X_arr.npy')
Carr = np.load('C_arr.npy')
T_H = np.linspace(-3,-11.2,42)
T_H = 10**T_H


# 

cmap = ListedColormap(['g', 'b', 'k', 'c', 'r'])
norm = BoundaryNorm([1e-8, 1e-7, 1e-6,1e-5,1e-4], cmap.N)

#
x = np.log10(T_H) 
y = np.log10(Xarr[:,0])
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Have to set the actual values used for colormapping separately.
lc = LineCollection(segments, cmap=cmap, norm=norm)
#lc.set_array(z)
lc.set_linewidth(3)

fig1 = plt.figure()
plt.gca().add_collection(lc)
plt.show()