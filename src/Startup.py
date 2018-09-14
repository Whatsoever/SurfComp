"""
Created on Mon Jul 30 16:55:33 2018

@author: Lützenkirchen, Heberling, Jara
"""

import sys
import os

# Getting path
d = os.getcwd()

file = os.path.join(os.getcwd(), os.listdir(os.getcwd())[0])
# Including path of src
sys.path.append(d)
r = os.listdir(d)

for i in range(0, len(r)):
    if os.path.isdir(r[i]) and r[i] != '__pycache__':
        sys.path.append(os.path.join(os.getcwd(),r[i]))