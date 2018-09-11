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
sys.path.append(file)