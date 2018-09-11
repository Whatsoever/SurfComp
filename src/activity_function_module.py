# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 13:24:35 2018

@author: DaniJ
"""
import numpy as np

def f_log_list(f_act_coef, ionic_strength, c, A, B, a, b, z):
    log_a = 0
    if f_act_coef == 'water':
        log_a = log_act_water(c)
    elif f_act_coef == 'Davis':
        log_a = log_Davies (ionic_strength, A, z)
    else:
        raise ValueError('[activity_function_module] Function not find it.')
    return log_a

def log_act_water (c):
    '''
        It calculates the log of the activity of water by applying the funtion of Garrels and Christ (1965)
    '''
      # This equation correspond to a_(H_2 O)=1-0.018*∑c_i 
    return numpy.log10(1 - 0.018*np.sum(c)) 

def log_Davies (ionic_strength, A, z):
    '''
        It applies the Davies equation
    '''
    sqrtI = np.sqrt(ionic_strength)
    brackets1 = sqrtI/(1+sqrtI)
    brackets2 = 0.3*ionic_strength
    brackets = brackets1 - brackets2
    
    return -A*z*z*brackets
    