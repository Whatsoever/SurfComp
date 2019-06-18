# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 22:11:31 2019
!python -m unittest AllTests
"""
import unittest
import numpy as np
from scipy import linalg
import four_layer_model_LNX as flm_lnX

class Test_FLM_LNX(unittest.TestCase):
    def test_mass_action_law(self):
        A = np.array([[2, 2, 2], [1, 1, 1],[2, 1, 3]])
        ln_X=np.array([1,2,3])
        ln_K=np.array([1,2,-1])
        result=flm_lnX.mass_action_law(ln_X, ln_K, A)
        Gold=np.array([13,8,12])
        np.testing.assert_array_equal(result, Gold)

        
    def test_u_componentvector(self):
        A = np.array([[2, 2, 2], [1, 1, 1],[2, 1, 3]])
        C = np.array([2,3,4])
        Gold=np.array([15,11,19])
        result=flm_lnX.u_componentvector(A,C)
        np.testing.assert_array_equal(result, Gold)
        
    def test_surface_charge_edgelayer_flm (self):
        C = 0.5
        psi_L0 = 4
        psi_L1 = 2
        result = flm_lnX.surface_charge_edgelayer_flm(C,psi_L0,psi_L1)
        Gold = 1
        self.assertEqual(result, Gold)
        
    def test_surface_charge_between_layer_flm(self):
        C_left = 0.5 
        C_right = 2 
        psi_mid = 4 
        psi_left = 2
        psi_right = 1
        Gold = 7
        result= flm_lnX.surface_charge_between_layer_flm(C_left, C_right, psi_mid, psi_left, psi_right)
        self.assertEqual(result, Gold)
        
    def test_charge_2_mol (self):
        charge = 2
        s = 4
        a = 5
        F = 5
        result = flm_lnX.charge_2_mol (charge, s, a, F)
        Gold = 8
        self.assertEqual(result, Gold)
        
    def test_surface_charge_diffusive_monovalentelectrolyte (self): 
        R = 8.314   # J/(k*mol)
        T = 298.15  # K
        epsilon = 8.854e-12 #F/m
        epsilon_0 =  80.2
        ionic_strength = 1e-6  # mol/L
        F = 96485.3329                     # the units are: s A / mol
        psi_d = 3   #V
        result = flm_lnX.surface_charge_diffusive_monovalentelectrolyte (R, T, epsilon, epsilon_0, ionic_strength, F, psi_d)
        Gold = -1.348819083717946e+21
        self.assertEqual(result, Gold)
        
        