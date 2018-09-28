# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 11:58:38 2018

@author: Lützenkirchen, Heberling, Jara
"""
import activity_function_module as afm

class Species:
    # Constructor
    def __init__(self, name):
        self.name = name  
# Aquous species (Aq_species) is a daughter class of Species
class Aq_Species (Species):
    # Constructor
    def __init__(self, name):
        '''Sets the inputted parameter name (e.g 'NaCl', 'Mg+2', 'Cl-') into the property name.'''
        Species.__init__(self, name)
        self.f_act_coef = 'Davis'
    # Setters
    def set_charge(self, charge):
        '''Sets the inputted parameter charge (e.g 0, 2, -1) into the property charge.'''
        self.charge = charge
    
    def set_gfw(self, molar_weight):
        '''Sets the inputted parameter molar_weight (e.g 58.44, 48.61, 35.4530) into the property gfw.
           It is assumed that the unit is g/mol.
        '''
        self.gfw = molar_weight
    
    
    def set_ionsizeparameter(self, a):
        '''
            The ion size paramete is set. The units of the input must be angstroms --> 1e-8 agnstrom = 1 cm
            and it will be converted to cm
        '''
        self.ionsizepar = a*1e-8
    
    def set_deviationionsizeparameter(self, b):
        '''
            The deviation ion size paramete is set.
            I am unsure about the units of this parameter. So far I leave it like that, it might change in the future
        '''
        self.devionsizepar = b
        
    def Set_f_activity_coefficient (self, String_name_func):
        '''
        Davis, Deby-Huckel
        '''
        
        self.f_act_coef = String_name_func
        
    #Getters
    
    # Special things
    def it_is_charge(self, boolean):
        self.charge_element = boolean
        
        
        ################### Calculations ############################
        
    def log_coefficient_activity(self, ionic_strength, c=0, A=0, B = 0, a = 0, b = 0, z = 0):
        if hasattr(self, 'charge'):
            z = self.charge
        if hasattr(self, 'ionsizepar'):
            a = self.ionsizepar
        if hasattr(self, 'devionsizepar'):
            b = self.devionsizepar
        log_act_coef = afm.f_log_list(self.f_act_coef, ionic_strength, c, A, B, a, b, z)
        return log_act_coef
    
    def dgamma_dionicstrength(self, actcoeff, ionic_strength, z=0, A=0):
        '''
            returns the derivative value of the activity regarding the inoic strength (dgamma/dionic_strength)
        '''
        if hasattr(self, 'charge'):
            z = self.charge
        return afm.dgamma_dionicstrength(self.f_act_coef, actcoeff, ionic_strength, z, A)