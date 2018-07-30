# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 11:58:38 2018

@author: Lützenkirchen, Heberling, Jara
"""

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
        
    # Setters
    def set_charge(self, charge):
        '''Sets the inputted parameter charge (e.g 0, 2, -1) into the property charge.'''
        self.charge = charge
    
    def set_gfw(self, molar_weight):
        '''Sets the inputted parameter molar_weight (e.g 58.44, 48.61, 35.4530) into the property gfw.
           It is assumed that the unit is g/mol.
        '''
        self.gfw = molar_weight
        
    #Getters