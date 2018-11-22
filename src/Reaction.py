# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:22:00 2018

@author: Lützenkirchen, Heberling, Jara
"""

class Reaction:
    # Constructor
    def __init__(self):
        pass
    # Setters
    def set_reaction(self, dic_reaction):
        '''
        In this current version, it is assumed that products have a positive stoichiometric
        and reactants a negative.
        e.g:
            Reaction ---> CO3-2 + 2 H+ = CO2 + H2O
            then:
                dic_reaction = {'CO3-2':-1,  'H+':-1, 'CO2':1, 'H2O':1}
        '''
        self.reaction = dic_reaction;
   
    def set_log_k(self, log_k):
        '''Sets the inputted parameter log_k (e.g -11.45) into the property log_k.'''
        '''
        This comment must be deleted once it is cleared changes in temperature or pressure.
        (It is assumed that the values is at 25°C and 1 atm like Phreeqc values). 
        Although, currently not formula to change the log_k that value of k regarding the temperature or pressure has been implemented.
        '''
        self.log_k = log_k
    
    def is_species_in_reaction (self, species):
        '''
            It looks if a given species is in the reaction. The species is given as a string.
        '''
        list_keys = [*self.reaction]
        return list_keys.count(species) == 1