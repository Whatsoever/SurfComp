# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:08:12 2018

@author: Lützenkirchen, Heberling, Jara
"""

class ChemSys:
    # Constructor
    def __init__(self, list_prim, list_val, DatabaseC):
        '''
            list_prim = ['Cl-', 'H2O', 'H+']        It is supossed to be the primary species in the system, it should be in accordance with the database.
            
            list_val  = [5 4 20]                    It is supposed to be the concentration values of the primary species stated in list_prim, namely Cl- = 5. So far units are mol/L          
        
            DatabaseC                               It is supposed to be the database that is going to be used in order to built the stochiometric matrix of these system, and the pertinant secondary species, log_k values, etc
        '''
        # Getting index of list_prim in the database
        indices = self.Index_InputPrimaryspeciesinDatabase ( list_prim, DatabaseC.name_primary_species)
        # Sorting list_prim and list_val according to the indices, AFTER THIS STEP THE LIST BECOMES TUPLE !!!!!!!!!!!
        indices, list_prim, list_val = zip(*sorted(zip(indices, list_prim, list_val)))
        # name_prymary_species is like list_prim but order following the order of the Database
        self.name_primary_species = list(list_prim)                        
        # It is the list of species = [Species1, Species2, Species3] related to name_primary_species
        self.list_primary_species = [DatabaseC.list_species[i] for i in indices]
        
    # Setters
    
    # Searching
    
    def Index_InputPrimaryspeciesinDatabase (self, list1, list2):
        '''
            Note_1: list1 <= list2
            Note_2: list1 is completely include in list2. Otherwise an error occurs
            Note_3: The function returns a list of indices of the position of list1 in list2. --> E.g. list1 =[a c], list2 = [a b c d] function returns listindices = [1,3]
        '''
        list_indices = []  
        for i in list1:
            # appends the index of the list2 that coincide with list1.
            list_indices.append(list2.index(i))
            
        return list_indices