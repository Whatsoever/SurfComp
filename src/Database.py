# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:55:33 2018

@author: Lützenkirchen, Heberling, Jara
"""
import numpy as np
import sympy as sp

class Database:
    # Constructor
    def __init__(self):
        pass
    # Setters
    def set_species_list(self, list_species):
        '''
            list_species is a list containing species classes. 
            Currently, e.g.:
                list_species = [Species1, Species2, ..., Species_n]   
            where each Species in the list is a Species class.
            
            
            The n_species and the list of species is set up.
        '''
        self.list_species = list_species
        self.n_species = len(list_species)
        
        
    def set_primary_species(self, names_primary_species):
        '''Sets the input parameter names_primary_species (e.g n = ['Cl-', 'H2O', 'H+']) in the property name_primary_species'''
        self.name_primary_species = names_primary_species
    
    def set_secondary_species(self, names_secondary_species):
        '''Sets the input parameter names_secondary_species (e.g n = ['Cl-', 'H2O', 'H+']) in the property name_secondary_species'''
        self.name_secondary_species = names_secondary_species
        
    def set_reaction_list(self, list_reactions):
        '''
        list_reactions is a list containing reaction classes. 
            Currently, e.g.:
                list_reactions = [Reaction1, Reaction2, ..., Reaction_m]   
            where each Reaction in the list is a Reaction class.
            
            
            The n_reactions and the list of reactions is set up.
        '''
        self.list_reactions = list_reactions
        self.n_reactions = len(list_reactions)
    
    def create_log_k_vector (self):
        '''Creates a list of the log_k vectors'''
        self.log_k_vector = []
        for i in range(0, len(self.list_reactions)):
            self.log_k_vector.append(self.list_reactions[i].log_k)
    
    def create_charge_vector (self):
        '''Creates a list of the charge values'''
        self.charge_vector = []
        for i in range(0, len(self.list_species)):
            self.charge_vector.append(self.list_species[i].charge)
    
    def create_gfw_vector (self):
        '''Creates a list of the charge values'''
        self.gfw_vector = []
        for i in range(0, len(self.list_species)):
            self.gfw_vector.append(self.list_species[i].gfw)
    
    # Matrix_Creation_From_Database
    def create_S (self):
        self.check_consistency_species()
        # Columns and rows definition of the S matrix
        self.S_columns = self.name_primary_species + self.name_secondary_species
        self.S_rows = range( self.n_reactions)
        
        # Instantiating S matrix
        self.S = np.zeros((len(self.S_rows), len(self.S_columns)))
        
        # The stoichiometric matrix must be fulled with the values of the Reaction classes stored in the list_reactions
        for i in range(0, self.n_reactions):
            d = [*self.list_reactions[i].reaction]
            for j in d:
                self.S[i,self.S_columns.index(j)] = self.list_reactions[i].reaction[j]
        qew,rew = np.linalg.qr(self.S)
        # Check that the S matrix is linear independent and if not output an erorr and point where is the problem
        self.check_Srows_linear_independent()
        self.Remove_Columns_S_zero()
        
    # Check that the database is consistent
    def check_consistency_species(self):
        # Check that n_species == len(list_species) == len(names_primary_species) + len(names_secondary_species)
        # n_species must be equal to len(list_species), because of set_species_list
        assert self.n_species == (len(self.name_primary_species) + len(self.name_secondary_species)), \
        "The length of the n_species on the sum of list and the names_primary_species and secondary species is not equal."
        # set in the list of primary and secondary species is different
        assert len(set(self.name_primary_species).intersection(set(self.name_secondary_species)))== 0, "The names on the primary and the secondary species are repeated."        
        # look that all the species in the list_spicies belong to primary_species or to secondary_species
        for i in range(0, self.n_species):
            assert self.list_species[i].name in self.name_primary_species or self.list_species[i].name in self.name_secondary_species, \
            "Species %r is not in the primary or secondary list" %self.list_species[i].name
            
    # Check that the S matrix is linear independent and if not output an erorr and point where is the problem
    def check_Srows_linear_independent(self):
        rank_S = np.linalg.matrix_rank(self.S)
        if rank_S != len(self.S_rows):
            self.S = None
            raise ValueError('One of your reaction equations is not linear dependent.')
            
    # Remove the columns that are zero from the S matrix and modifies therefore the S matrix and the s_columns
    def Remove_Columns_S_zero(self):
        col_pos = (~self.S.any(axis=0))
        for i in range(0, col_pos.size):
            if col_pos [i] == True:
                self.S = np.delete(self.S, i, 1)
                self.S_columns.pop(i)
    