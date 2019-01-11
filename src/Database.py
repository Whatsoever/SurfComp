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
   # def set_species_list(self, list_species):
        '''
            list_species is a list containing species classes. 
            Currently, e.g.:
                list_species = [Species1, Species2, ..., Species_n]   
            where each Species in the list is a Species class.
            
            
            The n_species and the list of species is set up.
            
            properties:
                names_aq_pri_sp                 e.g. All the properties that start with names... bla bla are supossed to be of type list such as r = ['H2O', 'Oh', ...]
                length_aq_pri_sp                e.g. All the properties with length qre supossed to be the number with the length of the list explained in the line above.
                names_aq_sec_sp
                length_aq_sec_sp
            methods:
            
        '''
        #self.list_species = list_species
        #self.n_species = len(list_species)

           
    def set_aq_list_pri_class (self, list_aq_pri_sp):
        self.list_aq_pri_sp = list_aq_pri_sp
        
    def set_aq_list_sec_class (self, list_aq_sec_sp):
        self.list_aq_sec_sp = list_aq_sec_sp
        
    def set_names_aq_primary_species (self, names_primary_species):
        '''Sets the input parameter names_primary_species (e.g n = ['Cl-', 'H2O', 'H+']) in the property name_primary_species'''
        self.names_aq_pri_sp = names_primary_species
        self.length_aq_pri_sp = len(names_primary_species)
    
    def set_names_aq_secondary_species (self, names_secondary_species):
        '''Sets the input parameter names_secondary_species (e.g n = ['Cl-', 'H2O', 'H+']) in the property name_secondary_species'''
        self.names_aq_sec_sp = names_secondary_species
        self.length_aq_sec_sp = len(names_secondary_species)
        
    def  set_aq_reactions_list(self, list_reactions):
        '''
        list_reactions is a list containing reaction classes. 
            Currently, e.g.:
                list_reactions = [Reaction1, Reaction2, ..., Reaction_m]   
            where each Reaction in the list is a Reaction class.
            
            
            The n_reactions and the list of reactions is set up.
        '''
        self.list_aq_reactions = list_reactions
#        self.n_reactions = len(list_reactions)
    
    def create_log_k_vector (self):
        '''Creates a list of the log_k vectors'''
        self.log_k_vector = []
        for i in range(0, len(self.list_aq_reactions)):
            self.log_k_vector.append(self.list_aq_reactions[i].log_k)
    
    def create_charge_vector (self):
        '''Creates a list of the charge values'''
        self.charge_vector = []
        for i in range(0, len(self.list_aq_pri_sp)):
            self.charge_vector.append(self.list_aq_pri_sp[i].charge)
        for i in range(0, len(self.list_aq_sec_sp)):
            self.charge_vector.append(self.list_aq_sec_sp[i].charge)
    
    
    def create_gfw_vector (self):
        '''Creates a list of the gfw values'''
        self.gfw_vector = []
        for i in range(0, len(self.list_aq_pri_sp)):
            self.charge_vector.append(self.list_aq_pri_sp[i].gwf)
        for i in range(0, len(self.list_aq_sec_sp)):
            self.charge_vector.append(self.list_aq_sec_sp[i].gwf)
            
    
    # Matrix_Creation_From_Database
    def create_S (self):
        self.check_consistency_species()
        # Columns and rows definition of the S matrix
        self.S_names_columns = self.names_aq_pri_sp  + self.names_aq_sec_sp 
        self.S_length_rows =  len(self.list_aq_reactions)
        
        # Instantiating S matrix
        self.S = np.zeros([self.S_length_rows , len(self.S_names_columns)])
        
        # The stoichiometric matrix must be fulled with the values of the Reaction classes stored in the list_reactions
        for i in range(0, len(self.list_aq_reactions)):
            d = [*self.list_aq_reactions[i].reaction]
            for j in d:
                self.S[i,self.S_names_columns.index(j)] = self.list_aq_reactions[i].reaction[j]
        qew,rew = np.linalg.qr(self.S)
        # Check that the S matrix is linear independent and if not output an erorr and point where is the problem
        self.check_Srows_linear_independent()

        
    # Check that the database is consistent
    def check_consistency_species(self):
        # Check that n_species == len(list_species) == len(names_primary_species) + len(names_secondary_species)
        # n_species must be equal to len(list_species), because of set_species_list
        assert self.length_aq_pri_sp + self.length_aq_sec_sp == (len(self.list_aq_pri_sp) + len(self.list_aq_sec_sp)), \
        "The length of the n_species on the sum of list and the names_primary_species and secondary species is not equal."
        # set in the list of primary and secondary species is different
        assert len(set(self.names_aq_pri_sp).intersection(set(self.names_aq_sec_sp)))== 0, "The names on the primary and the secondary species are repeated."        
        # look that all the species in the list_spicies belong to primary_species or to secondary_species and that there are in the correct order
    
        for i in range(0, self.length_aq_pri_sp):
            assert self.list_aq_pri_sp[i].name in self.names_aq_pri_sp, \
            "Species %r is not in the primary list" %self.list_aq_pri_sp[i].name 
        for i in range(0, self.length_aq_sec_sp):    
            assert self.list_aq_sec_sp[i].name  in self.names_aq_sec_sp, "Species %r is not in secondary list" %self.list_aq_sec_sp[i].name   #Actually if this conditions is fullfill the above condition can be removed, I leave it since it check twice (good redundancy?)
    # Check that the S matrix is linear independent and if not output an erorr and point where is the problem
    def check_Srows_linear_independent(self):
        rank_S = np.linalg.matrix_rank(self.S)
        if rank_S != self.S_length_rows:
            self.S = None
            raise ValueError('One of your reaction equations is not linear dependent.')
            

    