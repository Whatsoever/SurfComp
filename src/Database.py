# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:55:33 2018

@author: Lützenkirchen, Heberling, Jara
"""

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
        
    def ser_reaction_list(self, list_reactions):
        '''
        list_reactions is a list containing reaction classes. 
            Currently, e.g.:
                list_reactions = [Reaction1, Reaction2, ..., Reaction_m]   
            where each Reaction in the list is a Reaction class.
            
            
            The n_reactions and the list of reactions is set up.
        '''
        self.list_reactions = list_reactions
        self.n_reactions = len(list_reactions)
    
    # Matrix_Creation_From_Database
    def create_S (self):
        self.check_consistency_species()
    
    # Check that the database is consistent
    def check_consistency_species(self):
        # Check that n_species == len(list_species) == len(names_primary_species) + len(names_secondary_species)
        # n_species must be equal to len(list_species), because of set_species_list
        assert self.n_species == (len(self.names_primary_species) + len(self.name_secondary_species)), \
        "The length of the n_species on the sum of list and the names_primary_species and secondary species is not equal."
        # set in the list of primary and secondary species is different
        assert len(set(self.name_primary_species).intersection(set(self.name_secondary_species)))== 0, "The names on the primary and the secondary species are repeated."        
        # look that all the species in the list_spicies belong to primary_species or to secondary_species
        for i in range(0, self.n_species):
            assert self.list_species[i].name in self.name_primary_species or self.list_species[i].name in self.name_secondary_species, \
            "Species %r is not in the primary or secondary list" %self.list_species[i].name