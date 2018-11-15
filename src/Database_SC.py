# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 09:40:56 2018

@author: DaniJ
"""
from Database import Database as Parent

class Database_SC (Parent):
    '''
        Database_SC:         #Note for myself and other contributors, if you add or delete properties or methods of the class, documeted it here. Otherwise, it is a little caos (regarding my own experience)
            properties:
                names_aq_pri_sp                 e.g. All the properties that start with names... bla bla are supossed to be of type list such as r = ['H2O', 'Oh', ...]
                length_aq_pri_sp                e.g. All the properties with length qre supossed to be the number with the length of the list explained in the line above.
                names_aq_sec_sp
                length_aq_sec_sp
                names_sorpt_pri_sp
                length_sorpt_pri_sp
                names_sorpt_sec_sp
                lengt_sorpt_sec_sp
                list_aq_pri_sp                  e.g A list of 'Aq_Species' classes with the aqueous primary species
                list_aq_sec_sp
                list_sorpt_pri_sp
                list_sorpt_sec_sp
                list_aq_reactions               e.g. List of aqueous reactions
                list_sorpt_reactions
            methods:
                set_names_aq_primary_species (names_aq_pri_sp):
                set_names_aq_secondary_species (names_aq_sec_sp):
                set_names_sorpt_primary_species (names_sorpt_pri_sp):
                set_names_sorpt_secondary_species (names_sorpt_sec_sp):
                set_aq_list_pri_class (list_aq_pri_sp):
                set_aq_list_sec_class (list_aq_sec_sp):
                set_sorpt_list_pri_class (list_sorpt_pri_sp):
                set_sorpt_list_sec_class (list_sorpt_sec_sp):
                set_aq_reactions_list (list_aq_reactions):
                set_sorpt_reactions_list (list_sorpt_reactions):
                    
    '''
    # Constructor
    def __init__(self):
        pass

    # Set main properties
    def set_names_aq_primary_species (self, names_aq_pri_sp):
        self.names_aq_pri_sp = names_aq_pri_sp
        self.length_aq_pri_sp = len(names_aq_pri_sp)
    
    def set_names_aq_secondary_species (self, names_aq_sec_sp):
        self.names_aq_sec_sp = names_aq_sec_sp
        self.length_aq_sec_sp = len(names_aq_sec_sp)
    
    def set_names_sorpt_primary_species (self, names_sorpt_pri_sp):
        self.names_sorpt_pri_sp = names_sorpt_pri_sp
        self.length_sorpt_pri_sp = len(names_sorpt_pri_sp)
        
    def set_names_sorpt_secondary_species (self, names_sorpt_sec_sp):
        self.names_sorpt_sec_sp = names_sorpt_sec_sp
        self.lengt_sorpt_sec_sp = len(names_sorpt_sec_sp)
    
    def set_aq_list_pri_class (self, list_aq_pri_sp):
        self.list_aq_pri_sp = list_aq_pri_sp
        
    def set_aq_list_sec_class (self, list_aq_sec_sp):
        self.list_aq_sec_sp = list_aq_sec_sp
        
    def set_sorpt_list_pri_class (self, list_sorpt_pri_sp):
        self.list_sorpt_pri_sp = list_sorpt_pri_sp
        
        
    def set_sorpt_list_sec_class (self, list_sorpt_sec_sp):
        self.list_sorpt_sec_sp = list_sorpt_sec_sp
    
    def set_aq_reactions_list (self, list_aq_reactions):
        self.list_aq_reactions = list_aq_reactions
        
    def set_sorpt_reactions_list (self, list_sorpt_reactions):
        self.list_sorpt_reactions = list_sorpt_reactions
        
        
    #  Matrix_Creation from Database
    
    def create_pseudo_S (self):
        '''
            create_pseudo_S it creates the pseudo matrix S, where the potential site is still not consider. This part will be added later.
        '''