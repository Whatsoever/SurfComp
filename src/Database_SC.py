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
            methods:
                
    '''
    # Constructor
    def __init__(self):
        pass

    # Set main properties
    def set_names_aq_primary_species (names_aq_pri_sp):
        self.names_aq_pri_sp = names_aq_pri_sp
        self.length_aq_pri_sp = len(names_aq_pri_sp)
    
    def set_names_aq_secondary_species (names_aq_sec_sp):
        self.names_aq_sec_sp = names_aq_sec_sp
        self.length_aq_sec_sp = len(names_aq_sec_sp)
    
    def set_names_sorpt_primary_species (names_sorpt_pri_sp):
        self.names_sorpt_pri_sp = names_sorpt_pri_sp
        self.length_sorpt_pri_sp = len(names_sorpt_pri_sp)
        
    def set_names_sorpt_secondary_species (names_sorpt_sec_sp):
        self.names_sorpt_sec_sp = names_sorpt_sec_sp
        self.lengt_sorpt_sec_sp = len(names_sorpt_sec_sp)