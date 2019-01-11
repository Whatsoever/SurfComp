# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 09:40:56 2018

@author: DaniJ
"""
from Database import Database as Parent
import numpy as np
import sympy as sp

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
                pseudoS_names_columns           e.g. A list of the species in the column of the pseudoS matrix such as (following the order): names_aq_pri_sp + names_sorpt_pri_sp + names_aq_sec_sp + names_sorpt_sec_sp
                pseudoS_length_columns          length of the columns of the pseudoS matrix that should be equal to the length of the species in the columns given by pseudoS_names_columns
                pseudoS_length_rows             length of the rows in the pseudoS matrix which should be equal to the number of secondary aqueous and sorption species or the number of aqueous + sorption reactions
                pseudoS                         e.g Stoichiometric matrix without the unknown potential (This part is keep apart because is problem dependent).
                                                     Positives values in the Stoichiometric matrix mean products, negative values mean reactants
                log_k_vector                    e.g. vector of log(K), the K is supossed to be define by K = products/reactants      In the case of sorptions, K is the one given by chemistry and not electrostatic forces.
                                                     The order is electrostatic aqueous reactions then sorption constant
                charge_vector                   e.g. vector of the charge of the aqueous species (z) following the order of first primary species then secondary
                gfw_vector                      e.g. vector of the values for gram formula weight of the aqueous species following the order of first primary species then secondary. Namely, mol/gram
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
                create_pseudo_S    
                check_Srows_linear_independent
                check_consistency_species
                forloop_check_names
                check_main_attributes_defined
                create_log_k_vector
                create_charge_vector
                create_gfw_vector
    '''
    # Constructor
    def __init__(self):
        pass

    # Set main properties
   # def set_names_aq_primary_species (self, names_aq_pri_sp):  ----> It has been created in Database, It can be removed I thing
    #    self.names_aq_pri_sp = names_aq_pri_sp
     #   self.length_aq_pri_sp = len(names_aq_pri_sp)
    
   # def set_names_aq_secondary_species (self, names_aq_sec_sp):  ----> It has been created in Database, It can be removed I thing
    #    self.names_aq_sec_sp = names_aq_sec_sp
     #   self.length_aq_sec_sp = len(names_aq_sec_sp)
    
    def set_names_sorpt_primary_species (self, names_sorpt_pri_sp):
        self.names_sorpt_pri_sp = names_sorpt_pri_sp
        self.length_sorpt_pri_sp = len(names_sorpt_pri_sp)
        
    def set_names_sorpt_secondary_species (self, names_sorpt_sec_sp):
        self.names_sorpt_sec_sp = names_sorpt_sec_sp
        self.length_sorpt_sec_sp = len(names_sorpt_sec_sp)
    
 #   def set_aq_list_pri_class (self, list_aq_pri_sp):  ----> It has been created in Database, It can be removed I thing
  #      self.list_aq_pri_sp = list_aq_pri_sp
        
   # def set_aq_list_sec_class (self, list_aq_sec_sp):   ----> It has been created in Database, It can be removed I thing
    #    self.list_aq_sec_sp = list_aq_sec_sp
        
    def set_sorpt_list_pri_class (self, list_sorpt_pri_sp):
        self.list_sorpt_pri_sp = list_sorpt_pri_sp
        
        
    def set_sorpt_list_sec_class (self, list_sorpt_sec_sp):
        self.list_sorpt_sec_sp = list_sorpt_sec_sp
    
    #def set_aq_reactions_list (self, list_aq_reactions):  ----> It has been created in Database, It can be removed I thing
     #   self.list_aq_reactions = list_aq_reactions
       
    def set_sorpt_reactions_list (self, list_sorpt_reactions):
        self.list_sorpt_reactions = list_sorpt_reactions
        
        
    #  Matrix_Creation from Database
    
    def create_pseudo_S (self):
        '''
            create_pseudo_S it creates the pseudo matrix S, where the potential site is still not consider. This part will be added later.
        '''
        self.check_consistency_species()
        # defining length and names of columns
        self.pseudoS_names_columns = self.names_aq_pri_sp + self.names_sorpt_pri_sp + self.names_aq_sec_sp + self.names_sorpt_sec_sp
        self.pseudoS_length_columns = len(self.pseudoS_names_columns)
        # defining length of rows
        self.pseudoS_length_rows = len(self.list_aq_reactions) + len(self.list_sorpt_reactions)
    
        # Instantiating S matrix
        self.pseudoS = np.zeros((self.pseudoS_length_rows, self.pseudoS_length_columns))
        LT = self.list_aq_reactions + self.list_sorpt_reactions
        # Assigning the stoichiometric values
        for i in range(0, self.pseudoS_length_rows):
            d = [*LT[i].reaction]
            for j in d:
                self.pseudoS[i,self.pseudoS_names_columns.index(j)] = LT[i].reaction[j]
        qew,rew = np.linalg.qr(self.pseudoS)
        
        self.check_Srows_linear_independent()
        
    # Check that the S matrix is linear independent and if not output an erorr and point where is the problem
    def check_Srows_linear_independent(self):
        rank_S = np.linalg.matrix_rank(self.pseudoS)
        if rank_S != self.pseudoS_length_rows:
            self.pseudoS = None
            raise ValueError('One of your reaction equations is not linear dependent.')
            
    def check_consistency_species(self):
        # Check that the necessary attirbutes have been defined
        self.check_main_attributes_defined()
        
        # It is necessary to check length of some attributes
        assert len(self.list_aq_pri_sp) == self.length_aq_pri_sp, \
        "The length of the aqueous primary species and the list of aqueous primary species is not equal."
        assert len(self.list_aq_sec_sp) == self.length_aq_sec_sp, \
        "The length of the aqueous secondary species and the list of aqueous secondary species is not equal."
        assert len(self.list_sorpt_pri_sp) == self.length_sorpt_pri_sp, \
        "The length of the aqueous primary species and the list of aqueous primary species is not equal."
        assert len(self.list_sorpt_sec_sp) == self.length_sorpt_sec_sp, \
        "The length of the aqueous primary species and the list of aqueous primary species is not equal."
        # By convenction/logic or at least by the way how I defining things, the length of the reactions list must be equal to the number of secondary species
        assert len(self.list_aq_reactions) == self.length_aq_sec_sp, \
        "The length of the aqueous reactions list and the aqueous secondary species is not equal."
        assert len(self.list_sorpt_reactions) == self.length_sorpt_sec_sp, \
        "The length of the sorption reactions list and the sorption secondary species is not equal."
        
        # The order of the same attributes must also be check [This checking is somehow redundant - Since this one eliminates the previous one, but the previous one sould be a fast checking in comparison to this one]
        self.forloop_check_names(self.names_aq_pri_sp, self.list_aq_pri_sp)
        self.forloop_check_names(self.names_aq_sec_sp, self.list_aq_sec_sp)
        self.forloop_check_names(self.names_sorpt_pri_sp, self.list_sorpt_pri_sp)
        self.forloop_check_names(self.names_sorpt_sec_sp, self.list_sorpt_sec_sp)
        
        # Matrix independe must also be checked but not at this point
        return None
    
    def forloop_check_names(self, list_str, LObj):
        for i in range(0, len(LObj)):
            assert list_str[i] == LObj[i].name, 'The order of the list of string must be the same that the one of classes for primary/secondary aqueous and sorption species.'
        return None
    
    def check_main_attributes_defined(self):
        assert hasattr(self, 'names_aq_pri_sp'), 'names_aq_pri_sp has not been defined.'
        assert hasattr(self, 'length_aq_pri_sp'), 'length_aq_pri_sp has not been defined.'
        assert hasattr(self, 'names_aq_sec_sp'), 'names_aq_sec_sp has not been defined.'
        assert hasattr(self, 'length_aq_sec_sp'), 'length_aq_sec_sp has not been defined.'
        assert hasattr(self, 'names_sorpt_pri_sp'), 'names_sorpt_pri_sp has not been defined.'
        assert hasattr(self, 'length_sorpt_pri_sp'), 'length_sorpt_pri_sp has not been defined.'
        assert hasattr(self, 'names_sorpt_sec_sp'), 'names_sorpt_sec_sp has not been defined.'
        assert hasattr(self, 'length_sorpt_sec_sp'), 'length_sorpt_sec_sp has not been defined.'
        assert hasattr(self, 'list_aq_pri_sp'), 'list_aq_pri_sp has not been defined.'
        assert hasattr(self, 'list_aq_sec_sp'), 'list_aq_sec_sp has not been defined.'
        assert hasattr(self, 'list_sorpt_pri_sp'), 'list_sorpt_pri_sp has not been defined.'
        assert hasattr(self, 'list_sorpt_sec_sp'), 'list_sorpt_sec_sp has not been defined.'
        assert hasattr(self, 'list_aq_reactions'), 'list_aq_reactions has not been defined.'
        assert hasattr(self, 'list_sorpt_reactions'), 'list_sorpt_reactions has not been defined.'
        return None
    
    # Other things that can be created, and can be interesting or useful
    def create_log_k_vector (self):
        '''Creates a list of the log_k vectors'''
        self.log_k_vector = []
        for i in range(0, len(self.list_aq_reactions)):
            self.log_k_vector.append(self.list_aq_reactions[i].log_k)
        for i in range(0, len(self.list_sorpt_reactions)):
            self.log_k_vector.append(self.list_sorpt_reactions[i].log_k) 
    
    #def create_charge_vector (self):    ----> It has been created in Database, It can be removed I thing
     #   '''Creates a list of the charge values'''
      #  self.charge_vector = []
       # for i in range(0, len(self.list_aq_pri_sp)):
        #    self.charge_vector.append(self.list_aq_pri_sp[i].charge)
        #for i in range(0, len(self.list_aq_sec_sp)):
         #   self.charge_vector.append(self.list_aq_sec_sp[i].charge)
    
    def create_gfw_vector (self):
        '''Creates a list of the charge values'''
        self.gfw_vector = []
        for i in range(0, len(self.list_aq_pri_sp)):
            self.gfw_vector.append(self.list_aq_pri_sp[i].gfw)
        for i in range(0, len(self.list_aq_sec_sp)):
            self.gfw_vector.append(self.list_aq_sec_sp[i].gfw)