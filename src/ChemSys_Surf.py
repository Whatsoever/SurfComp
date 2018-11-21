# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:00:16 2018

@author: DaniJ

This module is supossed to contain the algorithms and information of Chemical speciation plus sorption.
It is a daughter of Database_SC but it can be used without a database. 
[If feasible (question of time), I will keep it apart]
"""

from Database_SC import Database_SC
import numpy as np

class ChemSys_Surf (Database_SC):
    '''
        ChemSys is a daughter class from Database_SC which is a daughter class of Database. Hence, they depend on these parameters.        
            #Note for myself and other contributors, if you add or delete properties or methods of the class, documeted it here. Otherwise, it is a little caos (regarding my own experience)
            properties:
                Faraday_constant
                temperature
                dielectric_constant
                A_activitypar
                B_activitypar
            methods:
                set_Faraday_constant
                set_temperature
                set_dielectric_constant
                calculate_dielectric_constant
                calculate_A_activitypar
                calculate_B_activitypar
                
        NOTE: Remark that ChemSys_Surf is a daughter class from Database_SC. Therefore, in order to create the pseudo S matrix (The stoichiometric matrix that does not contain the surface potential as unknown). Methods like ...
                ... set_names_aq_primary_species (names_aq_pri_sp), set_names_aq_secondary_species (names_aq_sec_sp), set_names_sorpt_primary_species (names_sorpt_pri_sp), set_names_sorpt_secondary_species (names_sorpt_sec_sp), set_aq_list_pri_class (list_aq_pri_sp), ...
                ... set_aq_list_sec_class (list_aq_sec_sp) can be used and must be used. However, it has to be check that the input given is in accordance with the own system, that can be done by ???????????
    '''
    # Constructor
    def __init__(self):
        self.Faraday_constant = 96485.3328959 # C/mol
        self.temperature = (273.15+25)   # It assumed that initially we are at T=25°C and we assume atmospheric pressure for dielectric and other constants
        self.universal_gas_constant = 8.314472  # J/(K*mol)
        self.calculate_dielectric_constant()
        pass
    
    # Instantiation of main attributes
    def define_system_from_input_and_database (self, database, n_aq_prim, list_aq_val, name_sorpt_pri, List_pri_sorpt_class = None):
        '''
            Given a database, the list of aqueous primary species, the list of aqueous values for the components associated to the primary species, the list of sorption of primary species 
            The system is defined.
            As extra List_pri_sorpt_class is given to update some species. list_sorpt_pri == list_pri_sorpt_class[i].name for i in length.
        '''
        # check that list_sorpt_pri is coherent with List_pri_sorpt_class
        assert len(n_aq_prim) == len(list_aq_val), \
        "The length of the aqueous primary species and the aqueous component values is not equal."
        if List_pri_sorpt_class is not None:
            assert len(name_sorpt_pri) == len(List_pri_sorpt_class), \
                "The length of the sorption primary species and the sorption list classes is not equal."
            for i in range(0, len(name_sorpt_pri)):
                assert i == name_sorpt_pri.index(List_pri_sorpt_class[i].name), 'The name or order of the list of names of sorption primary species and the list of classes of sorption primary species is not coherent.'
        # Instantiation of main attributes (Although not necessary, it is useful to keep sense)
        names_aq_pri_sp = n_aq_prim
        names_aq_sec_sp = []
        list_aq_pri_sp = []
        list_aq_sec_sp = []
        list_aq_reactions = []
        names_sorpt_pri_sp = name_sorpt_pri
        names_sorpt_sec_sp = []
        if List_pri_sorpt_class is not None:
            list_sorpt_pri_sp = List_pri_sorpt_class
        else:
            list_sorpt_pri_sp = []
        list_sorpt_sec_sp = []
        list_sorpt_reactions = []
    
        #  Drawn the list_aq_pri_sp & list_sorpt_pri_sp(if necessary) from Database
        index_list_pri_aq = self.search_index_list_classlist (names_aq_pri_sp, database.names_aq_pri_sp)
        for i in index_list_pri_aq:
            list_aq_pri_sp.append(database.list_aq_pri_sp[i])
        if List_pri_sorpt_class is None:
            index_list_sorpt = self.search_index_classlist_list (names_sorpt_pri_sp, database.names_sorpt_pri_sp)
            for i in index_list_sorpt:
                list_sorpt_pri_sp.append(database.list_sorpt_pri_sp[i])
        
        # Obtain list_aq_reactions, list_aq_sec_sp and names_aq_sec_sp  from names_aq_pri_sp
        index_aq_reactions, names_aq_sec_sp = self.search_index_list_listdictionaryreactions (names_aq_pri_sp, database.list_aq_reactions)
        
        index_list_sec_aq = self.search_index_list_classlist (names_aq_sec_sp, database.names_aq_sec_sp)
        for i in index_list_sec_aq:
            list_aq_sec_sp.append(database.list_aq_sec_sp[i])
        
        for i in index_aq_reactions:
            list_aq_reactions.append(database.list_aq_reactions[i])
        
        # Obtain list_sorpt_reactions, list_sorpt_sec_sp and names_sorpt_sec_sp  from names_aq_pri_sp + names_aq_sec_sp + names_sorpt_pri_sp
        index_sorpt_reactions, names_sorpt_sec_sp = self.search_index_list_listdictionaryreactions (names_aq_pri_sp + names_aq_sec_sp + names_sorpt_pri_sp, database.list_sorpt_reactions)
        
        index_list_sec_sorpt = self.search_index_list_classlist (names_sorpt_sec_sp, database.names_sorpt_sec_sp)
        for i in index_list_sec_sorpt:
            list_sorpt_sec_sp.append(database.list_sorpt_sec_sp[i])
        
        for i in index_sorpt_reactions:
            list_sorpt_reactions.append(database.list_sorpt_reactions[i])
            
        # Instantiation of main variables, hence definition of system to study
        self.set_names_aq_primary_species (names_aq_pri_sp)
        self.set_names_aq_secondary_species (names_aq_sec_sp)
        self.set_names_sorpt_primary_species ( names_sorpt_pri_sp)
        self.set_names_sorpt_secondary_species (names_sorpt_sec_sp)
        self.set_aq_list_pri_class (list_aq_pri_sp)
        self.set_aq_list_sec_class (list_aq_sec_sp)        
        self.set_sorpt_list_pri_class (list_sorpt_pri_sp)
        self.set_sorpt_list_sec_class (list_sorpt_sec_sp)
        self.set_aq_reactions_list (list_aq_reactions)
        self.set_sorpt_reactions_list (list_sorpt_reactions)

        self.set_vector_aqueous_component_value(list_aq_val)    
            
            
            
            
            
    def search_index_list_classlist (self, list1, list2):
        '''
            The function returns a list of indices of the position of list1 in list2. --> E.g. list1 =[a c], list2 = [a b c d] function returns listindices = [1,3]
            Precondition1: list1 <= list2
            Precondition2: list1 is completely include in list2. Otherwise an error occurs
        '''
        assert len(list1) <= len(list2), "List of species in the chemical system must be equal or smaller than the list os primary species on the database"
        list_indices = []  
        for i in list1:
            # appends the index of the list2 that coincide with list1.
            list_indices.append(list2.index(i))
        return list_indices
            
    def search_index_list_listdictionaryreactions (self, list1, list_dictionaries):
        '''
            The function returns two list. One with the indices of the reactions that occur in the ChemSys_Surf according to the inputed dictionary, and the other the secondary species in each reaction. 
            Both, list are in agremment. e.g. l_ind_reaction = [0, 4, 6, 9], l_secondary_species = ['A' 'B' 'C' 'F']  From reaction 0 of the database the secondary species obtained is A, from 6 is C, and so on.
        '''
        index_reactions = []
        name_aq_sec_sp  = []      
        
        for i in range(0, len(list_dictionaries)):
            temp_dict = list_dictionaries[i]
            temp_dict_list_keys = list(temp_dict.reaction.keys())
            n_s = 0
            for j in temp_dict_list_keys:
                count = list1.count(j)
                if count != 1 and count != 0:
                    raise ValueError('[ChemSys class, method Index_ReactionsinDatabase] It seems that the name_primary_species property is wrong.')
                elif count == 0:
                    n_s += 1
                    n_s_name = j
            if n_s == 1:
                index_reactions.append(i)
                name_aq_sec_sp.append(n_s_name)
                
        return   index_reactions, name_aq_sec_sp
            
    # Creating first pseudoS
    
    #Setters
    # aqueous component vector
    def set_vector_aqueous_component_value(self, list_aq_val):
        '''
            The value of vector
        '''
        self.aq_u_vector = list_aq_val
    
    
    
    
    # Faraday constant    
    def set_Faraday_constant (self, new_value):
        '''
            The Faraday constant is instantiated with the class. The Faraday constant has the value 96485.33289(59) C mol−1 [Obtained from WIKI: https://en.wikipedia.org/wiki/Faraday_constant]
            The constant is the relationship between the elementary charge or the magnitude of the charge of an electron ['e'] and the Avogrado constant (The number of particles in a mol) [NA]
            F = e * NA
            
            e ≈ 1.60217662×10−19 C
            NA ≈ 6.02214086×1023 mol−1
            
            Note of one of the authors: I do not think that it should be modified but maybe someone what to play with the value
        '''
        
        self.Faraday_constant = new_value
    
    # Temperature
    def set_temperature(self, new_T):
        '''
            Temperature is supposed to be given in kelvins.
        '''
        
        self.temperature = new_T
     
    # Universal gas constant
    def set_universal_gas_constant (self, r_value):
        '''
            Set the universal gas constant
        '''
        self.universal_gas_constant = r_value
        
     # dielectric constant
    def set_dielectric_constant (self, e_c):
        '''
            Set the dielectric constant of water
        '''
        self.dielectric_constant = e_c
        
        
        
        
        
        
        
    # Calculations

    # Dielectric constant of water    
    def calculate_dielectric_constant(self):
        '''
            Calculates the dielectric constant
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
        '''
        self.dielectric_constant = 2727.586 + 0.6224107*self.temperature - 466.9151*np.log(self.temperature) - (52000.87/self.temperature)
        
    def calculate_A_activitypar (self):
        '''
            Calculates the parameter A of the Debye Hueckel equation
            The units are supossed to be kg^(1/2)/mol^(1/2)
            Actually if you want the L/mol is possible to divide by the square of the density to obtain such value
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
        '''
        A = 1.82483e6*np.sqrt(self.waterdensity)
        B = (self.temperature*self.dielectric_constant)**(3/2)
        self.A_activitypar = A/B
        
    def calculate_B_activitypar (self):
        '''
            Calculates the parameter A of the Debye Hueckel equation
            The units are supossed to be kg^(1/2)/mol^(1/2)*cm
            Actually if you want the L/mol is possible to divide by the square of the density to obtain such value
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
            Here the equation is a bit different than that given in the book. The book takes the equation from 
            Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures; II, Debye Huckel parameters for activity coefficients and relative partial molal properties
            The differences is 10-8 and is related to the fact that they uses angstroms instead of cm
        '''
        A = 50.29158649e8*np.sqrt(self.waterdensity)
        B = np.sqrt(self.temperature*self.dielectric_constant)
        self.B_activitypar = A/B