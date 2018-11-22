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
            
        
        
    # Matrix_Creation_From_Database
    def create_S (self):
        # First we create the pseudoS matrix (if it does not exist) which has the following structure:
        #                                        Number_aqueous_primary_sp   Number_sorption_primary_sp    Number_aqueous_secondary_sp   Number_sorption_secondary_sp    
        #                         n_aqueousR1   |                                                                                                                   |
        #               pseudoS = nRn           |                                                                                                                   |
        #               	        n_sorptionR1  |                                            Stoichiometric values                                                  |
        #                         nRm           |                                                                                                                   |
        #
        #
        #  Remark: pseudoS is a matrix that is almost the sorption stoichiometric matrix. 
        #  The order of the columns is given by the Number_aqueous_primary_sp + Number_sorption_primary_sp + Number_aqueous_secondary_sp + Number_sorption_secondary_sp
        #  The order of the rows is first number of aqueous reactions followed by the number of the sorption reactions.
        if not hasattr(self, 'pseudoS'):
            self.create_pseudo_S()
        
        # Now the electrostatic variables must be added. These variables are treated as chemical species. They will be introduced between Number_sorption_primary_sp  and  Number_aqueous_secondary_sp.
        #
        # Each primary sorption class should have an attribute called type_sorption. The attribute will determine the number of surface potential variables that must be added to the stoichiometric matrix.
        # -CCM will add only one.
        #
        #
        # for the number of rows. Reactions that are aqueous have 0 has stoichiometric value. The stoichiometric values for the added surface potential species is obtained by the type of sorption and b the stoichiometric_value and the charge.
        if not hasattr(self, 'pseudoS') or hasattr(self, 'pseudoS'):
            self.create_electro_sorption_stoichiometric_M ()
        
        # defining length and names of columns
        self.S_names_columns = self.names_aq_pri_sp + self.names_sorpt_pri_sp + self.names_elec_sorpt + self.names_aq_sec_sp + self.names_sorpt_sec_sp
        self.S_length_columns = len(self.pseudoS_names_columns)
        # defining length of rows
        self.S_length_rows = len(self.list_aq_reactions) + len(self.list_sorpt_reactions)
        
        pseudo_S = copy(self.pseudoS)
        S_electro = copy(self.S_electro)
        pos_1 = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        
        S = np.concatenate((np.concatenate ((pseudo_S[:,:pos_1], S_electro), axis = 1), pseudo_s[:,pos_1:]), axis = 1)
        
        assert self.S_length_rows == S.shape[0]
        assert self.S_length_columns == S.shape[1]
        
        self.S = S
    
    # The stoichiometric matrix derived from sorption species.
    def create_electro_sorption_stoichiometric_M (self):
        '''
            The function assumes that some variables are already defined
        '''
        # create list of new boltzman surface potential variables from sorption species
        self.names_elec_sorpt = []
        for i in self.length_sorpt_pri_sp:
            if isinstance(self.list_sorpt_pri_sp[i].names_Boltz_psi, str):
                self.names_elec_sorpt.append(self.list_sorpt_pri_sp[i].names_Boltz_psi)
        self.length_names_elec_sorpt = len(self.names_elec_sorpt)
        # Block
        if not hasattr(self, 'pseudoS_length_rows'):
            # self.pseudoS_length_rows = len(self.list_aq_reactions) + len(self.list_sorpt_reactions)
            self.pseudoS_length_rows = self.length_aq_sec_sp + self.length_sorpt_sec_sp
        S_electro = np.zeros((self.pseudoS_length_rows, self.length_names_elec_sorpt))
        col_position = 0
        for i in self.length_sorpt_pri_sp:
            sub_B = self.create_stoichiometric_surfacepotential (self.names_sorpt_pri_sp, self.list_sorpt_pri_sp[i].type_sorption)
            if len(sub_B.shape) == 1:
                S_electro[:, col_position] = sub_B
                col_position += 1 
            elif len(sub_B.shape) == 2:
                old_col_position = col_position
                col_position = col_position + sub_B.shape[1]
                S_electro[:, old_col_position:col_position] = sub_B
        self.S_electro = S_electro
        
    # creates stoichiometric blocks    
    def create_stoichiometric_surfacepotential (self, name_pri_sp, type_sorpt):
        '''
        '''
        if type_sorpt == 'CCM':
            d = np.zeros((self.length_aq_sec_sp + self.length_sorpt_sec_sp))
            for i in range(0, self.length_sorpt_sec_sp):
                if self.list_sorpt_reactions[i].is_species_in_reaction (name_pri_sp):
                    names_species_in_reaction = [*self.list_sorpt_reactions[i].reaction]
                    summ_charges_times_stoichiometric = 0
                    for j in names_species_in_reaction:
                        if j in self.name_primary_species:
                            z = self.list_aq_pri_sp[self.name_primary_species.index(j)].charge
                            n = self.list_sorpt_reactions[i].reaction[j]
                            summ_charges_times_stoichiometric = summ_charges_times_stoichiometric  + (n*z)
                        elif j in self.name_secondary_species:   
                            z = self.list_aq_sec_sp[self.name_secondary_species.index(j)].charge
                            n = self.list_sorpt_reactions[i].reaction[j]
                            summ_charges_times_stoichiometric = summ_charges_times_stoichiometric  + (n*z)
                    d[self.length_aq_sec_sp + i] = summ_charges_times_stoichiometric
        return d
    
    
    
            
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
    
    # set names_electrostatic_variables
    def set_names_electrostatic_variables (self, names_elsctrostatic_var):
        '''
            The name of the electrostatic potentials that must be taken into account.
            Preferible define them using create_electro_sorption_stoichiometric_M
            Since the names_elsctrotatic_var and the amount in general should be related to a surface
        '''
        self.names_elec_sorpt = names_elsctrostatic_var
        self.length_names_elec_sorpt = len(self.names_elec_sorpt)
    
    # set the stoichiometric matrix given by
    def set_electro_sorption_stoichiometric_M (self, S_electro):
        '''
            The S matrix defined having as columns the surface variable potentials and as rows the reactions.
            Preferible define them using create_electro_sorption_stoichiometric_M
        '''
        self.S_electro =  S_electro
    
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