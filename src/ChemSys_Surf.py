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
from scipy import linalg
import scipy.integrate as integrate
from scipy import optimize
#import scipy as sp
class ChemSys_Surf (Database_SC):
    '''
        ChemSys is a daughter class from Database_SC which is a daughter class of Database. Hence, they depend on these parameters.        
            #Note for myself and other contributors, if you add or delete properties or methods of the class, documeted it here. Otherwise, it is a little caos (regarding my own experience)
            properties:
                Faraday_constant
                temperature
                dielectric_constant
                permittivity_free_space
                A_activitypar
                B_activitypar
                universal_gas_constant 
                ionic_strength_constant
                fix_ionic_strength 
                S
                S_electro
                names_elec_sorpt
                length_names_elec_sorpt
                U
                A_Borkovec
                B_Borkovec
                A_Borkovec_columns
                A_Borkovec_rows
                aq_u_vector
                waterdensity
                index_related_sorpt_pri
                
                
            methods:
                set_S
                set_vector_aqueous_component_value
                set_names_electrostatic_variables
                set_electro_sorption_stoichiometric_M
                set_universal_gas_constant
                set_Faraday_constant
                set_temperature
                set_dielectric_constant
                set_constant_ionic_strength
                set_permittivity_free_space
                calculate_dielectric_constant
                calculate_A_activitypar
                calculate_B_activitypar
                calculate_ionic_strength
                calculate_waterdensity
                calculate_u_electro 
                define_system_from_input_and_database
                create_S
                create_U
                remove_electro_mass_from_U
                separte_S_into_S1_and_S2
                create_electro_sorption_stoichiometric_M
                create_stoichiometric_surfacepotential
                search_index_list_classlist
                search_index_list_listdictionaryreactions
                instantiation_step
                speciation_Westall1980_CCM                                  #   NOTE --> probably speciation_Westall1980_CCM, speciation_Westall1980_TLM  can be unified in one algorithm, so far it is kept separated.
                speciation_Westall1980_TLM                                  #
                create_sorpt_vec 
                Boltzman_factor_2_psi 
                Jacobian_Speciation_Westall1980
                print_speciation
                speciation_Borkovec_1983_DLM
                get_z_vector
                calculate_log_activity_coefficient_aq_pri_species
                calculate_log_activity_coefficient_aq_sec_species
                
        NOTE: Remark that ChemSys_Surf is a daughter class from Database_SC. Therefore, in order to create the pseudo S matrix (The stoichiometric matrix that does not contain the surface potential as unknown). Methods like ...
                ... set_names_aq_primary_species (names_aq_pri_sp), set_names_aq_secondary_species (names_aq_sec_sp), set_names_sorpt_primary_species (names_sorpt_pri_sp), set_names_sorpt_secondary_species (names_sorpt_sec_sp), set_aq_list_pri_class (list_aq_pri_sp), ...
                ... set_aq_list_sec_class (list_aq_sec_sp) can be used and must be used. However, it has to be check that the input given is in accordance with the own system, that can be done by ???????????
    '''
    # Constructor
    def __init__(self):
        self.Faraday_constant = 96485.3328959 # C/mol
        self.temperature = (273.15+25)   # It assumed that initially we are at T=25°C and we assume atmospheric pressure for dielectric and other constants
        self.universal_gas_constant = 8.314472  # J/(K*mol)
        self.permittivity_free_space = 8.854187871e-12## Farrads = F --> F/m = C^2/(J*m) ALSO called vacuum permittivity, electri constant or distributed capacitance of the vacuum
        self.calculate_dielectric_constant()
        self.calculate_waterdensity()
        self.calculate_A_activitypar()
        self.calculate_B_activitypar()
        self.ionic_strength_constant = False
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
            
        
    def set_constant_ionic_strength (self, givenvalue):
        '''
            set the ionic_strength to a given value
        '''
        self.ionic_strength_constant = True
        self.fix_ionic_strength = givenvalue
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
        if not hasattr(self, 'S_electro') or not hasattr(self, 'pseudoS_length_rows'):
            self.create_electro_sorption_stoichiometric_M ()
        
        # defining length and names of columns
        self.S_names_columns = self.names_aq_pri_sp + self.names_sorpt_pri_sp + self.names_elec_sorpt + self.names_aq_sec_sp + self.names_sorpt_sec_sp
        self.S_length_columns = len(self.pseudoS_names_columns) + len(self.names_elec_sorpt)
        # defining length of rows
        self.S_length_rows = len(self.list_aq_reactions) + len(self.list_sorpt_reactions)
        
        pseudo_S = self.pseudoS.copy()
        S_electro = self.S_electro.copy()
        pos_1 = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        
        S = np.concatenate((np.concatenate ((pseudo_S[:,:pos_1], S_electro), axis = 1), pseudo_S[:,pos_1:]), axis = 1)
        
        assert self.S_length_rows == S.shape[0]
        assert self.S_length_columns == S.shape[1]
        
        self.S = S
        
    # Creation of the  Component matrix, [Westall does not really make a difference between stoichiometric matrix and U matrix, since somehow they are related]
    def create_U (self):
        if not hasattr(self, 'S'):
            self.create_S ()
            
        S1, S2 = self.separte_S_into_S1_and_S2()
        npri = self.length_aq_pri_sp +self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        I = np.identity(npri)
        Stop=-np.matmul(S1.transpose(), linalg.inv(S2.transpose()))
        U = np.concatenate((I, Stop), axis=1)
        U = self.remove_electro_mass_from_U (U)
        self.U = U
    
    
    # remove_electro_mass_from_U ()
    def remove_electro_mass_from_U (self, U):
        '''
            This methods should be used only in create_U not outside it.
        '''
        npri = self.length_aq_pri_sp +self.length_sorpt_pri_sp
        for i in range(0, self.length_names_elec_sorpt):
            U[npri, npri] = 0
            npri +=  1
        return U
    
    # Separate matrix from Primary and Secondary species
    def separte_S_into_S1_and_S2 (self):
        '''
            Separates primary and Secondary species matrices.
            e.g.:
                            Sp1  Sp1  Sp2                      
                     R1 ||  x11  x12  x13 ||                    || x11 x12 ||           || x11 ||
                S =  R2 ||  x21  x22  x23 ||        in to  S1 = || x21 x22 ||   and S2= || x21 ||
                     R3 ||  x31  x32  x33 ||                    || x31 x32 ||           || x32 ||
        ''' 
        np = self.length_aq_pri_sp +self.length_sorpt_pri_sp + len(self.names_elec_sorpt)
        S1 = self.S[:, 0:np].copy() 
        S2 = self.S[:, np:].copy()
        return S1, S2
    
    
    
    # The stoichiometric matrix derived from sorption species.
    def create_electro_sorption_stoichiometric_M (self):
        '''
            The function assumes that some variables are already defined
        '''
        # create list of new boltzman surface potential variables from sorption species
        self.names_elec_sorpt = []
        self.index_related_sorpt_pri = []
        for i in range(0,self.length_sorpt_pri_sp):
            if hasattr(self.list_sorpt_pri_sp[i], 'type_relation'):                         # related species should be defined in the list_sorpt_pri_sp after the leading species.
                self.index_related_sorpt_pri.append(self.names_sorpt_pri_sp.index(self.list_sorpt_pri_sp[i].type_relation))
            elif isinstance(self.list_sorpt_pri_sp[i].names_Boltz_psi, str):
                self.names_elec_sorpt.append(self.list_sorpt_pri_sp[i].names_Boltz_psi)
            elif isinstance(self.list_sorpt_pri_sp[i].names_Boltz_psi, list):
                for j in range(0, len(self.list_sorpt_pri_sp[i].names_Boltz_psi)):
                        self.names_elec_sorpt.append(self.list_sorpt_pri_sp[i].names_Boltz_psi[j])
        self.length_names_elec_sorpt = len(self.names_elec_sorpt)
        # Block
        if not hasattr(self, 'pseudoS_length_rows'):
            # self.pseudoS_length_rows = len(self.list_aq_reactions) + len(self.list_sorpt_reactions)
            self.pseudoS_length_rows = self.length_aq_sec_sp + self.length_sorpt_sec_sp
        S_electro = np.zeros((self.pseudoS_length_rows, self.length_names_elec_sorpt))
        col_position = 0
        track_dict = {}
        counter = 0
        for i in range(0, self.length_sorpt_pri_sp):
            if hasattr(self.list_sorpt_pri_sp[i], 'type_relation'):                 # related species should be defined in the list_sorpt_pri_sp after the leading species.
                sub_B = self.create_stoichiometric_surfacepotential (self.names_sorpt_pri_sp[i], self.list_sorpt_pri_sp[self.index_related_sorpt_pri[counter]].type_sorption)
                ind_start = track_dict['start_'+ self.names_sorpt_pri_sp[self.index_related_sorpt_pri[counter]]]
                ind_end =track_dict['end_'+ self.names_sorpt_pri_sp[self.index_related_sorpt_pri[counter]]]
                if len(sub_B.shape) == 1:
                    S_electro[:, ind_start:ind_end] = S_electro[:, ind_start:ind_end] + sub_B.reshape(sub_B.shape[0],1) 
                else:
                    S_electro[:, ind_start:ind_end] = S_electro[:, ind_start:ind_end] + sub_B 
                counter += 1
            else:       
                sub_B = self.create_stoichiometric_surfacepotential (self.names_sorpt_pri_sp[i], self.list_sorpt_pri_sp[i].type_sorption)
                if len(sub_B.shape) == 1:
                    S_electro[:, col_position] = sub_B
                    track_dict['start_'+self.names_sorpt_pri_sp[i]] = col_position
                    col_position += 1 
                    track_dict['end_'+self.names_sorpt_pri_sp[i]] = col_position
                elif len(sub_B.shape) == 2:
                    old_col_position = col_position
                    col_position = col_position + sub_B.shape[1]
                    S_electro[:, old_col_position:col_position] = sub_B
                    track_dict['start_'+self.names_sorpt_pri_sp[i]] = old_col_position
                    track_dict['end_'+self.names_sorpt_pri_sp[i]] = col_position
        self.S_electro = S_electro

        
               
        
        
    # creates stoichiometric blocks    
    def create_stoichiometric_surfacepotential (self, name_pri_sp, type_sorpt):
        '''
        '''
        if type_sorpt == 'CCM' or type_sorpt == 'DLM':
            d = np.zeros((self.length_aq_sec_sp + self.length_sorpt_sec_sp))
            for i in range(0, self.length_sorpt_sec_sp):
                if self.list_sorpt_reactions[i].is_species_in_reaction (name_pri_sp):
                    names_species_in_reaction = [*self.list_sorpt_reactions[i].reaction]
                    summ_charges_times_stoichiometric = 0
                    for j in names_species_in_reaction:
                        if j in self.names_aq_pri_sp:
                            z = self.list_aq_pri_sp[self.names_aq_pri_sp.index(j)].charge
                            n = self.list_sorpt_reactions[i].reaction[j]
                            summ_charges_times_stoichiometric = summ_charges_times_stoichiometric  + (n*z)
                        elif j in self.names_aq_sec_sp:   
                            z = self.list_aq_sec_sp[self.names_aq_sec_sp.index(j)].charge
                            n = self.list_sorpt_reactions[i].reaction[j]
                            summ_charges_times_stoichiometric = summ_charges_times_stoichiometric  + (n*z)
                    d[self.length_aq_sec_sp + i] = summ_charges_times_stoichiometric
        elif type_sorpt == 'TLM':
            d = np.zeros(((self.length_aq_sec_sp + self.length_sorpt_sec_sp), 3))
            for i in range(0, self.length_sorpt_sec_sp):
                if self.list_sorpt_reactions[i].is_species_in_reaction (name_pri_sp):
                    names_species_in_reaction = [*self.list_sorpt_reactions[i].reaction]
                    summ_charges_times_stoichiometric_o = 0
                    summ_charges_times_stoichiometric_b = 0
                    for j in names_species_in_reaction:
                        if j in self.names_aq_pri_sp:
                            z = self.list_aq_pri_sp[self.names_aq_pri_sp.index(j)].charge
                            n = self.list_sorpt_reactions[i].reaction[j]
                            if j =='H+' or j == 'OH-':
                                summ_charges_times_stoichiometric_o = summ_charges_times_stoichiometric_o + (n*z) 
                            else:
                                summ_charges_times_stoichiometric_b = summ_charges_times_stoichiometric_b + (n*z) 
                        elif j in self.names_aq_sec_sp: 
                            z = self.list_aq_sec_sp[self.names_aq_sec_sp.index(j)].charge
                            n = self.list_sorpt_reactions[i].reaction[j]
                            if j =='H+' or j == 'OH-':
                                summ_charges_times_stoichiometric_o = summ_charges_times_stoichiometric_o + (n*z) 
                            else:
                                summ_charges_times_stoichiometric_b = summ_charges_times_stoichiometric_b + (n*z) 
                    d[self.length_aq_sec_sp + i, 0] = summ_charges_times_stoichiometric_o
                    d[self.length_aq_sec_sp + i, 1] = summ_charges_times_stoichiometric_b
            
        
        return d
    
    
    def get_z_vector(self):
        z =[]
        for i in range(0, self.length_aq_pri_sp): 
            # if type(self.list_aq_pri_sp[i]) == Aq_Species:
            z.append(self.list_aq_pri_sp[i].charge)
            
        for i in range(0, self.length_aq_sec_sp):
            z.append(self.list_aq_sec_sp[i].charge)
        return z
            
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
    
    # set stoichiometric Matrix
    def set_S (self, S, names_species_columns):
        self.S = S
        self.S_length_rows = S.shape[0]
        self.S_length_columns = S.shape[1]
        self.S_names_columns = names_species_columns
        assert len(names_species_columns) == self.S_length_columns, 'The columns must have the same size that the list of strings containing the name of the species.'
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
        
    def set_permittivity_free_space (self, eo):
        '''
            Set permittivity of the free space, or distributed capacitance of the vacuum or vacuum permittivity etc
            Not recommended to be used. Unless sure of what are you doing
        '''
        self.permittivity_free_space = eo
        
        
        
        
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
    
    def calculate_waterdensity (self):
        '''
            Calculates the density of the water
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
        '''
        Tc = self.temperature - 273.15
        A = (Tc-3.9863)**2
        B = Tc + 288.9414
        C = Tc + 68.12963
        D = (A*B)/(508929.2*C)
        E = 0.011445*np.exp(-374.3/Tc)
        self.waterdensity = 1 - D + E  
        
        
        ############################################################################
        #####      instantiation_step ()
        #####
        #############################################################################
    def instantiation_step (self, type_I=1):
        '''
            
        '''
        if type_I == 1:
            c_ini = np.ones(self.S_length_columns)*1e-3
            
        return c_ini
        
        
        
        
        
        ############################################################################################################################################################
        ################# Speciation and related algorithms ########################################################################################################
        ############################################################################################################################################################
        
        #
    def speciation_Westall1980_CCM (self, tolerance = 1e-6, max_iterations = 100, c_guess = None):
        '''
            Implementation of the algorithm given in "Chemical Equilibrium Including Adsorption on Charged Surfaces" Westall, 1980
            ages 37-to-39
        '''
        # instantiation of unknowns
        if np.any(c_guess == None):
            c_guess = self.instantiation_step (type_I = 1)
        c_n =c_guess
        pos_start_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        pos_end_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        S1, S2 = self.separte_S_into_S1_and_S2()
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        # instantiation variables for loop
        counter_iterations = 0
        err = tolerance + 1
        while err>tolerance and counter_iterations < max_iterations:
            # Calculate U vector [If I am not wrong T_sigma must be calculated at every step, since it depends somehow in the surface potential, and it is unknown]
            u_electro = self.calculate_u_electro(c_n[pos_start_elec:pos_end_elec], c_n)
            T = np.concatenate ((T_chem, u_electro))
            # Calculate f or better said in this specific case Y
            Y = self.U.dot(c_n) - T
            # Calculate Z
            Z = self.Jacobian_Speciation_Westall1980(c_n, pos_start_elec, pos_end_elec)
            # Calculating the diff, Delta_X
            # In the paper Delta_X is X_old - X_new or as they called X_original - X_improved.
            # I am writing X_new- X-old, hence I use -Y instead of Y.
            delta_X = linalg.solve(Z,-Y)
        
            # The error will be equal to the maximum increment
            err = max(abs(delta_X))
        
            # Relaxation factor borrow from Craig M.Bethke to avoid negative values
            max_1 = 1
            max_2 =np.amax(-2*np.multiply(delta_X, 1/c_n[0:pos_end_elec]))
            Max_f = np.amax([max_1, max_2])
            Del_mul = 1/Max_f
        
        
            # Update
            c_n[0:pos_end_elec] = c_n[0:pos_end_elec] + Del_mul*delta_X   # Update primary species
            log_c2 = np.matmul(linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(c_n[0:pos_end_elec])))      # Update secondary
            c_n[pos_end_elec:] =10**log_c2
            counter_iterations += 1
        if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.')
        self.c = c_n
        return c_n
    
    
    def speciation_Westall1980_CCM_v2 (self, tolerance = 1e-6, max_iterations = 100, x = None):
        '''
            Implementation of the algorithm given in "Chemical Equilibrium Including Adsorption on Charged Surfaces" Westall, 1980
            ages 37-to-39
        '''
        # scipy.optimize.newton(func, x0, fprime=None, args=(), tol=1.48e-08, maxiter=50, fprime2=None, x1=None, rtol=0.0, full_output=False, disp=True)[source]
        S1, S2 = self.separte_S_into_S1_and_S2()
        pos_start_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        pos_end_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        
        
        #c_pri = optimize.newton(self.func_newton, x, args = (T_chem, pos_start_elec, pos_end_elec, S1, S2), fprime = self.Jacobian_Speciation_Westall1980_func)
        c_pri = optimize.fsolve(self.func_newton, x, args = (T_chem, pos_start_elec, pos_end_elec, S1, S2), fprime = self.Jacobian_Speciation_Westall1980_func)
        
        
        log_c2 = np.matmul(linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(c_pri)))      # Update secondary
        c2 =10**log_c2
        c_n = np.concatenate ((c_pri, c2))
        
        self.c = c_n
        
        return c_n
     
    def func_newton (self, x, T_chem, pos_start_elec, pos_end_elec, S1, S2):
        '''
            x is the vector of primary species
        '''
        log_c2 = np.matmul(linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(x)))      # Update secondary
        c2 =10**log_c2
        c_n = np.concatenate ((x, c2))
        u_electro = self.calculate_u_electro(x[pos_start_elec:pos_end_elec], c_n)
        T = np.concatenate ((T_chem, u_electro))
        Y = self.U.dot(c_n) - T
        return Y
    def Jacobian_Speciation_Westall1980_func (self, x, T_chem, pos_start_elec, pos_end_elec, S1, S2):
        log_c2 = np.matmul(linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(x)))      # Update secondary
        c2 =10**log_c2
        c_n = np.concatenate ((x, c2))
        return self.Jacobian_Speciation_Westall1980(c_n, pos_start_elec, pos_end_elec)
    
    def speciation_Westall1980_v3 (self, tolerance = 1e-6, max_iterations = 100, Ln_x = None, activity_b = False):
        '''
            Implementation of the algorithm given in "Chemical Equilibrium Including Adsorption on Charged Surfaces" Westall, 1980
            ages 37-to-39.
            That is the third version, here we will try to work with ln(X) as primary species instead of X. Such thing have an effect in the formulation.
            Specifically, the Newton-Rapshon jacobian of the system should become symetric (I am not taking into account activity, not sure if using activity and its derivatives the matrix is still symetric)
            The activity_b is just a boolean that if true, the speciaiton of the secondary species in  is done by substitution of
        '''
        # scipy.optimize.newton(func, x0, fprime=None, args=(), tol=1.48e-08, maxiter=50, fprime2=None, x1=None, rtol=0.0, full_output=False, disp=True)[source]
        S1, S2 = self.separte_S_into_S1_and_S2()
        pos_start_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        pos_end_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        
        lnK = self.log_k_vector/np.log10(np.e)      # Changing the base from log_10 to ln (log_e)
        
        #c_pri = optimize.newton(self.func_newton, x, args = (T_chem, pos_start_elec, pos_end_elec, S1, S2), fprime = self.Jacobian_Speciation_Westall1980_func)
        ln_c_pri = optimize.fsolve(self.residual_fun_v3, Ln_x, args = (lnK, T_chem, pos_start_elec, pos_end_elec, S1, S2, activity_b), fprime = self.Jacobian_Residual_fun_v3)
        
        
        ln_c2 = np.matmul(linalg.inv(S2), lnK - np.matmul(S1, ln_c_pri))
        
        
        c1 = np.exp(ln_c_pri)
        c2 = np.exp(ln_c2)
        c_n = np.concatenate ((c1, c2))
        
        self.c = c_n
        
        return c_n
    
    def residual_fun_v3 (self, x, lnK, T_chem, pos_start_elec, pos_end_elec, S1, S2, activity_b):
        '''
            This functions is not the 3rd version of an old function but it is related to the speciation_Westall1980_v3.
            The main algorithm uses the algorithms and formulas that can be found on the Westall paper but for the unknown variables it relies on ln X variables instead of just X variables.
            The function that I must bild is still Y = U*c -T
            what changes is how the c parameters are obtained. Before we assumed that our indepent variable was a sort of concentration, now the variable is exactly the lnX of the sort of concentration
            Hence the ecuation for c is translated into:
                c = exp(lnKi+sum(aik*lnX))
                but since we are using the stoichiometric matrix the relationship will be
                lnC2 = inv(S2)*lnk - inv(S2)*S1*lnX
                and c is the concatenation of c = exp(lnX) and exp(lnC2)
        '''
        if activity_b == False:
            c_n = self.speciation_no_activity_v3 (lnK, S1, S2, x)
        elif activity_b == True:
            c_n = self.speciation_activity_v3 (lnK, S1, S2, x)
        c1 = np.exp(x)    
        u_electro = self.calculate_u_electro(c1[pos_start_elec:pos_end_elec], c_n)
        T = np.concatenate ((T_chem, u_electro))
        Y = self.U.dot(c_n) - T
        return Y
    
    def Jacobian_Residual_fun_v3 (self, x, lnK, T_chem, pos_start_elec, pos_end_elec, S1, S2, activity_b):
        '''
            This functions is not the 3rd version of an old function but it is related to the speciation_Westall1980_v3.
        '''
        if activity_b == False:
            c_n = self.speciation_no_activity_v3 (lnK, S1, S2, x)
        elif activity_b == True:
            c_n = self.speciation_activity_v3 (lnK, S1, S2, x)
            
        return self.Jacobian_Speciation_Westall1980_modification_lnX (c_n, pos_start_elec, pos_end_elec)
    
    def speciation_no_activity_v3 (self, lnK, S1, S2, x):
        ln_c2 = np.matmul(linalg.inv(S2), lnK - np.matmul(S1, x))
        c1 = np.exp(x)
        c2 = np.exp(ln_c2)
        c_n = np.concatenate ((c1, c2))
        return c_n
    
    def speciation_activity_v3 (self, lnK, S1, S2, x):
        c_1 = np.exp(x)
        c_2 = np.zeros(S2.shape[1])
        
        c_2 = self.subfunction_of_speciation_activity_v3  (c_2, c_1, lnK, S1, S2)
        c_2 = optimize.fixed_point(self.subfunction_of_speciation_activity_v3, c_2, args = (c_1, lnK, S1, S2))
        
        #
      #  tolerance = 1e-8
       # n_max_iterations = 100
        #error = 1
        # I need to implement some sort of Picard method
        #c_1 = np.exp(x)

        #c_2 = np.zeros(S2.shape[1])
        
       # c_k = self.subfunction_of_speciation_activity_v3  (c_2, c_1, lnK, S1, S2)
        #counter = 0
       
        #while error > tolerance and counter < n_max_iterations:
         #   c_k1 = self.subfunction_of_speciation_activity_v3  (c_k, c_1, lnK, S1, S2)
          #  error = max(abs(c_k1-c_k))
           # print(error)
            #c_k = c_k1.copy()
            #counter += 1
        #if counter >= n_max_iterations:
         #   raise ValueError('Max number of iterations surpassed in speciation_activity_v3 (self, lnK, S1, S2, x.')
        c_n = np.concatenate((c_1, c_2))
        return c_n
    
    def subfunction_of_speciation_activity_v3 (self, c_2, c_1, lnK, S1, S2):
        c_a_pri = c_1[:self.length_aq_pri_sp]
        c_a_sec = c_2[:self.length_aq_sec_sp]
        ionic_strength = self.calculate_ionic_strength (np.concatenate((c_a_pri, c_a_sec)))
        log_a_coeff_aq_pri_sp = self.calculate_log_activity_coefficient_aq_pri_species (ionic_strength)
        a_coeff_aq_pri_sp = 10**(log_a_coeff_aq_pri_sp)
        log_a_coeff_aq_sec_sp = self.calculate_log_activity_coefficient_aq_sec_species (ionic_strength) 
        a_coeff_aq_sec_sp = 10**(log_a_coeff_aq_sec_sp)
        if 'H2O' in self.names_aq_pri_sp:
            ind = self.names_aq_pri_sp.index('H2O')
            c_a_pri_t = np.delte(c_a_pri, ind)
            a_coeff_aq_pri_sp [ind] = 1-(0.018*np.sum(np.concatenate ((c_a_pri_t, c_a_sec))))
        elif 'H2O' in self.names_aq_sec_sp:
            ind = self.names_aq_sec_sp.index('H2O')
            c_a_sec_t = np.delte(c_a_sec, ind)
            a_coeff_aq_sec_sp [ind] = 1-(0.018*np.sum(np.concatenate ((c_a_pri, c_a_sec_t))))
        
        c_1[:self.length_aq_pri_sp] = c_1[:self.length_aq_pri_sp]*a_coeff_aq_pri_sp
        ln_c1_a1 = np.log(c_1)
        ln_c2_a2 = np.matmul(linalg.inv(S2), lnK - np.matmul(S1, ln_c1_a1))
        ln_c2_a2[:self.length_aq_sec_sp] = ln_c2_a2[:self.length_aq_sec_sp] - np.log(a_coeff_aq_sec_sp)
        c_2 = np.exp(ln_c2_a2)
        print(c_2)
        return c_2
    
    def Jacobian_Speciation_Westall1980_modification_lnX (self, C, n_aq_plus_n_sorpt, n_primaryspecies):
        '''
            The jacobian matrix following an implementation based on the algorithm of  Westall (1980) 
            "Chemical equilibrium Including Adsorption on Charged Surfaces"
            Pages 37-to-39
            It is assumed that C is order first with the primary species and then with the secondary species such as C = [C1 C2]
            This function is identical to Jacobian_Speciation_Westall1980 but it has been modified considering the lnX the unknown variable.
            That means that the derivation of the residual function for the Newton-Raphson process is done by lnC1 (or LnX) and not C1 (or X)
            
            primary function:
                zjk = sum(aij*aik*Ci/Xk)  becomes now zjk = sum(aij*aik*Ci)
            For CCM:
                z_psipsi =  sum(aij*aipsi*Ci/Xpsi) + (s*a*C*R*T)/(F*F*Xpsi)
                becomes now
                z_psipsi =  sum(aij*aipsi*Ci) + (s*a*C*R*T)/(F*F)
            For TLM:
                
                
        '''
        # The first part treats all terms as it was a normal speciation
        Z = np.zeros((n_primaryspecies, n_primaryspecies))
        for i in range(0, n_primaryspecies):
            for j in range(0, n_primaryspecies):
                Z[i,j]= np.matmul(np.multiply(self.U[i,:], self.U[j,:]), C)
    
        # According to the point 2 of Table III of Westall the term C*sa/F*RT/Funknwon must be added to the electrostatic part
        # I am supposing here that all the sorption phases are CCM
        for i in range(0, self.length_sorpt_pri_sp):
            pos_unknown_vector = n_aq_plus_n_sorpt
            # I am supposing here that the sorption phases are CCM
            if self.list_sorpt_pri_sp[i].type_sorption == 'CCM':
                D1 = self.universal_gas_constant*self.temperature
                D2 = self.Faraday_constant
                F = ((self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant)
                Z[pos_unknown_vector,pos_unknown_vector] = Z[pos_unknown_vector, pos_unknown_vector] + (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
                pos_unknown_vector += 1
            # I am supposing here that the sorption phases are TLM
            elif self.list_sorpt_pri_sp[i].type_sorption == 'TLM':
                
                D1 = self.universal_gas_constant*self.temperature
                D2 = self.Faraday_constant
                D3 = self.Faraday_constant
                D4 = self.Faraday_constant
                F = ((self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant)
                # O-plane
                # plane 0 - 0
                Z[pos_unknown_vector,pos_unknown_vector] = Z[pos_unknown_vector, pos_unknown_vector] + (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
                # plane 0 - b
                Z[pos_unknown_vector,pos_unknown_vector+1] = Z[pos_unknown_vector, pos_unknown_vector+1] - (self.list_sorpt_pri_sp[i].C1*F)*(D1/D3)
                # plane 0 - d
                # plane b - 0
                Z[pos_unknown_vector + 1,pos_unknown_vector] = Z[pos_unknown_vector + 1,pos_unknown_vector] - (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
                # plane b - b
                Z[pos_unknown_vector + 1,pos_unknown_vector + 1] = Z[pos_unknown_vector + 1,pos_unknown_vector + 1] + ((self.list_sorpt_pri_sp[i].C1+self.list_sorpt_pri_sp[i].C2)*F)*(D1/D3)
                # plane b - d
                Z[pos_unknown_vector + 1,pos_unknown_vector + 2] = Z[pos_unknown_vector + 1,pos_unknown_vector + 2] - (self.list_sorpt_pri_sp[i].C2*F)*(D1/D4)
                # plane d - 0
                # plane d - b
                Z[pos_unknown_vector + 2,pos_unknown_vector + 1] = Z[pos_unknown_vector + 2,pos_unknown_vector + 1] - (self.list_sorpt_pri_sp[i].C2*F)*(D1/D3)
                # plane d - d
                # The part below correspond to the paper, which is wrong and must be deleted, once all part agree.
                    # A = -F/(2*R*T)
                #param = self.Faraday_constant/(2*(self.universal_gas_constant*self.temperature))
                #A = -param
                #
                #pos_C = self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt
                #C_aq = np.concatenate((C[:self.length_aq_pri_sp], C[pos_C : (pos_C + self.length_aq_sec_sp)]))
                #
                #I = self.calculate_ionic_strength(C_aq)
                #B = np.sqrt(8*self.permittivity_free_space*self.dielectric_constant*self.universal_gas_constant*self.temperature*I)
                #psi_d = self.Boltzman_factor_2_psi(C[pos_unknown_vector+2])
                #par_C = param*psi_d
                #C = np.cosh(par_C)
                #F_d = A*B*C
                #Z[pos_unknown_vector + 2,pos_unknown_vector + 2] = F_d + (self.list_sorpt_pri_sp[i].C2*F)*(D1/D4)
                pos_C = self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt
                C_aq = np.concatenate((C[:self.length_aq_pri_sp], C[pos_C : (pos_C + self.length_aq_sec_sp)]))
                I = self.calculate_ionic_strength(C_aq)
                B = np.sqrt(8*self.permittivity_free_space*self.dielectric_constant*self.universal_gas_constant*self.temperature*I)
                B_half = B/2
                C = np.cosh(-np.log(C[pos_unknown_vector+2])/2)
                F_d = C*B_half
                Z[pos_unknown_vector + 2,pos_unknown_vector + 2] = F_d + (self.list_sorpt_pri_sp[i].C2*F)*(D1/D4)
                
                pos_unknown_vector +=3
                
        return Z
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def speciation_Westall1980_TLM (self, tolerance = 1e-6, max_iterations = 100, c_guess = None):
        '''
            Implementation of the algorithm given in "Chemical Equilibrium Including Adsorption on Charged Surfaces" Westall, 1980
            ages 37-to-39
        '''
        # instantiation of unknowns
        if np.any(c_guess == None):
            c_guess = self.instantiation_step (type_I = 1)
        c_n = c_guess
        pos_start_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        pos_end_elec = self.length_aq_pri_sp + self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        
        S1, S2 = self.separte_S_into_S1_and_S2()
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        # instantation variables loop
        counter_iterations = 0
        err = tolerance + 1
        while err>tolerance and counter_iterations < max_iterations:
            # Calculate U vector [If I am not wrong T_sigma must be calculated at every step, since it depends somehow in the surface potential, and it is unknown]
            u_electro = self.calculate_u_electro(c_n[pos_start_elec:pos_end_elec], c_n)
            T = np.concatenate ((T_chem, u_electro))
            # Calculate f or better said in this specific case Y
            Y = self.U.dot(c_n) - T
            # Calculate Z
            Z = self.Jacobian_Speciation_Westall1980(c_n, pos_start_elec, pos_end_elec)
            # Calculating the diff, Delta_X
            # In the paper Delta_X is X_old - X_new or as they called X_original - X_improved.
            # I am writing X_new- X-old, hence I use -Y instead of Y.
            delta_X = linalg.solve(Z,-Y)
            #delta_X = sp.sparse.linalg.gmres(Z,-Y)
            #delta_X = delta_X[0]
            #print(delta_X)
            # The error will be equal to the maximum increment
            err = max(abs(delta_X))
            print(err)
            # Relaxation factor borrow from Craig M.Bethke to avoid negative values
            max_1 = 1
            max_2 =np.amax(-2*np.multiply(delta_X, 1/c_n[0:pos_end_elec]))
            Max_f = np.amax([max_1, max_2])
            Del_mul = 1/Max_f
        
        
            # Update
            c_n[0:pos_end_elec] = c_n[0:pos_end_elec] + Del_mul*delta_X   # Update primary species
            log_c2 = np.matmul(linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(c_n[0:pos_end_elec])))      # Update secondary
            c_n[pos_end_elec:] =10**log_c2
            counter_iterations += 1
        if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.')
        self.c = c_n
        return c_n







    def speciation_Borkovec_1983_DLM (self, tolerance = 1e-6, max_iterations = 100, c_guess = None, A_Borkovec = None, names_col = None, names_row = None ):      
        '''
            Implementation of the algorithm given in "Solution of the poisson-boltzman equation for surface excesses of ions in the diffuse layer at the oxide-electrolyte interface" Borkovec 1983
            There are some parts of this algorithm that are not clear for me, hence I will try to implement it as it is given in the paper.
        '''
        # modified matrices must be given:
        if A_Borkovec == None and not hasattr(self, 'A_Borkovec'):
            self.create_A_Borkovec()
        A = self.A_Borkovec
        if names_col == None:
            name_col = self.A_Borkovec_columns
        if names_row == None:
            name_row = self.A_Borkovec_rows
        
        #  The upper part can be expanded to add more outside inputs (Maybe later)
        # for equaiton 20, I need the right K
        S1, S2 = self.separte_S_into_S1_and_S2()
        l_k_comp = np.matmul(linalg.inv(S2),self.log_k_vector)
        
        
        K_eqn20_bulk = np.concatenate((np.zeros(self.length_aq_pri_sp), l_k_comp[:self.length_aq_sec_sp]))
        K_eqn20_surface = np.concatenate((np.zeros(self.length_sorpt_pri_sp), l_k_comp[self.length_aq_sec_sp:]))
        K_eqn20 = np.concatenate((K_eqn20_bulk, K_eqn20_surface))
        # Borkovec_1983- QUOTE (pag. 333) : To circumvent these difficulties one can use an iterative procedure consisting of an initial step to establish electroneutrality in the bulk, and then alternately (i) recomputing g with 
        #                                   the electroneutrality condition fulfilled, and ii) using the constant values of g in solving the equilibrium problems.
        # instantation variables loop
        counter_iterations = 0
        err = tolerance + 1
        ''' 
            Borkovec_1983 - QUOTE (page. 334) --> The initial step is made by treating the asymmetric electrolyte as a symmetric electrolyte of the same ionic strength, using eqn. (16) to evaluate g, and solving the equilibrium problem defined by eqns. (20)
            and (22). There is of course no requirement for electroneutrality when evaluating g by eqn. (16). Once the equilibrium problem is solved, albeit with only approximate values of g, the electroneutrality condition in the bulk is fulfilled, and corrected values
            of g can be evaluated from eqn. (11)
        '''
        # The values of the vector X must be instantiated
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        T = np.concatenate((T_chem, np.zeros(self.length_names_elec_sorpt)))
        # concentration of the components_ initial instantiation
        Xd = 1.1;   # <-- This part might be changed by an instantiation function
        X = np.concatenate((T_chem, np.array([Xd])))
        # initial guess, concentration species
        c = K_eqn20 + np.matmul(A,np.log10(X))      # This part here is equation 20
        c = 10**c
        z_vec = self.get_z_vector()
        
        ''' First part according to Borkovec 1983 - page 333-334, solving assuming symmetric electrolyte '''
        
        while err>tolerance and counter_iterations < max_iterations:
            I = self.calculate_ionic_strength(c[:self.length_aq_pri_sp + self.length_aq_sec_sp]) 
            # g must be calculated to create the matrix B
            g_vec = self.calculate_g_vec_Borkovec_1983_eqn_16(I, X[-1])
            # Now that vector g (assuming symmetrical electrolyte) --> I can build the B matrix and find Y
            B = self.create_B_Borkovec(A, g_vec)
            # Calculating Y. The Y is given in equation 22 in Borkovec(1983)
            Y = np.matmul(B.transpose(), c) - T
            # Now the jacobian must be created
            Z = self.create_jacobian_Borkovec_1983_symm(A, B, c, X, I, z_vec ,g_vec)
            delta_X = linalg.solve(Z,-Y)
            #print(delta_X)
            # The error will be equal to the maximum increment
            err = max(abs(delta_X))
            # Relaxation factor borrow from Craig M.Bethke to avoid negative values
            max_1 = 1
            max_2 =np.amax(-2*np.multiply(delta_X, 1/X))
            Max_f = np.amax([max_1, max_2])
            Del_mul = 1/Max_f
            
            # Update
            X = X + Del_mul*delta_X   # Update primary species
            c = K_eqn20 + np.matmul(A,np.log10(X))      # This part here is equation 20
            c = 10**c
            counter_iterations += 1
        if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.')
        
        X_o = X.copy()
        c_o = c.copy()
        
        ''' Second part, assuming no symmetric electrolyte
            This part today (14/12/2018) is to me not completely clear. Hence, I will see how these approach works.
            I notice that calculating g_vec_o is not equal to the old g_vec value:
                DISCUSS IT with Heberling und Luetzenkirchen
        '''
        g_vec_o = self.calculate_g_vec_Borkovec_1983_eqn_11(z_vec, c_o[:self.length_aq_pri_sp + self.length_aq_sec_sp], X_o[-1])  # Necessary for equation 36 of Borkovec 1983
        g_vec_o = np.array(g_vec_o)
        dg_dXd_vec_o = self.dg_dXd_vec_eqn_11(z_vec, c_o[:self.length_aq_pri_sp + self.length_aq_sec_sp], X_o[-1])         # Necessary for equation 36 of Borkovec 1983
        dg_dXd_vec_o = np.array(dg_dXd_vec_o)
        # instantation variables loop
        counter_iterations = 0
        err = tolerance + 1  
        while err>tolerance and counter_iterations < max_iterations:
            # g must be calculated to create the matrix B
            g_vec = g_vec_o + dg_dXd_vec_o*(X[-1]-X_o[-1])
            # Now that vector g (assuming asymmetrical electrolyte) --> I can build the B matrix and find Y
            B = self.create_B_Borkovec(A, g_vec)
            # Calculating Y. The Y is given in equation 22 in Borkovec(1983)
            Y = np.matmul(B.transpose(), c) - T
            # Now the jacobian must be created
            Z = self.create_jacobian_Borkovec_1983_asymm( A, B, c, X, z_vec, g_vec)
            delta_X = linalg.solve(Z,-Y)
            # The error will be equal to the maximum increment
            err = max(abs(delta_X))
            # Relaxation factor borrow from Craig M.Bethke to avoid negative values
            max_1 = 1
            max_2 =np.amax(-2*np.multiply(delta_X, 1/X))
            Max_f = np.amax([max_1, max_2])
            Del_mul = 1/Max_f
            
            # Update
            X = X + Del_mul*delta_X   # Update primary species
            c = K_eqn20 + np.matmul(A,np.log10(X))      # This part here is equation 20
            c = 10**c
            counter_iterations += 1
        if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.')    
        self.c_Borkovec = c
        return c
    
    def dg_dXd_vec_eqn_11(self, z_vec, cb, Xd):
        '''
            In eqn 36 of Borkovec there is a term evaluated at c_o and Xd_o which is necessary for the calculation of the g factors.
            The variable to integrate is in the integrand limit values of the equation.
            The way to develop it can be found here: https://math.stackexchange.com/questions/716596/derivative-of-definite-integral
            basically means:
                int[a(x),b(x)] f(t) dt = F(a(x), b(x))
                dF/dx = (dF/da)*(da/dx) - (dF/db)*(db/dx) = f(a(x))*(da/dx) - f(b(x))*(db/dx)
                so for our specific case:
                    a(Xd) = Xd  --> da/dXd = 1
                    b(Xd) = 1 --> db/dXd = 0
                    and f(t) will be the integrand of equation 11
        '''
        dg_dXd = []
        if Xd-1 >= 0 :
            b = 1
        else:
            b = -1
        sa_F = self.list_sorpt_pri_sp[0].sp_surf_area*(self.list_sorpt_pri_sp[0].solid_concentration_or_grams/self.Faraday_constant)
        alpha = self.alpha_Borkovec_1983() 
        partA = sa_F*b*alpha
        for i in range(0, len(z_vec)):
            zi = z_vec[i]
            partB = self.integrand_fun_Borkovec_1983_eqn_11(Xd, zi, z_vec, cb)
            dg_dXd.append(partA*partB)
        return dg_dXd
    

    
    def calculate_g_vec_Borkovec_1983_eqn_11 (self, z_vec, cb, Xd):
        '''
            This function should give a result value to the equation 11 stated in Borkovec 1983, If the parameters given are correct.
        '''
        g = []
        tol = 1e-4
        if Xd-1 >= 0 :
            b = 1
        else:
            b = -1
        sa_F = self.list_sorpt_pri_sp[0].sp_surf_area*(self.list_sorpt_pri_sp[0].solid_concentration_or_grams/self.Faraday_constant)
        alpha = self.alpha_Borkovec_1983() 
        partA = sa_F*b*alpha
        
        for i in range(0, len(z_vec)):
            zi = z_vec[i]
            partB = integrate.quad(self.integrand_fun_Borkovec_1983_eqn_11, 1, Xd, args = (zi,z_vec, cb))
            if partB[1] > tol:
                raise ValueError('equation 11, integration of integrand high numerical error')
            g.append(partA*partB[0])
        return g
    
    
    #def integrand_fun_Borkovec_1983_eqn_11 (self, x, zi,z_vec, cb):
     #   a = (x**zi)-1
      #  b= 0
       # for i in range(0, len(z_vec)):
        #    b = b + cb[i]*((x**z_vec[i])-1)
        #b = x*x*b
        #return a/b
        #https://scicomp.stackexchange.com/questions/30715/how-to-cope-with-the-following-singularity?noredirect=1#comment56672_30715
    def integrand_fun_Borkovec_1983_eqn_11 (self, x, zi,z_vec, cb):
        '''
            External help has been provided, here the link for this help. Maybe it should be asked to a mathematician working in this area.
            https://scicomp.stackexchange.com/questions/30715/how-to-cope-with-the-following-singularity?noredirect=1#comment56672_30715
            Actually the guy who provided the answer is : Lutz Lehmann from Humboldt-Universität zu Berlin (He is a mathematician). So I assume it is ok.
        '''
        a = self.term_integrand_fun_Borkovec_1983_eqn_11 (x, zi)
        b= 0
        for i in range(0, len(z_vec)):
            b = b + cb[i]*self.term_integrand_fun_Borkovec_1983_eqn_11(x, z_vec[i])
        b = x*x*b
        #b = (1e-20+max(0,b))**0.5
        b = abs(b)**0.5
        #print(b)
        #print(x)
        return a/b
    
    def term_integrand_fun_Borkovec_1983_eqn_11(self, X, z):
        if abs(X-1)>1e-8:
            return X**z-1       # If I am not close to zero I return (X^z)-1
        return z*(X-1)*(1+(z-1)*(X-1)/2.0*(1+(z-2)*(X-1)/3.0))
    
    def create_jacobian_Borkovec_1983_asymm (self, A, B, c, X, z_vector, g_vector):
        '''
            In the appendix (Borkovec_1983) is written the following, I quote:
                "(ii) For the case of the asymmetric electrolyte with k != d, we need only the first term of eqn.(A2), since in the iteration procedure we define the gi's to be function of Xd only."
            That means that the gi used is the one of equation 36, and hence gi is only function of Xd. Then the quoute continues with:
                "For k = d the derivative needed is simply the integrand of eqn. (11) evaluated at Xd"
        '''
        assert len(z_vector)==len(g_vector), " [create_jacobian_Borkovec_1983_symm] vectors of charge and vector of factor g are not equal. Something must be wrong."
        Nx = len(X)
        Ns = len(c)
        n_iprime = len(g_vector)
        Z = np.zeros((Nx, Nx))
        # There is a term that is repeated in all part of the matrix, also when k = d
        #
        # Sum(bij*aik* ci/Xk)
        # Such term will be the first to be calculated.
        for j in range(0, Nx):
            for k in range(0, Nx):
                for i in range(0, Ns):  #Sum(bij*aik* ci/Xk)
                        Z[j, k] = Z[j, k] + B[i,j]*A[i,k]*(c[i]/X[k])
                if k == (Nx-1):
                    Z[j, k] = Z[j, k] + self.term_A4_Borkovec_asym(n_iprime, z_vector, X[k], c[:n_iprime])
        return Z 
    
    def create_jacobian_Borkovec_1983_symm (self, A, B, c, X, I, z_vector, g_vector):
        '''
            Creating the jacobian for the Newton-Rapshon procedure. The algorithm is given in Borkovec(1983), you need to apply the info of the appendix, plus de info of the paper. 
            Some parameter are slightly tricky, but it seems that everything is ok, except the alpha parameter that I do no trust.
            This jacobian is treated as a symmetric electrolyte, namely equations (A.1 -A.4) of the appendix
            
            dY/dX = dYj/dXk
            
        '''
        assert len(z_vector)==len(g_vector), " [create_jacobian_Borkovec_1983_symm] vectors of charge and vector of factor g are not equal. Something must be wrong."
        Nx = len(X)
        Ns = len(c)
        n_iprime = len(g_vector)
        Z = np.zeros((Nx, Nx))
        # There is a term that is repeated in all part of the matrix, also when k = d
        #
        # Sum(bij*aik* ci/Xk)
        # Such term will be the first to be calculated.
        for j in range(0, Nx):
            for k in range(0, Nx):
                for i in range(0, Ns):  #Sum(bij*aik* ci/Xk)
                        Z[j, k] = Z[j, k] + B[i,j]*A[i,k]*(c[i]/X[k])
                if k != (Nx-1):
                    Z[j, k] = Z[j, k] + self.term_A2_and_A3_Borkovec(n_iprime, j, k, A, c, X,g_vector,z_vector, I) 
                elif k == (Nx-1): #There is one term for all K, except k = d and one for all 
                    Z[j, k] = Z[j, k] + self.term_A4_Borkovec_sym(n_iprime,I, z_vector, X[k], c)
                    
        return Z 
        
    def term_A4_Borkovec_asym(self, n_iprime, z_vector, Xd, c):  
        dg_dXd_vec = self.dg_dXd_vec_eqn_11(z_vector, c, Xd)
        b = sum(cb*z*dg for cb,z,dg in zip(c[:n_iprime],z_vector,dg_dXd_vec))
        return b
        
    def term_A2_and_A3_Borkovec(self, n_iprime, j, k, A, c, X, g_vector,z_vector, I):
        v = 0
        R = 0
        for iprime in range(0, n_iprime):
            v = v + ((z_vector[iprime]**2)/2)*A[iprime, k]*(c[iprime]/X[k])
        for iprime in range(0, n_iprime):
            R = R + c[iprime]*A[iprime, j]*(-g_vector[iprime]/(2*I))*v
        return R
    
    def term_A4_Borkovec_sym(self,  n_iprime, I, z_vector, X_d, c):
        R = 0
        alpha = self.alpha_Borkovec_1983()  
        for iprime in range(0, n_iprime):
            dgiprime_dXd = self.calculate_dg_dXd_Borkovec_1983_eqn_16 (I, alpha, X_d, z_vector[iprime])
            R = R + c[iprime]*z_vector[iprime]*dgiprime_dXd
        return R
    
    def alpha_Borkovec_1983 (self):
        '''
        I THINK THERE IS A TYPO HERE (parameter alpha); I AM USING EQUATION 13 BUT I THINK THE EQUATION IS WRONG: SO I USE A MODIFIED ONE; I MUST ASK THE AUTHORS
        '''
        return np.sqrt((self.dielectric_constant*self.permittivity_free_space)/(2*self.universal_gas_constant*self.temperature)) 
    
    def calculate_g_vec_Borkovec_1983_eqn_16 (self, I, X_d):    
        '''
            It calculates the g factors of the paper of Borkovec (1983) using equation 16. 
            Precondition: The concentration given is order: First primary aqueous species, in the same order that the list of the class. Then secondary species, in the same order as they are saved in the class.
        '''
        g = []
        
        alpha = self.alpha_Borkovec_1983()         
        
        
        for i in range(0, self.length_aq_pri_sp): 
            # if type(self.list_aq_pri_sp[i]) == Aq_Species:
            z = self.list_aq_pri_sp[i].charge
            g.append(self.calculate_g_Borkovec_1983_eqn_16 ( I, alpha, X_d, z))
        for i in range(0, self.length_aq_sec_sp):
        # if type(self.list_aq_sec_sp[i]) == Aq_Species:
            z = self.list_aq_sec_sp[i].charge
            g.append(self.calculate_g_Borkovec_1983_eqn_16 ( I, alpha, X_d, z))
        
        return g
            
    def  calculate_g_Borkovec_1983_eqn_16 (self, I, alpha, X_d, z):
        g = 2*alpha*(1/np.sqrt(I))*((X_d**(z/2))-1)*self.list_sorpt_pri_sp[0].sp_surf_area*(self.list_sorpt_pri_sp[0].solid_concentration_or_grams/self.Faraday_constant)
        return g
    
    def calculate_dg_dXd_Borkovec_1983_eqn_16 (self, I, alpha, X_d, z):
        dg_dXd = 2*alpha*(1/np.sqrt(I))*(z/2)*(X_d**((z/2)-1))*self.list_sorpt_pri_sp[0].sp_surf_area*(self.list_sorpt_pri_sp[0].solid_concentration_or_grams/self.Faraday_constant)
        return dg_dXd
    
    def create_A_Borkovec (self):
        if not hasattr(self, 'U'):
            self.create_U ()
        # HERE THE COLUMNS OF U are defined in the following way: Aqueous primary species (components) + Sorption primary species (components) + Electro components + Aqueous secondary species + Sorption Secondary species
        # The rows are the components which are formulated in the order: Aqueous primary species (components) + Sorption primary species (components) + Electro components
        #
        # Two steps are necessary (Assuming that what is written about U is true):
        #                       1) U must be transpose
        A_temp = self.U.transpose()
        
        # Now A_temp is almost A: In the columns it has: Aqueous primary species (components) + Sorption primary species (components) + Electro components
        # but the rows are in the following order: Aqueous primary species (components) + Sorption primary species (components) + Electro components + Aqueous secondary species + Sorption Secondary species
        # The second step:
        #                   2) The order of the rows must be modified to be: Bulk part, basically aqueous part + Surface part, basically surface species
        #                      Therefore it is decided to reorder in the folllowing way: Bulk: [Aqueous primary species (components)+ Aqueous secondary species] + Surface : [ Sorption primary species (components) +  Sorption Secondary species]
        # Furthermore, the row regarding the electrostatical potential that can be found in A_temp (A row made up of 0s) must be removed.
        n_comp = self.length_aq_pri_sp + self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        ABulk = np.concatenate((A_temp[:self.length_aq_pri_sp, :], A_temp[n_comp : n_comp + self.length_aq_sec_sp, :]))
        ASurface = np.concatenate ((A_temp[self.length_aq_pri_sp: self.length_aq_pri_sp + self.length_sorpt_pri_sp, :], A_temp[n_comp + self.length_aq_sec_sp :, :]))
        
        self.A_Borkovec = np.concatenate((ABulk, ASurface))
        self.A_Borkovec_columns = self.names_aq_pri_sp + self.names_sorpt_pri_sp + self.names_elec_sorpt
        self.A_Borkovec_rows = self.names_aq_pri_sp + self.names_aq_sec_sp + self.names_sorpt_pri_sp + self.names_sorpt_sec_sp
        
    def create_B_Borkovec (self, A, g):
        '''
            In Borkovec (1983), in table 2 is describe how the modified stoichiometry matrix B, must be build using A as model.
            Precondition: A is order according to g, g is order according to first the aqueous primary species followed by the secondary aqueous species.
        '''
        Nsb = self.length_aq_pri_sp + self.length_aq_sec_sp
        Ncb = self.length_aq_pri_sp
        
        B = A.copy()
        
        # Part A
        DG = np.diag(g) + np.identity(Nsb)
        ADG = np.matmul(DG,A[:Nsb, : Ncb])
        B [:Nsb, :Ncb] = ADG
    
        # Part B
        count=0
        for i in range(0, self.length_aq_pri_sp):
            z = self.list_aq_pri_sp[i].charge
            B[i,-1] = z*g[count]
            count += 1
        for i in range(0, self.length_aq_sec_sp): 
            z = self.list_aq_sec_sp[i].charge
            B[count,-1] = z*g[count]
            count += 1
            
        return B
        
        
    def calculate_u_electro (self, unknonw_boltzman_vect, C):
        '''
            T_depends in the surface sorption type somehow
        '''
        T_sigma = []
        pos_point_electro_unknown = 0
        for i in range(0, self.length_sorpt_pri_sp):
            if self.list_sorpt_pri_sp[i].type_sorption == 'CCM':
                x = unknonw_boltzman_vect [pos_point_electro_unknown]
                psi = self.Boltzman_factor_2_psi(x)
                charge_surface = self.list_sorpt_pri_sp[i].C1*psi
                T = charge_surface*((self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant)
                T_sigma.append(T)
                pos_point_electro_unknown += 1
            elif self.list_sorpt_pri_sp[i].type_sorption == 'TLM':
                x = unknonw_boltzman_vect [pos_point_electro_unknown : (pos_point_electro_unknown+3)]
                psi = self.Boltzman_factor_2_psi(x)
                charge_surface_0 = self.list_sorpt_pri_sp[i].C1*(psi[0]-psi[1])
                charge_surface_b = self.list_sorpt_pri_sp[i].C1*(psi[1]-psi[0]) + self.list_sorpt_pri_sp[i].C2*(psi[1]-psi[2])
                charge_surface_d = self.list_sorpt_pri_sp[i].C2*(psi[2]-psi[1])
                #print(charge_surface_0 +charge_surface_b+charge_surface_d)  Check that the sum of charges equals 0
                #charge_surface_d = self.list_sorpt_pri_sp[i].C2*(psi[2]-psi[0])
                D = (self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant
                T_0 = charge_surface_0*D
                T_b = charge_surface_b*D
                T_d = charge_surface_d*D
                # In T_d, it is assigned Y_d equation 14 from Westall
                pos_C = self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt
                C_aq = np.concatenate((C[:self.length_aq_pri_sp], C[pos_C : (pos_C + self.length_aq_sec_sp)]))
                I = self.calculate_ionic_strength(C_aq)
                B = np.sqrt(8*self.permittivity_free_space*self.dielectric_constant*self.universal_gas_constant*self.temperature*I)
                E = np.sinh((self.Faraday_constant*psi[2])/(2*(self.universal_gas_constant*self.temperature)))
                Y = B*E
                #print(Y-T_d)
                # print(Y+T_d)   I have an existencial doubt about these part. 
                #print(charge_surface_d+Y)
                #
                T_sigma.append(T_0); T_sigma.append(T_b); 
                T_sigma.append(Y+T_d)
                #T_sigma.append( charge_surface_d+T_d)
                pos_point_electro_unknown += 3
        #print([T_sigma])    
        return np.array(T_sigma)
        
        
    def Boltzman_factor_2_psi (self, x):
        D = self.universal_gas_constant*self.temperature
        psi = - np.log(x)*(D/self.Faraday_constant)
        return psi

    def  create_sorpt_vec (self):            
        T_sorpt = []
        for i in range(0, self.length_sorpt_pri_sp):
            T_sorpt.append(self.list_sorpt_pri_sp[i].T_solid)
        return T_sorpt
        
    def Jacobian_Speciation_Westall1980 (self, C, n_aq_plus_n_sorpt, n_primaryspecies):
        '''
            The jacobian matrix following an implementation based on the algorithm of  Westall (1980) 
            "Chemical equilibrium Including Adsorption on Charged Surfaces"
            Pages 37-to-39
            It is assumed that C is order first with the primary species and then with the secondary species such as C = [C1 C2]
        '''
        # The first part treats all terms as it was a normal speciation
        Z = np.zeros((n_primaryspecies, n_primaryspecies))
        for i in range(0, n_primaryspecies):
            for j in range(0, n_primaryspecies):
                Z[i,j]= np.matmul(np.multiply(self.U[i,:], self.U[j,:]), (C/C[j]))
    
        # According to the point 2 of Table III of Westall the term C*sa/F*RT/Funknwon must be added to the electrostatic part
        # I am supposing here that all the sorption phases are CCM
        for i in range(0, self.length_sorpt_pri_sp):
            pos_unknown_vector = n_aq_plus_n_sorpt
            # I am supposing here that the sorption phases are CCM
            if self.list_sorpt_pri_sp[i].type_sorption == 'CCM':
                D1 = self.universal_gas_constant*self.temperature
                D2 = self.Faraday_constant*C[pos_unknown_vector]
                F = ((self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant)
                Z[pos_unknown_vector,pos_unknown_vector] = Z[pos_unknown_vector, pos_unknown_vector] + (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
                pos_unknown_vector += 1
            # I am supposing here that the sorption phases are TLM
            elif self.list_sorpt_pri_sp[i].type_sorption == 'TLM':
                
                D1 = self.universal_gas_constant*self.temperature
                D2 = self.Faraday_constant*C[pos_unknown_vector]
                D3 = self.Faraday_constant*C[pos_unknown_vector+1]
                D4 = self.Faraday_constant*C[pos_unknown_vector+2]
                F = ((self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant)
                # O-plane
                # plane 0 - 0
                Z[pos_unknown_vector,pos_unknown_vector] = Z[pos_unknown_vector, pos_unknown_vector] + (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
                # plane 0 - b
                Z[pos_unknown_vector,pos_unknown_vector+1] = Z[pos_unknown_vector, pos_unknown_vector+1] - (self.list_sorpt_pri_sp[i].C1*F)*(D1/D3)
                # plane 0 - d
                # plane b - 0
                Z[pos_unknown_vector + 1,pos_unknown_vector] = Z[pos_unknown_vector + 1,pos_unknown_vector] - (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
                # plane b - b
                Z[pos_unknown_vector + 1,pos_unknown_vector + 1] = Z[pos_unknown_vector + 1,pos_unknown_vector + 1] + ((self.list_sorpt_pri_sp[i].C1+self.list_sorpt_pri_sp[i].C2)*F)*(D1/D3)
                # plane b - d
                Z[pos_unknown_vector + 1,pos_unknown_vector + 2] = Z[pos_unknown_vector + 1,pos_unknown_vector + 2] - (self.list_sorpt_pri_sp[i].C2*F)*(D1/D4)
                # plane d - 0
                # plane d - b
                Z[pos_unknown_vector + 2,pos_unknown_vector + 1] = Z[pos_unknown_vector + 2,pos_unknown_vector + 1] - (self.list_sorpt_pri_sp[i].C2*F)*(D1/D3)
                # plane d - d
                
                ###### ---This part below is what is written in the paper of Westall
                              # A = -F/(2*R*T)
                #param = self.Faraday_constant/(2*(self.universal_gas_constant*self.temperature))
                #A = -param
                #
                #pos_C = self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt
                #C_aq = np.concatenate((C[:self.length_aq_pri_sp], C[pos_C : (pos_C + self.length_aq_sec_sp)]))
                #
                #I = self.calculate_ionic_strength(C_aq)
                #B = np.sqrt(8*self.permittivity_free_space*self.dielectric_constant*self.universal_gas_constant*self.temperature*I)
                #psi_d = self.Boltzman_factor_2_psi(C[pos_unknown_vector+2])
                #par_C = param*psi_d
                #C = np.cosh(par_C)
                #F_d = A*B*C
                ########## This part below is my own assumption, since I think that the equation given by the paper is wrong derivated.
                pos_C = self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt
                C_aq = np.concatenate((C[:self.length_aq_pri_sp], C[pos_C : (pos_C + self.length_aq_sec_sp)]))
                I = self.calculate_ionic_strength(C_aq)
                B = np.sqrt(8*self.permittivity_free_space*self.dielectric_constant*self.universal_gas_constant*self.temperature*I)
                in_cosh = -np.log(C[pos_unknown_vector+2])/2
                F_d = (B/2)*np.cosh(in_cosh)*(1/C[pos_unknown_vector+2])
                
                Z[pos_unknown_vector + 2,pos_unknown_vector + 2] = F_d + (self.list_sorpt_pri_sp[i].C2*F)*(D1/D4)
                
                pos_unknown_vector +=3
                
        return Z

        
    def calculate_ionic_strength (self,c):
        '''
            Calculate the ion strength: The vector C is supossed to be a vector of concentrations that contains first the aqueous primary species followed by the aqueous secondary species. 
            Both primary and secondary species are supossed to be order in the same order that the one of the class, namely self.
        '''
        if self.ionic_strength_constant:
            return self.fix_ionic_strength
        Ionic_s=0
        count = 0
        for i in range(0, self.length_aq_pri_sp): 
            # if type(self.list_aq_pri_sp[i]) == Aq_Species:
            z = self.list_aq_pri_sp[i].charge
            Ionic_s = Ionic_s + c[count]*z*z
            count += 1
        for i in range(0, self.length_aq_sec_sp):
        # if type(self.list_aq_sec_sp[i]) == Aq_Species:
            z = self.list_aq_sec_sp[i].charge
            Ionic_s = Ionic_s + c[count]*z*z
            count += 1
         
        Ionic_s = 0.5*Ionic_s
        return Ionic_s
    
    
    def calculate_log_activity_coefficient_aq_pri_species (self, ionic_strength):
        log_coef_a=np.zeros(self.length_aq_pri_sp)
        for i in range(0, self.length_aq_pri_sp):
            if self.list_aq_pri_sp[i].name == 'H2O':
                # water has not coefficient activity (or it is 0). For water the activity is calculated directly with Garrels and Christ (1965) forumla
                log_coef_a[i] = 0
            else:
                log_coef_a[i] = self.list_aq_pri_sp[i].log_coefficient_activity(ionic_strength, A=self.A_activitypar, B = self.B_activitypar )
        return log_coef_a
    
    def calculate_log_activity_coefficient_aq_sec_species (self, ionic_strength):
        log_coef_a=np.zeros(self.length_aq_sec_sp)
        for i in range(0, self.length_aq_sec_sp):
            if self.list_aq_sec_sp[i].name == 'H2O':
                # water has not coefficient activity (or it is 0). For water the activity is calculated directly with Garrels and Christ (1965) forumla
                log_coef_a[i] = 0
            else:
                log_coef_a[i] = self.list_aq_sec_sp[i].log_coefficient_activity(ionic_strength, A=self.A_activitypar, B = self.B_activitypar )
        return log_coef_a
    
    
    def Bethke_algorithm (self, tolerance = 1e-6, max_n_iterations = 100, tolerance_psi = 1e-6, max_n_iterations_psi = 800, tolerance_big_loop = 1e-6, max_n_iterations_big_loop = 100):
        '''
            These algortihm implementation is based on Geochemical and Biogeochemical reaction modeling from Craig M.Bethke
            section_10.3
        '''
        # Check that water is on the database as primary species and has the first possition
        # So far, for simplicity I leave the thing with the H2O like that but it can be changed.
        ind = self.names_aq_pri_sp.index('H2O')
        if not (ind == 0):
            raise ValueError('[ChemSys/bethke_algorithm] -->To use this algortihm water must be on the first position of the primary species. \n')
        
        '''
        Separates primary and Secondary species matrices.
            e.g.:
                            H2O  sp_i     sp_p    sp_elec      sp_j    sp_q                
                     R1 ||  x11  x1ni     x1np    x1nelec      x1nj    x1nq ||                                    || x11 x12 ||           || x11 ||
                S =  R2 ||  x21  x2ni     x2np    x2nelec      x2nj    x2nq ||                        in to  S1 = || x21 x22 ||   and S2= || x21 ||
                     R3 ||  x31  x3ni     x3np    x3nelec      x3nj    x3nq ||                                    || x31 x32 ||           || x32 ||
        
        where rows R are reactions, and columns are H2O (water), sp_i (aqueous primary species), sp_p (sorption primary species - "uncomplexed" so-labelled by Bethe),
                                        sp_elec (Electrostatic part, Boltzman factor of the sorption charge), sp_j (aqueous secondary species), sp_q (sorption secondary species)
        These part can be separated in to S1 and S2:
            
             || x11  x1ni     x1np    x1nelec ||           || x1nj    x1nq ||
        S1 = || x21  x2ni     x2np    x2nelec ||   and S2= || x2nj    x2nq ||
             || x31  x3ni     x3np    x3nelec ||           || x3nj    x3nq ||
        
        # Or in S1, S2, S3 different ---> These S, separation will be more clear once the algorithm is done.
        
        The U is defined:
                                H2O  sp_i  sp_p  sp_elec   sp_j     sp_q
                  ∑ H2O     ||   1    0     0       0      v_wj     v_wq   ||
        U = V =   ∑ sp_i    ||   0    I     0       0      v_ij     v_iq   ||
                  ∑ sp_p    ||   0    0     I       0        0      v_pq   ||
                  ∑ sp_elec ||   0    0     0       0        0       z_q   ||
                  
         I call U = V because the nomenclature used by the Bethe
         
         As before, the algorithm can be divided in different U parts. For instance the last row, using the algorithm provided by Bethe must be decoupled of the matrix.
        ''' 
        
        # Instantiation of first guesses
        nw = 1                # So we are supossing that the initial amount of water is 1, actually it must change
        mi = (0.9*np.array(self.aq_u_vector[1:]))*nw
        Mp= self.create_sorpt_vec()
        mp = (0.9*np.array(Mp))
        Boltzfactor =  np.ones(self.length_names_elec_sorpt)              # Boltzfactor ==> exp(-psi*F/RT)  Later I will need psi but not yet. Now I only need the boltzman factor for mj and mp guesses
        
        S1, S2 = self.separte_S_into_S1_and_S2()
        S_prima = -np.matmul(linalg.inv(S2),S1)
        log_K_prima = np.matmul(linalg.inv(S2), self.log_k_vector)
        
        ionic_strength = 0  # self.calculate_ionic_strength (c_aqueouspecies)
        log_a_water = np.log10(1-(0.018*np.sum(mi)))              # calculating the log activity of water (water has not coefficient)
        log_a_coeff_aq_pri_sp = self.calculate_log_activity_coefficient_aq_pri_species (ionic_strength)
        log_a_coeff_aq_sec_sp = self.calculate_log_activity_coefficient_aq_sec_species (ionic_strength)
        
        mj_and_mq = self.log_speciation_secondaryspecies_Bethke (log_a_water, log_a_coeff_aq_pri_sp, log_a_coeff_aq_sec_sp, mi, mp, Boltzfactor,S_prima, log_K_prima)
        mj_and_mq = 10**mj_and_mq
        # separation
        mj = mj_and_mq[:self.length_aq_sec_sp]
        mq = mj_and_mq[self.length_aq_sec_sp:]
        
        ## Other parameters that must be calculated and are constant during the loops
        # length values
        length_aq_sorpt_pri = self.length_aq_pri_sp + self.length_sorpt_pri_sp
        length_pri = self.length_aq_pri_sp + self.length_sorpt_pri_sp + self.length_names_elec_sorpt
        
        # matrix that are keep constant through the loops
        U2 = self.U[:, length_pri:]
        
        M = np.concatenate((self.aq_u_vector, Mp))                      # The given component value for aqueous and surface species
        
        WV_and_WP= np.multiply(U2[0,:], U2[1:length_aq_sorpt_pri,:])    # The matrix WV and WP contain the terms v_wj*v_ij, v_wq*v_iq and the terms v_wq*v_pq
        
        
        I = np.identity(length_aq_sorpt_pri-1)                            # The delta of equation (10.33) of Craig M. Bethke's book
        
        
        
        
        Area_v = self.calculate_A_sf_Bethke()
        charge_background_solute = 1
        
        c_minus1 = np.zeros(length_pri + self.length_aq_sec_sp + self.length_sorpt_sec_sp)
    
        ## I have 2 loops: 1) It is a Newton-Raphson methods that must be solved. Once solved, the values are used to calculate a new value of the surface potential.
        # So this 2 loops are contained in other loop
        
        err_big_loop = 1
        counter_iterations_big_loop = 0
        while err_big_loop> tolerance_big_loop and counter_iterations_big_loop < max_n_iterations_big_loop:
            # Ini error parameter
            err = 1 
            counter_iterations = 0;

            # First loop, Newton-Raphson
            while err>tolerance and counter_iterations < max_n_iterations:
                #### Residual vector ####
                # water ##### 
                Jww = 55.5 + np.dot(U2[0,:],mj_and_mq)
                rw = nw*Jww
                
                # aqueous primary species ####
                Jiw = mi + np.matmul(U2[1:self.length_aq_pri_sp,:], mj_and_mq)
                ri = nw*Jiw
                
                # sorption primary species ####
                Jpw = mp + np.matmul(U2[self.length_aq_pri_sp:length_aq_sorpt_pri,:], mj_and_mq)    # Actually it should be only ∑_q v_pq*mq but the terms of v_pj are 0 (At least theoretically, I hope). So the equaiton is ok.
                rp =  nw*Jpw
                
                # assamble
                r = np.concatenate(([rw], ri, rp))
                # R functions evaluated
                R = r - M
                print(R)
                ####### Jacobian matrix #########
                # parameters Jww, Jiw, and Jpw already calculated
                # Jwp and Jwq are calculated together due to the U matrix is defined by using WV_and_WP*mj_and_mq
                jwp_and_jwq = np.matmul(WV_and_WP,mj_and_mq) 
                mi_and_mp = np.concatenate((mi,mp))
                Jwp_and_Jwq = np.multiply((nw/mi_and_mp), jwp_and_jwq)
            
                # If my intuition do not fool me, it should be possible to calculate the part of the Jacobian matrix [equation (10.34) of Craig books'] that comprises the terms Jii', Jip, Jpi, and Jpp'
                # In the same way that Jii when having only speciaiton (section 4 of the book) or in the same way that Jwp and Jwq was possible to be calculated together.
                
                Jii_Jip_and_Jpi_Jpp = nw*I + nw*self.Js_partB_calculation(mi_and_mp, mj_and_mq, U2[1:length_aq_sorpt_pri,:])
                
                # Assembling
                Jw = np.concatenate(([Jww],Jiw,Jpw))
                Jip = np.vstack((Jwp_and_Jwq, Jii_Jip_and_Jpi_Jpp))
                J = np.c_[Jw,Jip]
                
                # Solution Newthon-Raphson
                delta_c = linalg.solve(J,-R)
                err = max(abs(delta_c))
                #print(err)
                # relaxation factor
                max_1 = 1;
                max_2 =(-2*delta_c[0])/nw
                max_3 = np.amax(-2*np.multiply(delta_c[1:self.length_aq_pri_sp], 1/mi))
                max_4 = np.amax(-2*np.multiply(delta_c[self.length_aq_pri_sp:], 1/mp))
                Max_f = np.amax([max_1, max_2, max_3, max_4])
                Del_mul = 1/Max_f
            
                # Update guesses
                nw = nw + Del_mul*delta_c[0]
                mi = mi + Del_mul*delta_c[1:self.length_aq_pri_sp]
                mp = mp + Del_mul*delta_c[self.length_aq_pri_sp:]
                
                # Update secondaries
                
                ionic_strength = self.calculate_ionic_strength (np.concatenate(([55.5], mi, mj)))
                log_a_water = np.log10(1-(0.018*(np.sum(mi)+np.sum(mj))))              # calculating the log activity of water (water has not coefficient)
                log_a_coeff_aq_pri_sp = self.calculate_log_activity_coefficient_aq_pri_species (ionic_strength)
                log_a_coeff_aq_sec_sp = self.calculate_log_activity_coefficient_aq_sec_species (ionic_strength) 
                mj_and_mq = self.log_speciation_secondaryspecies_Bethke (log_a_water, log_a_coeff_aq_pri_sp, log_a_coeff_aq_sec_sp, mi, mp, Boltzfactor,S_prima, log_K_prima)
                mj_and_mq = 10**mj_and_mq
                
                mj = mj_and_mq[:self.length_aq_sec_sp]
                mq = mj_and_mq[self.length_aq_sec_sp:]
                
                counter_iterations += 1
            if counter_iterations >= max_n_iterations:
                raise ValueError('Max number of iterations in chemistry part surpassed.')
            # First loop terminated. Chemistry values establish
            # Second loop, loop of values of the psi must be started
            #### SECOND ITERATION LOOP #####
            # Newton approach for the psi potential
            # Parameter before loop
            a =  np.matmul(U2[length_aq_sorpt_pri,self.length_aq_sec_sp:], mq)
            da_dpsi = (self.Faraday_constant/(self.universal_gas_constant*self.temperature))*np.matmul(np.power(U2[length_aq_sorpt_pri,self.length_aq_sec_sp:],2), mq)
            # Ini error parameter
            err_psis = 1 
            counter_iterations_psi = 0;
            psi = self.Boltzman_factor_2_psi (Boltzfactor)
            while err_psis>tolerance_psi and counter_iterations_psi < max_n_iterations_psi:
                # calculate f and df
                f, df = self.calculate_f_df_psi_equation_10_37_38_Bethke(Area_v,psi, ionic_strength, nw, charge_background_solute)
                # R and DR
                R = f-a
                print(R)
                dR = df + da_dpsi
                # solution
                delta_psis = -R/dR
                # print(delta_psis) 

                # error 
                err_psis = max(abs(delta_psis))
                #print(err_psis)
                ## relaxation factor
                #max_1 = 1;
                #max_2 =(-2*delta_psis)/psi
                #Max_f = np.amax([max_1, max_2])
                #Del_mul = 1/Max_f
                
                # Update guesses
               # psi = psi + Del_mul*delta_psis
                psi = psi +delta_psis
                
                counter_iterations_psi += 1
            if counter_iterations_psi >= max_n_iterations_psi:
                raise ValueError('Max number of psi iterations in Surface potential part surpassed.')
            # New Boltzfactor
            Boltzfactor =  np.exp((-self.Faraday_constant/(self.universal_gas_constant*self.temperature))*psi)
            
            #### Update new values to redo loop 1 and 2
            mj_and_mq = self.log_speciation_secondaryspecies_Bethke (log_a_water, log_a_coeff_aq_pri_sp, log_a_coeff_aq_sec_sp, mi, mp, Boltzfactor,S_prima, log_K_prima)
            mj_and_mq = 10**mj_and_mq
                
            mj = mj_and_mq[:self.length_aq_sec_sp]
            mq = mj_and_mq[self.length_aq_sec_sp:]
            
            c_n = np.concatenate(([55.5087], mi, mp, Boltzfactor, mj, mq))
            
            err_big_loop = max(abs(c_n - c_minus1))                 # Absolute error
            
            
        self.c = c_n
        self.mass_water = nw
        return c_n

    
    
    def Js_partB_calculation(self, mi_and_mp, mj_and_mq, U):
        nc = len(mi_and_mp)
        Js = np.identity(nc)
        for i in range(0, nc):
            for j in range(0, nc):
                Js[i,j]= np.matmul(np.multiply(U[i,:], U[j,:]), (mj_and_mq/mi_and_mp[j]))
        return Js
    
    def log_speciation_secondaryspecies_Bethke (self, log_a_water, log_a_coeff_aq_pri_sp, log_a_coeff_aq_sec_sp, mi, mp, Boltzfactor,S_prima, log_K_prima):
        '''
            Speciation to find the log of aqueous secondary species (mj) and sorption secondary species (mq). 
            log(mj_and_mq) = log_K_prima +  S_prima*[log(mi_and_mp_and_Boltzfactor) + log_activity_coefficient_mi] - log_activity_coefficient_mj
            
            S_prima = np.matmul(np.linalg.inv(S2), self.log_k_vector)
            log_k_prima = np.matmul(np.linalg.inv(S2), self.log_k_vector)
            
        '''
        PartA = log_K_prima
        # making [log(mi_and_mp_and_Boltzfactor) + log_activity_coefficient_mi]
        # THIS algorithm uses Bethke approach hence water is at first position
        log_a_mi = np.log10(mi) + log_a_coeff_aq_pri_sp[1:]       # the log_a_coeff_aq_pri_sp has a 0 enctrance in the first position due to water, it is removed therefore
        log_mp = np.log10(mp)
        log_Boltzfactor = np.log10(Boltzfactor)
        
        log_pri = np.concatenate(([log_a_water], log_a_mi, log_mp, log_Boltzfactor))
        PartB = np.matmul(S_prima, log_pri)
        
        PartC = np.concatenate((log_a_coeff_aq_sec_sp, np.zeros(self.length_sorpt_sec_sp)))
        
        return PartA + PartB - PartC
        
    def calculate_A_sf_Bethke(self):
        '''
            return the area of the different surface.
            It is assumed that if primary sorption surfaces share an electrical potential the total area will be the sum of both.
        '''
        A_sf = []
        counter_related_sp = 0
        for i in range(0, self.length_sorpt_pri_sp):
            # It is assumed that related species have been defined after
            if hasattr(self.list_sorpt_pri_sp[i], 'type_relation'):
                A = self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams
                A_sf[self.index_related_sorpt_pri[counter_related_sp]] = A_sf[self.index_related_sorpt_pri[counter_related_sp]] + A
            else:
                A_sf.append(self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)
        return np.array(A_sf)
    
    def calculate_f_df_psi_equation_10_37_38_Bethke(self, Area_v,psi, ionic_strength, nw, charge_background_solute):
        '''
            It calculates a part of equation 10.37 and 10.38 of 'Goechemical and Biogeochemical reaction modeling' Craig M. Bethke
        '''
        #square_val = np.sqrt(8*self.dielectric_constant*self.permittivity_free_space*self.universal_gas_constant*self.temperature*1000*ionic_strength)
        square_val = np.sqrt(8*self.dielectric_constant*self.permittivity_free_space*self.universal_gas_constant*self.temperature*ionic_strength)
        inner_sincos = ((charge_background_solute*self.Faraday_constant)/(2*self.universal_gas_constant*self.temperature))*psi
        sin_in = np.sinh(inner_sincos)
        cos_in = np.cosh(inner_sincos)
        part_f = (1/(nw*self.Faraday_constant))*Area_v
        part_df = (charge_background_solute/(2*self.universal_gas_constant*self.temperature*nw))*Area_v
        
        f = part_f*square_val*sin_in
        df = part_df*square_val*cos_in
        
        return f,df
        
 #####################################################################################################################################################################################
###################################### Since I think I understand the algorithm but is not working with my implementation, I will re-do it. ##########################################
 #####################################################################################################################################################################################

    def speciation_Westall1980_TLMb (self, tolerance = 1e-6, max_n_iterations = 100, X_guess = None):
        '''
            My first Westall1980 algortihmm did not work. I try a similar implementation to the work of Westall to see if like that it works.
            Implementation of the algorithm given in "Chemical Equilibrium Including Adsorption on Charged Surfaces" Westall, 1980
            ages 37-to-39
        '''
        S1, S2 = self.separte_S_into_S1_and_S2()
        S_prima = -np.matmul(linalg.inv(S2),S1)
        ## Since I am copy the algortihm presented. I must started defining the same variables that they use.
        A = self.TMLb_obtain_matrix_A()
        log_kp = np.matmul(linalg.inv(S2), self.log_k_vector) 
        LOG_K =  np.concatenate((np.zeros(self.length_aq_pri_sp+self.length_sorpt_pri_sp), log_kp))
        
                                                       # The vector contains first 0, related to primary species and then the constants of sec aq reactions and sec sorption reactions
        
        if np.any(X_guess == None):
            X = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
            X = np.concatenate((X, np.ones(self.length_names_elec_sorpt)))         # X is the primary species vector, suposse to be order aq.prim, sorpt.prim, elect.species
        else:
            X = X_guess
        # mass-law  action
        log_C = LOG_K + np.matmul(A,np.log10(X))   # C must contain aq.prim, sorpt.prim, aq.sec, sorpt.sec
        C = 10**log_C
        
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        
        # parm loop
        err = 1 
        counter_iterations = 0;
        while err>tolerance and counter_iterations < max_n_iterations:
            # The T elec should be added
            Tele = self.TMLb_Telec(X[self.length_aq_pri_sp+self.length_sorpt_pri_sp:self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt])
            T = np.concatenate([T_chem, Tele])
            # Mass-Balance Equation
            Y = np.matmul(A.transpose(),C) - T
            y_end = self.TMLb_diffuseregionvalue(C, X[-1])
            Y[-1] = Y[-1] - y_end
            # I have Y. Now I need Z
            Z = self.TMLb_Jacobian(A, C, X)
            
            # solving
            delta_X = linalg.solve(Z,-Y)
            # The error will be equal to the maximum increment
            err = max(abs(delta_X))
            print(err)
            # Relaxation factor borrow from Craig M.Bethke to avoid negative values
            max_1 = 1
            max_2 =np.amax(-2*np.multiply(delta_X, 1/X))
            Max_f = np.amax([max_1, max_2])
            Del_mul = 1/Max_f
            # Update
            X = X + Del_mul*delta_X   # Update primary species
            log_C = LOG_K + np.matmul(A,np.log10(X))   # C must contain aq.prim, sorpt.prim, aq.sec, sorpt.sec
            C = 10**log_C
            counter_iterations += 1                        
        if counter_iterations >= max_n_iterations:
            raise ValueError('Max number of iterations surpassed.')
        self.X = X
        self.C = C
        return C
    
    def TMLb_obtain_matrix_A(self):
        '''
            From the U developed by the code. An "A" matrix similar to the one exposed by the formulation of Westall is created
        '''
        A = self.U
        #                           Species
        #   U = Sum of components[          ]   
        #
        A = np.delete(A, np.s_[self.length_aq_pri_sp+self.length_sorpt_pri_sp : self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt], axis=1)
        
        #
        # Now, I return the transpose matrix of A. Such transpose matrix is composed in the columns by the sum of the components: in principle in the following order aq.prim, sorpt.prim, boltz.factor
        # and the rows are the species, the boltz factor species are not included in the rows. In princile in the following order: aq.prim, sorpt.prim, aq.sec, sorpt.sec
        #
        return A.transpose()
    
    def TMLb_Telec(self,X):
        '''
            If the algorithm work, probably this part will be modified. So far, I am assuming here only one surface.
        '''
        psi = -np.log(X)*((self.temperature*self.universal_gas_constant)/self.Faraday_constant)
        sa_F = (self.list_sorpt_pri_sp[0].sp_surf_area*self.list_sorpt_pri_sp[0].solid_concentration_or_grams)/self.Faraday_constant
        
        T_sigma0 = sa_F*self.list_sorpt_pri_sp[0].C1*(psi[0]-psi[1])
        T_sigmabeta = sa_F*self.list_sorpt_pri_sp[0].C1*(psi[1]-psi[0])+sa_F*self.list_sorpt_pri_sp[0].C2*(psi[1]-psi[2])
        T_sigmad = sa_F*self.list_sorpt_pri_sp[0].C2*(psi[2]-psi[1])
        
        return [T_sigma0, T_sigmabeta, T_sigmad]
    
    def TMLb_diffuseregionvalue(self, C, Xd):
        '''
        '''
        # C contains 
        c_aq = np.concatenate((C[:self.length_aq_pri_sp],C[self.length_aq_pri_sp+self.length_sorpt_pri_sp:self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_aq_sec_sp]))
        ionic_strength = self.calculate_ionic_strength (c_aq)
        partA = np.sqrt(8*self.dielectric_constant*self.permittivity_free_space*self.universal_gas_constant*self.temperature*ionic_strength)
        psid = -np.log(Xd)*((self.temperature*self.universal_gas_constant)/self.Faraday_constant)
        partB = np.sinh((self.Faraday_constant*psid)/(2*self.universal_gas_constant*self.temperature))
        
        return partA*partB
    
    def TMLb_Jacobian(self, A, C, X):
        '''
         The jacobian must be calculated
        '''
        Nx = len(X)
        Ns = len(C)
        Z = np.zeros([Nx, Nx])
        # There is a step that takes place for all the entries in the matrix
        # There is a term that is repeated in all part of the matrix, also when k = d
        #
        # Sum(aij*aik* ci/Xk)
        # Such term will be the first to be calculated.
        for j in range(0, Nx):
            for k in range(0, Nx):
                for i in range(0, Ns):
                    Z[j, k] = Z[j, k] + A[i,j]*A[i,k]*(C[i]/X[k])
        
        # Now specific points are treated
        n_ele_plane0 = self.length_aq_pri_sp+self.length_sorpt_pri_sp
        sa_F = (self.list_sorpt_pri_sp[0].sp_surf_area*self.list_sorpt_pri_sp[0].solid_concentration_or_grams)/self.Faraday_constant
        C1 = self.list_sorpt_pri_sp[0].C1
        C2 = self.list_sorpt_pri_sp[0].C2
        
        RT = self.universal_gas_constant*self.temperature
        ## plane 00
        Z[n_ele_plane0, n_ele_plane0] = Z[n_ele_plane0, n_ele_plane0] + C1*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0]))
        ## plane 0b
        Z[n_ele_plane0, n_ele_plane0+1] = Z[n_ele_plane0, n_ele_plane0+1] - C1*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0+1]))
        ## plane 0d
        Z[n_ele_plane0, n_ele_plane0+2] = 0
        ## plane b0
        Z[n_ele_plane0+1, n_ele_plane0] = Z[n_ele_plane0+1, n_ele_plane0] - C1*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0]))
        ## plane bb
        Z[n_ele_plane0+1, n_ele_plane0+1] = Z[n_ele_plane0+1, n_ele_plane0+1] + (C1+C2)*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0+1]))
        ## plane bd
        Z[n_ele_plane0+1, n_ele_plane0+2] = Z[n_ele_plane0+1, n_ele_plane0+2] - C2*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0+2]))
        ## plane d0
        Z[n_ele_plane0+2, n_ele_plane0] = 0
        ## plane db
        Z[n_ele_plane0+2, n_ele_plane0+1] = Z[n_ele_plane0+2, n_ele_plane0+1] - C2*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0+1]))
        ## plane dd
        ### This part is what is written on the paper of Westall, which is 95% sure wrong.
       # partB =  C2*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0+2]))
        #A = -self.Faraday_constant/(2*RT)
        #c_aq = np.concatenate((C[:self.length_aq_pri_sp],C[self.length_aq_pri_sp+self.length_sorpt_pri_sp:self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_aq_sec_sp]))
        #ionic_strength = self.calculate_ionic_strength (c_aq)
        #B = np.sqrt(8*self.dielectric_constant*self.permittivity_free_space*self.universal_gas_constant*self.temperature*ionic_strength) 
        #psid = -np.log(X[-1])*(RT/self.Faraday_constant)
        #C = np.cosh((self.Faraday_constant*psid)/(2*RT))
        #partA = A*B*C
        #Z[n_ele_plane0+2, n_ele_plane0+2] = partA + partB
       
        ### This is what I think it should be and so far has proven to give the proper results.
        partB =  C2*sa_F*(RT/(self.Faraday_constant*X[n_ele_plane0+2]))
        c_aq = np.concatenate((C[:self.length_aq_pri_sp],C[self.length_aq_pri_sp+self.length_sorpt_pri_sp:self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_aq_sec_sp]))
        ionic_strength = self.calculate_ionic_strength (c_aq)
        B = np.sqrt(8*self.dielectric_constant*self.permittivity_free_space*self.universal_gas_constant*self.temperature*ionic_strength) 
        A =(1/(2*X[n_ele_plane0+2]))
        C = np.cosh((-np.log(X[n_ele_plane0+2])/2))
        partA = A*B*C
        Z[n_ele_plane0+2, n_ele_plane0+2] = partA + partB
        
        return Z
            
##################################################
###############################
########################### Postprocessing
###################################################
        
    def print_speciation (self):
        #ionic_strength = self.calculate_ionic_strength (self.c)
        #log_activity_coefficient = self.calculate_log_activity_coefficient (ionic_strength, self.c)
        #v_activity = self.c*(10**log_activity_coefficient)
        for i in range(0, self.S_length_columns):
            print(self.S_names_columns[i] + ' : ' + str(self.c[i]) + '\n')
            #   print(self.list_species[i].name + ' : ' + str(v_activity[i]) + '\n')
    
    def print_speciation_Borkovec (self):
        #ionic_strength = self.calculate_ionic_strength (self.c)
        #log_activity_coefficient = self.calculate_log_activity_coefficient (ionic_strength, self.c)
        #v_activity = self.c*(10**log_activity_coefficient)
        for i in range(0, len(self.A_Borkovec_rows)):
            print(self.A_Borkovec_rows[i] + ' : ' + str(self.c_Borkovec[i]) + '\n')
            #   print(self.list_species[i].name + ' : ' + str(v_activity[i]) + '\n')
        
        
        
        