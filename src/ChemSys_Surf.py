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
import scipy.integrate as integrate

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
                
            methods:
                set_S
                set_vector_aqueous_component_value
                aq_u_vector
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
        Stop=-np.matmul(S1.transpose(),np.linalg.inv(S2.transpose()))
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
        for i in range(0,self.length_sorpt_pri_sp):
            if isinstance(self.list_sorpt_pri_sp[i].names_Boltz_psi, str):
                self.names_elec_sorpt.append(self.list_sorpt_pri_sp[i].names_Boltz_psi)
            if isinstance(self.list_sorpt_pri_sp[i].names_Boltz_psi, list):
                for j in range(0, len(self.list_sorpt_pri_sp[i].names_Boltz_psi)):
                    self.names_elec_sorpt.append(self.list_sorpt_pri_sp[i].names_Boltz_psi[j])
        self.length_names_elec_sorpt = len(self.names_elec_sorpt)
        # Block
        if not hasattr(self, 'pseudoS_length_rows'):
            # self.pseudoS_length_rows = len(self.list_aq_reactions) + len(self.list_sorpt_reactions)
            self.pseudoS_length_rows = self.length_aq_sec_sp + self.length_sorpt_sec_sp
        S_electro = np.zeros((self.pseudoS_length_rows, self.length_names_elec_sorpt))
        col_position = 0
        for i in range(0, self.length_sorpt_pri_sp):
            sub_B = self.create_stoichiometric_surfacepotential (self.names_sorpt_pri_sp[i], self.list_sorpt_pri_sp[i].type_sorption)
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
            delta_X = np.linalg.solve(Z,-Y)
        
            # The error will be equal to the maximum increment
            err = max(abs(delta_X))
        
            # Relaxation factor borrow from Craig M.Bethke to avoid negative values
            max_1 = 1
            max_2 =np.amax(-2*np.multiply(delta_X, 1/c_n[0:pos_end_elec]))
            Max_f = np.amax([max_1, max_2])
            Del_mul = 1/Max_f
        
        
            # Update
            c_n[0:pos_end_elec] = c_n[0:pos_end_elec] + Del_mul*delta_X   # Update primary species
            log_c2 = np.matmul(np.linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(c_n[0:pos_end_elec])))      # Update secondary
            c_n[pos_end_elec:] =10**log_c2
            counter_iterations += 1
        if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.')
        self.c = c_n
        return c_n
     
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
            delta_X = np.linalg.solve(Z,-Y)
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
            log_c2 = np.matmul(np.linalg.inv(S2), self.log_k_vector - np.matmul(S1, np.log10(c_n[0:pos_end_elec])))      # Update secondary
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
        l_k_comp = np.matmul(np.linalg.inv(S2),self.log_k_vector)
        
        
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
        T = np.concatenate((T_chem, np.zeros(self.length_sorpt_pri_sp)))
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
            # Now that vector g (assuming symmetrical electrolyte is build I can build the B matrix and find Y)
            B = self.create_B_Borkovec(A, g_vec)
            # Calculating Y. The Y is given in equation 22 in Borkovec(1983)
            Y = np.matmul(B.transpose(), c) - T
            # Now the jacobian must be created
            Z = self.create_jacobian_Borkovec_1983_symm(A, B, c, X, I, z_vec ,g_vec)
            delta_X = np.linalg.solve(Z,-Y)
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
        dg_dXd_vec_o = self.dg_dXd_vec_eqn_11(z_vec, c_o[:self.length_aq_pri_sp + self.length_aq_sec_sp], X_o[-1])         # Necessary for equation 36 of Borkovec 1983
        # instantation variables loop
        counter_iterations = 0
        err = tolerance + 1  
        while err>tolerance and counter_iterations < max_iterations:
            # g must be calculated to create the matrix B
            g_vec = g_vec_o + dg_dXd_vec_o*(X[-1]-X_o[-1])
            # Now that vector g (assuming symmetrical electrolyte is build I can build the B matrix and find Y)
            
            
            
            
            
            
            
            counter_iterations += 1
        
            
        self.c_Borkovec = c
        return c
    
    def dg_dXd_eqn_11(self, z_vec, cb, Xd):
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
            dg_dXd.append(partQ*partB)
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
            partB = integrate.quad(self.integrand_fun_Borkovec_1983_eqn_11, 1, Xd, args = (zi,z_vec, cb), points = [0.9,1.1])
            if partB[1] > tol:
                raise ValueError('equation 11, integration of integrand high numerical error')
            g.append(partA*partB[0])
        return g
    
    def integrand_fun_Borkovec_1983_eqn_11 (self, x, zi,z_vec, cb):
        a = (x**zi)-1
        b= 0
        for i in range(0, len(z_vec)):
            b = b + cb[i]*((x**z_vec[i])-1)
        b = x*x*b
        return a/b
    
    
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
                    Z[j, k] = Z[j, k] + self.term_A4_Borkovec(n_iprime,I, z_vector, X[k], c)
                    
        return Z 
        
        
        
    def term_A2_and_A3_Borkovec(self, n_iprime, j, k, A, c, X, g_vector,z_vector, I):
        v = 0
        R = 0
        for iprime in range(0, n_iprime):
            v = v + ((z_vector[iprime]**2)/2)*A[iprime, k]*(c[iprime]/X[k])
        for iprime in range(0, n_iprime):
            R = R + c[iprime]*A[iprime, j]*(-g_vector[iprime]/(2*I))*v
        return R
    
    def term_A4_Borkovec(self,  n_iprime, I, z_vector, X_d, c):
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
                #
                T_sigma.append(T_0); T_sigma.append(T_b); T_sigma.append(Y+T_d)
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
                    # A = -F/(2*R*T)
                param = self.Faraday_constant/(2*(self.universal_gas_constant*self.temperature))
                A = -param
                #
                pos_C = self.length_aq_pri_sp+self.length_sorpt_pri_sp+self.length_names_elec_sorpt
                C_aq = np.concatenate((C[:self.length_aq_pri_sp], C[pos_C : (pos_C + self.length_aq_sec_sp)]))
                #
                I = self.calculate_ionic_strength(C_aq)
                B = np.sqrt(8*self.permittivity_free_space*self.dielectric_constant*self.universal_gas_constant*self.temperature*I)
                psi_d = self.Boltzman_factor_2_psi(C[pos_unknown_vector+2])
                par_C = param*psi_d
                C = np.cosh(par_C)
                F_d = A*B*C
                Z[pos_unknown_vector + 2,pos_unknown_vector + 2] = F_d + (self.list_sorpt_pri_sp[i].C2*F)*(D1/D4)
                
                pos_unknown_vector +=3
                
        return Z

        
    def calculate_ionic_strength (self,c):
        '''
            Calculate the ion strength: The vector C is supossed to be a vectorof concentrations that contains first the aqueous primary species followed by the aqueous secondary species. 
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
        
        ###################################################
        ################################
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
        
        
        
        