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
        npri = self.length_aq_pri_sp +self.length_sorpt_pri_sp + len(self.names_elec_sorpt)
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
        for i in range(0, len(self.names_elec_sorpt)):
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
        if type_sorpt == 'CCM':
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
        # instantiation variables for loop
        counter_iterations = 0
        err = tolerance + 1
        sorpt_u_vector = self.create_sorpt_vec()
        T_chem = np.concatenate ((self.aq_u_vector, sorpt_u_vector))
        while err>tolerance and counter_iterations < max_iterations:
            # Calculate U vector [If I am not wrong T_sigma must be calculated at every step, since it depends somehow in the surface potential, and it is unknown]
            u_electro = self.calculate_u_electro(c_n[pos_start_elec:pos_end_elec])
            T = np.concatenate ((T_chem, u_electro))
            # Calculate f or better said in this specific case Y
            Y = self.U.dot(c_n) - T
            # Calculate Z
            Z = self.Jacobian_Speciation_CCM_Westall1980(c_n, pos_start_elec, pos_end_elec)
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
        if counter_iterations >= max_iterations:
            raise ValueError('Max number of iterations surpassed.')
        self.c = c_n
        return c_n
                
                
    def calculate_u_electro (self, unknonw_boltzman_vect):
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
            
        return np.array(T_sigma)
        
        
    def Boltzman_factor_2_psi(self, x):
        D = self.universal_gas_constant*self.temperature
        psi = - np.log(x)*(D/self.Faraday_constant)
        return psi

    def  create_sorpt_vec(self):            
        T_sorpt = []
        for i in range(0, self.length_sorpt_pri_sp):
            T_sorpt.append(self.list_sorpt_pri_sp[i].T_solid)
        return T_sorpt
        
    def Jacobian_Speciation_CCM_Westall1980(self, C, n_aq_plus_n_sorpt, n_primaryspecies):
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
            if self.list_sorpt_pri_sp[i].type_sorption == 'CCM':
                p = i + n_aq_plus_n_sorpt
                D1 = self.universal_gas_constant*self.temperature
                D2 = self.Faraday_constant*C[i]
                F = ((self.list_sorpt_pri_sp[i].sp_surf_area*self.list_sorpt_pri_sp[i].solid_concentration_or_grams)/self.Faraday_constant)
                Z[p,p] = Z[p,p] + (self.list_sorpt_pri_sp[i].C1*F)*(D1/D2)
        return Z

        
        
        
        ###################################################
        ################################
        ########################### Postprocessing
        ###################################################
        
    def print_speciation(self):
        #ionic_strength = self.calculate_ionic_strength (self.c)
        #log_activity_coefficient = self.calculate_log_activity_coefficient (ionic_strength, self.c)
        #v_activity = self.c*(10**log_activity_coefficient)
        for i in range(0, self.S_length_columns):
            print(self.S_names_columns[i] + ' : ' + str(self.c[i]) + '\n')
            #   print(self.list_species[i].name + ' : ' + str(v_activity[i]) + '\n')
        
        
        
        