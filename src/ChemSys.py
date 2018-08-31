# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:08:12 2018

@author: Lützenkirchen, Heberling, Jara
"""
from Database import Database
import numpy as np

class ChemSys (Database):
    # Constructor
    def __init__(self, list_prim, list_val, DatabaseC):
        '''
            list_prim = ['Cl-', 'H2O', 'H+']        It is supossed to be the primary species associated to a component.
            
            list_val  = [5 4 20]                    It is supposed to be the concentration values of the components associated to primary species that can be found in list_prim, namely Cl- = 5. So far units are mol/L          
        
            DatabaseC                               It is supposed to be the database that is going to be used in order to built the stochiometric matrix of these system, and the pertinant secondary species, log_k values, etc
            
            precondition: The DatabaseC class should have run the create_S or check_consistency_species, previously. Otherwise, The method will fail.
        '''
        DatabaseC.check_consistency_species()       # This method is run to assure that what is written below this line has sense.
        # Getting index of list_prim in the database
        indices_ps = self.index_InputPrimaryspeciesinDatabase ( list_prim, DatabaseC.name_primary_species)
        # Sorting list_prim and list_val according to the indices, AFTER THIS STEP THE LIST BECOMES TUPLE !!!!!!!!!!!
        indices_ps, list_prim, list_val = zip(*sorted(zip(indices_ps, list_prim, list_val)))
        
        self.u_comp_vec = np.array(list_val)
        # name_prymary_species is like list_prim but order following the order of the Database
        self.name_primary_species = list(list_prim)                        
        # It is the list of species = [Species1, Species2, Species3] related to name_primary_species
        self.list_species = [DatabaseC.list_species[i] for i in indices_ps]
        # Check which Reactions in the database occur in the ChemSys --> from there the secondary species will be obtained.
        indices_r_inDatabase, self.name_secondary_species = self.index_reactionsinDatabase ( DatabaseC.list_reactions, DatabaseC.n_reactions)
        # Add secondary species in the list of species
        self.add_secondary_species (DatabaseC)
        # Get the list of reactions from the Database
        self.list_reactions = [DatabaseC.list_reactions[i] for i in indices_r_inDatabase]
        self.n_reactions = len(self.list_reactions)
        self.n_species = len(self.list_species)
        # For calculations --> This functions are linked to the parents
        self.create_S()
        self.create_log_k_vector()
        self.create_charge_vector()
        
    # Setters
    
    # Initializations functions
    # Searching
    # Returns the indice of the chemical species in the database
    def index_InputPrimaryspeciesinDatabase (self, list1, list2):
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
    
    # Returns the index of the reactions that can occur according to the database
    def index_reactionsinDatabase (self, l_dic_react, len_l_dic_react):
        '''
            The function returns two list. One with the indices of the reactions that occur in the ChemSys according to the inputed Database, and the other the secondary species in each reaction. 
            Both, list are in agremment. e.g. l_ind_reaction = [0, 4, 6, 9], l_secondary_species = ['A' 'B' 'C' 'F']  From reaction 0 of the database the secondary species obtained is A, from 6 is C, and so on.
        '''
        list_indices_r=[]
        list_secondary_species_r=[]
        for i in range(0, len_l_dic_react):
            temp_d_r = l_dic_react[i].reaction
            list_temp_d_r = list(temp_d_r.keys())
            n_s = 0
            for d in list_temp_d_r:
                c = self.name_primary_species.count(d)
                if c != 1 and c!= 0:
                    raise ValueError('[ChemSys class, method Index_ReactionsinDatabase] It seems that the name_primary_species property is wrong.')
                elif c == 0:
                    n_s += 1
                    n_s_name = d
            if n_s == 1:
                list_indices_r.append(i)
                list_secondary_species_r.append(n_s_name)       
        return list_indices_r, list_secondary_species_r
    
    # Adds secondary species
    def add_secondary_species (self, DatabaseC):
        '''
            Adds the secondary species in the list of species
        '''
        l = DatabaseC.name_primary_species + DatabaseC.name_secondary_species
        for i in self.name_secondary_species:
            pos = l.index(i)
            self.list_species.append(DatabaseC.list_species[pos])
            

    ###################### Calculations #######################################
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
        S1 = self.S[:, 0:len(self.name_primary_species)].copy() 
        S2 = self.S[:, len(self.name_primary_species):].copy()
        return S1, S2
    
    # Component matrix functions
    
    # In This function the component matrix U is calculated as steated by Steefel and MacQuarrie, 1996 or as it is the same stated by Saaltink et al., 1998  (equation 24) 
    # It is supposed to be obtained by means of a Gauss-Jordan elimination
    def calculate_U_f1 (self):
        '''
            Calculates U from the Stoichiometric matrix by using "Gauss-Jordan elimination". Saaltink et al., 1998  (equation 24) 
            Actually the Gaus elimination is applied to the transpose of St
        '''
        S1, S2 = self.separte_S_into_S1_and_S2()
        # Create Identity
        I = np.identity (self.n_species - self.n_reactions)
        # Use  Saaltink et al., 1998  (equation 24) 
        St = np.matmul(-S1.transpose(), np.linalg.inv(S2.transpose()))  # -S1t*inv(S2t)
        self.U = np.concatenate((I, St), 1)
        self.check_U_consistancy()
    
    def check_U_consistancy(self):
        '''
        multiplies U*S_transpose and checks that all the values are zero, which implies that U is the kernel of the S_transpose, necessary for the use of components
        '''
        R = np.matmul(self.U, self.S.transpose())
        b = not np.any(R)
        if b == False:
            raise ValueError('[ChemSys Class/ check_U_consistancy] Apparently the U matrix is not the Kernel of the transpose of the stoichiometric matrix.')
    
    # Calculations
    # Speciation calculations
    def caculate_speciation(self, algorithm_num):
        if algorithm == 1:
            self.speciation_algorithm1()
        else:
            raise ValueError('Not algorithm for speciation with this number.')
    
    def speciation_algorithm1(self, tolerance = 1e-8, max_n_iterations = 100):
        # Tolerance
        #tolerance = 1e-4
        # Initial guess c1 (c2 must also be inizialitated)
        c = self.Instatiation_step(1)
        nps = length(self.name_primary_species)
        c1_n = c[:nps].copy()
        c2_n0 = c[nps:].copy()
        delta_c = 0;
        # start loop
        b = False        # True will mean that no solution has been found
        counter_iterations = 0;
        #max_n_iteration = 100;
        
        while not b and counter_iterations < max_n_iterations:
            #update cn
            cn = cn + delta_c
            # Compute c2, in order to do such think, compute I, activity, and then c2.
            ionic_strength = self.calculate_ionic_strength (c1_n, c2_n0)
            activity_coefficient = self.calculate_activity (ionic_strength)
            c2_n1 = self.speciation_secondaryspecies (c1_n, activity_coefficient)
            #compute f
            f = self.calculate_NR_function (c1_n, c2_n1)
            # compute df/dc1
            df_dc1 = self.Jacobian_NR_function(c1_n,c2_n1)
            # delta_c = c_n+1 - c_n
            delta_c = -self.division_vector(f, df_dc1)
            # Converge?
            b = delta_c < tolerance
            if b == True:
                break
            else:
                counter_iterations += 1
                
        
        
        
    def c_ini = Instantiation_step (self, type_I=1):
            if type_I == 0:
                c_ini = np.ones(self.n_species)*1e-10
                c_ini = c_ini.transpose()
            elif type_I == 1:
                c_ini = self.NewtonRaphson_noactivitycoefficient()
            else:
                raise ValueError('Not algorithm for instantiationwith these number.')
        
    def NewtonRapshon_noactivitycoefficient(self, tolerance = 1e-10, max_n_iterations = 100):
        c_guess = Instantiation_step (self, type_I=0)
        c_n = c_guess
        
        b = False        # True will mean that no solution has been found
        counter_iterations = 0;
        #max_n_iteration = 100;
        
        while not b and counter_iterations < max_n_iterations:
            # Calculate F; F = [U*c-u; S*log(c)-logK];
            Upper_Part = self.U.dot(c_n) - self.u_comp_vec
            Lower_Part = self.S.dot(np.log10(c_n)) - np.array(self.log_k_vector)'
    