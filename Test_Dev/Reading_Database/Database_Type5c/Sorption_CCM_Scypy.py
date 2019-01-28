# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:06:19 2018

@author: DaniJ
"""

import numpy as np
from Species import Aq_Species
from Reaction import Reaction
from Database import Database
from Database_SC import Database_SC
from ChemSys import ChemSys
from ChemSys_Surf import ChemSys_Surf
from creating_databaseObject_from_text_type1 import *
from getting_Information_from_databaseSC_file_v1 import *
from getting_informationSC_input_file_v1 import *
from Read_input_file_type1 import *

######################################################################################################################################################################
#Applying alorithm

# Read database
database_file = 'Type5_Database_SC.txt'
# Instantiating database
n_aq_sp_pri, n_aq_sp_sec, n_sorpt_sp_pri, n_sorpt_sp_sec, Aq_sp_list_pri, Aq_sp_list_sec, Sorp_sp_list_pri, Sorp_sp_list_sec, Aq_list_react, Sorp_list_react = getting_Information_from_databaseSC_file_v1 (database_file)
#
DS = Database_SC()
DS.set_names_aq_primary_species ( n_aq_sp_pri)
DS.set_names_aq_secondary_species ( n_aq_sp_sec)
DS.set_names_sorpt_primary_species (n_sorpt_sp_pri)
DS.set_names_sorpt_secondary_species (n_sorpt_sp_sec)
DS.set_aq_list_pri_class (Aq_sp_list_pri)
DS.set_aq_list_sec_class (Aq_sp_list_sec)
DS.set_sorpt_list_pri_class (Sorp_sp_list_pri)
DS.set_sorpt_list_sec_class (Sorp_sp_list_sec)
DS.set_aq_reactions_list (Aq_list_react)
DS.set_sorpt_reactions_list (Sorp_list_react)

DS.create_pseudo_S()

#print(DS.pseudoS)
# In this case, unless change on the database, this matrix should be the Pseudo S matrix of the database clase.
DSPseudoS = np.array([[1, 0, 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0], [0,  0, -1,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0], [-3,  0,  0, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0], [-1,  0,  0, -1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0], \
                    [-2,  0,  0, -1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0], [-1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0], [1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0], \
                    [-2,  0,  0, -1, -1,  1,  0,  0,  0, 0,  0,  0,  1,  0,  0], [-1,  0,  0, -1, -1,  1,  0,  0,  0,  0,  0, 0,  0,  1,  0], [0, 0,  0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])

print(np.array_equal(DS.pseudoS, DSPseudoS))



# Reading input
infile = 'Type5_Input_SC_CCM.txt'
# Instantiating input
list_aq_component, list_aq_value, names_pri_sorpt, list_sorption_comp  = getting_informationSC_input_file_v1 (infile)

#
CS = ChemSys_Surf()
CS.define_system_from_input_and_database ( DS, list_aq_component, list_aq_value, names_pri_sorpt, List_pri_sorpt_class = list_sorption_comp)
CS.create_pseudo_S()
#print(CS.pseudoS) 


# In this case, unless change on the database, this matrix should be the Pseudo S matrix of the database clase.
CSPseudoS = np.array([[1, 0, 0,  1,  0,  0,  0,  0,  0,  0,  0,  0], [-3, -1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0], [-1,  -1,  0,  0,  0,   1,  0,  0,  0,  0,  0,  0], \
                    [-2, -1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0], [-1, 0, -1,  0,  0,  0,  0,   1,  0,  0,  0,  0], [1,  0, -1,  0,  0,  0,  0,  0,  1,  0,  0,  0], \
                    [-2, -1, -1,  1,  0,  0,  0, 0,  0,   1,  0,  0], [-1, -1, -1,  1,  0,  0,  0,  0,  0, 0,   1,  0], [0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0,  1]])

print(np.array_equal(CS.pseudoS, CSPseudoS))


CS.create_S ()
#print(CS.S) 

CS_S = np.array([[1, 0, 0, 0, 1,  0,  0,  0,  0,  0,  0,  0,  0], [-3, -1,  0, 0, 0,  1,  0,  0,  0,  0,  0,  0,  0], [-1,  -1,  0,  0, 0,  0,   1,  0,  0,  0,  0,  0,  0], \
                    [-2, -1,  0, 0, 0,  0,  0,  1,  0,  0,  0,  0,  0], [-1, 0, -1, -1, 0,  0,  0,  0,   1,  0,  0,  0,  0], [1,  0, -1, 1, 0,  0,  0,  0,  0,  1,  0,  0,  0], \
                    [-2, -1, -1,  0, 1,  0,  0,  0, 0,  0,   1,  0,  0], [-1, -1, -1, 1,  1,  0,  0,  0,  0,  0, 0,   1,  0], [0, -1, -1, 3,  0, 0, 0, 0, 0, 0, 0, 0,  1]])

print(np.array_equal(CS.S, CS_S))


CS.create_U ()
print(CS.U)
CS.create_log_k_vector()

        #if type_I == 5:
            # Case example 5
            # c_ = [H+, AsO4-3, SurfOH, Boltzman_SurfOH, OH-, H3AsO4, HAsO4-2, H2AsO4-, SurfOH2+, SurfO-, SurfH2AsO4, SurfHAsO4-, SurfOHAsO4-3]
#c_guess = np.array([7.864e-09, 2.626e-28, 5.373e-02, 9.975e-01, 1.315e-06, 5.634e-32, 6.126e-25, 4.211e-26, 8.133e-03, 8.133e-03, 1.219e-10, 2.488e-08, 4.910e-19])
c_guess = np.array([7.864e-09, 2.626e-28, 5.373e-02, 9.975e-01, 1.315e-06, 5.634e-32, 6.126e-25, 4.211e-26, 8.133e-03, 8.133e-03, 1.219e-10, 2.488e-08, 4.910e-19])
x_guess =  np.array([7.864e-09, 2.626e-28, 5.373e-02, 9.975e-01])
CS.speciation_Westall1980_CCM (c_guess = c_guess)
CS.print_speciation()
#CS.speciation_Westall1980_CCM_v2 (x= x_guess)
#CS.print_speciation()

