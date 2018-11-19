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
from creating_databaseObject_from_text_type1 import *
from getting_Information_from_DatabaseSC_file_v1 import *
from Read_input_file_type1 import *
# Read database
database = 'Type5_Database_JustSolution.txt'
# Instantiating database
Species_list, Reaction_list, primary_species_list, secondary_species_list = creating_databaseObject_from_text_type1 (database)

D = Database()
D.set_species_list(Species_list)
D.set_reaction_list(Reaction_list)
D.set_primary_species(primary_species_list)
D.set_secondary_species(secondary_species_list)

# Check S
D.create_log_k_vector()
D.create_charge_vector()
D.create_S()

#print(D.log_k_vector == [14.0, 10.329, 2.98, 3.224, 11.399, 11.435, 16.681, -12.78, -11.44])


#
#
# Reading input
infile = 'Type5_Input_JustSolution.txt'
# Instantiating input
list_prim, list_val = Read_input_file_type1 (infile)
C = ChemSys( list_prim, list_val, D)
C.calculate_U_f1()



c = C.speciation_noactivitycoefficient(tolerance = 1e-12)

c2 = C.speciation_noactivitycoefficient_Westall1980(tolerance = 1e-12)

#C.print_speciation()

print(c)
print(c2)


#######################################   FROM HERE STARTS THE SORPTION PART ##########################################################################################

# Read database
database_file = 'Type5_Database_SC.txt'
# Instantiating database
n_aq_sp_pri, n_aq_sp_sec, n_sorpt_sp_pri, n_sorpt_sp_sec, Aq_sp_list_pri, Aq_sp_list_sec, Sorp_sp_list_pri, Sorp_sp_list_sec, Aq_list_react, Sorp_list_react = getting_Information_from_DatabaseSC_file_v1 (database_file)
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
PseudoS = np.array([[1, 0, 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0], [0,  0, -1,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0], [-3,  0,  0, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0], [-1,  0,  0, -1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0], \
                    [-2,  0,  0, -1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0], [-1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0], [1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0], \
                    [-2,  0,  0, -1, -1,  1,  0,  0,  0, 0,  0,  0,  1,  0,  0], [-1,  0,  0, -1, -1,  1,  0,  0,  0,  0,  0, 0,  0,  1,  0], [0, 0,  0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])

print(np.array_equal(DS.pseudoS, PseudoS))