# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 10:25:11 2018

@author: DaniJ
"""
import numpy as np
from Species import Aq_Species
from Reaction import Reaction
from Database import Database
from ChemSys import ChemSys
from creating_databaseObject_from_text_type1 import *
from Read_input_file_type1 import *
# Read database
database = 'Type3_Database.txt'
# Instantiating database
names_aqueous_primary_species, names_aqueous_secondary_species, Aq_Species_list_pri, Aq_Species_list_sec, Aq_Reaction_list = creating_databaseObject_from_text_type1 (database)

D = Database()
D.set_names_aq_primary_species(names_aqueous_primary_species)
D.set_aq_list_pri_class(Aq_Species_list_pri)
D.set_names_aq_secondary_species (names_aqueous_secondary_species)
D.set_aq_list_sec_class( Aq_Species_list_sec)
D.set_aq_reactions_list(Aq_Reaction_list)

# Check S
D.create_log_k_vector()
D.create_charge_vector()
D.create_S()

#print(D.log_k_vector == [14.0, 10.329, 2.98, 3.224, 11.399, 11.435, 16.681, -12.78, -11.44])


#
#
# Reading input
infile = 'Type3_Input.txt'
# Instantiating input
list_prim, list_val = Read_input_file_type1 (infile)
C = ChemSys( list_prim, list_val, D)
C.calculate_U_f1()
# Calculating Speciation
#C.calculate_speciation (1)
#print(np.matmul(C.U, C.S.transpose()))
# Pre-processing Results
#c = C.speciation_noactivitycoefficient( )
#c = C.speciation_noactivitycoefficient( tolerance = 1e-7, max_n_iterations = 150)
c = C.speciation_algorithm1( tolerance = 1e-12)
C.print_speciation()

C.print_ionic_strength()
C.print_activity()