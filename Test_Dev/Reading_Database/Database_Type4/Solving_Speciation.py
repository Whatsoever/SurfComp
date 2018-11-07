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
database = 'Type4_Database.txt'
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
infile = 'Type4_Input.txt'
# Instantiating input
list_prim, list_val = Read_input_file_type1 (infile)
C = ChemSys( list_prim, list_val, D)
C.calculate_U_f1()
# Calculating Speciation
#C.calculate_speciation (1)
#print(np.matmul(C.U, C.S.transpose()))
# Pre-processing Results

c = C.speciation_algorithm2_reduced_problem()
C.print_speciation()
C.calculate_ionic_strength(c)