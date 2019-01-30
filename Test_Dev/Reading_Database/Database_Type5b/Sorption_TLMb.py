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

# Read database
database_file = 'Type5_Database_SCb.txt'
# Instantiating database
n_aq_sp_pri, n_aq_sp_sec, n_sorpt_sp_pri, n_sorpt_sp_sec, Aq_sp_list_pri, Aq_sp_list_sec, Sorp_sp_list_pri, Sorp_sp_list_sec, Aq_list_react, Sorp_list_react = getting_Information_from_databaseSC_file_v1 (database_file)
#
DS2 = Database_SC()
DS2.set_names_aq_primary_species ( n_aq_sp_pri)
DS2.set_names_aq_secondary_species ( n_aq_sp_sec)
DS2.set_names_sorpt_primary_species (n_sorpt_sp_pri)
DS2.set_names_sorpt_secondary_species (n_sorpt_sp_sec)
DS2.set_aq_list_pri_class (Aq_sp_list_pri)
DS2.set_aq_list_sec_class (Aq_sp_list_sec)
DS2.set_sorpt_list_pri_class (Sorp_sp_list_pri)
DS2.set_sorpt_list_sec_class (Sorp_sp_list_sec)
DS2.set_aq_reactions_list (Aq_list_react)
DS2.set_sorpt_reactions_list (Sorp_list_react)
# Reading input
infile = 'Type5_Input_SC_TLMb.txt'
# Instantiating input
list_aq_component, list_aq_value, names_pri_sorpt, list_sorption_comp  = getting_informationSC_input_file_v1 (infile)
CS1 = ChemSys_Surf()
CS1.define_system_from_input_and_database ( DS2, list_aq_component, list_aq_value, names_pri_sorpt, List_pri_sorpt_class = list_sorption_comp)
CS1.create_pseudo_S()
CS1.create_S ()
CS1_S = np.array([[ 1.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.], [-3., -1.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,0.,  0.], [-1., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0., 0.,  0.], \
                    [-2., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0., 0.,  0.], [-1.,  0., -1., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0., 0.,  0.], [1.,  0., -1.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0., 0.,  0], \
                    [-2., -1., -1., -3.,  3.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  1., 0.,  0.], [-1., -1., -1., -2.,  3.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0., 1.,  0.], [0., -1., -1.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,0.,  1]])

print(np.array_equal(CS1.S, CS1_S))


CS1.create_U ()
CS1.create_log_k_vector()
# c_ = ['H+', 'Na+', 'Cl-', 'SOH', 'SOH_Boltzf_psi_0', 'SOH_Boltzf_psi_beta', 'SOH_Boltzf_psi_diffuselayer', 'OH-', 'SOH2+', 'SO-', 'SOH2Cl', 'SONa']
c_guess = np.array([1.602e-08, 7.412e-24, 2.503e-05, 5.371e-02, 4.872e-01, 4.902e-01 , 6.325e-01, 6.366e-07, 8.143e-03, 8.117e-03, 8.945e-35, 2.500e-05])
# X = ['H+', 'Na+', 'Cl-', 'SOH', 'SOH_Boltzf_psi_0', 'SOH_Boltzf_psi_beta', 'SOH_Boltzf_psi_diffuselayer']
X_g = np.array([1.602e-08, 7.412e-24, 2.503e-05, 5.371e-02, 4.872e-01, 4.902e-01 , 6.325e-01])
#CS1.set_constant_ionic_strength (7.293e-06)
#CS1.speciation_Westall1980_TLM (tolerance = 1e-8, max_iterations = 100, c_guess = c_guess)
#CS1.speciation_Westall1980_TLMb (tolerance = 1e-8, max_n_iterations = 100, X_guess = X_g)
#CS1.speciation_Westall1980_CCM_v2 (x= X_g)
Ln_X = np.log(X_g)
#CS1.speciation_Westall1980_v3 ( tolerance = 1e-8, max_iterations = 100, Ln_x = Ln_X)
#CS1.speciation_Westall1980_v3 ( tolerance = 1e-6, max_iterations = 100, Ln_x = Ln_X, activity_b = True)
CS1.print_speciation()



