# -*- coding: utf-8 -*-
"""
Created on 02/01/2019

@author: DaniJ
"""
#
#               THIS BENCHMARK IS BASED ON THE WORK OF Craig. Bethke "Geochemical and Biogeochemical Reaction Modeling"
#

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

# Bethe book
# Read database
database_file = 'Type7_Database_SC.txt'
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




# Reading input
infile = 'Type7_Input_SC_DL.txt'
# Instantiating input
list_aq_component, list_aq_value, names_pri_sorpt, list_sorption_comp  = getting_informationSC_input_file_v1 (infile)

#
CS = ChemSys_Surf()
CS.define_system_from_input_and_database ( DS, list_aq_component, list_aq_value, names_pri_sorpt, List_pri_sorpt_class = list_sorption_comp)


CS.create_pseudo_S()

#####
#### Following the example of Craig Bethke, section 10.4. The reactions having Hg have not been added because where not found in the llnl database
#### pseudo S shoudl have as columns: H2O, H+, Na+, Cl-, Pb+2, SO4-2, Hfo_sOH, Hfo_wOH, OH-, NaCl, Hfo_sOH2+, Hfo_sO-, Hfo_wOH2+, Hfo_wO-, Hfo_sOPb+, Hfo_wOPb+, Hfo_wSO4-, Hfo_wOHSO4-2
### And as rows, it should have the 2 reactions of the aqueous secondary species (OH-, NaCl) and 8 of the secondary surfaces species (Hfo_wOH2+, Hfo_wO-, Hfo_sOPb+, Hfo_wOPb+, Hfo_wSO4-, Hfo_wOHSO4-2)
### if nothing has been changed, it should be like:
pseudoS = np.array([[-1.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.], \
                    [ 0.,  0., -1., -1.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.], \
                    [ 0., -1.,  0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  1.,  0.,  0., 0.,  0.,  0.,  0.,  0.], \
                    [ 0.,  1.,  0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  1.,  0., 0.,  0.,  0.,  0.,  0.], \
                    [ 0., -1.,  0.,  0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  1., 0.,  0.,  0.,  0.,  0.], \
                    [ 0.,  1.,  0.,  0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  0., 1.,  0.,  0.,  0.,  0.], \
                    [ 0.,  1.,  0.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  1.,  0.,  0.,  0.], \
                    [ 0.,  1.,  0.,  0., -1.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  1.,  0.,  0.],  \
                    [ 1., -1.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  1.,  0.], \
                    [ 0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  1.]])

print(np.array_equal(CS.pseudoS, pseudoS))


CS.create_S ()

#### Following the example of Craig Bethke, section 10.4. The reactions having Hg have not been added because where not found in the llnl database
#### pseudo S shoudl have as columns: H2O, H+, Na+, Cl-, Pb+2, SO4-2, Hfo_sOH, Hfo_wOH, Hfo_sOH_Boltzf_psi_diffuselayer, OH-, NaCl, Hfo_sOH2+, Hfo_sO-, Hfo_wOH2+, Hfo_wO-, Hfo_sOPb+, Hfo_wOPb+, Hfo_wSO4-, Hfo_wOHSO4-2
#### Note:Hfo_sOH_Boltzf_psi_diffuselayer is the boltzman value of the psi_d of the diffuse layer
#### ALTHOUGH it is called Hfo_sOH_Boltzf_psi_diffuselayer it is actually the surface of Hfo_sOH and Hfo_wOH since they share the same surface and therefore the same electrostatic potential. The model just assume to sites weak and strong
### And as rows, it should have the 2 reactions of the aqueous secondary species (OH-, NaCl) and 8 of the secondary surfaces species (Hfo_wOH2+, Hfo_wO-, Hfo_sOPb+, Hfo_wOPb+, Hfo_wSO4-, Hfo_wOHSO4-2)
### if nothing has been changed, it should be like:
S = np.array([[-1.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.], \
              [ 0.,  0., -1., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.], \
              [ 0., -1.,  0.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0.,  1.,  0., 0.,  0.,  0.,  0.,  0.,  0.], \
              [ 0.,  1.,  0.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  1., 0.,  0.,  0.,  0.,  0.,  0.], \
              [ 0., -1.,  0.,  0.,  0.,  0.,  0., -1., -1.,  0.,  0.,  0.,  0., 1.,  0.,  0.,  0.,  0.,  0.], \
              [ 0.,  1.,  0.,  0.,  0.,  0.,  0., -1.,  1.,  0.,  0.,  0.,  0., 0.,  1.,  0.,  0.,  0.,  0.], \
              [ 0.,  1.,  0.,  0., -1.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0., 0.,  0.,  1.,  0.,  0.,  0.], \
              [ 0.,  1.,  0.,  0., -1.,  0.,  0., -1., -1.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  1.,  0.,  0.], \
              [ 1., -1.,  0.,  0.,  0., -1.,  0., -1.,  1.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  1.,  0.], \
              [ 0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  2.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  1.]])

print(np.array_equal(CS.S, S))


CS.create_U ()

#### Following the example of Craig Bethke, section 10.4. The reactions having Hg have not been added because where not found in the llnl database and therefore also the Hg component
#### U shoudl have as columns: H2O, H+, Na+, Cl-, Pb+2, SO4-2, Hfo_sOH, Hfo_wOH, Hfo_sOH_Boltzf_psi_diffuselayer, OH-, NaCl, Hfo_sOH2+, Hfo_sO-, Hfo_wOH2+, Hfo_wO-, Hfo_sOPb+, Hfo_wOPb+, Hfo_wSO4-, Hfo_wOHSO4-2
### And as tows, it should have the component made up of the following primary species =>  H2O, H+, Na+, Cl-, Pb+2, SO4-2, Hfo_sOH, Hfo_wOH, Hfo_sOH_Boltzf_psi_diffuselayer
U = np.array([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1., -0.,  0., -0., -0., -0., -0., -0., -1., -0.], \
              [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -0.,  1., -1., 1., -1., -1., -1.,  1., -0.], \
              [ 0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0., -0.,  1., -0., -0., -0., -0., -0., -0., -0., -0.], \
              [ 0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0., -0.,  1., -0., -0., -0., -0., -0., -0., -0., -0.], \
              [ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0., -0., -0., -0., -0., -0., -0.,  1.,  1., -0., -0.], \
              [ 0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0., -0., -0., -0., -0., -0., -0., -0., -0.,  1.,  1.], \
              [ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0., -0., -0.,  1.,  1., -0., -0.,  1., -0., -0., -0.], \
              [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0., -0., -0., -0., -0., 1.,  1., -0.,  1.,  1.,  1.], \
              [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -0., -0.,  1., -1., 1., -1.,  1.,  1., -1., -2.]])

print(np.array_equal(CS.U, U))

CS.create_log_k_vector()

# Case Craig Bethke, section 10.4.
# c_ = [H2O, H+, Na+, Cl-, Pb+2, SO4-2, Hfo_sOH, Hfo_wOH, Hfo_sOH_Boltzf_psi_diffuselayer, OH-, NaCl, Hfo_sOH2+, Hfo_sO-, Hfo_wOH2+, Hfo_wO-, Hfo_sOPb+, Hfo_wOPb+, Hfo_wSO4-, Hfo_wOHSO4-2]
#c_guess = np.array([])
#tolerance = 1e-6, max_n_iterations = 100
CS.Bethke_algorithm ()

#CS.speciation_Borkovec_1983_DLM(tolerance = 1e-8) --> give similar results als Phreeqc

#CS.print_speciation()
CS.print_speciation_Borkovec()

