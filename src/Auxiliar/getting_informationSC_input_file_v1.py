# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:34:45 2018

@author: DaniJ
"""

import compilation_functions_reading_file as cfrf


def getting_informationSC_input_file_v1 (text):
    #text = 'Type1.txt' Reading text
    f = open(text, 'r')
    lines = f.readlines()
    f.close()
    
    # Defining output variables # The creation of these below variables has not much sense, since it can be created directly from the function call, but it helps in keeping a coherence.
    list_aq_component = []          # A list containing the name of the (primary species) components in the system such as ['A'  'CO3' ...]
    list_aq_value = []              # The value of the components defined in the list above [5  6e-6 ...] 
    list_sorption_comp = []         # List of primary surface species
    names_pri_sorpt = []            # list_names such as ['SOH'  'S_wOH' ... etc]
    # First the startindex and the endindex of Solution and Surface is found
    [s_index_sol, end_index_sol] = cfrf.find_Indices_Solution(lines)
    [s_index_sorpt, end_index_sorpt] = cfrf.find_Indices_Sorption(lines)
    
    # filling 
    list_aq_component, list_aq_value = cfrf.read_block_Solution (lines[s_index_sol:end_index_sol])   
    names_pri_sorpt, list_sorption_comp = cfrf.read_block_Sorption (lines[s_index_sorpt:end_index_sorpt])
    return list_aq_component, list_aq_value, names_pri_sorpt,list_sorption_comp