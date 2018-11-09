# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 12:06:04 2018

@author: DaniJ
"""
from Species import *


# Searching functions

def find_Indices_Aqueous_PrimaryBlock(list_text):
    string_startindex = 'Primary_Species\n'
    list_p_end = ['Secondary_Species\n', 'Surface_Primary\n', 'Surface_Secondary\n']
    s_index = list_text.index(string_startindex)
    e_index = find_last_index(list_text, s_index, list_p_end)
    return s_index, e_index

def find_Indices_Aqueous_SecondaryBlock(list_text):
    string_startindex='Secondary_Species\n'
    list_p_end = ['Primary_Species\n', 'Surface_Primary\n', 'Surface_Secondary\n']
    s_index = list_text.index(string_startindex)
    e_index = find_last_index(list_text, s_index, list_p_end)
    return s_index, e_index

def find_Indices_Sorption_PrimaryBlock(list_text):
    string_startindex='Surface_Primary\n'
    list_p_end = ['Primary_Species\n', 'Secondary_Species\n', 'Surface_Secondary\n']
    s_index = list_text.index(string_startindex)
    e_index = find_last_index(list_text, s_index, list_p_end)
    return s_index, e_index

def find_Indices_Sorption_SecondaryBlock(list_text):
    string_startindex='Surface_Secondary\n'
    list_p_end = ['Primary_Species\n', 'Secondary_Species\n', 'Surface_Primary\n']
    s_index = list_text.index(string_startindex)
    e_index = find_last_index(list_text, s_index, list_p_end)
    return s_index, e_index

def find_last_index(list_text, min_point, possible_ends_list):
    ind_last = len(list_text)
    for i in possible_ends_list:
        temp_index = list_text.index(i)
        if temp_index > min_point and temp_index < ind_last:
            ind_last = temp_index
    return ind_last

# Functions for reading blocks (i.e. Primary Species, Secondary Species, etc) and outputing parameters

def read_block_PrimarySpecies (list_datablock_primary):
    # initialization of variables to return
    names_primary_species = []
    List_Aq_classes_primary_species = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_datablock_primary):
        temp =list_datablock_primary[line_counter].strip() 
        if temp == '' or temp[0] == '#' or temp == 'Primary_Species':
            line_counter += 1
        else:
            words_line = temp.split()
            # appending primary Species
            names_primary_species.append(words_line[0])
            # Creating an aqueous species
            S = Aq_Species(words_line[0])
            # Taking into account the exception such as water or charge
            if words_line[0] == 'e-':
                S.it_is_charge(True)
            elif words_line[0] == 'H2O':
                S.Set_f_activity_coefficient ('water')
            S.set_charge (int(words_line[1]))
            S.set_gfw (float(words_line[2]))    
            if len(words_line) > 3:
                if words_line[3][0] != '#':
                    S.set_ionsizeparameter(float(words_line[3]))
            if len(words_line) > 4:        
                if words_line[4][0] != '#':
                    S.set_deviationionsizeparameter(float(words_line[4]))
            List_Aq_classes_primary_species.append(S)
            line_counter += 1
            
    return names_primary_species, List_Aq_classes_primary_species
    

