# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 12:06:04 2018

@author: DaniJ
"""
from Species import *
from Reaction import *

# Searching functions

# for database type
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
        try:
            temp_index = list_text.index(i)
        except ValueError:
            pass
        if temp_index > min_point and temp_index < ind_last:
            ind_last = temp_index
    return ind_last


# for Input type

def find_Indices_Solution (list_text):
    string_startindex='Solution\n'
    list_p_end = ['Sorption\n']
    s_index = list_text.index(string_startindex)
    e_index = find_last_index(list_text, s_index, list_p_end)
    return s_index, e_index

def find_Indices_Sorption (list_text):
    string_startindex='Sorption\n'
    list_p_end = ['Solution\n']
    s_index = list_text.index(string_startindex)
    e_index = find_last_index(list_text, s_index, list_p_end)
    return s_index, e_index




# Functions for reading blocks (i.e. Primary Species, Secondary Species, etc) and outputing parameters
    
# for database type

def read_block_Primary_Species (list_datablock_primary):
    # initialization of variables to return
    names_primary_species = []
    List_Aq_classes_primary_species = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_datablock_primary):
        temp = list_datablock_primary[line_counter].strip() 
        if temp == '' or temp[0] == '#' or temp == 'Primary_Species':
            line_counter += 1
        else:
            words_line = temp.split()
            words_line = remove_coments(words_line)             # Remove the parts that contain comments "#"  
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

def read_block_Secondary_Species (list_datablock_secondary): 
    # initialization of variables to return
    names_secondary_species = []
    List_Aq_classes_secondary_species = []
    List_aq_reactions = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_datablock_secondary):
        temp = list_datablock_secondary[line_counter].strip() 
        if temp == '' or temp[0] == '#' or temp == 'Secondary_Species':
            line_counter += 1
        else:
            words_line = temp.split()
            words_line = remove_coments(words_line)             # Remove the parts that contain comments "#"                                         
            if words_line[0] == '-log_k':
                R.set_log_k(float(words_line[1]))
                List_aq_reactions.append(R)
                line_counter +=1
            elif words_line[0] == '-a':
                S.set_ionsizeparameter(float(words_line[1]))
                line_counter += 1
            elif words_line[0] == '-b':
                S.set_deviationionsizeparameter(float(words_line[1]))
                line_counter += 1
            else:
                R = Reaction()
                dic_reaction = {}
                for i in range(0, len(words_line), 2):
                    if i == 0:
                        dic_reaction[words_line[i]] = int(1)
                    else:
                        dic_reaction[words_line[i]] = int(words_line[i+1])
                R.set_reaction(dic_reaction)
                S = Aq_Species(words_line[0])
                names_secondary_species.append(words_line[0])
                S.set_charge (int(words_line[1]))
                List_Aq_classes_secondary_species.append(S)
                line_counter += 1
    return names_secondary_species,List_Aq_classes_secondary_species, List_aq_reactions 


def read_block_Surface_Primary (list_datablock_surfpri): 
    # initialization of variables to return
    names_primary_sorptspecies = []
    List_Sorpt_classes_primary_species = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_datablock_surfpri):
       temp = list_datablock_surfpri[line_counter].strip() 
       if temp == '' or temp[0] == '#' or temp == 'Surface_Primary':
           line_counter += 1 
       else:
           words_line = temp.split()
           words_line = remove_coments(words_line)             # Remove the parts that contain comments "#"  
           S = Surf_species(words_line[0])
           names_primary_sorptspecies.append(words_line[0])
           List_Sorpt_classes_primary_species.append(S)
           line_counter += 1 
        
    return names_primary_sorptspecies, List_Sorpt_classes_primary_species


def read_block_Surface_Secondary (list_datablock_surfsec):    
    # initialization of variables to return
    names_sorption_secondary_species = []
    List_Sorption_secondary_species = []
    List_sorption_reactions = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_datablock_surfsec):
        temp = list_datablock_surfsec[line_counter].strip() 
        if temp == '' or temp[0] == '#' or temp == 'Surface_Secondary':
            line_counter += 1
        else:
            words_line = temp.split()
            words_line = remove_coments(words_line)             # Remove the parts that contain comments "#"  
            if words_line[0] == '-log_k':
                R.set_log_k(float(words_line[1]))
                List_sorption_reactions.append(R)
                line_counter +=1
            else:
                R = Reaction()
                dic_reaction = {}
                for i in range(0, len(words_line), 2):
                    if i == 0:
                        dic_reaction[words_line[i]] = int(1)
                    else:
                        dic_reaction[words_line[i]] = int(words_line[i+1])
                R.set_reaction(dic_reaction)
                S = Surf_species(words_line[0])
                S.set_surface_charge (int(words_line[1]))
                names_sorption_secondary_species.append(words_line[0])
                List_Sorption_secondary_species.append(S)
                line_counter += 1    
      
    return names_sorption_secondary_species, List_Sorption_secondary_species, List_sorption_reactions


# for input type
def read_block_Solution (list_solution):
    # initialization of variables to return
    names_aq_component = []
    values_aq_comp = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_solution):
        temp = list_solution[line_counter].strip() 
        if temp == '' or temp[0] == '#' or temp == 'Solution':
            line_counter += 1
        else:
            words_line = temp.split()   
            words_line = remove_coments(words_line)             # Remove the parts that contain comments "#"  
            names_aq_component.append(words_line[0])
            values_aq_comp.append(float(words_line[1]))
            line_counter += 1
    return names_aq_component, values_aq_comp

def read_block_Sorption (list_sorption):
    # initialization of variables to return
    list_sorption_class = []
    names_pri_sorpt = []
    # variables for loop
    line_counter = 0
    while line_counter < len(list_sorption):
        temp = list_sorption[line_counter].strip() 
        if temp == '' or temp[0] == '#' or temp == 'Sorption':
            line_counter += 1
        else:
            words_line = temp.split()
            words_line = remove_coments(words_line)             # Remove the parts that contain comments "#"  
            if words_line[0] == '-type':
                if words_line[1] == 'related':
                    Sorpt_pri_sp.set_type_sorption(words_line[1])
                    Sorpt_pri_sp.set_type_relation(words_line[2])
                else:
                    Sorpt_pri_sp.set_type_sorption(words_line[1])
                line_counter += 1
            elif words_line[0] == '-C1':
                Sorpt_pri_sp.set_capacitance_1(words_line[1])
                line_counter += 1
            elif words_line[0] == '-C2':
                Sorpt_pri_sp.set_capacitance_2(words_line[1])
                line_counter += 1
            else:
                Sorpt_pri_sp = Surf_species(words_line[0])
                names_pri_sorpt.append(words_line[0])
                Sorpt_pri_sp.moles_component_solid(words_line[1])
                Sorpt_pri_sp.specific_surface_area(words_line[2])
                Sorpt_pri_sp.solid_concentration(words_line[3])
                list_sorption_class.append (Sorpt_pri_sp)
                line_counter += 1
    return names_pri_sorpt, list_sorption_class





################## Auxiliar from auxiliar functions ####################
def remove_coments (list_of_words):
    # Remove the parts that contain comments "#"
    if list_of_words.count('#') == 1 :
        remover_value = list_of_words.index('#')
        list_of_words = list_of_words[:remover_value]
    return list_of_words