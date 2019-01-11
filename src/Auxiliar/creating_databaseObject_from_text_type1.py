# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:39 2018

@author: DaniJ
"""

''' Auxiliar function to read database Type1'''

import compilation_functions_reading_file as cfrf


'''
    The FOLLOWING function is quite big, it can be (in the future), slice in smaller sections.
'''




def creating_databaseObject_from_text_type1 (text):
#text = 'Type1.txt'
    f = open(text, 'r')
    lines = f.readlines()
    f.close()    
    [s_index_aq_pri, end_index_aq_pri] = cfrf.find_Indices_Aqueous_PrimaryBlock(lines)
    [s_index_aq_sec, end_index_aq_sec] = cfrf.find_Indices_Aqueous_SecondaryBlock(lines)
    
    names_aqueous_primary_species, Aq_Species_list_pri= cfrf.read_block_Primary_Species (lines[s_index_aq_pri:end_index_aq_pri])   
    names_aqueous_secondary_species, Aq_Species_list_sec, Aq_Reaction_list = cfrf.read_block_Secondary_Species (lines[s_index_aq_sec:end_index_aq_sec]) 
    
    
    
    return names_aqueous_primary_species, names_aqueous_secondary_species, Aq_Species_list_pri, Aq_Species_list_sec, Aq_Reaction_list
 