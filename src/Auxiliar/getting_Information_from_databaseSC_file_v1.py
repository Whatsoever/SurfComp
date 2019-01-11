# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 11:44:06 2018

@author: DaniJ
"""
import compilation_functions_reading_file as cfrf


def getting_Information_from_databaseSC_file_v1 (text):
    #text = 'Type1.txt' Reading text
    f = open(text, 'r')
    lines = f.readlines()
    f.close()
    
    # Defining output variables # The creation of these below variables has not much sense, since it can be created directly from the function call, but it helps in keeping a coherence.
    
    Aq_Reaction_list=[]         # Each aqueous reaction is defined in an reaction class, that is given here.
    Sorpt_Reaction_list = []    # Each sorption species is defined in an reaction class, that is given here.
                                # Note: the 4 list defined above can be put together, it is keep separated for a better handle
    names_aqueous_primary_species = []      
    names_aqueous_secondary_species = []
    names_sorption_primary_species = []
    names_sorption_secondary_species = []
    
    # This file altough pretty similar to creating_databaseObject_from_text_type1 is going to proceed differently.
    # First the startindex and the endindex of Primary_Species, Secondary_Species, Surface_Primary, and Surface_Secondary will be find and given.
    [s_index_aq_pri, end_index_aq_pri] = cfrf.find_Indices_Aqueous_PrimaryBlock(lines)
    [s_index_aq_sec, end_index_aq_sec] = cfrf.find_Indices_Aqueous_SecondaryBlock(lines)
    [s_index_sorpt_pri, end_index_sorpt_pri] = cfrf.find_Indices_Sorption_PrimaryBlock(lines)
    [s_index_sorpt_sec, end_index_sorpt_sec] = cfrf.find_Indices_Sorption_SecondaryBlock(lines)
    
    # filling 
    names_aqueous_primary_species, Aq_Species_list_pri= cfrf.read_block_Primary_Species (lines[s_index_aq_pri:end_index_aq_pri])   
    names_aqueous_secondary_species, Aq_Species_list_sec, Aq_Reaction_list = cfrf.read_block_Secondary_Species (lines[s_index_aq_sec:end_index_aq_sec]) 
    names_sorption_primary_species, Sorp_Species_list_pri = cfrf.read_block_Surface_Primary(lines[s_index_sorpt_pri:end_index_sorpt_pri])
    names_sorption_secondary_species, Sorpt_Species_list_sec, Sorpt_Reaction_list = cfrf.read_block_Surface_Secondary (lines[s_index_sorpt_sec:end_index_sorpt_sec]) 
    
    return names_aqueous_primary_species, names_aqueous_secondary_species, names_sorption_primary_species, names_sorption_secondary_species, Aq_Species_list_pri, Aq_Species_list_sec, Sorp_Species_list_pri, Sorpt_Species_list_sec, Aq_Reaction_list, Sorpt_Reaction_list
    
