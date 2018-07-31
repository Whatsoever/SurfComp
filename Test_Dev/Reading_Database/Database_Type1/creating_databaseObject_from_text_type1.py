# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:39 2018

@author: DaniJ
"""

''' Auxiliar function to read database Type1'''

#text = 'Type1.txt'

from Species import Aq_Species
from Reaction import Reaction
#def creating_databaseObject_from_text_type1 (text):
text = 'Type1.txt'
f = open(text, 'r')
lines = f.readlines()
f.close()    
Species_list =[];
Reaction_list=[];
line_counter = 0
while line_counter < len(lines):
    if lines[line_counter].strip()[0] == '#' or lines[line_counter].strip() == '':
        line_counter +=1
    elif lines[line_counter].strip() == 'Primary_Species':
        line_counter +=1
        primary_species_list =[]
       #while NotBlock(lines[i].strip()):
        while lines[line_counter].strip != 'Secondary_Species' and line_counter < len(lines):
            temp =lines[line_counter].strip()
            if temp == '#' or temp == '':
                line_counter +=1
            else:
                words_line = temp.split()
                primary_species_list.append(words_line[0])
                S = Aq_Species(words_line[0])
                S.set_charge (int(words_line[1]))
                S.set_gfw (float(words_line[2]))
                Species_list.append(S)
                line_counter +=1
                
    elif lines[line_counter].strip != 'Secondary_Species':
        line_counter +=1
        secondary_species_list =[]
        while lines[line_counter].strip != 'Primary_Species' and line_counter < len(lines):
            temp =lines[line_counter].strip()
            if temp == '#' or temp == '':
                line_counter +=1
            else:
                words_line = temp.split()
                if words_line[0] == '-log_k':
                    R.set_log_k(words_line[1])
                    Reaction_list.append(R)
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
                    secondary_species_list.append(words_line[0])
                    S.set_charge (int(words_line[1]))
                    Species_list.append(S)
                   # S.set_gfw (int(words_line[1]))  <-- function must be created
  #  return
 
#def NotBlock (instring):
   # b = true
    #if instring == 'Primary_Species':
    #    b = false
   # elif instring == 'Secondary_Species':
   #     b = false
   # return b      