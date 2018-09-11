# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 18:34:31 2018

@author: DaniJ
"""

def Read_input_file_type1 (text):
#text = 'Type1.txt'
    f = open(text, 'r')
    lines = f.readlines()
    f.close()   
    list_names = []
    list_values = []
    line_counter = 0
    while line_counter < len(lines):
        if lines[line_counter].strip() == '' or lines[line_counter].strip()[0] == '#':
            line_counter += 1
        elif lines[line_counter].strip() == 'Solution':
            line_counter += 1
            while line_counter < len(lines):
                temp =lines[line_counter].strip()
                if temp == '' or temp[0] == '#':
                    line_counter += 1
                else:
                    words_line = temp.split()
                    list_names.append(words_line[0])
                    list_values.append(float(words_line[1]))
                    line_counter += 1
    return list_names, list_values