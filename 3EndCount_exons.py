#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 11:02:04 2017

@author: yinzh
"""
import sys
from yz_seq_functions import line_2_list

try:
    input_bed = sys.argv[1]
    input_annotation = sys.argv[2]
    output_file_name = sys.argv[3]
except Exception:
    print("Usage: 3EndCount_exons.py input.bed annotation outfile_name")
    sys.exit()
output = open(output_file_name + "_count.txt", 'w+')
count_dict={}

with open (input_annotation) as counts:
    for line in counts:
        line_list = line_2_list(line, "\t")
        if line_list[0] in count_dict:
            count_dict[line_list[0]].update ({int(line_list[1]):0})
        else:
            count_dict.update({line_list[0]:{}})


with open (input_bed) as in_data:
    for line in in_data:
        new_value = 0
        line_list = line_2_list(line, "\t")
        if line_list[0] in count_dict:
            if (int(line_list[2]) in count_dict[line_list[0]]) & (line_list[5]=='+') :
                new_value = count_dict[line_list[0]][int(line_list[2])] + int(line_list[4])
                count_dict[line_list[0]][int(line_list[2])] = new_value
            elif (int(line_list[1]) in count_dict[line_list[0]]) & (line_list[5]=='-') :
                new_value = count_dict[line_list[0]][int(line_list[1])] + int(line_list[4])
                count_dict[line_list[0]][int(line_list[1])] = new_value
            else:
                continue
        else:
            continue
       
for chrom in sorted(count_dict):
    for base in sorted(count_dict[chrom]):
        output.write(chrom + '_' + str(base) + '\t' + str(count_dict[chrom][base]) + '\n')
        
