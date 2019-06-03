#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:49:39 2017

@author: yinzh
"""
def sequence_mismatch (input_seq, input_dict, input_linker_length, j):
    if input_seq[(input_linker_length+j):(input_linker_length+j+6):] in input_dict[input_seq[0:j:]]:
        return 1
    else:
        return 0



# Main
import itertools
import sys
try:
    input_r1 = sys.argv[1]
    input_r2 = sys.argv[2]
    linker_seq = sys.argv[3]  #linker sequence
    output_file_name = sys.argv[4]
    randomer_length = sys.argv[5]
except:
    print("Usage: ./2OMePCRDupRm.py read1.fastq read2.fastq linker_sequence output_prefix randomer_length")
    sys.exit(1)

output_read1 = open(output_file_name+"_R1.fastq", 'w')
output_read2 = open(output_file_name+"_R2.fastq", 'w')
scanned_dict = {}
linker_length = len(linker_seq)
i = int(randomer_length)
R2_counter=0
R2_list={}
R1_counter=0
inclusion_R1=0
inclusion_R2=0
with open (input_r2, 'r') as read2:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read2]*4):
        if line2[0:i:] not in scanned_dict:
            scanned_dict.update({line2[0:i:]: {line2[(linker_length + i):(linker_length + i + 6):]:line2[(linker_length + i):(linker_length + i + 6):]}})
            output_read2.write(line1)
            output_read2.write(line2)
            output_read2.write(line3)
            output_read2.write(line4)
            R2_list.update({str(R2_counter):'0'})
            inclusion_R2+=1
        else:
            comp = sequence_mismatch (line2, scanned_dict, linker_length, i)
            if comp == 1:
                R2_counter+=1
                continue
            else:
                scanned_dict[line2[0:i:]].update({line2[(linker_length + i):(linker_length + i + 6):]:line2[(linker_length + i):(linker_length + i + 6):]})
                output_read2.write(line1)
                output_read2.write(line2)
                output_read2.write(line3)
                output_read2.write(line4)
                R2_list.update({str(R2_counter):'0'})
                inclusion_R2+=1
        R2_counter+=1
output_read2.close()
with open (input_r1, 'r') as read1:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read1]*4):
        if str(R1_counter) in R2_list:
            output_read1.write(line1)
            output_read1.write(line2)
            output_read1.write(line3)
            output_read1.write(line4)
            inclusion_R1+=1
        R1_counter += 1
output_read1.close()

print("total ",R2_counter," read2s, wrote:",inclusion_R2,"reads.")
print("total ",R1_counter," read1s, wrote:",inclusion_R1,"reads.")
