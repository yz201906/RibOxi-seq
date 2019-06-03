#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:49:39 2017
Updated on Wed Dec 05 13:30:55 2018

@author: yinzh
"""
def line_2_list (input_line, separator):
    line_list = (input_line.rstrip("\n")).split(separator)
    return line_list

# Main
import itertools
import sys
import os
try:
    input_r1 = sys.argv[1]
    input_r2 = sys.argv[2]  
    output_file_name = sys.argv[3]
    randomer_length = sys.argv[4]
except:
    print("Usage: ./2OMePCRDupRm.py read1.fastq read2.fastq linker_sequence output_prefix randomer_length")
    sys.exit(1)

output_read1 = open(output_file_name+"_R1.fastq", 'w+')
output_read2 = open(output_file_name+"_R2.fastq", 'w+')
R2_tmp = open("read_2_tmp", 'w+')
R1_tmp = open("read_1_tmp", 'w+')
scanned_dict = {}
i = int(randomer_length)
print("Randomer sizes=",i)
inclusion_R1=0
inclusion_R2=0

with open (input_r2, 'r') as read2:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read2]*4):
    	if len(line2)>20:
		R2_tmp.write(line2.rstrip('\n')[6:6+i] + ',' + line2.rstrip('\n')[6+i:6+i+4]+"\n")
        else:
		R2_tmp.write('AAAA' + ',' + 'AAAA'+"\n")
R2_tmp.close()

R1_list=[]
with open (input_r1, 'r') as read1:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read1]*4):
        R1_list.append(line2[0:4])

R1_counter=0
with open ("read_2_tmp", 'r') as tmp:
	for line in tmp:
		R1_tmp.write(line.rstrip('\n')+','+R1_list[R1_counter]+'\n')
		R1_counter+=1
R1_tmp.close()

R2_counter=0
read_number={}
with open ("read_1_tmp", 'r') as comparison:
	for line in comparison:
		line_list = line_2_list(line, ',')
		if line_list[0] not in scanned_dict:
			scanned_dict.update({line_list[0]:{line_list[1]:{line_list[2]:''}}})
			read_number.update({str(R2_counter):''})
		elif line_list[1] not in scanned_dict[line_list[0]]:
			scanned_dict[line_list[0]].update({line_list[1]:{line_list[2]:''}})
			read_number.update({str(R2_counter):''})
		elif line_list[2] not in scanned_dict[line_list[0]][line_list[1]]:
			scanned_dict[line_list[0]][line_list[1]].update({line_list[2]:''})
			read_number.update({str(R2_counter):''})
		R2_counter+=1

R2_counter=0
with open (input_r2, 'r') as read2:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read2]*4):
    	if str(R2_counter) in read_number:
	    output_read2.write(line1)
            output_read2.write(line2)
            output_read2.write(line3)
            output_read2.write(line4)
            inclusion_R2+=1
        R2_counter+=1

R1_counter=0
with open (input_r1, 'r') as read2:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read2]*4):
    	if str(R1_counter) in read_number:
    	    output_read1.write(line1)
            output_read1.write(line2)
            output_read1.write(line3)
            output_read1.write(line4)
            inclusion_R1+=1
        R1_counter+=1
print("Input: ",R2_counter," read2s, wrote:",inclusion_R2,"reads.")
print("Input: ",R1_counter," read1s, wrote:",inclusion_R1,"reads.")
os.remove("read_1_tmp")
os.remove("read_2_tmp")
