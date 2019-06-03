#!/usr/bin/python

from yz_seq_functions import line_2_list, expand_gtf
import sys

try:
    input_gtf = sys.argv[1]
    input_bed = sys.argv[2]
    output_name = sys.argv[3]

except Exception:
    print("Usage:base_annotation.py input_gtf Oxidized.bed output_file")
    sys.exit()
output = open(output_name, 'w+')

gtf_dict = {}
with open(input_gtf, 'r') as gtf:
    for line in gtf:
        line_list = line_2_list(line, "\t")
        if line_list[2] == "CDS":
            if line_list[0] in gtf_dict:
            	if line_list[3] not in gtf_dict[line_list[0]]:
            		gtf_dict[line_list[0]].update({line_list[3]:[line_list[3], line_list[4]]})
                else:
                    continue
            else:
                gtf_dict.update({line_list[0]:{line_list[3]:[line_list[3], line_list[4]]}})

with open(input_bed, 'r') as bed:
    for line in bed:
        line_list = line_2_list(line, "\t")
        if line_list[0] in gtf_dict:
            if int(line_list[4]) >= 1:
                for base in sorted(gtf_dict[line_list[0]]):
            	    base_list=gtf_dict[line_list[0]][base]
                    if int(base_list[0]) <= int(line_list[2]) <= int(base_list[1]):
                        expand_gtf(line_list[0], base_list[0], base_list[1], output)
                        del gtf_dict[line_list[0]][base]
		    else:
			continue

            else:
                continue
        else:
            continue

