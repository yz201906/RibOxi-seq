# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 18:45:48 2017

@author: yinzh
"""

def line_2_list (input_line, separator):
    line_list = (input_line.rstrip("\n")).split(separator)
    return line_list

def expand_gtf (chrom, start, end, output_file):
    for x in range (int(start), int(end) + 1):
        output_file.write(chrom + '\t' + str(x) + '\n')
