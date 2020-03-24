#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:19:19 2019

@author: yinzh
"""
# Main
import itertools
import sys
from riboxi_functions import *

try:
    input_read = sys.argv[1]
    umi_length = sys.argv[2]
    output_file_name = sys.argv[3]
    in_line_barcode = sys.argv[4]
    adapter_sequence = sys.argv[5]
except:
    print("Usage: ./move_umi.py ")
    sys.exit(1)

output_read = open(output_file_name + ".fastq", 'w')
umiDict = {}
with open(input_read, 'r') as read2:
    for line1, line2, line3, line4 in itertools.zip_longest(*[read2] * 4):
        if spot_mispriming(line2[len(line2) - len(adapter_sequence) - 1:len(line2):], in_line_barcode) == 1:
            continue
        processed_reads = move_umi(line1, line2, int(umi_length))
        output_read.write(processed_reads[0])
        output_read.write(processed_reads[1])
        output_read.write(line3)
        output_read.write(trim_quality(line4, int(umi_length)))
output_read.close()
