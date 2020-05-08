#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:19:19 2019

@author: yinzh
"""

import argparse
import itertools
import os
import sys

from riboxi_functions import *

# Main
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--input_read', required=True, help="File name of assembled fastqs.")
parser.add_argument('-l', "--umi_length", required=True, help="Length of the umi as integer.")
parser.add_argument('-o', '--output_file_name', required=True, help="Desired name for output file.")
parser.add_argument('-i', '--in_line_barcode', required=True, help="Sequence of in-line barcode.")
parser.add_argument('-a', '--adapter_sequence', required=True, help="3'-adapter sequence including in-line barcode.")
parser.add_argument('-m', '--mode', required=False)
args = parser.parse_args()

try:  # If parameters are not from pipeline, check input.
    if not args.mode or args.mode != "pipeline":
        if os.path.isfile(os.curdir + "/" + args.input_read) is False:
            print(args.species)
            raise CannotOpenFile

        length = num_str_to_int(args.umi_length)
        if length == 0:
            print("UMI length")
            raise IsNotInt

        for bases in args.in_line_barcode:
            if bases not in "ATCGN":
                print("Adapter sequence")
                raise IsNotDNABaseWithN

        for bases in args.adapter_sequence:
            if bases not in "ATCGN":
                print("Adapter sequence")
                raise IsNotDNABaseWithN

    output_read = open(args.output_file_name + ".fastq", 'w')
    umiDict = {}
    reads_in = 0
    reads_out = 0
    with open(args.input_read, 'r') as read2:
        # noinspection PyUnresolvedReferences
        for line1, line2, line3, line4 in itertools.zip_longest(*[read2] * 4):
            reads_in += 1
            if sequence_compare(line2[len(line2) - len(args.adapter_sequence) - 1:len(line2):],
                                args.in_line_barcode) == 2:
                continue
            processed_reads = move_umi(line1, line2, int(args.umi_length))
            output_read.write(processed_reads[0])
            output_read.write(processed_reads[1])
            output_read.write(line3)
            output_read.write(trim_quality(line4, int(args.umi_length)))
            reads_out += 1
    output_read.close()
    print("Number of merged reads as input:" + str(reads_in))
    print("Number of non-mispriming reads:" + str(reads_out))
    print("% passed reads:" + str((reads_out / reads_out)*100) + "%")
except IsNotInt:
    print("can only be a non-zero integer.")
    sys.exit(1)
except CannotOpenFile:
    print("cannot be located, please ensure the path is correct.")
    sys.exit(1)
except IsNotDNABaseWithN:
    print("can only contain A, T, C, G and N.")
    sys.exit(1)
