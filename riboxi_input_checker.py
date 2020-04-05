#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 19:12:31 2020

@author: yinzh
"""
import argparse
import os
import sys

from house_keeping_functions import *

parser = argparse.ArgumentParser()
parser.add_argument('-l', "--umi_length", required=True, help="Length of the umi as integer.")
parser.add_argument('-g', '--genome_path', required=True, help="Path to genome.")
parser.add_argument('-b', '--genome_basename', required=True, help="Genome basename.")
parser.add_argument('-s', '--species', required=True, help="Specify name of species file in the current dir.")
parser.add_argument('-a', '--adapter_sequence', required=True, help="3'-adapter sequence including in-line barcode.")
parser.add_argument('-i', '--in_line_barcode', required=True, help="Sequence of in-line barcode.")
parser.add_argument('-c', '--cpu_threads', required=False)
args = parser.parse_args()

try:
    length = num_str_to_int(args.umi_length)
    if length == 0:
        print('UMI length')
        raise IsNotInt

    if os.path.isdir(args.genome_path) is False:
        print(args.genome_path)
        raise CannotOpenFile

    if os.path.isfile(os.curdir + "/" + args.species) is False:
        print(args.species)
        raise CannotOpenFile

    for bases in args.adapter_sequence:
        if bases not in "ATCGN":
            print("Adapter sequence")
            raise IsNotDNABaseWithN

    for bases in args.in_line_barcode:
        if bases not in "ATCG":
            print("In-line barcode")
            raise IsNotDNABase

    threads = num_str_to_int(args.cpu_threads)
    if threads == 0:
        print('CPU threads')
        raise IsNotInt

except IsNotInt:
    print("can only be a non-zero integer.")
    sys.exit(1)
except CannotOpenFile:
    print("cannot be located, please ensure the path is correct.")
    sys.exit(1)
except IsNotDNABaseWithN:
    print("can only contain A, T, C, G and N.")
    sys.exit(1)
except IsNotDNABase:
    print("can only contain A, T, C, G.")
    sys.exit(1)
