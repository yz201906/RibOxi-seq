#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 08 21:09:55 2019

@author: yinzh
"""
import argparse
import copy
import glob
import multiprocessing as mp
import os
import sys

import pandas as pd

from riboxi_functions import *

# Main
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bed_files', required=True, help="Names of the files separated with commas.")
parser.add_argument('-s', '--species', required=True, help="Specify species file in the current dir.")
parser.add_argument('-g', '--gtf', required=True, help="GTF file is required (path to gtf file).")
parser.add_argument('-2b', "--genome_2_bit", required=True, help="Path to the 2bit genome.")
parser.add_argument('-c', '--cpu_threads', required=False)
parser.add_argument('-m', '--mode', required=False)
args = parser.parse_args()

try:  # If parameters are not from pipeline, check input.
    if not args.mode or args.mode != "pipeline":
        if len((args.bed_files.lstrip(',')).rstrip(',')) == 0:
            print("Sample list")
            raise IsEmpty
        for samples in args.bed_files:
            if os.path.isfile(os.curdir + "/" + samples + '.bed') is False:
                print(samples + '.bed')
                raise CannotOpenFile
        if os.path.isfile(os.curdir + "/" + args.species) is False:
            print(args.species)
            raise CannotOpenFile
        if os.path.isfile(args.gtf) is False:
            print(args.gtf)
            raise CannotOpenFile
        if os.path.isfile(args.genome_2_bit) is False:
            print(args.genome_2_bit)
            raise CannotOpenFile
        if args.cpu_threads:
            threads = num_str_to_int(args.cpu_threads)
            if threads == 0:
                print('CPU threads')
                raise IsNotInt

except IsEmpty:
    print("cannot be empty.")
    sys.exit(1)
except CannotOpenFile:
    print("cannot be located, please ensure the path is correct.")
    sys.exit(1)

threads = 1
if args.cpu_threads:
    threads = num_str_to_int(args.cpu_threads)
    if threads >= mp.cpu_count():
        threads = mp.cpu_count() - 1

input_samples = args.bed_files.strip(' ')
if ' ' in args.bed_files:
    file_list = line_2_list(input_samples, ' ')
else:
    file_list = [input_samples]
print(file_list)

my_species = open(args.species, 'r')
all_chromosomes = line_2_list(my_species.readline(), ',')
print(all_chromosomes)

for csv in glob.glob("*.csv"):
    os.remove(csv)
for tsv in glob.glob("*.tsv"):
    os.remove(tsv)
coding_dict = {}
non_coding_dict = {}
with open(args.gtf+'_coding.gtf', 'r') as annotation:
    for line in annotation:
        line_list = line_2_list(line, '\t')
        gtf_gene_info = line_2_list((line_list.pop(8)), ";")
        for naming in gtf_gene_info:
            if "gene_name" in naming:                              #Parse for gene names
                gtf_gene_name = line_2_list(naming, " ")[2]
        if line_list[0] not in all_chromosomes:
            continue
        get_gtf_dict(line_list[0], int(line_list[3]), int(line_list[4]), gtf_gene_name.strip('\"'), coding_dict)
with open(args.gtf+'_non_coding.gtf', 'r') as annotation:
    for line in annotation:
        line_list = line_2_list(line, '\t')
        gtf_gene_info = line_2_list((line_list.pop(8)), ";")
        for naming in gtf_gene_info:
            if "gene_name" in naming:                              #Parse for gene names
                gtf_gene_name = line_2_list(naming, " ")[2]
        if line_list[0] not in all_chromosomes:
            continue
        get_gtf_dict(line_list[0], int(line_list[3]), int(line_list[4]), gtf_gene_name.strip('\"'), non_coding_dict)

pos_dict = {}
for sample in file_list:
    with open(sample + '.bed', 'r') as sample_count:
        for line in sample_count:
            line_list = line_2_list(line, '\t')
            if line_list[0] not in all_chromosomes:
                continue
            chr_base_dictionary(line_list[0], int(line_list[1]) + 1, int(line_list[2]), line_list[5], pos_dict, coding_dict, non_coding_dict)

for sample in file_list:
    tmp_dict = copy.deepcopy(
        pos_dict)  # making a deep copy of the dict for each sample, so counts are recorded separately
    with open(sample + '.bed', 'r') as sample_count:
        for line in sample_count:
            line_list = line_2_list(line, '\t')
            if line_list[0] not in all_chromosomes:
                continue
            update_tmp_dictionary(line_list[0], int(line_list[1]) + 1, int(line_list[2]), line_list[5], tmp_dict)
    write_to_csv(tmp_dict, sample)
pos_dict = {}
gtf_list_dict = {}

df = pd.concat([pd.read_csv(f, dtype={"chr": "string"}) for f in glob.glob('*.csv')], axis=1, ignore_index=False)
head = df.head()
print(head)
for csv in glob.glob("*.csv"):
    os.remove(csv)
counts_df = df.loc[:, ~df.columns.duplicated()]
head = counts_df.head()
print(head)
counts_df.to_csv('all_counts.tsv', sep='\t', index=False)

final_file = open('final_counts.tsv', 'w+')
all_counts = open('all_counts.tsv', 'r')
header_list = line_2_list(all_counts.readline(), '_\t')
header = header_list[0] + '\t' + header_list[2] + '\t' + header_list[3] + '\t' + '\n'
#header = header_list[0] + '\t' + header_list[2] + '\t' + 'seq' + '\t' + header_list[3] + '\t' + '\n'
final_file.write(header)
print("The header is: " + header)
all_counts.close()

pool = mp.Pool(threads)
print("Using " + str(threads) + " cpu threads to generate table annotation...")
with open('all_counts.tsv', 'r') as counts:
    next(counts)
    # results = [pool.apply(get_annotated_table, args=(line, args.genome_2_bit)) for line in counts]
    # for row in results:
    #     final_file.write(row)
        
    results = [pool.apply(get_annotated_table, args=(line,)) for line in counts] #without seq
    for row in results:
        final_file.write(row)
    
    # for line in counts:                                                    #without multiprocess
    #     final_file.write(get_annotated_table(line, args.genome_2_bit))
    
pool.close()

final_file.close()
counts.close()
