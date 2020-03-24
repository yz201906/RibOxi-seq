#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 08 21:09:55 2019

@author: yinzh
"""
import sys
import pandas as pd
import argparse
import os
import copy
import glob
from riboxi_functions import *

# Main
species_list = ['mouse', 'human']
mouse = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
human = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',
         'chrUn_gl000220']
parser = argparse.ArgumentParser()
parser.add_argument("bed_files", help="Names of the files separated with ',':")
parser.add_argument("species", help="human or mouse?")
parser.add_argument("gtf", help="GTF file is required.")
parser.add_argument("genomeTwoBit", help="reference genome fasta file is required.")
args = parser.parse_args()

file_list = line_2_list(args.bed_files, ',')

if args.species not in species_list:
    sys.exit("Error: Species not supported")
else:
    all_chromosomes = args.species
    print(all_chromosomes)

for csv in glob.glob("*.csv"):
    os.remove(csv)
for tsv in glob.glob("*.tsv"):
    os.remove(tsv)
gtf_dict = {}
gtf_list_dict = {}
with open(args.gtf, 'r') as annotation:
    for line in annotation:
        line_list = line_2_list(line, '\t')
        get_gtf_dict(line_list[0], int(line_list[3]), int(line_list[4]), (line_2_list(line, ' '))[1].strip('\"'),
                     gtf_dict)

for chromosome in sorted(gtf_dict):
    for gene in gtf_dict[chromosome]:
        if chromosome not in gtf_list_dict:
            gtf_list_dict.update({chromosome: [[gtf_dict[chromosome][gene][0], gtf_dict[chromosome][gene][1], gene]]})
        else:
            gtf_list_dict[chromosome].append([gtf_dict[chromosome][gene][0], gtf_dict[chromosome][gene][1], gene])

for chromosome in gtf_list_dict:
    gtf_list_dict[chromosome] = sorted(gtf_list_dict[chromosome])

pos_dict = {}
for sample in file_list:
    with open(sample + '.bed', 'r') as sample_count:
        for line in sample_count:
            line_list = line_2_list(line, '\t')
            if line_list[0] not in all_chromosomes:
                continue
            chr_base_dictionary(line_list[0], int(line_list[1]) + 1, int(line_list[2]), line_list[5], pos_dict,
                                gtf_list_dict)

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

df = pd.concat([pd.read_csv(f) for f in glob.glob('*.csv')], axis=1, ignore_index=False)
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
header = header_list[0] + '\t' + header_list[2] + '\t' + 'seq' + '\t' + '\t' + header_list[
    3] + '\t' + '\n'
final_file.write(header)
with open('all_counts.tsv', 'r') as counts:
    next(counts)
    for line in counts:
        seq = ''
        line_list = line_2_list(line, '\t')
        seq = get_rna_seq(line_list[0], int(line_list[1]) - 16, int(line_list[1]) + 15, args.genomeTwoBit)
        if line_list[2] == '-':
            seqRC = reverse_compliment(seq)
            seq = seqRC
        final_file.write(line_list[0] + '\t' + line_list[1] + '\t' + line_list[3].rstrip(
            '_') + '\t' + seq + '\t' + '\t' + line_2_list(line, '_\t')[1] + '\n')
