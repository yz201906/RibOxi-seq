#!/usr/bin/python
import argparse
import copy
import glob
import os
import sys
import pandas as pd

from yz_seq_functions import *

# Main
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
if args.species == 'mouse':
    all_chrs = mouse
elif args.species == 'human':
    all_chrs = human
else:
    sys.exit("Error: Please specify either human or mouse.")
print(all_chrs)

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

for chr in sorted(gtf_dict):
    for gene in gtf_dict[chr]:
        if chr not in gtf_list_dict:
            gtf_list_dict.update({chr: [[gtf_dict[chr][gene][0], gtf_dict[chr][gene][1], gene]]})
        else:
            gtf_list_dict[chr].append([gtf_dict[chr][gene][0], gtf_dict[chr][gene][1], gene])

for chr in gtf_list_dict:
    gtf_list_dict[chr] = sorted(gtf_list_dict[chr])

pos_dict = {}
for sample in file_list:
    with open(sample + '.bed', 'r') as sample_count:
        for line in sample_count:
            line_list = line_2_list(line, '\t')
            if line_list[0] not in all_chrs:
                continue
            chrBaseDictionary(line_list[0], int(line_list[1]) + 1, int(line_list[2]), line_list[5], pos_dict,
                              gtf_list_dict)

for sample in file_list:
    tmp_dict = copy.deepcopy(
        pos_dict)  # making a deep copy of the dict for each sample, so counts are recorded separately
    with open(sample + '.bed', 'r') as sample_count:
        for line in sample_count:
            line_list = line_2_list(line, '\t')
            if line_list[0] not in all_chrs:
                continue
            updateTmpDictionary(line_list[0], int(line_list[1]) + 1, int(line_list[2]), line_list[5], tmp_dict)
    writeToCSV(tmp_dict, sample)
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
counts_df.to_csv('allCounts.tsv', sep='\t', index=False)

final_file = open('finalCounts.tsv', 'w+')
allCounts = open('allCounts.tsv', 'r')
header_list = line_2_list(allCounts.readline(), '_\t')
header = header_list[0] + '\t' + header_list[2] + '\t' + 'seq' + '\t' + 'mispriming' + '\t' + header_list[
    3] + '\t' + '\n'
final_file.write(header)
with open('allCounts.tsv', 'r') as counts:
    next(counts)
    for line in counts:
        seq = ''
        line_list = line_2_list(line, '\t')
        seq = getRNASeq(line_list[0], int(line_list[1]) - 16, int(line_list[1]) + 15, args.genomeTwoBit)
        if line_list[2] == '-':
            seqRC = reverseCompliment(seq)
            seq = seqRC
        mispriming = str(spotMisPriming(seq))
        final_file.write(line_list[0] + '\t' + line_list[1] + '\t' + line_list[3].rstrip(
            '_') + '\t' + seq + '\t' + mispriming + '\t' + line_2_list(line, '_\t')[1] + '\n')
