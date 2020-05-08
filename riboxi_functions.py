#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 18:21:10 2018

@author: yinzh
"""
from os import popen

from house_keeping_functions import *


def move_umi(read_name, seq, umi_length):
    """Takes in read names and actual read sequence and length umi, removes first n bases (n = umi_length)
    from sequence and append to read name line"""
    umi = seq[0:umi_length]
    seq_umi_removed = (seq[umi_length::]).rstrip("\n")
    name_list = line_2_list(read_name, " ")
    return name_list[0] + "_" + umi + " " + name_list[1] + "\n", seq_umi_removed + "\n"


def trim_quality(read_quality, umi_length):
    """Takes in read quality line of a fastq record and removes the first n characters (n=umi_length)"""
    return read_quality[umi_length::]


def get_gtf_dict(chromosome, start, end, gene_id, gtf_dict):
    """Parses GTF annotation lines and generate gtf_dictionary that has chromosomal information for each gene."""
    if chromosome not in gtf_dict:
        gtf_dict.update({chromosome: {gene_id: [start, end]}})
    else:
        if gene_id not in gtf_dict[chromosome]:
            gtf_dict[chromosome].update({gene_id: [start, end]})
        else:
            if start < gtf_dict[chromosome][gene_id][0]:
                gtf_dict[chromosome][gene_id][0] = start
            if end > gtf_dict[chromosome][gene_id][1]:
                gtf_dict[chromosome][gene_id][1] = end


def find_gene_id(base_list, base, left, right):
    """Takes a list as input (a list of [start position, end position, gene id] *pre-sorted*). Also takes in the
    arguments of base position, 0, length of list, which are used to search the base position in the list and returns
    the corresponding gene id. Implemented with binary search"""
    while left <= right:
        mid = left + (right - left) // 2
        if base_list[mid][0] <= base <= base_list[mid][1]:
            return base_list[mid][2]
        elif base_list[mid][0] > base:
            right = mid - 1
        elif base_list[mid][1] < base:
            left = mid + 1
    return -1


def chr_base_dictionary(chromosome, start, end, strand, pos_dict, gtf_list_dict):
    """Takes an empty or non-empty dictionary and an annotation dictionary (gtf_list_dict). Create a list inside
    dictionary inside dictionary nested structure in dict {Chr1:{3'end position:[count, gene id, strandedness]}}. """
    if strand == '+':
        if chromosome in pos_dict:
            if end not in pos_dict[chromosome]:
                gene_id = find_gene_id(gtf_list_dict[chromosome], end, 0, len(gtf_list_dict[chromosome]) - 1)
                if gene_id != -1:
                    pos_dict[chromosome].update({end: [0, gene_id, strand]})
        else:
            gene_id = find_gene_id(gtf_list_dict[chromosome], end, 0, len(gtf_list_dict[chromosome]) - 1)
            if gene_id != -1:
                pos_dict.update({chromosome: {end: [0, gene_id, strand]}})
    else:
        if chromosome in pos_dict:
            if start not in pos_dict[chromosome]:
                gene_id = find_gene_id(gtf_list_dict[chromosome], start, 0, len(gtf_list_dict[chromosome]) - 1)
                if gene_id != -1:
                    pos_dict[chromosome].update({start: [0, gene_id, strand]})
        else:
            gene_id = find_gene_id(gtf_list_dict[chromosome], start, 0, len(gtf_list_dict[chromosome]) - 1)
            if gene_id != -1:
                pos_dict.update({chromosome: {start: [0, gene_id, strand]}})


def update_tmp_dictionary(chromosome, start, end, strand, tmp_dict):
    """Takes a nested dictionary object along with chromosome number, start, end positions and strandedness as input.
    Count the ends and update the dictionary. """
    if strand == "+":
        if end in tmp_dict[chromosome]:
            tmp_dict[chromosome][end][0] += 1
    else:
        if start in tmp_dict[chromosome]:
            tmp_dict[chromosome][start][0] += 1


def write_to_csv(tmp_dict, sample_name):
    """Takes a dictionary, sort it by chromosome, and write to csv with a file name given by sample_name argument"""
    tmp_file = open(sample_name + '.csv', 'w+')
    tmp_file.write('chr,base_,strand_,gene_,' + sample_name + '\n')
    for chromosome in tmp_dict:
        for base in sorted(tmp_dict[chromosome]):
            gene = line_2_list(tmp_dict[chromosome][base][1], '_')[0]
            tmp_file.write(chromosome + ',' + str(base) + ',' + tmp_dict[chromosome][base][2] + ',' + gene + '_,' + str(
                tmp_dict[chromosome][base][0]) + '\n')
    tmp_file.close()


def get_rna_seq(chromosome, start, end, two_bit):
    """Takes chromosome, sequence start and end positions and a binary 2bit genome. Returns the sequence string."""
    process = popen(
        'twoBitToFa -seq=' + chromosome + ' -start=' + str(start) + ' -end=' + str(
            end) + ' ' + two_bit + ' stdout').read()
    return line_2_list(process.upper(), '\n')[1]

def get_annotated_table(table_row, genome_2_bit):
    line_list = line_2_list(table_row, '\t')
    seq = get_rna_seq(line_list[0], int(line_list[1]) - 16, int(line_list[1]) + 15, genome_2_bit)
    if line_list[2] == '-':
        seq_rc = reverse_compliment(seq)
        seq = seq_rc
    out_string = line_list[0] + '\t' + line_list[1] + '\t' + line_list[3].rstrip('_') + '\t' + seq + '\t' + '\t' + \
                 line_2_list(table_row, '_\t')[1] + '\n'
    return out_string