#!/usr/bin/python
# from subprocess import check_output
from os import popen


def line_2_list(input_line, separator):
    line_list = (input_line.rstrip("\n")).split(separator)
    return line_list


def get_gtf_dict(chr, start, end, gene_id, dict):
    if chr not in dict:
        dict.update({chr: {gene_id: [start, end]}})
    else:
        if gene_id not in dict[chr]:
            dict[chr].update({gene_id: [start, end]})
        else:
            if start < dict[chr][gene_id][0]:
                dict[chr][gene_id][0] = start
            if end > dict[chr][gene_id][1]:
                dict[chr][gene_id][1] = end


def findGeneId(base_list, base, l, r):
    while l <= r:
        mid = l + (r - l) // 2
        if base_list[mid][0] <= base <= base_list[mid][1]:
            return base_list[mid][2]
        elif base_list[mid][0] > base:
            r = mid - 1
        elif base_list[mid][1] < base:
            l = mid + 1
    return -1


def chrBaseDictionary(chr, start, end, strand, dict, gtf_dict):
    if strand == '+':
        if chr in dict:
            if end not in dict[chr]:
                gene_id = findGeneId(gtf_dict[chr], end, 0, len(gtf_dict[chr]) - 1)
                if gene_id != -1:
                    # gene_id='NOT_ANNOTATED'
                    dict[chr].update({end: [0, gene_id, strand]})
        else:
            gene_id = findGeneId(gtf_dict[chr], end, 0, len(gtf_dict[chr]) - 1)
            if gene_id != -1:
                # gene_id='NOT_ANNOTATED'
                dict.update({chr: {end: [0, gene_id, strand]}})
    else:
        if chr in dict:
            if start not in dict[chr]:
                gene_id = findGeneId(gtf_dict[chr], start, 0, len(gtf_dict[chr]) - 1)
                if gene_id != -1:
                    # gene_id='NOT_ANNOTATED'
                    dict[chr].update({start: [0, gene_id, strand]})
        else:
            gene_id = findGeneId(gtf_dict[chr], start, 0, len(gtf_dict[chr]) - 1)
            if gene_id != -1:
                # gene_id='NOT_ANNOTATED'
                dict.update({chr: {start: [0, gene_id, strand]}})


def updateTmpDictionary(chr, start, end, strand, dict):
    if strand == "+":
        if end in dict[chr]:
            dict[chr][end][0] += 1
    else:
        if start in dict[chr]:
            dict[chr][start][0] += 1


def writeToCSV(tmp_dict, sampleName):
    tmp_file = open(sampleName + '.csv', 'w+')
    tmp_file.write('chr,base_,strand_,gene_,' + sampleName + '\n')
    for chr in tmp_dict:
        for base in sorted(tmp_dict[chr]):
            gene = line_2_list(tmp_dict[chr][base][1], '_')[0]
            tmp_file.write(chr + ',' + str(base) + ',' + tmp_dict[chr][base][2] + ',' + gene + '_,' + str(
                tmp_dict[chr][base][0]) + '\n')
    tmp_file.close()


def getRNASeq(chr, start, end, twoBit):
    process = popen(
        'twoBitToFa -seq=' + chr + ' -start=' + str(start) + ' -end=' + str(end) + ' ' + twoBit + ' stdout').read()
    return line_2_list(process, '\n')[1]


def reverseCompliment(seq):
    seqRC = ''
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == "A":
            seqRC += "T"
        elif seq[i] == "T":
            seqRC += "A"
        elif seq[i] == "C":
            seqRC += "G"
        else:
            seqRC += "C"
    return seqRC


def spotMisPriming(seq):
    # Takes string of 31 bases and compare the 17th-22nd to CTGTAGG, returns True or False
    mp = "CTGTAGG"
    match_count = 0
    for i in range(1, 7):
        if mp[i] == seq[i + 16]:
            match_count += 1
    if mp[0] == seq[16] and match_count >= 4:
        return True
    else:
        return False
