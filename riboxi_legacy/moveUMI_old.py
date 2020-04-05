#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:19:19 2019

@author: yinzh
"""


def moveUMI(readName, seq):
    umi = seq[-7:-1].rstrip("\n")
    seqUMIRemoved = seq[:-7]
    nameList = line_2_list(readName, " ")
    return (nameList[0] + "_" + umi + " " + nameList[1] + "\n", seqUMIRemoved + "\n")


def trimQuality(readQuality):
    return (readQuality[:-7] + "\n")


# Main
import itertools
import sys

from yz_seq_functions import line_2_list

try:
    inputRead = sys.argv[1]
    outputFileName = sys.argv[2]
except:
    print("Usage: ./2OMePCRDupRm.py ")
    sys.exit(1)

outputRead = open(outputFileName + ".fastq", 'w')
umiDict = {}
with open(inputRead, 'r') as read2:
    for line1, line2, line3, line4 in itertools.izip_longest(*[read2] * 4):
        outputRead.write(moveUMI(line1, line2)[0])
        outputRead.write(moveUMI(line1, line2)[1])
        outputRead.write(line3)
        outputRead.write(trimQuality(line4))
outputRead.close()
