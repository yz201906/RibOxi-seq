#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 00:23:32 2020

@author: yinzh
"""
from pylcs import lcs


def num_str_to_int(num_str):
    for char in num_str:
        if char not in ['1', '2', '3', '4', '5', '6', '7', '8', '9']:
            return 0
    return int(num_str)


def line_2_list(input_line, separator):
    line_list = (input_line.rstrip("\n")).split(separator)
    return line_list


def reverse_compliment(seq):
    """Takes a sequence string and returns a complimentary sequence string."""
    seq_rc = ''
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == "A":
            seq_rc += "T"
        elif seq[i] == "T":
            seq_rc += "A"
        elif seq[i] == "C":
            seq_rc += "G"
        else:
            seq_rc += "C"
    return seq_rc


def sequence_compare(seq_1, seq_2):
    if lcs(seq_1, seq_2) >= len(seq_1) - 1:
        return 1
    else:
        return 2


class Error(Exception):
    """Base class for other exceptions"""
    pass


class CannotOpenFile(Error):
    """Raised when cannot open a file"""
    pass


class IsNotString(Error):
    """Raised when a string input is required"""
    pass


class IsNotInt(Error):
    """Raised when a int input is required"""
    pass


class IsEmpty(Error):
    """Raised when empty input"""
    pass


class IsNotDNABase(Error):
    """Raised when input contains characters other than ATCG"""
    pass


class IsNotDNABaseWithN(Error):
    """Raised when input contains characters other than ATCGN"""
    pass
