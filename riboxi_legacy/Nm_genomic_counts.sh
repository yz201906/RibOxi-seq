#!/bin/bash
# -*- coding: utf-8 -*-
#"""
#Created on Thu Jul 22 12:30:01 2019

#@author: yinzh
#"""
while IFS='' read -r line || [[ -n "$line" ]]; do
  samplelist="$samplelist $line"
done <"$1"
species=$2
gtf=$3
fa=$4
list=$(echo $samplelist | tr ' ' ,)
echo $list
echo $species
Nm_bed_parsing_new.py $list $species $gtf $fa
