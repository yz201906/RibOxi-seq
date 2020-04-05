#!/bin/bash
# -*- coding: utf-8 -*-
#"""
#Created on Thu Dec 05 19:59:37 2018

#@author: yinzh
#"""
set -e
while IFS='' read -r line || [[ -n "$line" ]]; do
  samplelist="$samplelist $line"
done <"$1"
umi_length="$2"
genome_path="$3"
genome="$4"
species="$5"
adaptor_sequence="$6"
in_line_barcode="$7"
cpu_threads="$8"
if [[ -z "$cpu_threads" ]]; then
  cpu_threads=1
fi

riboxi_input_checker.py -l "$umi_length" -g "$genome_path" -b "$genome" -s "$species" -a "$adaptor_sequence" -i "$in_line_barcode" -c "$cpu_threads"

i=1
while [[ $i -le $umi_length ]]; do
  umi_N_bases="${umi_N_bases}N"
  ((i = i + 1))
done
echo "Randomer: ""$umi_N_bases"

for read in $samplelist; do
  if [ ! -f 'trimmed_'"$read"'.fastq' ]; then
    echo "Trimming $read reads...

    "
    #Remove read-through adapter sequences. The fixed sequence under -A option corresponds to immumina universal 5' PCR adapter.
    cutadapt -j $cpu_threads \
      --minimum-length 20:20 \
      -a "$adaptor_sequence" \
      -A "$umi_N_bases"'AGATCGGAAGAGCGGTTCAG' \
      -e 0.1 \
      -o 'dt_'"$read"'_R1.fastq' \
      -p 'dt_'"$read"'_R2.fastq' "$read"'_R1.fastq' "$read"'_R2.fastq' >"$read"'.report'
    #Use PEAR to merge Read1s and Read2s into a single read
    echo "


Merging PE reads with PEAR..."
    echo "Effective command: -f 'dt_'$read'_R1.fastq' -r 'dt_'$read'_R2.fastq' -o $read -n 20 -j $cpu_threads" >>"$read"'.report'
    pear -f 'dt_'"$read"'_R1.fastq' -r 'dt_'"$read"'_R2.fastq' -o "$read" -n 20 -j $cpu_threads >>"$read"'.report'
    #Move the UMI from reads to read identifier line (first line of each fastq record). Also discard reads that do not contain in-line barcode
    echo "


Moving UMIs to read names..."
    echo "Effective command: -r $read'.assembled.fastq' -l $umi_length -o 'umiRemoved_'$read -i $in_line_barcode -a $adaptor_sequence" >>"$read"'.report'
    move_umi.py -r "$read"'.assembled.fastq' -l "$umi_length" -o 'umiRemoved_'"$read" -i "$in_line_barcode" -a "$adaptor_sequence" -m 'pipeline'
    #Remove 3' adapter sequence which also contains the inline barcode for mis-priming mitigation
    echo "


Triming 3' adapter...

    "
    cutadapt -j $cpu_threads \
      --discard-untrimmed \
      --minimum-length 20 \
      -a "$adaptor_sequence" \
      -e 0.2 \
      -o 'trimmed_'"$read"'.fastq' 'umiRemoved_'"$read"'.fastq' >>"$read"'.report'
  else
    echo "


$read reads already processed.
"
  fi
done

#Alignment of processed reads
for sample in $samplelist; do
  R1="trimmed_""$sample"".fastq"
  if [ ! -d "$sample" ]; then
    mkdir -p "$sample"
  fi
  if [ ! -f "$sample""/Aligned.sortedByCoord.out.bam" ]; then
    echo "


Aligning ""$sample"" with STAR..."
    STAR --runThreadN $cpu_threads \
      --outMultimapperOrder Random \
      --outFilterScoreMinOverLread 0.4 \
      --alignEndsType EndToEnd \
      --outFilterMatchNminOverLread 0.4 \
      --outFilterMismatchNoverLmax 0.05 \
      --limitOutSJcollapsed 20000000 \
      --limitIObufferSize 550000000 \
      --limitBAMsortRAM 2160722171 \
      --genomeDir "$genome_path" \
      --readFilesIn "$R1" \
      --outFileNamePrefix ./"$sample"/ \
      --outSAMtype BAM SortedByCoordinate
  else
    echo "


$sample"" already aligned using STAR."
  fi
done

mkdir -p "bed_files"
cp "$species" "bed_files"
mkdir -p "genomecov"
for sample in $samplelist; do
  cd "$sample"
  echo "

Processing ""$sample"" STAR alignment files..."
  samtools index "Aligned.sortedByCoord.out.bam"
  echo "

Deduplicating aligned reads..."
  umi_tools dedup -I "Aligned.sortedByCoord.out.bam" -S "$sample""_dedup.bam" >>"$read"".report"
  echo "

Converting BAM to BED..."
  bedtools bamtobed -i "$sample""_dedup.bam" >"../bed_files/""$sample"".bed"
  echo "

Generating genome coverage..."
  bedtools genomecov \
    -trackline \
    -trackopts "name=""$sample" \
    -bga -3 -ibam "$sample""_dedup.bam" >"../genomecov/""$sample"".genomecov"
  cd ../
done
# shellcheck disable=SC2002
cat "$genome_path""/""$genome"".gtf" | cut -d";" -f1 >"$genome_path""/""$genome""_cut.gtf" #Cut out everything after first ";" from every row of the GTF file for easier parsing
cd "bed_files"
echo "

Counting 3' ends and generating count table...


"
input_list=$(echo "$samplelist" | tr ' ' ,)
echo "
Input samples: $input_list
Effective command: '-b '$input_list' -s '$species' -g '$genome_path'/'$genome'_cut.gtf -2b '$genome_path'/'$genome'.2bit' -c $cpu_threads
"
riboxi_bed_parsing.py -b "$input_list" -s "$species" -g "$genome_path"/"$genome"'_cut.gtf' -2b "$genome_path"/"$genome"'.2bit' -m 'pipeline' -c $cpu_threads
cd ".."
#rm "$genome_path""/""$genome""_cut.gtf"
#rm "dt_*"
#rm "trimmed_*"
#rm_"umi_removed_"
#rm "assembled_*"
#rm "unassembled_*"
