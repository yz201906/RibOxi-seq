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
adapter_sequence="$6"
in_line_barcode="$7"
cpu_threads="$8"
paired="$9"
linker_sequence=`echo $adapter_sequence | sed "s/in_line_barcode//"`
if [[ -z "$cpu_threads" ]]; then
  cpu_threads=1
fi

riboxi_input_checker.py -l "$umi_length" -g "$genome_path" -b "$genome" -s "$species" -a "$adapter_sequence" -i "$in_line_barcode" -c "$cpu_threads"

i=1
while [[ $i -le $umi_length ]]; do
  umi_N_bases="${umi_N_bases}N"
  ((i = i + 1))
done
echo "Randomer: ""$umi_N_bases"
echo "Full adapter: ""$adapter_sequence"
echo "in_line_barcode: ""$in_line_barcode" 

for read in $samplelist; do
  if [ ! -f 'trimmed_'"$read"'.fastq' ]; then
    if [ $paired = "paired" ]; then
      echo "Trimming $read reads..."
      
      #Remove read-through adapter sequences. The fixed sequence under -A option corresponds to Illumina universal PCR primer.
      cutadapt -j $cpu_threads \
        --minimum-length 20:20 \
        -a "$adapter_sequence" \
        -A "$umi_N_bases"'AGATCGGAAGAGCGGTTCAG' \
        -e 0.2 \
        -o 'dt_'"$read"'_R1.fastq' \
        -p 'dt_'"$read"'_R2.fastq' "$read"'_R1.fastq' "$read"'_R2.fastq' >"$read"'.report'

      #Use PEAR to merge Read1s and Read2s into a single read
      echo "Merging PE reads with PEAR..."
      echo "Effective command: -f 'dt_'$read'_R1.fastq' -r 'dt_'$read'_R2.fastq' -o $read -n 20 -j $cpu_threads" >>"$read"'.report'
      pear -f 'dt_'"$read"'_R1.fastq' -r 'dt_'"$read"'_R2.fastq' -o "$read" -n 20 -j $cpu_threads >>"$read"'.report'

      #Move the UMI from reads to read identifier line (first line of each fastq record).
      echo "Moving UMIs to read names..."
      echo "Effective command: -r $read'.assembled.fastq' -l $umi_length -o 'umiRemoved_'$read -i $in_line_barcode -a $adapter_sequence" >>"$read"'.report'
      move_umi.py -r "$read"'.assembled.fastq' -l "$umi_length" -o 'umiRemoved_'"$read" -i "$in_line_barcode" -a "$adapter_sequence" -m 'pipeline' >>"$read"".report"
      echo "Triming 3' adapter..."
      
      #Remove 3'-adapter including 3' inline barcode
      cutadapt -j $cpu_threads \
        --discard-untrimmed \
        --minimum-length 20 \
        -a "$adapter_sequence" \
        -e 0.1 \
        -o 'trimmed_'"$read"'.fastq' 'umiRemoved_'"$read"'.fastq' >>"$read"'.report'
    else
      echo "Trimming $read reads..."
      
      # remove linker sequence excluding in-line barcode
      cutadapt -j $cpu_threads \
        --minimum-length 20 \
        -a "$linker_sequence" \
        -e 0.2 \
        -o 'dt_'"$read"'.fastq' "$read"'.fastq' >"$read"'.report'
      echo "Moving UMIs to read names..."
      
      #Move the UMI from reads to read identifier line (first line of each fastq record).
      echo "Effective command: -r 'dt_'$read'.fastq' -l $umi_length -o 'umiRemoved_'$read -i $in_line_barcode -a $adapter_sequence" >>"$read"'.report'
      move_umi.py -r 'dt_'$read'.fastq' -l "$umi_length" -o 'umiRemoved_'"$read" -i "$in_line_barcode" -a "$adapter_sequence" -m 'pipeline' >>"$read"".report"
      echo "Triming 3' adapter..."
      
      #Remove 3' inline barcode
      cutadapt -j $cpu_threads \
        --discard-untrimmed \
        --minimum-length 20 \
        -a "$in_line_barcode" \
        -e 0.01 \
        -o 'trimmed_'"$read"'.fastq' 'umiRemoved_'"$read"'.fastq' >>"$read"'.report'
    fi
  else
    echo "$read reads already processed."
  fi
done

#Alignment of processed reads

for sample in $samplelist; do
  R1="trimmed_""$sample"".fastq"
  if [ ! -d "$sample" ]; then
    mkdir -p "$sample"
  fi
  if [ ! -f "$sample"'/Aligned.sortedByCoord.out.bam' ]; then
    echo "Aligning $sample with STAR..."
    align_status=1
    STAR --runThreadN $cpu_threads \
         --outMultimapperOrder Random \
         --outFilterScoreMinOverLread 0.8 \
         --outFilterMatchNminOverLread 0.8 \
         --outFilterMismatchNoverReadLmax 0.1 \
         --outFilterMismatchNoverLmax 0.05 \
         --alignEndsType EndToEnd \
         --limitOutSJcollapsed 20000000 \
         --limitIObufferSize 550000000 \
         --limitBAMsortRAM 2160722171 \
         --genomeDir "$genome_path" \
         --readFilesIn "$R1" \
         --outFileNamePrefix ./"$sample"/ \
         --genomeLoad LoadAndKeep \
         --outSAMtype BAM SortedByCoordinate
  else
    echo "$sample already aligned using STAR."
    align_status=0
  fi
done

if [ $align_status == 1 ]; then
  STAR --genomeLoad Remove \
       --genomeDir "$genome_path"
fi 

mkdir -p "bed_files"
cp "$species" "bed_files"
mkdir -p "genomecov"
for sample in $samplelist; do
  cd "$sample"
  samtools view -h "Aligned.sortedByCoord.out.bam" | sed '/chrM/d' | sed '/KI27/d' | sed '/GL00/d' | samtools view -bSh > "clean.bam"
  samtools view -H clean.bam
  echo "Processing $sample STAR alignment files..."
  samtools index "clean.bam"
  echo "Deduplicating aligned reads..."
  umi_tools dedup -I "clean.bam" -S $sample'_dedup.bam' >>$sample'.report'
  echo "Converting BAM to BED..."
  bedtools bamtobed -i $sample'_dedup.bam' >'../bed_files/'$sample'.bed'
  echo "Generating genome coverage..."
  bedtools genomecov \
    -trackline \
    -trackopts 'name='$sample \
    -bga \
    -3 \
    -ibam $sample"_dedup.bam" > '../genomecov/'$sample'.genomecov'
  cd ../
done

# shellcheck disable=SC2002
grep $'\t'"gene"$'\t' $genome_path'/'$genome'.gtf' > $genome'_subset.gtf'
grep 'protein_coding' $genome'_subset.gtf' > $genome'_coding.gtf'
grep -v 'protein_coding' $genome'_subset.gtf' | grep -v "misc_RNA" > $genome'_non_coding.gtf'
cd "bed_files"
echo "Counting 3 prime ends and generating count table..."
echo "Input samples: $samplelist"
echo "Effective command: -b $samplelist -s $species -g '../'$genome'_subset.gtf' -2b $genome_path'/'$genome'.2bit' -c $cpu_threads"
riboxi_bed_parsing.py -b "$samplelist" -s "$species" -g '../'"$genome" -2b "$genome_path"'/'"$genome"'.2bit' -c "$cpu_threads" -m 'pipeline'
echo "Done."
rm all_counts.tsv
cd ../
rm $genome"_non_coding.gtf"
rm $genome"_coding.gtf"
rm $genome"_subset.gtf"
#rm dt_*
#rm trimmed_*
#rm umi_removed_
#rm assembled_*
#rm unassembled_*
