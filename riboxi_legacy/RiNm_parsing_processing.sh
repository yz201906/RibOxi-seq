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
aligner="$2"
genome_path="$3"
genome="$4"
species="$5"
fa="$6"
cpu_threads="$7"

for read in $samplelist; do
  if [ ! -f "trimmed_"$read".fastq" ]; then
    echo "Trimming $read reads..."
    cutadapt -j $cpu_threads \
      --minimum-length 20:20 \
      -a CTGTAGGCACCATCAATNNNNNN \
      -A AGATCGGAAGAGCGGTTCAG \
      -e 0.1 \
      -o "dt_"$read"_R1.fastq" \
      -p "dt_"$read"_R2.fastq" $read"_R1.fastq" $read"_R2.fastq" >$read".report"
    cutadapt -j $cpu_threads \
      --minimum-length 20:20 \
      -a CTGTAGGCACCAT \
      -A AGATCGGAAG \
      -e 0.1 \
      -o "dt2_"$read"_R1.fastq" \
      -p "dt2_"$read"_R2.fastq" "dt_"$read"_R1.fastq" "dt_"$read"_R2.fastq" >>$read".report"
    pear -f "dt2_"$read"_R1.fastq" -r "dt2_"$read"_R2.fastq" -o $read -n 20 -j $cpu_threads >>$read".report"
    moveUMI_old.py $read".assembled.fastq" "umiRemoved_"$read
    cutadapt -j $cpu_threads \
      --discard-untrimmed \
      --minimum-length 20 \
      -a CTGTAGGCACCATCAAT \
      -e 0.2 \
      -o "trimmed_"$read".fastq" "umiRemoved_"$read".fastq" >>$read".report"
    rm 'dt'*'.fastq'
    rm *'assembled'*
    rm *'discarded'*
    rm 'umi'*
  else
    echo $read" reads already processed."
  fi
done

for sample in $samplelist; do
  R1="trimmed_"$sample".fastq"
  if [ ! -d "$sample" ]; then
    mkdir -p "$sample"
  fi
  if [[ "$aligner" == "STAR" ]]; then
    if [ ! -f $sample"/Aligned.sortedByCoord.out.bam" ]; then
      echo "Aligning "$sample" with STAR..."
      STAR --runThreadN $cpu_threads \
        --outMultimapperOrder Random \
        --outFilterScoreMinOverLread 0.4 \
        --alignEndsType EndToEnd \
        --outFilterMatchNminOverLread 0.4 \
        --outFilterMismatchNoverLmax 0.05 \
        --limitOutSJcollapsed 20000000 \
        --limitIObufferSize 550000000 \
        --limitBAMsortRAM 2160722171 \
        --genomeDir $genome_path \
        --readFilesIn $R1 \
        --outFileNamePrefix ./$sample/ \
        --outSAMtype BAM SortedByCoordinate
    else
      echo $sample" already aligned using STAR."
    fi
  elif [[ "$aligner" == "tophat2" ]]; then
    if [ ! -f $sample"/accepted_hits.bam" ]; then
      echo "Aligning "$sample" with Tophat2..."
      tophat -p $cpu_threads \
        --b2-very-sensitive \
        --library-type fr-secondstrand \
        -o $sample \
        $genome_path/$genome $R1
    else
      echo $sample" already aligned using Tophat2."
    fi
  else
    echo "Aligner has to be either STAR or tophat2."
    exit 1
  fi
done

mkdir -p "bed_files"
mkdir -p "genomecov"
for sample in $samplelist; do
  cd $sample
  if [[ "$aligner" == "STAR" ]]; then
    align_output="Aligned.sortedByCoord.out.bam"
  else
    align_output="accepted_hits.bam"
  fi
  echo "Processing "$sample" STAR alignment files..."
  samtools index $align_output
  echo "Deduplicating aligned reads..."
  umi_tools dedup -I $align_output -S $sample"_dedup.bam" >>$read".report"
  echo "Converting BAM to BED..."
  bedtools bamtobed -i $sample"_dedup.bam" >"../bed_files/"$sample".bed"
  echo "Generating genome coverage..."
  bedtools genomecov \
    -trackline \
    -trackopts "name="$sample \
    -bga -3 -ibam $sample"_dedup.bam" >"../genomecov/"$sample".genomecov"
  cd ../
done
cat $genome_path'/'$genome".gtf" | cut -d';' -f1 >$genome_path'/'$genome"_cut.gtf"
cd "bed_files"
echo "Counting 3' ends..."
Nm_genomic_counts.sh "samplelist" $species $genome_path'/'$genome"_cut.gtf" $genome_path'/'$fa
cd ".."
rm $genome_path'/'$genome"_cut.gtf"
