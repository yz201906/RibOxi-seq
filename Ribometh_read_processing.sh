#!/bin/bash
# -*- coding: utf-8 -*-
#"""
#Created on Sun Dec 09 14:36:31 2018
#@author: yinzh
#"""
while IFS='' read -r line || [[ -n "$line" ]]; do
	samplelist="$samplelist $line"
done < "$1"
read -p "STAR or tophat2? (Case sensitive): " aligner
read -p "Path/to/reference/genome/ and base name? (For example: home/hg19/ hg19. Must have already been built into index with same base name.)" genome_path genome
read -p "Adapter sequence:" adapter
read -p "Output sample postfix (Used for distinguishing data aligned to different genomes/indeces):" postfix

if stat --printf='' *.fastq 2>/dev/null; then
	echo "Already decompressed."
else
	gunzip *.fastq.gz
fi

for sample in $samplelist; do
	if (stat --printf='' $sample"_R1.fastq" 2>/dev/null); then
        	echo "Reads already concatenated."
	else
		echo "Concatenating..."
		cat *$sample*"R1"*".fastq" > $sample"_R1.fastq"
	fi
done

for sample in $samplelist; do
	if [ ! -f "trimmed_"$sample"_R1.fastq" ]; then
    	cutadapt -j 6 --minimum-length 20 --discard-untrimmed -a $adapter -o "trimmed_"$sample"_R1.fastq" $sample"_R1.fastq" > $sample"_processing_report.txt"
	else
		echo "Reads already processed."
	fi
done

for sample in $samplelist; do
	R1="trimmed_"$sample"_R1.fastq"
	if [[ "$aligner" == "STAR" ]]; then
		if [ ! -f $sample$postfix"/Aligned.sortedByCoord.out.bam" ]; then
			STAR --runThreadN 5 --sjdbOverhang 30 --genomeSAindexNbases 5 --outMultimapperOrder Random --genomeDir $genome_path --readFilesIn $R1 --outFileNamePrefix ./$sample$postfix/ --outSAMtype BAM SortedByCoordinate
		else
			echo "Reads alread aligned using STAR."
		fi
	elif [[ "$aligner" == "tophat2" ]]; then
		if [ ! -f $sample$postfix"/accepted_hits.bam" ]; then
			tophat -g 1 -p 6 --b2-very-sensitive --no-coverage-search -o $sample$postfix $genome_path/$genome $R1 >> $sample"_processing_report.txt"
		else
        	echo "Reads alread aligned using Tophat2."
        fi
	else
		echo "Aligner has to be either STAR or tophat2."
		exit 1
	fi
done

for sample in $samplelist; do
	cd $sample$postfix
	if [[ "$aligner" == "STAR" ]]; then
		echo "Processing STAR alignment files."
		samtools view -F 256 "Aligned.sortedByCoord.out.bam" | samtools view -b | samtools sort -o "sorted_"$sample$postfix".bam"
                bedtools bamtobed -i "sorted_"$sample$postfix".bam" > $sample$postfix".bed"
	else
		echo "Processing Tophat2 alignment files." 
		samtools view -F 256 "accepted_hits.bam" -b -o "unique.bam"
		samtools sort "unique.bam" -o "sorted_"$sample$postfix".bam"
		bedtools bamtobed -i "sorted_"$sample$postfix".bam" > $sample$postfix".bed"
	fi
	cd ../
done

