#!/bin/bash
# -*- coding: utf-8 -*-
#"""
#Created on Thu Dec 05 19:59:37 2018
#Updated on Sun Dec 09 12:39:05 2018
#@author: yinzh
#"""
while IFS='' read -r line || [[ -n "$line" ]]; do
	samplelist="$samplelist $line"
done < "$1"
read -p "STAR or tophat2? (Case sensitive): " aligner
read -p "Path/to/reference/genome/ and base name? (For example: home/hg19/ hg19. Must have already been built into index and have GTF in the same directory with same base name.)" genome_path genome

if stat --printf='' *.fastq 2>/dev/null; then
	echo "Already decompressed."
else
	gunzip *.fastq.gz
fi

for sample in $samplelist; do
	if (stat --printf='' $sample"_R1.fastq" 2>/dev/null) && (stat --printf='' $sample"_R2.fastq" 2>/dev/null); then
        	echo "Reads already concatenated."
	else
		cat $sample*"R1"*".fastq" > $sample"_R1.fastq"
		cat $sample*"R2"*".fastq" > $sample"_R2.fastq"
	fi
done

for read in $samplelist; do
if [ ! -f "trimmed_"$read"_R2.fastq" ]; then
    	cutadapt -j 5 --minimum-length 20:20 --match-read-wildcards -a NNNNATCACGCTGTAGGCACCATCAATGAC -A NNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -e 0.1 -o "dt_"$read"_R1.fastq" -p "dt_"$read"_R2.fastq" $read"_R1.fastq" $read"_R2.fastq" > $read"_processing_report.txt"
    	#cutadapt -j 5 --minimum-length 20:20 --discard-untrimmed -G GTCATTGATGGTGCCTACAG -e 0.1 -o "dt2_"$read"_R1.fastq" -p "dt2_"$read"_R2.fastq" "dt_"$read"_R1.fastq" "dt_"$read"_R2.fastq" >> $read"_processing_report.txt"
    	#PE_RibOxiPCRDupRm.py "dt2_"$read"_R1.fastq" "dt2_"$read"_R2.fastq" "collapsed_"$read 4 >> $read"_processing_report.txt"
    	#cutadapt -j 5 --minimum-length 20:20 --discard-untrimmed -g NNNN -G ^CGTGATNNNN -e 0.1 -o "trimmed_"$read"_R1.fastq" -p "trimmed_"$read"_R2.fastq" "collapsed_"$read"_R1.fastq" "collapsed_"$read"_R2.fastq" >> $read"_processing_report.txt"

	cutadapt -j 5 --minimum-length 20:20 --match-read-wildcards --discard-untrimmed -g ^NNNN -G GTCATTGATGGTGCCTACAGCGTGATNNNN -e 0.2 -o "trimmed_"$read"_R1.fastq" -p "trimmed_"$read"_R2.fastq" "dt_"$read"_R1.fastq" "dt_"$read"_R2.fastq" >> $read"_processing_report.txt"

else
	echo "Reads already processed."
fi
done

for sample in $samplelist; do
	R1="trimmed_"$sample"_R1.fastq"
        R2="trimmed_"$sample"_R2.fastq"
	if [[ "$aligner" == "STAR" ]]; then
		if [ ! -f $sample"/Aligned.sortedByCoord.out.bam" ]; then
			STAR --runThreadN 6 --outMultimapperOrder Random --genomeDir $genome_path --readFilesIn $R1 $R2 --outFileNamePrefix ./$sample/ --outSAMtype BAM SortedByCoordinate
		else
			echo "Reads alread aligned using STAR."
		fi
	elif [[ "$aligner" == "tophat2" ]]; then
		if [ ! -f $sample"/accepted_hits.bam" ]; then
			tophat -p 6 -r 2 --mate-std-dev 25 --b2-very-sensitive --library-type fr-secondstrand --no-mixed --no-discordant -o $sample $genome_path/$genome $R1 $R2 >> $sample"_processing_report.txt"
		else
                        echo "Reads alread alignedusing Tophat2."
                fi
	else
		echo "Aligner has to be either STAR or tophat2."
		exit 1
	fi
done

for sample in $samplelist; do
	cd $sample
	if [[ "$aligner" == "STAR" ]]; then
		echo "Processing STAR alignment files."
		samtools view -f 255 -b "Aligned.sortedByCoord.out.bam" -o "unique.bam"
		samtools view -F 0x10 -h -b "unique.bam" -o $sample"_forward.bam"
                samtools view -f 0x10 -h -b "unique.bam" -o $sample"_reverse.bam"
                samtools sort $sample"_forward.bam" -o $sample"_forward_sorted.bam"
                samtools sort $sample"_reverse.bam" -o $sample"_reverse_sorted.bam"
                bedtools bamtobed -i $sample"_forward.bam" > $sample"_forward.bed"
                bedtools bamtobed -i $sample"_reverse.bam" > $sample"_reverse.bed"
                3EndCountFromBED_hg19_rRNA.py $sample"_forward.bed" ../$sample"_forward_count.txt"
                3EndCountFromBED_hg19_rRNA.py $sample"_reverse.bed" ../$sample"_reverse_count.txt"
	else
		echo "Processing Tophat2 alignment files."
		#samtools view -F 256 -b "accepted_hits.bam" -o "unique.bam"
		#samtools view -F 0x10 -h -b "unique.bam" -o $sample"_forward.bam"
		#samtools view -f 0x10 -h -b "unique.bam" -o $sample"_reverse.bam"
		#samtools sort $sample"_forward.bam" -o $sample"_forward_sorted.bam"
		#samtools sort $sample"_reverse.bam" -o $sample"_reverse_sorted.bam"
                #bedtools bamtobed -i $sample"_forward.bam" > $sample"_forward.bed"
		#bedtools bamtobed -i $sample"_reverse.bam" > $sample"_reverse.bed"
		bedtools bamtobed -i "accepted_hits.bam" > $sample".bed"
		3EndCountFromBED_hg19_rRNA.py $sample".bed" ../$sample"_forward_count.txt"
		5EndCountFromBED_hg19_rRNA.py $sample".bed" ../$sample"_reverse_count.txt"
	fi
	cd ../
done

