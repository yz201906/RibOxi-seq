## RibOxi-Seq pipeline (Works with Illumina PE platforms)
Custom scripts for processing sequencing reads from fastqs to visualization and 2'-OMe count table

### Tested environment
Ubuntu server 18.04 LTS 

Python 3.6.9

### Required packages and python libraries
*cutadapt* (tested on version 2.8)

*pear* (tested on version 0.9.11)

*STAR* (tested on version 2.7.1)

*samtools* (tested on version 1.8)

*umi_tools* (package also requires: numpy, scipy, cython, pysam, future, regex and matplotlib)

*bedtools* (tested on version 2.28.0)

*pandas*

*twoBitToFa*, *faToTwoBit* (From UCSC https://genome.ucsc.edu/goldenpath/help/twoBit.html)


### Setup
Clone the repository or download scripts to desired location. Either ``export`` the directory containing the scripts to ``PATH`` or add a line such as this ``export PATH="path\to\directory:$PATH"`` to the end of .bashrc file.

### Usage
You can go from fastqs all the way to count tables and visualization using ``riboxi_pipeline.sh``, or the python scripts and listed packages can be run separately which might be helpful in first few passes to optimize some of the options.

**Note:** 
The pipeline and scripts assumes read structures demonstrated in the workflow pdf file, that is, randomer sequence at the 3'-end of read1 sequence and [in-line barcode+linker] at the 3'-end of read 2 seuqnece (while the in-line barcode is at the very end relative to the rest of the linker). 
Also, additional adjustments to options used in *cutadapt*, *pear*, *STAR* etc. can be done by modifying the ``riboxi_pipeline.sh`` script.
For count table generation, I currently only listed 2 common species, human and mouse, see begginning of ``riboxi_bed_parsing.py``:
```javascript
species_list = ['mouse', 'human']
mouse = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
human = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',
         'chrUn_gl000220']
```
For other species, one can modify the ``species_list`` and add additional lists containing chromosomal informations. For example, after verify exact naming on the corresponding GTF or GFF3 files, make the following changes: 
```javascript
species_list = ['mouse', 'human', 'arabidopsis']
mouse = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
human = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',
         'chrUn_gl000220']
arabidopsis = ['chr1', 'chr2',...]
```

**To run the entire pipeline**, first make sure to place a file with samples names (one sample per line) inthe same directory with fastq files; also make sure that you have already generated 2bit genome with the same basename to your reference genome and it is placed in the same directory with your index; then run:
```
riboxi_pipeline.sh [file_containing_sample_list] [umi_length_as_number] [genome_path] [genome_basename] [species_exactly_as_mentioned_above] [3'-adapter_sequence] [in-line_barcode_sequence]
```
In the directory where the pipeline is run, there will be intermediary files starting with "dt" and "trimmed" generated, which is from cutadapt read-through adapter trimming step. There are also assembled and unassembled read files from pear. These are for potential troubleshooting, and they can be automatically removed by uncommenting ``riboxi_pipeline.sh`` line 104-108. Trimming, read-merging and de-duplication statistics can be found in files named as "[sample_name].report".
There will also be new directories named "bed_files" and "genomecov" in addition to one directory per sample for respective alignment output. The "bed_files" directory houses bed files and final count table, while the genomecov directory houses .genomecov files which can be directly uploaded to UCSC genome browser for visualization.

**To use scripts separately**:

``move_umi.py``: make sure read-through adapters are trimmed and R1 and R2 are merged into a single read. Run:
```
move_umi.py [input_fastq] [length_of_umi_as_number] [output_name] [sequence_of_in-line_barcode] [3'-linker_sequence]
```
To generate the count table: 

Place your reference_genome.2bit in the same directory with reference genome of interest and make sure they have same base_name.
Prepare the genome annotation file (GTF used here) by cutting out the portion that will not be used:
```javascript
cat "path/to/genome/genome.gtf" | cut -d";" -f1 > "path/to/genome/genome_cut.gtf"
```
Then navigate to where bed files are and run:
```javascript
riboxi_genomic_counts.sh [file_containing_sample_list] [species] [path/to/genome/genome_cut.gtf] [genome_name.2bit]
```
The sh script will call ``riboxi_bed_parsing.py``

## Sample data availability
This is the current working pipeline that we are working with internally, of which sequencing data is not yet publicly available.
Sequencing data from original publications: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102516  and  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96999 can be analyzed using scripts from the legacy_pipeline directory. The differences are attibuted to changes in umi and linker designs.

### Legacy scripts usage
