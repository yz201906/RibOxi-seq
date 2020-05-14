## RibOxi-Seq pipeline (Works with Illumina PE platforms)
Custom scripts for processing sequencing reads from fastqs to visualization and 2'-OMe count table

**For overal workflow and visual demonstration, please refer to the PDF file**

### Tested environment
Ubuntu server 18.04 LTS 

Python 3.6.9 and 2.7.11

R 3.6.3

### Required tools and python libraries
*cutadapt* (tested on version 2.8)

*pear* (tested on version 0.9.11)

*STAR* (tested on version 2.7.1)

*samtools* (tested on version 1.8)

*umi_tools* (package also requires: numpy, scipy, cython, pysam, future, regex and matplotlib)

*bedtools* (tested on version 2.28.0)

*pandas*

*twoBitToFa*, *faToTwoBit* (From UCSC https://genome.ucsc.edu/goldenpath/help/twoBit.html)

pylcs (tested on version 0.0.6)

### Required R packages:
CRAN:

*readr, reshape2, tidyverse, shiny, DT*

Bioconductor:

*Gviz, GenomicRanges, BSgenome, BSgenome.Hsapiens.UCSC.hg19, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Mmusculus.UCSC.mm10*


### Setup
Clone the repository or download scripts to desired location. Either ``export`` the directory containing the scripts to ``PATH`` or add a line such as this ``export PATH="path\to\directory:$PATH"`` to the end of .bashrc file. Use command: `chmod +x directory/containing/pipeline/*` to make the scripts executable locally.

For the Shiny app, place the riboxi_shinyapp directory anywhere as long as you can run R command on it.

### Usage
You can go from fastqs all the way to count tables and visualization using ``riboxi_pipeline.sh``, or the python scripts and listed packages can be run separately which might be helpful in first few passes to optimize some of the options.

**Note:** 
The pipeline and scripts assumes read structures demonstrated in the workflow pdf file, that is, randomer sequence at the 3'-end of read1 sequence and [in-line barcode+linker] at the 3'-end of read 2 seuqnece (while the in-line barcode is at the very end relative to the rest of the linker). 
For species of interest, place a file containing chromosome names under the directory where the pipeline is run. The file should have exact same name as the input for species. The chromosomes should be in a single line separated by commas. e.g.:
`chr1,chr2,chr3,chr4,chr5,chr6,chr7 `
**To run the entire pipeline**, first make sure to place a file with samples names (one sample per line) inthe same directory with fastq files; also make sure that you have already generated 2bit genome with the same basename to your reference genome and it is placed in the same directory with your index; then run:
```
riboxi_pipeline.sh [file_containing_sample_list] [umi_length_as_number] [genome_path] [genome_basename] [species_exactly_as_mentioned_above] [3'-adapter_sequence] [in-line_barcode_sequence] [number_of_cpu_threads](optional, default: 1)
```
In the directory where the pipeline is run, there will be intermediary files starting with "dt" and "trimmed" generated, which is from cutadapt read-through adapter trimming step. There are also assembled and unassembled read files from pear. These are for potential troubleshooting, and they can be automatically removed by uncommenting ``riboxi_pipeline.sh`` line 155-160. Trimming, read-merging, msi-priming filtering and de-duplication statistics can be found in files named as "[sample_name].report" or stdout.
There will also be new directories named "bed_files" and "genomecov" in addition to one directory per sample for respective alignment output. The "bed_files" directory houses bed files and final count table, while the genomecov directory houses .genomecov files which can be directly uploaded to UCSC genome browser for visualization.

**ShinyApp**
The count data needs to be pre-processed for the shinyapp, to do this, run ``data_table_prep.R`` inside the directory where 'final_counts.tsv' is located with additional arguments for group parsing. For example, if one sample group is wild type and the other is knock down (if reference genome used for alignment is hg19), then: ``data_table_prep.R WT KD hg19``. An R session image will be saved with processed data table and neccessary variables that will be loaded by the shinyapp.
The app allows extensive data filtering and counts visualization with gene models, which can then be easily downloaded. More capabilities, such as overlaying multiple samples, are in development.


**To use scripts separately**:

``move_umi.py``: make sure read-through adapters are trimmed and R1 and R2 are merged into a single read. Run:
```
move_umi.py -r input_fastq -l length_of_umi_as_number -o output_name -i sequence_of_in-line_barcode -a 3'-adapter_sequence
```
To generate the count table: 

Place your reference_genome.2bit in the same directory with reference genome of interest and make sure they have same base_name.
Prepare the genome annotation file (GTF used here) by cutting out the portion that will not be used:
```javascript
cat "path/to/genome/genome.gtf" | cut -d";" -f1 > "path/to/genome/genome_cut.gtf"
```
Then navigate to where bed files are and run:
```javascript
riboxi_bed_parsing.py -b sample_names_separated_by_commas -s species(file in current dir) -g path/to/genome/genome_cut.gtf -2b genome_name.2bit/path -c [optional: number of cpu processes]
```

## Sample data availability
This is the current working pipeline that we are working with internally, of which sequencing data is not yet publicly available.
Sequencing data from original publications: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102516  and  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96999 can be analyzed using scripts from the legacy_pipeline directory. The differences are attibuted to changes in umi and linker designs.

I have prepared some sample data and prcessed it through the entire pipeline into the ``raw_data.RData``, which is included inside the app directory together with the data table prior to the processing.

### Legacy scripts usage
