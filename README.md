# RibOxi-Seq pipeline (Works with Illumins PE platforms)
Custom scripts for processing sequencing reads from fastq to visualization and 2'-OMe count table

### Tested environment
Ubuntu server 18.04 LTS, Python 3.6.9
### Required packages and python libraries
*cutadapt* (tested on version 2.8)
*pear* (tested on version 0.9.11)
*STAR* (tested on version 2.7.1)
*samtools* (tested on version 1.8)
*umi_tools*
*bedtools* (tested on version 2.28.0)
*pandas*
### Setup
Clone the repository or download scripts to desired location. Either 'export' the directory containing the scripts to 'PATH' or add a line such as this 'export PATH="path\to\directory:$PATH"' to the end of .bashrc
### Usage
