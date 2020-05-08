#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# Title     : Prepare data tables and store in a .RData file
# Objective : Generate .RData files for riboxi_shinyapp (count data and annotation table)
# Created by: yinzh
# Created on: 5/4/2020
require(readr)
require(reshape2)
require(tidyverse)
require(shiny)
final_counts <- read_delim("final_counts.tsv",
    "\t", escape_double = FALSE, trim_ws = TRUE)
final_counts <- final_counts[,colSums(is.na(final_counts))<nrow(final_counts)]
counts_melt<-melt(final_counts, id=c("gene","chr","base","seq"))
colnames(counts_melt)<-c("gene", "chr", "base", "seq", "sample_id", "counts")

head(counts_melt)
print("Saving counts to RData...")
saveRDS(counts_melt, file = "raw_data_melt.rds")
print("Done.")
