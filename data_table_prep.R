#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# Title     : Prepare data tables and store in a .RData file
# Objective : Generate .RData files for riboxi_shinyapp (count data and annotation table). Takes 2 commandline arguments each as basename of a grouping, first as control, second as treatment.
# Created by: yinzh
# Created on: 5/4/2020
library(readr)
library(tidyverse)

final_counts <- read_delim("final_counts.tsv",
    "\t", escape_double = FALSE, trim_ws = TRUE)
final_counts <- final_counts[,colSums(is.na(final_counts))<nrow(final_counts)]

Control_num<-ncol(dplyr::select(final_counts, starts_with(args[1])))
Treatment_num<-ncol(dplyr::select(final_counts, starts_with(args[2])))
raw_data<-cbind(dplyr::select(final_counts, chr, base, gene, seq,), dplyr::select(final_counts, starts_with(args[1])), dplyr::select(final_counts, starts_with(args[2])))
colnames(raw_data)<-c('chr', 'base', 'gene', 'seq', make.unique(rep(args[1], Control_num), sep = '.'), make.unique(rep(args[2], Treatment_num), sep = '.'))

counts_sum <- raw_data %>% select(starts_with(args[1])|starts_with((args[2])))%>% summarise_all(sum)
counts_means <- rowMeans(select(counts_sum, starts_with(args[1])))
counts_sum
norm_factors <- counts_sum %>% summarise_all(funs(./ counts_means))
norm_factors
sample_names <- colnames(counts_sum)
sample_names
for (samples in sample_names) {
  raw_data[[samples]] <- (raw_data[[samples]])/(norm_factors[[samples]])
}

rm(final_counts)
rm(Control_num)
rm(Treatment_num)
rm(counts_sum, norm_factors, counts_means, sample_names, samples)

head(raw_data)
print("Saving to RData...")
save.image('raw_data.RData')
print("Done.")
