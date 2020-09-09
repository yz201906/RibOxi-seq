#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# Title     : Prepare data tables and store in a .RData file
# Objective : Generate .RData files for riboxi_shinyapp (count data and annotation table). Takes 2 commandline arguments each as basename of a grouping, first as control, second as treatment.
# Created by: yinzh
# Created on: 5/4/2020
library(readr)
library(tidyverse)
library(DESeq2)
library(ggrepel)

final_counts <- read_delim("final_counts.tsv",
    "\t", escape_double = FALSE, trim_ws = TRUE)
final_counts <- final_counts[,colSums(is.na(final_counts))<nrow(final_counts)]

Control_num<-ncol(dplyr::select(final_counts, starts_with(args[1])))
Treatment_num<-ncol(dplyr::select(final_counts, starts_with(args[2])))

count_data <- final_counts%>%tidyr::unite("chr", chr, base, gene, sep="_")%>%dplyr::select(chr, starts_with(args[1]), starts_with(args[2])) %>% as.data.frame()
raw_data <- final_counts %>% dplyr::select(chr, base, gene, starts_with(args[1]), starts_with(args[2])) %>% as.data.frame()
rownames(count_data)<-count_data$chr
count_data$chr<-NULL

fill <- c(rep(args[1],Control_num),rep(args[2],Treatment_num))
col_data <- data.frame(matrix(data = fill,nrow=(Control_num+Treatment_num),1,dimnames = list(colnames(count_data),c("condition"))))
dds <- DESeqDataSetFromMatrix(countData = count_data,colData = col_data,design = ~ condition)
dds <- estimateSizeFactors(dds)
#print("Using the following factors to normalize counts (Median of ratios)...")
#sizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)

#raw_data<-cbind(dplyr::select(final_counts, chr, base, gene), normalized_counts)
colnames(raw_data)<-c('chr', 'base', 'gene', make.unique(rep(args[1], Control_num), sep = '.'), make.unique(rep(args[2], Treatment_num), sep = '.'))
#rownames(raw_data)<-NULL

print("Performing r log transformation")
rl_dist <- rlogTransformation(dds, fitType='local')
pca_plot_rl_dist <- DESeq2::plotPCA(rl_dist, intgroup="condition") + 
  geom_point() +
  geom_text_repel(aes(label=name)) + 
  theme(panel.background = element_rect(fill = 'transparent', color = 'black'), 
        panel.grid = element_blank(), legend.position = "none")

#raw_data<-cbind(dplyr::select(final_counts, chr, base, gene, seq), dplyr::select(as.data.frame(normalized_counts), starts_with(args[1]), starts_with(args[2])))
#colnames(raw_data)<-c('chr', 'base', 'gene', 'seq', make.unique(rep(args[1], Control_num), sep = '.'), make.unique(rep(args[2], Treatment_num), sep = '.'))

rm(rl_dist, fill)
rm(final_counts)
rm(Control_num)
rm(Treatment_num)
rm(normalized_counts) 
rm(count_data, dds, col_data)

head(raw_data)
print("Saving to RData...")
save.image('raw_data.RData')
print("Done.")
