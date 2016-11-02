#' **Script:** `1_mirna_abundance.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/5_mirna_abundance_variance/scripts`
#' 
#' **Date:**  03/15/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Input File(s):** 
#' 
#' 1. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`
#' 
#' 2. `5_normalized_dge_object_and_voom_output.Rata` 
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/5_mirna_abundance_variance`
#' 
#' **Output File(s):** `1_mirna_abundance_normreadcounts_dgecpm_summary.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to determine the overall abundance of each of the 295 miRNA expression profiles in the 174 F2 MSUPRP population.
#' 
#' A comparison will be made between the voom output and the dge object output to see which is the more appropriate miRNA abundance to report. 
#' The counts used to determine total miRNA abundance will be from the dge object, filtered for expression 
#' (at least 1 cpm in at least 1/4 of the samples) and normalized relative to library size (calcNormFactors()).
#' A simple rowSums will calculate the total normalized read counts for each of the 295 miRNAs, and the variance of each miRNA can also be calculated from this data.
#' 
#' THIS ANALYSIS COMPLETED IN R/3.2.0


#' ## Install libraries
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/5_mirna_abundance_variance/scripts")
library(limma)
library(edgeR)

#' ## Load data
#' 
#' Voom-adjusted log-cpm from gpdata pheno object:
load("../../3_build_dge_object_for_eqtl/4_msuprp_mirna_gpdata_pheno_counts.Rdata")

#' Filtered miRNA expression profiles (not normalized read counts)
load("../../3_build_dge_object_for_eqtl/5_normalized_dge_object_and_voom_output.Rata")
ls()

#' ## Analysis
dim(final_MSUPRP_miRNA$pheno)
str(final_MSUPRP_miRNA$pheno)
final_MSUPRP_miRNA$pheno[1:10,1:10,1]

#' Calculate the total normalized abundance of each miRNA across all samples (colSums):
miRNAsum<-colSums(final_MSUPRP_miRNA$pheno)
head(miRNAsum)

#' Reorder them by most abundant to least abundant for easy identification:
ordvoommiRNAsum<-miRNAsum[order(miRNAsum, decreasing = TRUE),]
head(ordvoommiRNAsum)
tail(ordvoommiRNAsum)

#' This is the vector of the 30 most abundant miRNAs in the dataset (~16% of all miRNAs)
ordvoommiRNAsum[1:30]

#' Notice, many of the expected skeletal muscle miRNAs are in the most abundant miRNAs. 
#' 
names(ordvoommiRNAsum[1:30])


#' Look also at the most abundant miRNAs based on just cpm read counts (filtered; output from dge object):
#' 
#' First, take cpm of dge read counts (normalization factors included)
dge$samples[1:5,]
dge$counts[1:5,1:5]
dgcpm<-cpm(dge)
dgcpm[1:5,1:5]

#' Then, do a similar rowSums of the cpm-adjusted dge object (normalization factors included) to get the abundance there (see the difference in distribution)
dgcpmsums<-rowSums(dgcpm)
head(dgcpmsums)

orddgcpmsums<-dgcpmsums[order(dgcpmsums, decreasing = TRUE)]
head(orddgcpmsums)

orddgcpmsums[1:30]
names(orddgcpmsums[1:30])

sum(orddgcpmsums[1:5])
sum(orddgcpmsums)

sum(orddgcpmsums[1:5])/sum(orddgcpmsums)

#' Check that the top 30 miRNAs are the same: 
match(names(orddgcpmsums[1:30]), names(ordvoommiRNAsum[1:30]))
#' Notice slight shift in order of abundance between voom and normalized cpm, but essentially the same order.

#' ## Visualize
#' 
#' Distribution of reads output from voom adjustment (adjusted for statistical accuracy in GBLUP and GWA analysis)
plot(ordvoommiRNAsum/sum(ordvoommiRNAsum))

#' Distribution of normalized cpm for the read counts: (MUCH cooler to look at, and a more relatable picture of the output from sequencing)
#' 
#' Here you can actually see the contribution of each miRNA to the total number of reads in the dataset. Notice the high variability of a few miRNA being highly expressed while most have low expression
plot(orddgcpmsums/sum(orddgcpmsums))


#' ## Save data
save(orddgcpmsums, file="../1_mirna_abundance_normreadcounts_dgecpm_summary.Rdata")