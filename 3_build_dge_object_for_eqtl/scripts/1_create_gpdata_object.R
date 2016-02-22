#' **Script:** `1_create_gpdata_object.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`
#' 
#' **Date:**  2/19/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
#' 2. `/mnt/research/pigsnp/MSUPRP/carcass_quality/data/`
#' 3. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/`
#' 
#' **Input File(s):** 
#' 
#' 1. `1_mean_mature_mirna_expression_unfiltered.Rdata`
#' 2. `MSUPRP_meat.RData`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. `3_msuprp_mirna_gpdata.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)
#' 
#' ## Objectives
#' The objective of this script is to create the gpdata object needed for the miRNA eQTL analysis.
#' This gpdata object will be filtered from the MSUPRP_meat gpdata object in the following ways:
#' 
#' 1. The covariate data will be reduced to include only animals in this analysis (n=174), and will have the additional column of growth_group as a factor.  
#' 
#' 2. The geno data will be reduced to include only animals in this analysis (n=174), and filtered for removal of fixed SNPs, SNPs with maf < 0.10, and SNPs on sex chromosomes. 
#' 
#' ## This analysis conducted using R/3.2.0, not R/3.1.0
#' 
#' ## Install libraries
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts/")

library(synbreed)
library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' ## Load data
rm(list=ls())
#' Load the miRNA expression data:
load("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")
#' Load the MSUPRP_meat gpdata object:
load("/mnt/research/pigsnp/MSUPRP/carcass_quality/data/MSUPRP_meat.RData")

#' ## Analysis
#' ### 1. A matrix of gene expression, with dimensions genes x samples (335 x 174)
dim(no.zero.dfmeanrcround)
colnames(no.zero.dfmeanrcround)
head(rownames(no.zero.dfmeanrcround))

#' Create vector of animal IDs from column names of expression matrix:
pigid <- colnames(no.zero.dfmeanrcround)
length(pigid)
pigid

#' ### 2. A data frame of metadata (trait data) with dimensions samples x categories (174 x categories), including growth_group as a factor:
#' 
#' 
#' Discard the animals from the gpdata geno object not present in the miRNA expression matrix:
todisc <- rownames(MSUPRP_meat$geno)[!rownames(MSUPRP_meat$geno) %in% pigid]

redMSU <- discard.individuals(MSUPRP_meat, todisc)

#' Remove covariate data from animals not in this analysis:
todisc <- redMSU$covar$id[!redMSU$covar$id %in% rownames(redMSU$geno)]

redMSU <- discard.individuals(redMSU, todisc)

#' Create the growth_group column as a factor by combining selcrit and Status
redMSU$covar <- data.frame(redMSU$covar,
        growth_group=paste(redMSU$covar$selcrit, redMSU$covar$Status, sep="-"))

#' Change the levels of the factor growth_group to have lma-L the first level
redMSU$covar$growth_group<-relevel(redMSU$covar$growth_group, ref = "lma-L")

is.factor(redMSU$covar$growth_group)

head(redMSU$covar)

if (sum(redMSU$covar$id!=pigid)!=0) stop ("pigids of redMSU$covar not correct")

#' The dataframe of covariates is complete.
#' 
#' ### 3. A matrix of genotype data, filtered for removing fixed SNPs, SNPs with low maf, and SNPs on sex chromosomes. 

#' Extract the genotypes from the GPData object:
genomat <- redMSU$geno

head(genomat[,1:10])

dim(genomat)

if (sum(rownames(genomat)!=redMSU$covar$id) != 0) stop ("rows of geno not the same as ids of covar")
if (sum(colnames(genomat)!=colnames(MSUPRP_meat$geno))!=0) stop ("columns of genomat not the same as MSUPRP_meat$geno")

#' Filter the matrix of genotypes: 
#' 
#' Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):
sdv <- apply(genomat, 2, sd)
length(sdv)
sum(sdv == 0)

#' Remove fixed SNPs from genotype matrix:
genomat <- genomat[,sdv>0]

dim(genomat)

#' Filter SNPs for minor allele frequency (maf > 0.10 retained):
#'
#' First, get the allelic frequency:
af <- colMeans(genomat)/2
summary(af)

#' How many SNPs have a af < 0.10 and >0.90?
#' 
#' We filter this way since we don't 'know' the minor allele in each SNP based on allele frequency:
sum(af<0.10)
sum(af>0.90)
sum(c(af<0.10,af>0.90))

#' Retain all SNPs with af > 0.10 and <0.90 (maf < 0.10 discarded):
genomat <- genomat[ ,af>0.10 & af<0.90]

#' Dimensions of remaining SNP marker matrix:
dim(genomat)

if (sum((colMeans(genomat)/2)<=0.10) !=0) stop ("maf filtering did not work correctly")
if (sum((colMeans(genomat)/2)>=0.90) !=0) stop ("maf filtering did not work correctly")
if (sum(rownames(genomat)!=redMSU$covar$id)!=0) stop ("rownames of marker matrix not equal to trait data")

#' Eliminate markers on sex chromosomes:
table(redMSU$map$chr)
sexchr <- rownames(redMSU$map)[redMSU$map$chr == 19]
length(sexchr)

sum((colnames(genomat) %in% sexchr))

genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]

#' This object contains the markers that are not fixed, have a "maf" > 0.10, and are not mapped to the sex chromosomes:
dim(genomatfil)

if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(redMSU$covar$id != rownames(genomatfil)) != 0) stop ("rownames of trait data and genotype matrix not the same")
if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")

#' Create the todel vector, containing the names of the SNPs in the gpdata geno object NOT found in the filtered genotype matrix:
todel <- colnames(redMSU$geno)[!colnames(redMSU$geno) %in% colnames(genomatfil)]
length(todel)

#' That means that, in total, we are removing 7010 SNPs from the gpdata geno object. 
#' 
#' Using discard.markers allows us to remove both the markers we don't want, and the map information we don't want, all in one step.
final_MSUPRP_miRNA<-discard.markers(redMSU, todel)

#' The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 
summary_final_MSUPRP_miRNA<-summary(final_MSUPRP_miRNA)
summary_final_MSUPRP_miRNA$geno
summary_final_MSUPRP_miRNA$pedigree

#' 
#' ## Save data
#' 
#' Save the gpdata object with filtered genotypes, map information, and covariate information:
save(final_MSUPRP_miRNA, summary_final_MSUPRP_miRNA, file="../3_msuprp_mirna_gpdata.Rdata")