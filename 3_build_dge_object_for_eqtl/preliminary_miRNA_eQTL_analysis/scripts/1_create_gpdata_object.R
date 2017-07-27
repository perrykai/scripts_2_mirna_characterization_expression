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
#' The genotype data will be compiled as follows:
#' 
#' 1. First, extract the F2 population genotypes from the MSUPRP_meat$geno object and calculate MAFs
#' 
#' 2. Subset the genotype data to include the 174 F2 animals in this analysis
#' 
#' 3. Filter those genotypes for removal of fixed SNPs, those located on sex chromosomes, and MAF < 0.10 (being sure to take both ends of the distribution)
#' 
#' The covariate data will be reduced to include only animals in this analysis (n=174), and will have the additional column of growth_group as a factor.  
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
#' 
#' no.zero.dfmeanrcround is a matrix of gene expression, with dimensions genes x samples (335 x 174)
dim(no.zero.dfmeanrcround)
colnames(no.zero.dfmeanrcround)
head(rownames(no.zero.dfmeanrcround))

#' Create vector of animal IDs from the column names of the expression matrix for subsetting the gpdata object:
pigid <- colnames(no.zero.dfmeanrcround)
length(pigid)
pigid

#' ### A data frame of metadata (trait data) with dimensions samples x categories (174 x categories), including growth_group as a factor:
#' 
#' 
#' Remove covariate data from animals not in this analysis:
todisc <- MSUPRP_meat$covar$id[!MSUPRP_meat$covar$id %in% pigid]
length(todisc)
redMSU <- discard.individuals(MSUPRP_meat, todisc)

dim(redMSU$pedigree)
dim(redMSU$geno)

#' Create the growth_group column as a factor by combining selcrit and Status
redMSU$covar <- data.frame(redMSU$covar,
        growth_group=paste(redMSU$covar$selcrit, redMSU$covar$Status, sep="-"))

#' Change the levels of the factor growth_group to have lma-L the first level
redMSU$covar$growth_group<-relevel(redMSU$covar$growth_group, ref = "lma-L")

is.factor(redMSU$covar$growth_group)
dim(redMSU$covar)
head(redMSU$covar)

if (sum(redMSU$covar$id!=pigid)!=0) stop ("pigids of redMSU$covar not correct")

#' The dataframe of covariates is complete.
#' 

#' ### Genotype Filtering
#' 
#' 1. First, extract the F2 population genotypes from the MSUPRP_meat$geno object and calculate MAFs:
dim(MSUPRP_meat$geno)

#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:
geno_f2<-MSUPRP_meat$geno[!((as.numeric(rownames(MSUPRP_meat$geno))<=1000) | (as.numeric(rownames(MSUPRP_meat$geno))>=6000)),]
dim(geno_f2)

#' Calculate allele frequency for all the F2 animals:
allele_freq<-colMeans(geno_f2,na.rm=T)/2
length(allele_freq)
summary(allele_freq)

genomat<-redMSU$geno
dim(genomat)

#' 
#' 2. Filter the subset of 174 animals' genotypes for MAF < 0.10 (being sure to take both ends of the distribution), removal of fixed SNPs, and those on sex chromosomes.
#' 
#' Filter SNPs for minor allele frequency calculated using entire F2s (maf >= 0.10 retained):
#' How many SNPs have a af < 0.10 and >0.90?
#' 
#' We filter this way since we don't 'know' the minor allele in each SNP based on allele frequency:

length(which(allele_freq<0.10)) + length(which(allele_freq>0.90))

length(names(which(allele_freq>=0.10 & allele_freq<=0.90)))

length(allele_freq) - length(names(which(allele_freq>=0.10 & allele_freq<=0.90)))

#' Retain all SNPs with af >= 0.10 and <=0.90 (maf < 0.10 discarded):
genomat <- genomat[,allele_freq>=0.10 & allele_freq <=0.90]

#' Dimensions of remaining SNP marker matrix:
dim(genomat)

#' Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):
sdv <- apply(genomat, 2, sd)
length(sdv)

if(sum(names(sdv) != colnames(genomat)) != 0) stop ("SNP names not the same between redMSU$geno and sdv")

sum(sdv == 0)

#' To make sure the proper number of SNPs were deleted, calculate the difference between the two datasets:
ncol(genomat) - sum(sdv == 0)

#' Remove fixed SNPs from genotype matrix:
genomat <- genomat[,sdv>0]

dim(genomat)

#' Eliminate markers on sex chromosomes:
table(redMSU$map$chr)
sexchr <- rownames(redMSU$map)[redMSU$map$chr == 19]
length(sexchr)

sum((colnames(genomat) %in% sexchr))

genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]

#' This object contains the markers that are not fixed, have a "maf" > 0.10, and are not mapped to the sex chromosomes:
dim(genomatfil)

if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")

ncol(MSUPRP_meat$geno)
ncol(genomatfil)

ncol(MSUPRP_meat$geno) - ncol(genomatfil)

#' Create the todel vector, containing the names of the SNPs in the gpdata geno object NOT found in the filtered genotype matrix:
#' 
todel <- colnames(redMSU$geno)[!colnames(redMSU$geno) %in% colnames(genomatfil)]
length(todel)
#' That means that, in total, we have removed 7165 SNPs from the gpdata geno object. 
#' 
#' Using discard.markers allows us to remove both the markers we don't want, and the map information we don't want, all in one step.
final_MSUPRP_miRNA <- discard.markers(redMSU, todel)
dim(final_MSUPRP_miRNA$geno)

#' The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 
#' 
summary_final_MSUPRP_miRNA<-summary(final_MSUPRP_miRNA)
summary_final_MSUPRP_miRNA$geno
summary_final_MSUPRP_miRNA$pedigree
summary_final_MSUPRP_miRNA$covar
#' 
#' ## Save data
#' 
#' Save the gpdata object with filtered genotypes, map information, and covariate information:
save(final_MSUPRP_miRNA, summary_final_MSUPRP_miRNA, file="../3_msuprp_mirna_gpdata.Rdata")