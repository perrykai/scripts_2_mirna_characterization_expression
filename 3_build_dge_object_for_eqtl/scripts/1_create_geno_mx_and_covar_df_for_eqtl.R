#' **Script:** `1_create_metadata_object.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`
#' 
#' **Date:**  2/18/16
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
#' 3. `ssc.gff3`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. `3_covardata_for_eQTL_analysis.Rdata`
#' 2. `4_filtered_genotypes_for_G_matrix.Rdata`
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
#' The objective of this script is to create the other components of data needed for the miRNA eQTL analysis.
#' These data objects are needed for this analysis:
#' 
#' 1. A matrix of gene expression, with dimensions genes x samples (335 x 174)
#' 
#' 2. A data frame of metadata (trait data) with dimensions samples x categories (174 x categories)
#' 
#' 3. A matrix of genotype data, filtered for removing fixed SNPs, those SNPs with low maf, and those SNPs on sex chromosomes. 
#' 
#' ## This analysis conducted using R/3.2.0, not R/3.1.0
#' 
#' ## Install libraries

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
#' Load the gp data object:
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

#' ### 2. A data frame of metadata (trait data) with dimensions samples x categories (174 x categories)
#' 
#' Extract data frame of covariates (trait data) from GPData object:
covaMSU<-MSUPRP_meat$covar
rownames(covaMSU)<-covaMSU$id
colnames(covaMSU)
head(covaMSU)
dim(covaMSU)

if (sum(rownames(covaMSU) != MSUPRP_meat$covar$id) != 0) stop ("rownames not the same as MSUPRP_meat id column")

#' Subset the covaMSU dataframe, retaining only those animals in pigid:
redcovaMSU<-covaMSU[(rownames(covaMSU) %in% pigid),]
dim(redcovaMSU)
head(redcovaMSU)

if (sum(rownames(redcovaMSU) != pigid) != 0) stop ("rows are not correct")
if (sum(rownames(redcovaMSU) != redcovaMSU$id) != 0) stop ("rownames do not match id column")

#' Create the growth_group column as a factor by combining Status and selcrit
growth_group<-paste(redcovaMSU$selcrit, redcovaMSU$Status, sep = "-")
length(growth_group)
head(growth_group)

redcovaMSU$growth_group<-as.factor(growth_group)

#' Change the levels of the factor growth_group to have lma-L the first level
redcovaMSU$growth_group<-relevel(redcovaMSU$growth_group, ref = "lma-L")

is.factor(redcovaMSU$growth_group)

head(redcovaMSU)
tail(redcovaMSU)

dim(redcovaMSU)

#' The matrix of covariates is complete. This will be utilized in the GBLUP when combined with the matrix of normalized read counts.
#' 
#' ### 3. A matrix of genotype data, filtered for removing fixed SNPs, SNPs with low maf, and SNPs on sex chromosomes. 
#' 
#' Discard the animals from the GPData Geno object not present in the miRNA expression matrix:
todisc <- rownames(MSUPRP_meat$geno)[!rownames(MSUPRP_meat$geno) %in% pigid]

redgenoMSU <- discard.individuals(MSUPRP_meat, todisc)

#' Extract the genotypes from the GPData object:
genomat <- redgenoMSU$geno

head(genomat[,1:10])

dim(genomat)

if (sum(rownames(genomat)!=rownames(redcovaMSU)) != 0) stop ("rows of geno not the same as rows of pheno")

#' Filter the matrix of genotypes: 
#' 
#' Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):
sdv <- apply(genomat, 2, sd)
length(sdv)
sum(sdv == 0)

#' Remove fixed SNPs from genotype matrix:
genomat <- genomat[,sdv>0]

dim(genomat)

#' Filter SNPs for minor allele frequency (maf > 0.1 retained):
maf <- colMeans(genomat)/2

head(maf)

#' How many SNPs have a maf <= 0.1?
sum(maf<=0.1)

#' Retain all SNPs with maf > 0.1:
genomat <- genomat[,maf>0.1]

#' Dimensions of remaining SNP marker matrix:
dim(genomat)

if (sum((colMeans(genomat)/2)<=0.1) !=0) stop ("maf filtering did not work correctly")

if (sum(rownames(genomat)!=rownames(redcovaMSU))!=0) stop ("rownames of marker matrix not equal to trait data")

#' Eliminate markers on sex chromosomes:
sexchr <- rownames(redgenoMSU$map)[redgenoMSU$map$chr == 19]
length(sexchr)

sum((colnames(genomat) %in% sexchr))

genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]

dim(genomatfil)

if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(rownames(redcovaMSU) != rownames(genomatfil)) != 0) stop ("rownames of trait data and genotype matrix not the same")
if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")

#' The filtered matrix of SNPs is complete. This will be used in the GBLUP and GWAS as the G matrix (genomic relationship matrix) after standardizing the SNPs. 

#' 
#' ## Save data
#' 
#' Save the trait data (covariates):
save(redcovaMSU,file = "../3_covardata_for_eQTL_analysis.Rdata")

#' Save the filtered genotype matrix:
save(genomatfil,file = "../4_filtered_genotypes_for_G_matrix.Rdata")
