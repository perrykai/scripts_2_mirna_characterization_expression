#' **Script: ** `2_create_dge_g_objects_for_eqtl.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`
#' 
#' **Date:**  2/23/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
#' 
#' 3. `/mnt/research/pigsnp/MSUPRP/carcass_quality/data/MSUPRP_meat.RData`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_msuprp_mirna_gpdata.Rdata`
#' 
#' 2. `2_mature_mirna_annotation.Rdata`
#' 
#' 3. `2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata`
#' 
#' 4. `MSUPRP_meat.RData`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`
#' 
#' 2. `5_normalized_dge_object_and_voom_output.Rata`
#' 
#' 3. `6_Z_G_miRNA_gblup_gwas.Rdata`
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
#' The objective of this script is to create the dge object by combining the filtered, rounded read counts and the annotation file for the mature miRNA expression profiles.
#' 
#' This object will then be used in transforming the read counts to log cpm, and filtering miRNAs expressed at 1 cpm in less than 1/4 the samples (If 1cpm in >= 44 libraries, retain miRNA). 
#' The normalized counts from this dge object will be used in the eQTL analysis as the expression/phenotypic data, added to the gpdata pheno object.
#' 
#' Additionally, the matrix of SNP markers will be standardized in this script, and the standardized matrix (Z) will be multiplied by its transpose (Z') to form the G matrix. 
#' The G matrix will be used in the eQTL analysis; the Genomic Relationship Matrix. 
#' 
#' # This analysis conducted with R/3.2.0
#' 

#' ## Install libraries
library(synbreed)
library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)

rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts/")

#' Load Yeni's function for standardizing the Z matrix:
load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")
ls()

#' ## Load data
#' 
#' The mature miRNA annotation:
load("../2_mature_mirna_annotation.Rdata")
#' The rounded, filtered mean read counts:
load("../../2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")
#' Load the gpdata object created for the filtered markers:
load("../3_msuprp_mirna_gpdata.Rdata")
#' Load the MSUPRP_meat gpdata object for use in standardizing the SNP genotypes:
load("/mnt/research/pigsnp/MSUPRP/carcass_quality/data/MSUPRP_meat.RData")
ls()

#' ## Analysis
#' 
#' ### 1. Create the dge object and transform the read counts: log-cpm, then filter genes by expression (1 cpm in 44 or more samples retained):
#' 
#' Create the dge object:
dge<-DGEList(counts=no.zero.dfmeanrcround,genes=total.mature.annot2)
dim(dge)
dge[1:5,1:5]

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.dge<-cpm(dge)
dim(cpm.dge)
cpm.dge[1:5,1:5]

if (sum(rownames(no.zero.dfmeanrcround)!=rownames(cpm.dge))!=0) stop ("miRNAs not the same between read counts and cpm")
if (sum(colnames(no.zero.dfmeanrcround)!=colnames(cpm.dge))!=0) stop ("animal ids not the same between read counts and cpm")

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (174/4=43.5-->44)
filtercpm<-rowSums(cpm.dge>=1)>=44
sum(filtercpm)

#' We are removing 40 miRNA profiles from the analysis
#' 
#' So, keep the miRNA profiles in dge based on those retained in the cpm-filtering step:
#' 
#' This retains the rounded, filtered mean read counts, not the cpm.
dge<-dge[filtercpm,]
names(dge)
dge[1:5,1:5]
dim(dge$counts)

if (sum(colnames(dge)!=colnames(cpm.dge))!=0) stop ("colnames not the same between dge and cpm.dge")

#' Apply the TMM normalization:
dge<-calcNormFactors(dge)
head(dge$samples)
hist(dge$samples$norm.factors)

dge<-estimateCommonDisp(dge,verbose=TRUE)
dge$common.dispersion

#' This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.
#' Save this object before voom transformation so we can access the pseudo.counts if needed.

#' Run voom transformation:
v<-voom(dge,plot=TRUE)
names(v)
dim(v$weights)
dim(v$genes)
v$genes[1:5,1:5]
v$weights[1:5,1:5]
v$targets[1:5,]

#' Extract the voom precision weights and transform them for use in the GBLUP error term:
wt<-t(v$weights)
dim(wt)
wt[1:5,1:5]

#' Standardize the voom weights: [(1/sqrt(wt))/mean(1/sqrt(wt))]
wtsq<-1/sqrt(wt)
wtcen<-as.matrix(sweep(wtsq,2,FUN="/",STATS=colMeans(1/sqrt(wt))))
rownames(wtcen)<-rownames(t(v$E))
colnames(wtcen)<-colnames(t(v$E))
dim(wtcen)
wtcen[1:5,1:5]

#' Now, replace the gpdata "pheno" object with the voom-adjusted read counts:
#' 
#' The pheno object must be a data frame with rows of animals and columns of miRNAs.
#' So, the voom-transformed read count matrix (E) must be transposed to comply, and then be made into an array of 3 dimensions:
final_MSUPRP_miRNA$pheno <- array(t(v$E), dim=c(ncol(v$E),nrow(v$E),1))
rownames(final_MSUPRP_miRNA$pheno)<- colnames(v$E)
colnames(final_MSUPRP_miRNA$pheno)<- rownames(v$E)
dim(final_MSUPRP_miRNA$pheno)
final_MSUPRP_miRNA$pheno[1:5,1:5,1]

if (sum(final_MSUPRP_miRNA$pheno[,,1] != t(v$E)) != 0) stop ("pheno object not equal to voom counts")
if (sum(rownames(final_MSUPRP_miRNA$pheno) != rownames(final_MSUPRP_miRNA$geno)) != 0) stop ("pheno object rownames not equal to geno object rownames")
if (sum(rownames(final_MSUPRP_miRNA$pheno) != final_MSUPRP_miRNA$covar$id) != 0) stop ("pheno object rownames not equal to covar id column")
sum(rownames(no.zero.dfmeanrcround) %in% colnames(final_MSUPRP_miRNA$pheno))

#' Create the correct gpdata summary for the updated final_MSUPRP_miRNA:
summary_final_MSUPRP_miRNA<-summary(final_MSUPRP_miRNA)
summary_final_MSUPRP_miRNA

#' This object is ready to be saved, and used in GBLUP and GWA analyses.
#' 
#' ---

#' ### 2. Standardize the matrix of SNPs and create G matrix:
#' 
#' Using zstandard, the function created by Yeni Bernal Rubio:
#' 
#' This Z matrix is constructed using the allele frequencies of the entire F2 population and extracting the markers present in this analysis:
#'
#' First, extract the F2 population genotypes from the MSUPRP_meat$geno object:
geno<-MSUPRP_meat$geno
dim(geno)

#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:
geno_f2<-geno[!((as.numeric(rownames(geno))<=1000) | (as.numeric(rownames(geno))>=6000)),]
dim(geno_f2)

#' Calculate allele frequency for the F2s:
all_freq<-colMeans(geno_f2,na.rm=T)/2
length(all_freq)
summary(all_freq)

#' Subset the allele frequencies for the markers in this analysis (filtered for fixed SNPs, maf<0.10, and SNPs on sex chromosomes):
all_freq<-all_freq[colnames(final_MSUPRP_miRNA$geno)]
length(all_freq)

if (sum(length(all_freq)!=ncol(final_MSUPRP_miRNA$geno))!= 0) stop ("allele freq not the same length as number of markers")

#' Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:
Z<-zstandard(final_MSUPRP_miRNA$geno, alfreq=all_freq, procedure="heterogeneous")
dim(Z)

G<-Z%*%t(Z)
summary(diag(G))

IQR(diag(G))

#' ## Save data
#'
#' Save full gpdata object with updated pheno data:
save(final_MSUPRP_miRNA, summary_final_MSUPRP_miRNA, file="../4_msuprp_mirna_gpdata_pheno_counts.Rdata")
#' Save dge object and voom output with weights:
save(dge, v, wtcen, file = "../5_normalized_dge_object_and_voom_output.Rata")
#' Save the Z and G matrices:
save(Z, G, file = "../6_Z_G_miRNA_gblup_gwas.Rdata")
