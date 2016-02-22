#' **Script: ** `2_create_dge_g_objects_for_eqtl.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`
#' 
#' **Date:**  2/18/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_msuprp_mirna_gpdata.Rdata`
#' 
#' 2. `2_mature_mirna_annotation.Rdata`
#' 
#' 3. `2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata`
#' 
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Output File(s):** ``
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
#' The objective of this script is to create the dge object by combining the filtered, rounded read counts and the annotation file for the mature miRNA expression profiles.
#' 
#' This object will then be used in transforming the read counts to log cpm, and filtering miRNAs expressed at 1 cpm in less than 1/4 the samples (If 1cpm in >= 44 libraries, retain miRNA). 
#' The normalized counts from this dge object will be used in the eQTL analysis as the expression/phenotypic data.
#' 
#' Additionally, the matrix of SNP markers will be standardized in this script, and the standardized matrix (Z) will be multiplied by its transpose (Z') to form the G matrix. 
#' The G matrix will be used in the eQTL analysis, as the Genomic Relationship Matrix. 
#' 
#' Finally, the data frame containing covariate data will be joined with the expression data (voom output), needed for the GBLUP and GWA:
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
#' ### 1. Create the dge object and transform the read counts: log-cpm, then filter genes by expression (1 cpm in 87 or more samples retained):
#' 
#' Create the dge object:
dge<-DGEList(counts=no.zero.dfmeanrcround,genes=total.mature.annot2)
dim(dge)
dge[1:5,1:5]

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.dge<-cpm(dge)
dim(cpm.dge)
cpm.dge[1:5,1:5]

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (174/4=43.5-->44)
filtercpm<-rowSums(cpm.dge>=1)>=44
sum(filtercpm)

removeddge<-dge[!filtercpm,]
names(removeddge)
removeddge[1:5,1:5]
rowSums(removeddge$counts)
apply(removeddge$counts,1,table)
dim(removeddge$counts)

dge<-dge[filtercpm,]
names(dge)
dge[1:5,1:5]
dim(dge$counts)

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
voomwt<-t(v$weights)
dim(voomwt)
voomwt[1:5,1:5]

voomwtsq<-1/sqrt(voomwt)
voomwtcen<-as.matrix(sweep(voomwtsq,2,FUN="/",STATS=colMeans(1/sqrt(voomwt)))) #Need to ask Deborah/Juan what this step does
rownames(voomwtcen)<-rownames(t(v$E))
colnames(voomwtcen)<-colnames(t(v$E))
dim(voomwtcen)
voomwtcen[1:5,1:5]

#' Now, replace the gpdata "pheno" object with the voom-adjusted read counts:
#' 
#' The pheno object must be a data frame with rows of animals and columns of miRNAs.
#' So, the voom-transformed read count matrix must be transposed to comply, and then be made into an array of 3 dimensions:
final_MSUPRP_miRNA$pheno <- array(t(v$E), dim=c(ncol(v$E),nrow(v$E),1))
rownames(final_MSUPRP_miRNA$pheno)<- colnames(v$E)
colnames(final_MSUPRP_miRNA$pheno)<- rownames(v$E)

final_MSUPRP_miRNA$pheno[1:5,1:5,1]

if (sum(final_MSUPRP_miRNA$pheno[,,1] != t(v$E)) != 0) stop ("pheno object not equal to voom counts")
if (sum(rownames(final_MSUPRP_miRNA$pheno) != rownames(final_MSUPRP_miRNA$geno)) != 0) stop ("pheno object rownames not equal to geno object rownames")
if (sum(rownames(final_MSUPRP_miRNA$pheno) != final_MSUPRP_miRNA$covar$id) != 0) stop ("pheno object rownames not equal to covar id column")
sum(rownames(no.zero.dfmeanrcround) %in% colnames(final_MSUPRP_miRNA$pheno))





#' ### 2. Standardize the matrix of SNPs and create G matrix:
#' 
#' Standardize the FULL matrix of SNP markers (all genotyped animals in the MSUPRP), using Yeni's function:
#' 
#' Load Yeni's function for standardizing the Z matrix:
load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")

#Calculate allele frequency for the full population, and subset my markers:
pg<-colMeans(MSUPRP_meat$geno)/2
Mf2<-final_MSUPRP_miRNA$geno
pg<-pg[colnames(Mf2)]
length(pg) == ncol(final_MSUPRP_miRNA$geno)

Zhetero<-zstandard(final_MSUPRP_miRNA$geno, alfreq=pg, procedure="heterogeneous")
dim(Zhetero)
Zhetero[1:5,1:5]

#' Calculate the genomic relationship matrix (G) by multiplying Z * Z':
Ghetero<-Zhetero %*% t(Zhetero)
dim(Ghetero)
Ghetero[1:5,1:5]

summary(diag(Ghetero))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.6271  0.9067  0.9472  0.9432  0.9820  1.0710
IQR(diag(Ghetero))
# [1] 0.07528552




#' Create the design for the GBLUP model:
design <- c(~sex + growth_group)

##rst.gblup<-mcapply(rsp, function(x) lrt(gblup(rsp=x)))
# CHECK THIS WITH DEBORAH- MAKE SURE GWAS WILL INCLUDE ALL MIRNAS, NOT ONLY SIGNIFICANTLY HERITABLE ONES (SIG WITH LRT)


#' ## Save data
save(final_MSUPRP_miRNA, file = "../3_msuprp_mirna_gpdata.Rdata")
save(dge, file = "../4_normalized_dge_object_before_voom.Rata")
save(Z, G, file = "../5_Z_G_for_miRNA_eQTL.Rdata")
