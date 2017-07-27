#' The objective of this script is to compare multiple methods of calculating the G genomic relationship matrix using the function created by Yeni Bernal Rubio.
#' 
#' There are multiple options for this calculation:
#' 
#' 1. Using either the entire MSUPRP genotype markers for calculating the expected allele frequencies, or only the F2 generation of genotypes
#' 
#' 2. Using either the "heterogeneous" or "homogeneous" protocols when standardizing the Z matrix: difference here is the contents of the denominator,
#' 
#' 	ie how the markers are standardized. "Heterogeneous" makes G analogous to the numerator relationship matrix A. 
#' 
#' Additionally, testing of the filtering done on the cpm-transformed read counts prior to voom transformation will be accomplished in this script.

library(synbreed)
library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)

rm(list=ls())

#' The mature miRNA annotation:
load("../2_mature_mirna_annotation.Rdata")
#' The rounded, filtered mean read counts:
load("../../2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")
#' Load the gpdata object for the miRNA data:
load("../3_msuprp_mirna_gpdata.Rdata")
#' Load the MSUPRP_meat gpdata object:
load("/mnt/research/pigsnp/MSUPRP/carcass_quality/data/MSUPRP_meat.RData")
#' Load Yeni's function for standardizing the Z matrix:
load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")

#' ### 2. Standardize the matrix of SNPs and create G matrix:
#' 
#' Standardize the FULL matrix of SNP markers (all genotyped animals in the MSUPRP), using Yeni's function:
#' 
#' 
#' USE ZSTANDARD FUNCTION WITH ALLELE FREQ FROM ONLY F2 POPULATION:
geno<-MSUPRP_meat$geno
Mf2<-final_MSUPRP_miRNA$geno
dim(geno)
# [1]  1015 45329

#' Filter just the F2
#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1
#' Thus take those away 
geno_f2<-geno[!((as.numeric(rownames(geno))<=1000) | (as.numeric(rownames(geno))>=6000)),]
dim(geno_f2)
# [1]   940 45329

#' Calculate allele frequency for the F2s:
all_freq<-colMeans(geno_f2,na.rm=T)/2
length(all_freq)
# [1] 45329

#' Subset my markers:
all_freq<-all_freq[colnames(Mf2)]
length(all_freq)
# [1] 38319

Zf2hetero<-zstandard(final_MSUPRP_miRNA$geno, alfreq=all_freq, procedure="heterogeneous")
dim(Zf2hetero)
# [1]   174 38319

Gf2hetero<-Zf2hetero%*%t(Zf2hetero)
summary(diag(Gf2hetero))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.6311  0.9076  0.9508  0.9473  0.9882  1.0850
IQR(diag(Gf2hetero))
# [1] 0.08059075

#' FINAL CHECK: HOMOGENEOUS WITH ONLY THE 174 ANIMALS:
#' First, extract the 174 F2 animals from MSUPRP_meat$geno to calculate allele frequencies:
pigid<-rownames(final_MSUPRP_miRNA$geno)

todisc <- rownames(MSUPRP_meat$geno)[!rownames(MSUPRP_meat$geno) %in% pigid]
length(todisc)

redMSU <- discard.individuals(MSUPRP_meat, todisc)
dim(redMSU$geno)
sum(rownames(redMSU$geno) != rownames(final_MSUPRP_miRNA$geno))

#Calculate allele frequency for all markers in the 174 animals:
miRNA_freq<-colMeans(redMSU$geno)/2

#Then subset the markers in my dataset:
miRNA_freq<-miRNA_freq[colnames(final_MSUPRP_miRNA$geno)]
length(miRNA_freq)

ZmiRNAhomo<-zstandard(final_MSUPRP_miRNA$geno, alfreq=miRNA_freq, procedure="homogeneous")

GmiRNAhomo<- ZmiRNAhomo %*% t(ZmiRNAhomo)
summary(diag(GmiRNAhomo))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.6036  0.9002  0.9359  0.9337  0.9702  1.0640

plot(GmiRNAhomo, Gf2hetero)
abline(0,1,col="red")



#==========================
#' Also testing different thresholds of filtering for the cpm.dge step in the transformation of the read counts. 
#' There are some interesting expression profiles removed when I filter at 1cpm > 87 samples 
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

rowSums(cpm.dge)
length(rowSums(cpm.dge))

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (174/4=43.5-->44)
filtercpm<-rowSums(cpm.dge>=1)>=44
sum(filtercpm)

filtercpm2<-rowSums(cpm.dge>=1)>=87 #(1/2 of samples)
sum(filtercpm2)

filtercpm3<-rowSums(cpm.dge>=1)>=58 #(1/3 of samples)
sum(filtercpm3)

#' Instead of filtering the dge object, I'll filter the cpm.dge object to look at the cpm of the miRNAs we're eliminating.
removed44cpm.dge<-cpm.dge[!filtercpm,]
dim(removed44cpm.dge)
#' So, we're removing 40 miRNAs from the analysis at this threshold.
removed44cpm.dge[1:5,1:5]

#' This provides the number of non-zero cpm in the removed44cpm.dge object: (How many of the samples have a non-zero read count in the miRNAs we're removing?)
rowSums(removed44cpm.dge!=0)

#' This allows me to look at the cpm of the miRNAs we're removing from the analysis to decide if we're being too stringent with filtering:
apply(removed44cpm.dge,1,table)

plot(density(removed44cpm.dge[1,]),ylim=c(0,25),xlim=c(-1,10))
for (i in 2:nrow(removed44cpm.dge)) {
	lines(density(removed44cpm.dge[i,]))
}


#' Now do the same tests with a more stringent threshold: Are we removing any potentially interesting miRNAs by filtering at at least 1cpm for 87 samples?
removed87cpm.dge<-cpm.dge[!filtercpm2,]
#' So, we're removing 51 miRNAs from the analysis at this threshold.
dim(removed87cpm.dge)
removed87cpm.dge[1:5,1:5]
rowSums(removed87cpm.dge!=0)
apply(removed87cpm.dge,1,table)

removed58cpm.dge<-cpm.dge[!filtercpm3,]
dim(removed58cpm.dge)
#' So, we're removing 45 miRNAs from the analysis at this threshold.
rowSums(removed58cpm.dge)
apply(removed58cpm.dge,1,table)
