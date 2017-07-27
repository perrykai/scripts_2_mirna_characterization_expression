#' **Script:** `1_create_gpdata_dge_G_Z_objects.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`
#' 
#' **Date:**  2/19/16 UPDATED 5/31/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
#' 2. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL`
#' 3. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/`
#' 
#' **Input File(s):** 
#' 
#' 1. `1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata`
#'
#' 2. `gpData_PRKAG3.Rdata`
#' 
#' 3. `2_mature_mirna_annotation.Rdata`
#' 
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
#' 
#' The Objectives of this script are to: 
#' 
#' Create the gpdata object needed for the miRNA eQTL analysis.
#' 
#' Add in PRKAG3 SNPs to gpdata object (Jun 2017) 
#'
#' Filter the miRNA expression profiles by abundance across libraries and apply the voom transformation for use in the miRNA eQTL analysis.
#' 
#' Create the Z and G matrices for use in the miRNA eQTL analysis.
#' 
#' ---
#' 
#' The gpdata object will be filtered from the MSUPRP_meat gpdata object in the following ways:
#' 
#' The covariate data will be reduced to include only animals in this analysis (n=174), and will have the additional column of growth_group as a factor.  
#' 
#' The genotype data will be compiled as follows:
#' 
#' 1. First, extract the F2 population genotypes from the MSUPRP_meat$geno object and calculate AFs
#' 
#' 2. Subset the genotype data to include the 174 F2 animals in this analysis
#' 
#' 3. Filter those genotypes for MAF < 0.10 (being sure to take both ends of the distribution), removal of fixed SNPs, those located on sex chromosomes
#' 
#' 4. Remove these markers from both the SNP dataset and the map data.
#' 
#' The phenotype data in the gpdata object will consist of the filtered, normalized miRNA expression profiles
#' 
#' ## This analysis conducted using R/3.2.0, not R/3.1.0
#' 
#' ## Install libraries
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts/")

library(synbreed)
library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' ## Load data
rm(list=ls())

#' Load Yeni's function for standardizing the Z matrix:
load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")
ls()
#' Load the miRNA expression data:
load("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")
#' Load the MSUPRP_meat gpdata object:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/gpData_PRKAG3.Rdata")
#' Load the mature miRNA annotation:
load("../2_mature_mirna_annotation.Rdata")
ls()

rm(PRKAG3)

#' ## Analysis
#' 
#' no.zero.dfmeanrcround is a matrix of gene expression, with dimensions genes x samples (335 x 174)
dim(no.zero.dfmeanrcround)
head(colnames(no.zero.dfmeanrcround))
head(rownames(no.zero.dfmeanrcround))

#' Create vector of animal IDs from the column names of the expression matrix for subsetting the gpdata object:
pigid <- colnames(no.zero.dfmeanrcround)
length(pigid)
head(pigid)

#' ---
#' 
#' ### Covariate Data

#' Dimensions samples x categories (174 x categories), including growth_group as a factor:
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
#' ---
#' 
#' ### Genotype Filtering
#' 
#' First, extract the F2 population genotypes from the MSUPRP_meat$geno object and calculate MAFs:
dim(MSUPRP_meat$geno)

#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:
geno_f2<-MSUPRP_meat$geno[!((as.numeric(rownames(MSUPRP_meat$geno))<=1000) | (as.numeric(rownames(MSUPRP_meat$geno))>=6000)),]
dim(geno_f2)

#' Calculate allele frequency for all the F2 animals:
allele_freq<-colMeans(geno_f2,na.rm=T)/2
length(allele_freq)
summary(allele_freq)

#' Filter the subset of 174 animals' genotypes for MAF < 0.10 (being sure to take both ends of the distribution), removal of fixed SNPs, and those on sex chromosomes.
#' 
#' Extract the genotype object from the reduced gpdata object (the 174 animals and all the SNP markers):
genomat<-redMSU$geno
dim(genomat)

#' Filter SNPs for minor allele frequency calculated using entire F2s (maf >= 0.10 retained):
#' How many SNPs have a af < 0.10 and >0.90?
#' 
#' We filter this way since we don't know which allele is "minor" in each SNP based on allele frequency (if A or B is the minor allele):

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
#' 
#' JUNE 13, 2017: Error in synbreed package code; it saves the map object as a list when using "add.Markers"
#' function, which was used by DV to add the PRKAG3 markers.
#' 
#' So, need to extract the map object, convert to data.frame, and replace in the MSUPRP_meat gpdata object, prior to filtering out SNPs.
map <- data.frame(chr=redMSU$map[[1]],pos=redMSU$map[[2]])
rownames(map) <- colnames(redMSU$geno)
head(map)
redMSU$map <- map
str(redMSU$map)

#' Check to make sure the PRKAG3 markers successfully added to the map object
grep("PRKAG3", rownames(map))
map[grep("PRKAG3", rownames(map)),]

#' Check for NA's in the genotypes:
table(redMSU$geno[,grep("PRKAG3", colnames(redMSU$geno))])

table(redMSU$map$chr)

#' Identify SNPs on sex chr (chr == 19)
sexchr <- rownames(redMSU$map)[redMSU$map$chr == 19]
length(sexchr)

#' How many SNPs in my dataset are on the sex chr:
sum((colnames(genomat) %in% sexchr))

#' Filter out the SNPs on the sex chromosome
genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]

#' This object contains the markers that are not fixed, have a maf > 0.10, and are not mapped to the sex chromosomes:
dim(genomatfil)

if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")

ncol(MSUPRP_meat$geno)
ncol(genomatfil)

#' Check on how many SNPs are being deleted:
ncol(MSUPRP_meat$geno) - ncol(genomatfil)

#' Create the todel vector, containing the names of the SNPs in the gpdata geno object NOT found in the filtered genotype matrix:
#' 
todel <- colnames(redMSU$geno)[!colnames(redMSU$geno) %in% colnames(genomatfil)]
length(todel)
#' That means that, in total, we have removed 7165 SNPs from the gpdata geno object. 
#' 
#' Using discard.markers allows us to remove both the markers we don't want, and the map information we don't want, all in one step.
MSUPRP_miRNA <- discard.markers(redMSU, todel)
dim(MSUPRP_miRNA$geno)

#' The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 
#' 
summary_MSUPRP_miRNA<-summary(MSUPRP_miRNA)
summary_MSUPRP_miRNA$geno
summary_MSUPRP_miRNA$pedigree
summary_MSUPRP_miRNA$covar
str(summary_MSUPRP_miRNA$pheno)

#' ---
#' 
#' ### Phenotype Data (miRNA Expression Profile Preparation) 
#' 
#' Create the dge object and transform the read counts: log-cpm, then filter genes by expression (1 cpm in 44 or more samples retained):
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
MSUPRP_miRNA$pheno <- array(t(v$E), dim=c(ncol(v$E),nrow(v$E),1))
rownames(MSUPRP_miRNA$pheno)<- colnames(v$E)
colnames(MSUPRP_miRNA$pheno)<- rownames(v$E)
dim(MSUPRP_miRNA$pheno)
MSUPRP_miRNA$pheno[1:5,1:5,1]

if (sum(MSUPRP_miRNA$pheno[,,1] != t(v$E)) != 0) stop ("pheno object not equal to voom counts")
if (sum(rownames(MSUPRP_miRNA$pheno) != rownames(MSUPRP_miRNA$geno)) != 0) stop ("pheno object rownames not equal to geno object rownames")
if (sum(rownames(MSUPRP_miRNA$pheno) != MSUPRP_miRNA$covar$id) != 0) stop ("pheno object rownames not equal to covar id column")
sum(rownames(no.zero.dfmeanrcround) %in% colnames(MSUPRP_miRNA$pheno))

#' Create the correct gpdata summary for the updated final_MSUPRP_miRNA:
summary_MSUPRP_miRNA<-summary(MSUPRP_miRNA)
summary_MSUPRP_miRNA$pheno[,1:10]

#' This object is ready to be saved, and used in GBLUP and GWA analyses.
#' 
#' ---
#' 
#' ### Standardize the matrix of SNPs and create G matrix:
#' 
#' Using zstandard, the function created by Yeni Bernal Rubio:
#' 
#' This Z matrix is constructed using the allele frequencies of the entire F2 population and extracting the markers present in this analysis:
#'
#' Need to filter the allele_freq object to retain the correct SNP markers.
length(allele_freq)
head(colnames(MSUPRP_miRNA$geno))
allele_freq<-allele_freq[colnames(MSUPRP_miRNA$geno)]
length(allele_freq)

if(sum(names(allele_freq) != colnames(MSUPRP_miRNA$geno)) !=0) stop ("allele freq don't match SNP matrix")

#' Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:
Z<-zstandard(MSUPRP_miRNA$geno, alfreq=allele_freq, procedure="heterogeneous")
dim(Z)

G<-Z%*%t(Z)
dim(G)
summary(diag(G))

IQR(diag(G))

#' 
#' ## Save data
#' 
#' Save full gpdata object with filtered genotypes, map information, and covariate information, and updated pheno data:
save(MSUPRP_miRNA, summary_MSUPRP_miRNA, file="../3_msuprp_mirna_gpdata_pheno_counts.Rdata")
#' Save dge object and voom output with weights:
save(dge, v, wtcen, file = "../4_normalized_dge_object_and_voom_output.Rdata")
#' Save the Z and G matrices:
save(Z, G, file = "../5_Z_G_miRNA_gblup_gwas.Rdata")