#' **Script:** `1_miRNA_gblup_gwas.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts/`
#' 
#' **Date:**  2/24/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Input File(s):** 
#' 
#' 1. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`
#' 
#' 2. `5_normalized_dge_object_and_voom_output.Rata`
#' 
#' 3. `6_Z_G_miRNA_gblup_gwas.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`
#' 
#' **Output File(s):** `1_gblup_gwa_results_summary_full.Rdata`
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
#' The objective of this script is to conduct the first gblup and gwa scan for the 174 MSUPRP F2 pig miRNA expression profiles.
#' 
#' This analysis will utilize the gwaR package developed by our lab group: <https://github.com/steibelj/gwaR>
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
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts/")

#' ## Load data
#' 
#' The final miRNA gpdata object with the voom-adjusted counts as phenotype data:
load("../../3_build_dge_object_for_eqtl/4_msuprp_mirna_gpdata_pheno_counts.Rdata")
#' The normalized dge object and voom centered precision weights:
load("../../3_build_dge_object_for_eqtl/5_normalized_dge_object_and_voom_output.Rata")
#' The standardized Z and G matrices:
load("../../3_build_dge_object_for_eqtl/6_Z_G_miRNA_gblup_gwas.Rdata")
ls()

#' ## Analysis
#' 
#' Create the design for the GBLUP model:
design_1<-c(~sex + growth_group)

#' Create miRnames, the list of the miRNAs to be input into the gblup analysis:
miRnames<-colnames(final_MSUPRP_miRNA$pheno)
length(miRnames)

#' ## Run gblup:
system.time({
rst.gblup<-lapply(miRnames, function(x) gblup(rsp=x, data=final_MSUPRP_miRNA, design=design_1, G=G, vdata=NULL, wt=wtcen, pos=c(T,T)))

names(rst.gblup)<-miRnames
})


#' Check standard error, if NA eliminate from analysis:
std<-do.call(rbind, lapply(rst.gblup,function(x) x$coefm[6:7,2]))
#' Check how many NAs in standard error:
sum(is.na(rowSums(std)))
#' Retain only those miRNAs with a non-NA Standard Error
miRnames<-miRnames[!is.na(rowSums(std))]
length(miRnames)

#' Perform Likelihood Ratio Test on gblup results:
system.time({
like<-lapply(rst.gblup, lrt)
names(like)<-miRnames
})

#' Multiple test corrections: FDR:
#' 
#' Obtain p-values from LRT
gblup.pval<-unlist(lapply(like, function(x) x$pvalue))

#' Calculate q-values from LRT p-values:
system.time({
	gblup.qval<-qvalue(gblup.pval,lambda=0)$qvalue
})

#' Total number of miRNAs with significant h2 at FDR < 0.05
sum(gblup.qval<0.05)

#' ## Run GWA:
#' 
#' First, transpose standardized Z matrix:
Zt<-t(Z)
dim(Zt)

system.time({
	rst.gwa<-lapply(miRnames, run.gwa, data=final_MSUPRP_miRNA, design=design_1, G=G, vdata=NULL, wt=wtcen, x=Zt, LRT=F, threshold=0.01, returnz=T, pos=c(T,T))

	names(rst.gwa)<-miRnames
	})

#' Calculate pvalues from Z-scores (Gualdron Duarte 2014)
system.time({
	gwa.pval<-lapply(rst.gwa, getpvalue, log.p=F, is.z=T)
	})

#' Multiple test correction for the GWA: FDR
system.time({
	gwa.qval<-do.call(cbind, lapply(gwa.pval, function(x) qvalue(x)$qvalues))
	})

#' Check significance per gene at FDR < 0.01
threshold <- 0.01
sig <- colSums(apply(gwa.qval, 2, "<", threshold))
length(sig[sig!=0])
sum(sig)
sig


#' Check significance per gene at FDR < 0.05:
threshold5 <- 0.05
sig5 <- colSums(apply(gwa.qval, 2 , "<", threshold5))
length(sig5[sig5!=0])
sum(sig5)
sig5

#' Matrix of gblup results:
summary.rst.gblup<-do.call(rbind, 
	lapply(miRnames, function(x) cbind(t(like[[x]]$vars[1]),
	 h2=like[[x]]$vars[1,3],
		like[[x]]$llik,
		lrtpvalue=like[[x]]$pvalue)))

rownames(summary.rst.gblup)<-miRnames

head(summary.rst.gblup)

#' Matrix of standard errors of GBLUP:
std<-std[miRnames,]
colnames(std)<-paste("stdEr", colnames(std), sep="-")

#' Matrix of GWA Zscores:
rst.gwa<-do.call(cbind,rst.gwa)

#' Matrix of pvalues from GWA:
gwa.pval<-do.call(cbind, gwa.pval)

#' Merge eQTL results into one R object:
eqtl<-list(gblup=cbind(summary.rst.gblup[,1:3], std, summary.rst.gblup[,4:7], qvalue=gblup.qval), gwa=rst.gwa, gwa.pval=gwa.pval, gwa.qval=gwa.qval)

#' ## Save data
#' 
save(eqtl, file="../1_gblup_gwa_results_summary_full.Rdata")