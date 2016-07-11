#' **Script:** `2_miRNA_eqtl_summary.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`
#' 
#' **Date:**  2/25/16 UPDATED 6/29/16
#' 
#' **Input File Directory:**  
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Input File(s):** 
#'  
#' 1. `1_gblup_gwa_results_summary_full.Rdata`
#' 
#' 2. `5_normalized_dge_object_and_voom_output.Rata`
#' 
#' 3. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`
#' 
#' **Output File(s):** `2_miRNA_eqtl_summary.Rdata`
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
#' The objective of this script is to summarize the number of eQTL peaks per miRNA output from the first eQTL analysis of the 174 F2 MSUPRP pig miRNA expression profiles.
#' 
#' ## Install libraries
#'

rm(list=ls())

#' Load eqtl function Rdata containing stb function, which will summarize the eQTL results
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

library(limma)
library(edgeR)
library(gwaR)
library(plyr)

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")

#' ## Load data
#' 
#' Load the eQTL output:
load("../1_gblup_gwa_results_summary_full.Rdata")

#' Load the dge object to obtain the mature miRNA annotation:
load("../../3_build_dge_object_for_eqtl/5_normalized_dge_object_and_voom_output.Rata")

#' Load the final_MSUPRP_mirna gpdata object for the mapping information:
load("../../3_build_dge_object_for_eqtl/4_msuprp_mirna_gpdata_pheno_counts.Rdata")
ls()

#' ## Analysis
#' 
#' ### Summarize the heritability of the miRNAs, output from GBLUP:
#' 
#' The average heritability of all the miRNA expression profiles:
mean(eqtl$gblup$h2)
summary(eqtl$gblup$h2)

#' The average heritability of the miRNAs found to be significantly heritable at FDR < 0.05:
#'
#' How many significantly heritable miRNAs found at FDR < 0.05? (Should match from eQTL scan output md file)
sum(eqtl$gblup$qvalue<0.05)

# #' View the gblup results ordered by p-value
# eqtl$gblup[order(eqtl$gblup$lrtpvalue),]

#' Extract the significantly heritable miRNAs from the eqtl$gblup dataset and calculate mean h2:
sigh2<-eqtl$gblup[eqtl$gblup$qvalue<0.05,]
dim(sigh2)
mean(sigh2$h2)
summary(sigh2$h2)

#' Define the minimum p-value that is not significant (based on q-value < 0.05) as the threshold for plotting significant points
eqtl$gblup$lrtpvalue[eqtl$gblup$qvalue >= 0.05]
summary(eqtl$gblup$lrtpvalue[eqtl$gblup$qvalue >= 0.05])
sigthresh<-min(eqtl$gblup$lrtpvalue[eqtl$gblup$qvalue >= 0.05])
sigthresh

#' Plot h2 vs -log10(p-value) like before to determine trend in significance and h2:
plot(eqtl$gblup$h2, -log10(eqtl$gblup$lrtpvalue),
	xlab = expression("Heritability"~(h^{2})),
    ylab = "-log10(p-value)",
    main = "-log10(p-value) vs Heritability")
points(eqtl$gblup$h2[-log10(eqtl$gblup$lrtpvalue)>(-log10(sigthresh))],
       -log10(eqtl$gblup$lrtpvalue)[-log10(eqtl$gblup$lrtpvalue)>(-log10(sigthresh))],
       pch=19,col="red")
abline(a = -log10(sigthresh), b = 0, lty = 5)


#' ---
#' 
#' ### The Summary Table of GWAS Results:
#' 
#' 
#' Assign the correct names to the different objects:
map <- final_MSUPRP_miRNA$map
colnames(map)
annotation <- dge$genes
colnames(annotation)

#' Annotation must contain columns "chr", "start", "end", and "strand"
colnames(annotation)<-c("Name","chr","start","end","width","strand","type","Alias","Precursors")
head(annotation)

# #' This function (stb) returns a summary table of the eQTL peaks per chromosome/per gene at FDR < 0.01:
# rsumtb1 <- stb(qval=eqtl$gwa.qval, map=map, annotation=annotation, Z=eqtl$gwa, threshold=0.01, gene.name="Precursors", pergene=T)
# dim(rsumtb1)

# #' miR-eQTL peaks at FDR<0.01:
# rsumtb1

#' This function (stb) returns a summary table of the eQTL peaks per chromosome/per gene at FDR < 0.05:
rsumtb5 <- stb(qval=eqtl$gwa.qval, map=map, annotation=annotation, Z=eqtl$gwa, threshold=0.05, gene.name="Precursors", pergene=T)
dim(rsumtb5)

#' miR-eQTL peaks at FDR<0.05:
rsumtb5

# #' #### Summary of GWAS results at FDR < 0.01
# #' 
# #' Number of eQTL peaks per chromosome:
# table(rsumtb1$chr.snp)

# #' Names of associated miRNAs:
# table(rsumtb1$Gene)

# #' Names of associated markers:
# table(as.character(rsumtb1$SNP))

#' #### Summary of GWAS results at FDR < 0.05
#' 
#' Number of eQTL peaks per chromosome:
table(rsumtb5$chr.snp)

#' Names of associated miRNAs:
table(rsumtb5$Gene)

#' Names of associated markers:
table(as.character(rsumtb5$SNP))

#' ### Determining the ranges of associated SNPs per eQTL peak on SSC15 (for ISAG abstract):
#' 
#' First, create the summary table at FDR 5% again, this time with pergene=F to identify all markers associated with each eQTL peak:
fullrsumtb5 <- stb(qval=eqtl$gwa.qval, map=map, annotation=annotation, Z=eqtl$gwa, threshold=0.05, gene.name="Precursors", pergene=F)
dim(fullrsumtb5)
head(fullrsumtb5)

#' Summarize the number of SNPs associated with each miRNA eQTL (some have multiple peaks)
numsnps<-by(fullrsumtb5, as.character(fullrsumtb5$Gene), nrow)
numsnps<-ldply(numsnps, fun=NULL, id=names(numsnps))
colnames(numsnps) <- c("miRNA", "numsnps")
numsnps


#' Extract from the summary table the markers associated with SSC15
ssc15fullrsumtb5<-fullrsumtb5[fullrsumtb5$chr.snp == "15",]
dim(ssc15fullrsumtb5)
ssc15fullrsumtb5

#' Then, extract the min and max of each miRNA's associated SNPs to determine the range of the eQTL peaks, and check for overlap
ssc15min<-by(ssc15fullrsumtb5$pos.snp, as.character(ssc15fullrsumtb5$Gene), min)
ssc15max<-by(ssc15fullrsumtb5$pos.snp, as.character(ssc15fullrsumtb5$Gene), max)

sum(names(ssc15min) != unique(ssc15fullrsumtb5$Gene))
sum(names(ssc15max) != unique(ssc15fullrsumtb5$Gene))

#' Extract the minimum and maximum SNP positions from the summary table, then merge the two together:
ssc15min<-ldply(ssc15min, fun=NULL, id=names(ssc15min))
colnames(ssc15min) <- c("miRNA", "minpos")
ssc15min

ssc15max<-ldply(ssc15max, fun=NULL, id=names(ssc15max))
colnames(ssc15max) <- c("miRNA", "maxpos")
ssc15max

ssc15summary<-merge(ssc15min, ssc15max, by="miRNA")
ssc15summary

#' Extract the associated number of SNPs per miRNA for SSC15 and add that to the data.frame as well:
ssc15snps<-by(ssc15fullrsumtb5, as.character(ssc15fullrsumtb5$Gene), nrow)
ssc15snps<-ldply(ssc15snps, fun=NULL, id=names(ssc15snps))
colnames(ssc15snps) <- c("miRNA", "numsnps")
ssc15snps

ssc15summary<-merge(ssc15snps, ssc15summary, by = "miRNA")

ssc15summary$range<-ssc15summary$maxpos - ssc15summary$minpos
ssc15summary

#' ### Extract peak range data from all miRNA eQTL peaks
#' 
#' Using Deborah's function "peakrng", I can obtain the same information on the range of each peak based on miRNA:
peakrng
#' nmt = name transcript (in this case, miRNA); so, make a list of the miRNAs and loop through to get the peak range information for each miRNA with significant eQTL peaks
#' sumtb = output from the summary table function, with pergene = FALSE

sigmirnames <- unique(as.character(fullrsumtb5$Gene))
sigmirnames

mirpeakranges<-lapply(sigmirnames, function(x) peakrng(nmt=x, fullrsumtb5))
mirpeakranges

#' ---
#' 
#' ### Creating Manhattan plots of the six miRNA with the highest numbers of associated SNP markers (for ISAG poster)
#' 
#' First, convert the map positions to absolute map positions using Deborah's function "absmap"
absmap

#' Notice how the map object's positions are relative to chromosome:
head(map)

#' Use the absmap function to convert the chromosomal map positions to absolute map positions:
absposmap<-absmap(map)
head(absposmap)

#' manhpt function from Deborah Velez:
manhpt

#' Create a vector of the names of the six miRNAs with the highest number of associated SNPs for the ISAG poster:
mannmsig10<-sigmirnames[c(7,11,14)]
mannmsig6<-sigmirnames[c(8,9,13)]

#' Then, provide the absposmap object to the manhpt function (also from Deborah's func_eqtl.Rdata) and loop through the vector of miRNA names to create the Manhattan plots:
#' 
#' Manhattan Plots with y-axis 0 - 10:
#' 
#+ manhattan_plotssig10, fig.width=6, fig.height=12
par(oma=c(2,2,2,2),mfrow=c(3,1))
for (i in 1:length(mannmsig10)){
	manhpt(nm=mannmsig10[[i]], abspos=absposmap, rst=eqtl$gwa.qval[,mannmsig10[[i]]], map=map, fdr=0.05, annotation=NULL, pvalues=FALSE, arrow=TRUE, pQTL=FALSE, ylim=c(0,10))
}

#' Manhattan Plots with y-axis 0 - 6
#'
#+ manhattan_plotssig6, fig.width=6, fig.height=12
par(oma=c(2,2,2,2),mfrow=c(3,1))
for (i in 1:length(mannmsig6)){
	manhpt(nm=mannmsig6[[i]], abspos=absposmap, rst=eqtl$gwa.qval[,mannmsig6[[i]]], map=map, fdr=0.05, annotation=NULL, pvalues=FALSE, arrow=TRUE, pQTL=FALSE, ylim=c(0,6))
}

#' ## Save data
# save(rsumtb1, file = "../2_eqtl_summary_table_fdr1.Rdata")
save(rsumtb5, file = "../3_eqtl_summary_table_fdr5.Rdata")
save(fullrsumtb5, file = "../4_eqtl_summary_pergeneF_fdr5.Rdata")