#' **Script: ** `3_miRNA_eQTL_characterization.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`
#' 
#' **Date:**  `06/14/17`
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/eQTL/paper/output/`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_eqtl_summary_tables.Rdata`
#' 
#' 2. `MSUPRP_miRNA.Rdata`
#' 
#' 3. `funct_eqtl.Rdata`
#' 
#' 4. ``
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis`
#' 
#' **Output File(s):** 
#' 
#' 1. `4_miReQTL_local-distant-regulators.Rdata`
#' 
#' 2. `5_miReQTL_peaks-equal-chr.Rdata`
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
#' The objective of this script is to characterize the miRNA eQTL results. This includes:
#'
#' Using the summary file from the miRNA eQTL analysis determine the number of miRNA eQTL peaks
#' that have local regulation and distant regulation.
#' 
#' Looking only at the eQTL peaks on the same chromosome or overlapping with the mapped gene position check how many markers are before and after the gene position.
#' 
#' Create a data frame containing the results of the eQTL scan along with a column
#' specifying the eQTL regulator type, local (cis), distant (tran)
#' 
#' This analysis will be completed using functions created by DV in directory /mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z`
#' 
#' ## Install libraries

library(synbreed)
library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' ## Load data
#'
rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")

#' Load DV's functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

#' Load the gwa results:
load("../2_gwa_results_summary.Rdata")

#' Load the miRNA eQTLsummary tables:
load("../3_eqtl_summary_tables_maps.Rdata")

#' Load the dge object to obtain the mature miRNA annotation:
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")

#' Load the MSUPRP_miRNA gpdata object for the mapping information:
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")

#' Load the miRNA annotation files to obtain info for miRNA eQTL with multiple precursors:
load("../../3_build_dge_object_for_eqtl/1_precursor_mirna_annotation.Rdata")
load("../../3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")

#' Load the MSUPRP_meat gpdata object:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/gpData_PRKAG3.Rdata")
rm(PRKAG3)
ls()

#' ## Analysis

#' eQTL peaks
head(sum.eqtl)
str(sum.eqtl)

#' ---
#' 
#' Separate the first column for use in the Target Prediction analysis, substituting the "ssc" for "hsa":
hsamir<-data.frame(V1=unique(gsub("ssc-", "", sum.eqtl$miRNA)))
hsamir

#' Save that table for use in target prediction analysis:
write.table(hsamir, file="../../6_miRNA_eQT_target_prediction/1_mirna_names.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#' ---
#' 
#' eQTL peaks with unknown position
as.character(sum.eqtl$miRNA[is.na(sum.eqtl$chr.miR)])

#' Hotspots
tail(sort(table(as.character(fullsum.eqtl$SNP))))
mk <- table(as.character(fullsum.eqtl$SNP))

table(mk)
#' Names of hotspot SNP (SNP associated with at least 4 miRNAs)
nm <- names(mk[mk >= 4])
nm

#' Check map positions of hotspot SNP:
MSUPRP_miRNA$map[nm,]
#' All my hotspot SNP are on chromosome 15; LD?
abspos<-absposmap
abspos[nm]


#' Want to extract the miRNA information for each potential hotspot SNP
htsp.mir<-list()
for(i in nm){
htsp.mir[[i]]<-fullsum.eqtl[grep(i, as.character(fullsum.eqtl$SNP)), c("miRNA", "chr.miR", "precursors")]
}

htsp.mir

#' Build a matrix of qvalues for all the miRNAs.
#' 
#' This will be used in plotting the eQTL map and determining local vs distal-acting miRNA eQTL
sigmirqval<-do.call(cbind, lapply(unique(as.character(sum.eqtl$miRNA)), function(x) rst.gwa[rst.gwa$miRNA==x, "gwa.qval"]))

colnames(sigmirqval)<-unique(as.character(sum.eqtl$miRNA))
rownames(sigmirqval)<-rownames(MSUPRP_miRNA$map)
dim(sigmirqval)
head(sigmirqval)
str(sigmirqval)

#' Check for correctness in building this matrix:
names(rst.gwa)
head(rst.gwa)
#' Rownames didn't change between datasets:
if (sum(unique(rst.gwa$SNPid)!=rownames(sigmirqval)) !=0) stop("Rownames changed between datasets")
#' Print the sum of the q-values not equal between the gwa results and the subset matrix:
for (i in colnames(sigmirqval)){
	print(sum(rst.gwa[rst.gwa$miRNA==i,"gwa.qval"] != sigmirqval[,i]))
}
#' Looks good.

#' Extract the q-values for 17 significant miRNA associations with hotspot SNPs:
htsp.qv <- sigmirqval[nm,]
dim(htsp.qv)
head(htsp.qv)


#' ---
#'
#' Table of positions for significant miRNA-marker associations within a hotspot using gene-wise FDR
#' 
#' tbpos = table of positions; takes the qvalues and absolute positions of hotspot SNPs & miRNAs (htsp.qv), 
#' and puts them into a table for use in the eQTL map function later on.
threshold <- 0.05
#' ttp.htsp = table to plot. hotspot
ttp.htsp <- tbpos(qval=htsp.qv, abspos=abspos, threshold=threshold)
dim(ttp.htsp)
ttp.htsp

#' Table of positions for all gene-marker associations and for significant eQTL SNPs
#' 
#' Notice, it removed the miRNA that doesn't have assembly information, miR-140-5p.
#' 
#' It can't add it to the eQTL map if it doesn't have mapping information. 
ttp.all <- tbpos(qval=sigmirqval, abspos=abspos, threshold=threshold)
dim(ttp.all)
head(ttp.all)

#' Also lose 4 SNPs, which are the ones significantly associated to miR-140-5p only:
table(sigmirqval[,"ssc-miR-140-5p"]<threshold)
names(which(sigmirqval[,"ssc-miR-140-5p"]<threshold))
table(ttp.all$SNP %in% names(which(sigmirqval[,"ssc-miR-140-5p"]<threshold)))

#' Build a data.frame of the miRNA eQTL peak positions and their lengths
ttp.peaks <- data.frame(miRNA=sum.eqtl[,"miRNA"], SNP=sum.eqtl[,"SNP"],
	pos.snp=abspos[as.character(sum.eqtl$SNP)],
	pos.miR=abspos[as.character(sum.eqtl$miRNA)],
	diff=abs(abspos[as.character(sum.eqtl$miRNA)] - abspos[as.character(sum.eqtl$SNP)]),
	qvalue=sum.eqtl[,"qvalue"])
head(ttp.peaks)


#' Local eQTL (na.omit to remove the miRNA with no assembly data)
local <- na.omit(mirpeaks[mirpeaks$chr.miR == mirpeaks$chr.snp,])
dim(local)
local

#' Distant eQTL (na.omit to remove the miRNA with no assembly data)
distant <- na.omit(mirpeaks[mirpeaks$chr.miR != mirpeaks$chr.snp,])
dim(distant)
distant

#' Compute the distance between the mapped position of the gene expression and the position of the peak
#+ results='hide'
distancemir<-function(peaks){
    dist <- data.frame(start.min=peaks$start.miR - peaks$min.pos, max.end=peaks$max.pos - peaks$end.miR, diff=abs(peaks$pos.snp - ((peaks$end.miR + peaks$start.miR)/2)))
    return(dist)
}

#' Identify the local regulators:
distL <- distancemir(local)
dim(distL)
distL

#' Make sure the correct amount of miRNAs are in the "distant" category (distD not used downstream):
distD <- distancemir(distant)
dim(distD)
distD

#' miRNA expression mapped within the eQTL peak
cisI <- local[distL$start.min > 0 & distL$max.end > 0,]
nrow(cisI)
cisI

#' eQTL peak within the mapped position of the miRNA expression (only one marker within peak)
cisII <- local[distL$start.min < 0 & distL$max.end < 0,]
nrow(cisII)

#' miRNA expressions mapped in close proximity to the eQTL peak (less than 10MB)
# miRNA maps to the right hand side of its peak
cisIIIa <- local[distL$start.min > 0 & distL$max.end < 0,]
cisIIIa <- data.frame(cisIIIa, dist=abs(ifelse(cisIIIa$max.pos - cisIIIa$start.miR < 0,
	cisIIIa$max.pos - cisIIIa$start.miR, cisIIIa$pos.snp - cisIIIa$start.miR)),
	contained=ifelse(cisIIIa$max.pos - cisIIIa$start.miR < 0, "No", "Yes"))
cisIIIa<-data.frame(cisIIIa,position=rep("right",nrow(cisIIIa)))
nrow(cisIIIa)

# miRNA maps to the left hand side of its peak
cisIIIb <- local[distL$start.min < 0 & distL$max.end > 0,]
cisIIIb <- data.frame(cisIIIb, dist=abs(ifelse(cisIIIb$end.miR - cisIIIb$min.pos < 0,
	cisIIIb$end.miR - cisIIIb$min.pos, cisIIIb$end.miR - cisIIIb$pos.snp)),
	contained=ifelse(cisIIIb$end.miR - cisIIIb$min.pos < 0, "No", "Yes"))
cisIIIb<-data.frame(cisIIIb,position=rep("left",nrow(cisIIIb)))
nrow(cisIIIb)

#' miRNA overlapping peak region
cisIII <- rbind(cisIIIa,cisIIIb)
cisIII <- cisIII[cisIII$contained == "Yes",]
nrow(cisIII)

#' miRNAs on the same chromosome as peak
cisIV <- rbind(cisIIIa,cisIIIb)
cisIV <- cisIV[!cisIV$contained == "Yes",]
nrow(cisIV)

#' eQTL mapped less than 5Mb from miRNA
idx <- abs(cisIV$dist) < 5
nrow(cisIV[idx,])

#' eQTL mapped on same chromosome but at a distance greater than 5MB
idx <- abs(cisIV$dist) > 5
nrow(cisIV[idx,])
tranI <- cisIV[idx,]
cisIV <- cisIV[!idx,]
tranI

#' eQTL mapped on a different chromosome than the mapped position of the miRNA expression
tranII <- distant[distant$range.peak > 0 & distant$chr.miR %in% 1:19,]
nrow(tranII)
tranII

#' eQTL mapped on a different chromosome than the mapped position of the miRNA expression but with only one marker in peak
tranIII <- distant[distant$range.peak == 0 & distant$chr.miR %in% 1:19,]
nrow(tranIII)
tranIII

#' ---
#'  
#' Function to plot eQTL peaks
plot.GMA <-
function (ttp, abspos, ...) {
    plot(ttp[,"pos.snp"],ttp[,"pos.miR"],ylim=c(min(abspos),
         max(abspos)),xlim=c(min(abspos),max(abspos)), pch=19,
         cex.main=2,cex.lab=1.3,
         cex.axis=1.2, ...)
    abline(0,1, col="blue", lwd=3)
}

#' eQTL Plot
#' 
#' Change map positions to Mb for mapping
absposmb<-abspos/1e6
head(absposmb)
ttp.peaks$pos.snp<-ttp.peaks$pos.snp/1e6
ttp.peaks$pos.miR<-ttp.peaks$pos.miR/1e6
head(ttp.peaks)

#+ eQTLmap-regulators, fig.align='center'
par(oma=c(5,2,2,2))
plot.GMA(ttp=ttp.peaks, abspos=absposmb, xlab="SNP position Mb",
		ylab="Transcript position Mb", main="miRNA eQTL Map")

points(ttp.peaks[rownames(cisI),"pos.snp"],
	ttp.peaks[rownames(cisI),"pos.miR"],
	col="black", pch=20)

points(ttp.peaks[rownames(tranII),"pos.snp"], ttp.peaks[rownames(tranII),"pos.miR"],
	col="blue", pch=20)

points(ttp.peaks[rownames(tranIII),"pos.snp"], ttp.peaks[rownames(tranIII),"pos.miR"],
	col="forestgreen", pch=20)

points(ttp.peaks[ttp.peaks$SNP %in% nm[4],"pos.snp"], ttp.peaks[ttp.peaks$SNP %in% nm[4],"pos.miR"],
	col="red", pch=20)


add_legend(-0.75,-0.75, legend=c("Overlapping",
	"Different Chromosome", "Single SNP Different-Chr", "Putative Hotspot"),
	pch=c(19,19,19,19,19), bty="o",ncol=2, cex=1.2)
	add_legend(-0.75,-0.75, legend=c("Overlapping",
	"Different Chromosome", "Single SNP Different-Chr", "Putative Hotspot"),
	pch=c(20,20,20,20,20), col=c("black",
	"blue", "forestgreen", "red"), bty="n",ncol=2, cex=1.2)


#' Save maps for markers_around_peaks.R analysis
#' 
#' Complete map (unfiltered SNPs)
comp.snp.map <- data.frame(chr=MSUPRP_meat$map[[1]],pos=MSUPRP_meat$map[[2]])
rownames(comp.snp.map)<-colnames(MSUPRP_meat$geno)
dim(comp.snp.map)
head(comp.snp.map)
#' Reduced map (SNPs used in this eQTL analysis)
snp.map <- MSUPRP_miRNA$map
dim(snp.map)
head(snp.map)
#' Change the scale of the snp position to Mb
snp.map$pos <- snp.map$pos /1e6
comp.snp.map$pos <- comp.snp.map$pos /1e6

#' ## Save data
#' 
#' Save local and distant peaks, and SNP maps:
save(cisI, tranII, tranIII, snp.map, comp.snp.map,
	file="../4_prelim_miRNA_eQTL_local-distant-regulators.Rdata")

#' ----
#' 
#' Looking only at the eQTL peaks on the same chromosome or overlapping with the mapped gene position check how many markers are before and after the gene position.
#' 
#' Add rownames to peaks
rownames(cisIII) <- cisIII$gene

#' Number of snp per chromosome
# Complete map
nsnpall <- c(table(comp.snp.map[,"chr"]))
nsnpall
# Filtered map
nsnpred <- c(table(snp.map[,"chr"]),0)
names(nsnpred) <- 1:19
nsnpred
# Percent
perc <- round(nsnpred / nsnpall * 100)
perc

#' Plot number of snp per chromosome
#+ snpmap, fig.align='center', fig.width=16, fig.height=12
nsnp <- rbind(nsnpred, nsnpall=nsnpall - nsnpred)
bp <- barplot(nsnp, col=c("blue","black"), ylim=c(0,max(nsnpall)+10),
	main="SNP per Chromosome", xlab="chromosome", ylab="#", cex.main=2,
	cex.axis=1.5, cex.lab=1.5, legend=c("Retained", "Removed"))
text(bp, nsnp[1,]-130, labels=nsnp[1,], col="white", cex=1.5)
text(bp, colSums(nsnp)-130, labels=nsnp[2,], col="white", cex=1.5)


#' --------
#' 
#' Classify eQTL Peaks
#' 
#' Local eQTL
cis <- rbind(cisI, cisII, cisIII[cisIII$contained == "Yes",1:13])
cis$miRNA <- as.character(cis$miRNA)
cis$SNP <- as.character(cis$SNP)
nrow(cis)
cisIII <- cisIII[!cisIII$contained == "Yes",]

cis <- data.frame(cis, regulator=rep("cis", nrow(cis)))

#' Distant eQTL

# Different chromosome
tran <- rbind(tranII, tranIII)
tran$miRNA <- as.character(tran$miRNA)
tran$SNP <- as.character(tran$SNP)
nrow(tran)

tran <- data.frame(tran, regulator=rep("tran", nrow(tran)))
tran

regul <- rbind(cis, tran)
rownames(regul) <- NULL

#' Add eQTL peaks lacking assembly information:
regul <- rbind(regul, data.frame(mirpeaks[!mirpeaks$miRNA %in% regul$miRNA,],
	regulator=rep(NA, sum(!mirpeaks$miRNA %in% regul$miRNA))))

#' eQTL plot with updated classifications:
#+ miRNA_eQTL_local_distant, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
par(oma=c(5,2,2,2))
plot.GMA(ttp=ttp.peaks, abspos=absposmb, xlab="SNP position Mb",
		ylab="Transcript position Mb", main="miRNA eQTL Map")

points(ttp.peaks[rownames(cis),"pos.snp"],
	ttp.peaks[rownames(cis),"pos.miR"],
	col="black", pch=20)

points(ttp.peaks[rownames(tran),"pos.snp"], ttp.peaks[rownames(tran),"pos.miR"],
	col="forestgreen", pch=20)

points(ttp.peaks[ttp.peaks$SNP %in% nm[4],"pos.snp"], ttp.peaks[ttp.peaks$SNP %in% nm[4],"pos.miR"],
	col="red", pch=20)

add_legend(-0.75,-0.75, legend=c("Local",
	"Distant", "Putative Hotspot"),
	pch=c(19,19,19), bty="o",ncol=3, cex=1.2)
	add_legend(-0.75,-0.75, legend=c("Local",
	"Distant", "Putative Hotspot"),
	pch=c(20,20,20), col=c("black",
	"forestgreen", "red"), bty="n",ncol=3, cex=1.2)

#' 
#' ## Save data
#' 
#' Save the data frame containing all the regulators, defining them as cis or trans (local or distant)
save(regul, file="../5_miRNA_eQTL_local_distal_regulators.Rdata")