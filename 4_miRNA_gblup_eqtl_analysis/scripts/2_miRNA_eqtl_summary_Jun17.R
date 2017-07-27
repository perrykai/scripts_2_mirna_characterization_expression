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
#' 1. `1_gblup_results_summary.Rdata`
#'
#' 2. `2_gwa_results_summary.Rdata`
#' 
#' 3. `3_msuprp_mirna_gpdata_pheno_counts.Rdata`
#'
#' 4. `4_normalized_dge_object_and_voom_output.Rdata`
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
ls()

library(limma)
library(edgeR)
library(gwaR)
library(plyr)

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")

#' ## Load data
#'
#' Load the gblup and eQTL output:
load("../1_gblup_results_summary.Rdata")

load("../2_gwa_results_summary.Rdata")

#' Load the dge object to obtain the mature miRNA annotation:
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")

#' Load the MSUPRP_miRNA gpdata object for the mapping information:
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
ls()

#' ## Analysis
#'
#' ### Summarize the heritability of the miRNAs, output from GBLUP:
#'
#' The average heritability of all the miRNA expression profiles:
mean(summary.rst.gblup$h2)
summary(summary.rst.gblup$h2)

#' The average heritability of the miRNAs found to be significantly heritable at FDR < 0.05:
#'
#' How many significantly heritable miRNAs found at FDR < 0.05? (Should match from eQTL scan output md file)
sum(summary.rst.gblup$qvalue<0.05)

#' Extract the significantly heritable miRNAs from the summary.rst.gblup dataset and calculate mean h2:
sigh2<-summary.rst.gblup[summary.rst.gblup$qvalue<0.05,]
dim(sigh2)
mean(sigh2$h2)
summary(sigh2$h2)

#' Define the minimum p-value that is not significant (based on q-value < 0.05) as the threshold for plotting significant points
summary(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
sigthresh<-min(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
sigthresh

#' Plot h2 vs -log10(p-value) like before to determine trend in significance and h2:

#+ vol_sig_h2, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))

plot(summary.rst.gblup$h2, -log10(summary.rst.gblup$lrtpvalue),
        xlab = expression("Heritability"~(h^{2})),
    ylab = "-log10(p-value)",
    main = "Significance vs Heritability")
points(summary.rst.gblup$h2[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       -log10(summary.rst.gblup$lrtpvalue)[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       pch=19,col="red")
abline(a = -log10(sigthresh), b = 0, lty = 5)

#' ---
#'
#' ### The Summary Table of GWAS Results:
#'
#'
#' Assign the correct names to the different objects:
map <- MSUPRP_miRNA$map
colnames(map)
annotation <- dge$genes
head(annotation)
colnames(annotation)

#' Annotation must contain columns "chr", "start", "end", and "strand"
colnames(annotation)<-c("Name","chr","start","end","width","strand","type","Alias","Precursors")
annotation$mid.mir<-round(rowMeans(annotation[,c("start", "end")], na.rm=TRUE))
head(annotation)

#' Add the mid.mir to the map object for later use in manhattan plots (arrow where midpoint of transcript lies)
#' 
#' First, substitute the chromosome number for the "chr" column of the annotation file and to add to the chr column of the map:
annotation$chr<-gsub("X", "19", annotation$chr)
annotation$chr<-as.numeric(gsub("chr", "", annotation$chr))
table(annotation$chr)
head(annotation)

#' Extract the chromosome and the position of the miRNA and build a data.frame to add to the map object later:
mirpos<-data.frame(chr=annotation$chr,
	pos=annotation$mid.mir, row.names=annotation$Name)
str(mirpos)
head(mirpos)
dim(mirpos)

if(sum(mirpos$chr != annotation$chr, na.rm=TRUE) !=0) stop ("chr of miR did not add correctly")
if (sum(mirpos$pos != annotation$mid.mir, na.rm=TRUE) !=0) stop("mid-position of miRNA did not add correctly")

eqtlsum<-function(gwarst, map, annot, threshold=0.05, pergene=TRUE){
sigeqtl<-gwarst[gwarst$gwa.qval<threshold,]
mir<-sigeqtl$miRNA
head(mir)
snp<-sigeqtl$SNPid
head(snp)

eqtlrst<-data.frame(
	miRNA=mir, 
	chr.miR=annot[mir, "chr"],
	start.miR=annot[mir,"start"],
	end.miR=annot[mir,"end"],
	mid.miR=annot[mir,"mid.mir"],
	strand=annot[mir,"strand"],
	miRBase.ID=annot[mir,"Alias"],
	precursors=annot[mir,"Precursors"],
	SNP=snp,
	chr.snp=map[match(snp, rownames(map)),"chr"],
	pos.snp=map[match(snp, rownames(map)),"pos"],
	snp.sign=sigeqtl[match(snp,sigeqtl$SNPid),"SNP.sign"],
	qvalue=sigeqtl[match(snp,sigeqtl$SNPid),"gwa.qval"]
	)
if(pergene){
        id<-unique(eqtlrst[,"miRNA"])
        x<-list()

        for (i in id){
            a<-eqtlrst[eqtlrst[,"miRNA"]==i,]

            if (length(unique(a[,"chr.snp"]))==1){
                a<-a[order(a[,"qvalue"])[1],]

                } 
            else {

                b<-a[order(a[,"qvalue"]),]
                a<-list()

                    for (j in unique(b$chr.snp)){
                        a[[j]]<-b[b[,"chr.snp"]==j,][1,]
                }

            a<-do.call(rbind,a)
                }

        x[[i]]<-a

        }

    eqtlrst<-do.call(rbind,x)
    }
    rownames(eqtlrst)<-NULL
    return(eqtlrst)
}

#' Create the summary table, using pergene = TRUE to get gene-wise eQTL peaks:
sum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=TRUE)
dim(sum.eqtl)
head(sum.eqtl)
#' How many unique miRNAs have eQTL?
length(unique(sum.eqtl$miRNA))

#' miR-eQTL peaks at FDR<0.05:
sum.eqtl

#' #### Summary of GWAS results at FDR < 0.05
#'
#' Number of eQTL peaks per chromosome:
table(sum.eqtl$chr.snp)

#' Names of associated miRNAs:
table(sum.eqtl$miRNA)
#' Chromosomes of associated miRNAs:
table(sum.eqtl$chr.miR)

#' Names of associated peak markers:
table(as.character(sum.eqtl$SNP))

#' ---
#' 
#' ### Determining the ranges of associated SNPs per eQTL peak on SSC15 (for ISAG abstract):
#'
#' First, create the summary table at FDR 5% again, this time with pergene=F to identify all markers associated with each eQTL peak:
fullsum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=FALSE)
dim(fullsum.eqtl)

#' Summarize the number of SNPs associated with each miRNA eQTL (some have multiple peaks)
numsnps<-by(fullsum.eqtl, as.character(fullsum.eqtl$miRNA), nrow)
numsnps<-ldply(numsnps, fun=NULL, id=names(numsnps))
colnames(numsnps) <- c("miRNA", "numsnps")
numsnps
sum(numsnps$numsnps)

#' ---
#' 
#' ### Extract peak range data from all miRNA eQTL peaks
#'
#' I can obtain information on the range of each peak based on miRNA (adapted from Deborah's function "peakrng"):
peakrngmir<-function(nmt, sumtb) {

# Positions of Snps in eQTL peak
    nmt <- nmt
    map.CI <- sumtb[sumtb$miRNA == nmt,c("SNP","chr.snp","pos.snp","qvalue")]

# Number of associated markers by chromosome
    chr <- table(map.CI$chr.snp)
    cat("Number of markers per chromosomal peak for",nmt,":")
    print(chr)
    idx <- as.numeric(names(chr))

    min.pos <- unlist(lapply(idx, function(x) min(map.CI[map.CI$chr.snp == x,"pos.snp"])))
    max.pos <- unlist(lapply(idx, function(x) max(map.CI[map.CI$chr.snp == x,"pos.snp"])))

    start.miR <- rep(sumtb[sumtb$miRNA == nmt,"start.miR"][1], length(idx))
    end.miR <- rep(sumtb[sumtb$miRNA == nmt,"end.miR"][1], length(idx))

# Identify the position of marker extreams for each peak
    peaks <- data.frame(miRNA=rep(nmt,length(idx)),
                chr.miR=rep(sumtb[sumtb$miRNA == nmt,"chr.miR"][1], length(idx)),
                start.miR=start.miR, end.miR=end.miR, range.miR=end.miR-start.miR,
                miRBase.ID=rep(sumtb[sumtb$miRNA == nmt,"miRBase.ID"][1], length(idx)),
                chr.snp=idx, range.peak=max.pos-min.pos,
                min.pos=min.pos, max.pos=max.pos, num.snp=as.vector(chr))

    return(peaks)
}

#' nmt = name transcript (in this case, miRNA); so, make a list of the miRNAs and loop through to get the peak range information for each miRNA with significant eQTL peaks
#' sumtb = output from the summary table function, with pergene = FALSE

sigmirnames <- unique(as.character(sum.eqtl$miRNA))
sigmirnames

mirpeaks<-data.frame(do.call(rbind, lapply(sigmirnames, peakrngmir, fullsum.eqtl)), sum.eqtl[,c("SNP", "pos.snp")])
mirpeaks

#' ---
#' 
#' ### Creating Manhattan plots of the six miRNA with the highest numbers of associated SNP markers (for ISAG poster)
#'
#' First, convert the map positions to absolute map positions using Deborah's function "absmap"
#' 
#' Notice how the map object's positions are relative to chromosome:
head(map)
dim(map)

#' Add the map positions of the miRNA with eQTL
map.full<-rbind(map, mirpos[sigmirnames,])
dim(map.full)

if(nrow(map.full) - nrow(map) != length(sigmirnames)) stop ("miRNA map positions not added correctly")

#' Use the absmap function to convert the chromosomal map positions to absolute map positions:
absposmap<-absmap(map.full)
head(absposmap)
tail(absposmap)
#' Notice it didn't include the miRNA with no map position (ssc-miR-140-5p)
absposmap[which(names(absposmap) %in% sigmirnames)]


#' Divide by 1e6 to get absolute position in Mb (nicer x-axis for plots)
head(absposmap/1e6)

#' Use sigpval function (from DV) to calculate the significant p-value cutoff for plotting manhattan plots for each miRNA
#' 
#' Then, provide the vector of miRNA names and the absposmap object to the manhpt function (also from Deborah's func_eqtl.Rdata) and loop through the vector of miRNA names to create the Manhattan plots:

#' ## Visualize
#'
sigmirnames
#+ man_plot_pval, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames){
	pvalcutoff<-sigpval(i, pvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.pval"], qvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.qval"], fdr=0.05)
	pvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.pval"]
	names(pvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=pvals,map=map.full,annotation=NULL,pvalues=TRUE,cutoff=pvalcutoff, arrow=TRUE)

}

#' The Manhattan Plots of q-values should be put on equal y-axes for comparison on poster.
#' 
#' Three of the peaks are strong enough signals to be on a y-axis of 0-10 (-log10(qval))
#' 
#' The remaining 14 peaks are less strong, on a y-axis of 0-4
#' 
sigmirnames4<-sigmirnames[c(1:4,6,8,9,11,12,14:17)]
sigmirnames4
sigmirnames8<-sigmirnames[c(5,7,10)]
sigmirnames8
sigmirnames10<-sigmirnames[13]
sigmirnames10

#+ man_plot_qval_y4, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames4){
	qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
	names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,4))

}

#+ man_plot_qval_y8, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames8){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,8))

}

#+ man_plot_qval_y10, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames10){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,10))

}

#' ## Save data
save(sum.eqtl, fullsum.eqtl, absposmap, map.full, mirpeaks, file = "../3_eqtl_summary_tables_maps.Rdata")