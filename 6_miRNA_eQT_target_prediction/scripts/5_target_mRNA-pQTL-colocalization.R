#' **Script:** `5_target_mRNA-pQTL-colocalization.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/scripts`
#' 
#' **Date:**  6/29/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/`
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/`
#' 
#' 1. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/`
#' 
#' **Input File(s):** 
#' 
#' 1. `data.Rdata`, `pQTL_ALL.Rdata`, `inrange_function.Rdata`
#' 
#' 1. `10_mRNA_miRNA_correlation_output.Rdata`
#' 
#' 1. `gpData_PRKAG3.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/`
#' 
#' **Output File(s):** ``
#' 
#' 1. `12_target_mrna_colocalized_pqtl.Rdata`
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
#' The objective of this script is to co-localize the significant target mRNAs with pQTL identified in the MSUPRP dataset
#' 
#' ## Install libraries
#' 
rm(list=ls())
library(methods)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' Session Information
sessionInfo()
#' ## Load data
#' 
#' Load required R objects 
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/data.Rdata")
load("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/10_mRNA_miRNA_correlation_output.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/pQTL_ALL.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/gpData_PRKAG3.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/inrange_function.Rdata")
 
#' ## Analysis
#' 
#' ### Annotation of corrleated target genes
genes <- dge$genes
colnames(genes)[1] <- "chr"
annot <- do.call(rbind, lapply(names(sig.mrnaR), function(x) 
	data.frame(genes[sig.mrnaR[[x]],], miRNA=rep(x, length(sig.mrnaR[[x]])), 
		rst.corR[[x]][sig.mrnaR[[x]],c("cor", "pvalue", "qvalue")])))

#' ### Phenotypic QTL MSUPRP
#' q-values pQTL GWAS
qval <- pqtl$gwa.qval
dim(qval)

#' p-values pQTL GWAS
pval <- pqtl$gwa.pval
dim(pval)

#' Standardized SNP effects pQTL GWA
sdEff <- pqtl$gwa
dim(sdEff)

#' Marker Map
map <- data.frame(chr=MSUPRP_meat$map[[1]], pos=MSUPRP_meat$map[[2]])
rownames(map) <- colnames(MSUPRP_meat$geno)
dim(map)

#' Retain marker information for all pQTL
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(map[names(sig[[x]]),], 
	std.eff=sdEff[names(sig[[x]]),x], pvalue=pval[names(sig[[x]]),x], qvalue=sig[[x]]))
names(sig) <- colnames(qval)
sig <- sig[unlist(lapply(sig, nrow)) > 0]
length(sig)
names(sig)

#' If a pQTL contains more than one peak split each peak
idx <- lapply(sig, function(x) as.numeric(names(table(x$chr))))
idx <- idx[unlist(lapply(idx, function(x) length(x) > 1))]

mp <- lapply(1:length(idx), function(x) lapply(1:length(idx[[x]]), function(y) sig[[names(idx)[x]]][sig[[names(idx)[x]]]$chr == idx[[x]][y],]))
names(mp) <- names(idx)

for (i in 1:length(mp)){
names(mp[[i]]) <- paste(names(mp)[[i]], 1:length(mp[[i]]), sep=".")
}

qtl <- sig[!names(sig) %in% names(mp)]
for (i in 1:length(mp)){
qtl <- c(qtl, mp[[i]])
}

#' pQTL genomic regions
qtlP <- do.call(rbind, lapply(names(qtl), function(x) data.frame(pheno=strsplit(x, "[.]")[[1]][1], 
	chr=unique(qtl[[x]]$chr), start=min(qtl[[x]]$pos), end=max(qtl[[x]]$pos))))
rownames(qtlP) <- names(qtl)
qtlP <- qtlP[order(qtlP$chr),]
tmp <- list()
for(i in unique(qtlP$chr)){
	x <- qtlP[qtlP$chr == i,]
	tmp[[i]] <- x[order(x$start),]
}
qtlP <- do.call(rbind,tmp)
dim(qtlP)

#' ### Co-localization mRNA with pQTL 
#' 
#' Check all mRNA within a eQTL
win <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single=NULL, range=c(start="start", end="end")))
names(win) <- rownames(qtlP)


#' mRNA overlaping left side of pQTL
left <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="start", range=NULL))
names(left) <- rownames(qtlP)

#' mRNA overlaping right side of pQTL
right <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="end", range=NULL))
names(right) <- rownames(qtlP)

#' Merge all mRNA co-localizing with pQTL  
# Merge within and left side
coloc <- lapply(names(win), function(x) rbind(win[[x]], 
	left[[x]][!as.character(left[[x]]$geneID) %in% as.character(win[[x]]$geneID) | 
	!as.character(left[[x]]$miRNA) == as.character(win[[x]]$miRNA),]))
names(coloc) <- names(win)

#' When we create the merged coloc & right side object, 
#' we receive a warning. 
#' 
#' Check why we're getting the warning:
right[[36]][!rownames(right[[36]]) %in% rownames(coloc[[36]]) | !right[[36]]$miRNA == coloc[[36]]$miRNA,]
#' For the num_ribs phenotype (36), one gene overlaps the pQTL peak to the right that is not within the coloc object.
#' So, when building the merged object, the two objects will have different lengths, triggering the warning 
#' and adding the last line of the longer list to the merged object (hence the duplicate, below).

# Merge coloc and right side
coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]], 
	right[[x]][!as.character(right[[x]]$geneID) %in% as.character(coloc[[x]]$geneID) |
	!as.character(right[[x]]$miRNA) == as.character(coloc[[x]]$miRNA),]))
names(coloc) <- names(win)


#' Final list of mRNA targets significantly correlated with miRNA and co-localizing with a pQTL
coloc <- do.call(rbind, lapply(names(coloc), function(x) data.frame(coloc[[x]], 
	pheno=rep(strsplit(x, "[.]")[[1]][1], nrow(coloc[[x]])))))
rownames(coloc) <- NULL
dim(coloc)
head(coloc)

#' Investigate the occurences of multiple rows:
sort(table(as.character(coloc$genes)))

#' This is the only gene-miRNA-phenotype that is duplicated in the table; 
#' all other occurrences of multiple gene names differ either by phenotype, miRNA-associated, or XLOC gene ID.
coloc[coloc$genes=="VASH1",]

#' Extract the pertinent columns for the significant negatively-associated mRNAs overlapping pQTLs:
negcoloc<-coloc[coloc$cor < 0,c("chr", "geneID", "genes", "miRNA", "cor", "pheno")]
negcoloc
as.character(unique(negcoloc$genes))
table(negcoloc$miRNA)
negcoloc[negcoloc$miRNA=="ssc-miR-874",]
negcoloc[negcoloc$miRNA=="ssc-miR-140-5p",]
negcoloc[negcoloc$miRNA=="ssc-miR-6782-3p",]

#' ## Save data
save(coloc, negcoloc, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/12_target_mrna_colocalized_pqtl.Rdata")