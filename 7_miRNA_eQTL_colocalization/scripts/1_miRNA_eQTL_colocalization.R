#' **Script:** `1_miRNA_eQTL_colocalization.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/7_miRNA_eQTL_colocalization/scripts`
#' 
#' **Date:**  07/14/17 - 7/18/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/`
#'
#' 3. `/mnt/research/pigeqtl/analyses/microRNA/3_build_dge_object_for_eqtl/`
#' 
#' 4. `/mnt/research/pigeqtl/analyses/microRNA/4_miRNA_gblup_eqtl_analysis/`
#' 
#' **Input File(s):** 
#' 
#' 1. `pQTL_ALL.Rdata`, `inrange_function.Rdata`
#' 
#' 2. `gpData_PRKAG3.Rdata`
#' 
#' 3. `2_mature_mirna_annotation.Rdata`, `3_msuprp_mirna_gpdata_pheno_counts.Rdata`, `4_normalized_dge_object_and_voom_output.Rdata`, `5_Z_G_miRNA_gblup_gwas.Rdata`
#' 
#' 4. `1_gblup_results_summary.Rdata`, `2_gwa_results_summary.Rdata`, `5_miRNA_eQTL_local_distal_regulators.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/7_miRNA_eQTL_colocalization/`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_miRNA_eQTL_pQTL_colocalized_peaks.Rdata`
#' 
#' 2. `2_miRNA_eQTL_pQTL_colocalized_peaks.txt`
#' 
#' 3. `3_pQTL_summary.txt`
#' 
#' 4. `4_summary_phenotypes_174.txt`
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
#' The objective of this script is to investigate if any miRNA eQTL co-localize with pQTL and/or mRNA eQTL.
#' 
#' ## Install libraries
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/7_miRNA_eQTL_colocalization/scripts")

rm(list=ls())

library(methods)
library(limma)
library(edgeR)

#' ## Load data
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/pQTL_ALL.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/gpData_PRKAG3.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/inrange_function.Rdata")

#' Load DV functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

load("../../3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")
load("../../3_build_dge_object_for_eqtl/5_Z_G_miRNA_gblup_gwas.Rdata")
load("../../4_miRNA_gblup_eqtl_analysis/1_gblup_results_summary.Rdata")
load("../../4_miRNA_gblup_eqtl_analysis/2_gwa_results_summary.Rdata")
load("../../4_miRNA_gblup_eqtl_analysis/5_miRNA_eQTL_local_distal_regulators.Rdata")

ls()

#' ## Analysis
#' 
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

# Number of associated markers per QTL peak
length(qtl)
unlist(lapply(qtl, nrow))


#' Co-localization eQTL with pQTL
# Check all eQTL peaks within a pQTL
sig <- cbind(regul[,c(7:14)], regul[,c(1:6)])
# sig[,c(2:4,7,11:13)] <- sig[,c(2:4,7,11:13)] * 1e6

eqtl.pqtl <- lapply(qtl, function(x) inrange(chr=unique(x$chr), start=min(x$pos),
        end=max(x$pos), map=sig, single="pos.snp", range=NULL))
unlist(lapply(eqtl.pqtl, nrow))[unlist(lapply(eqtl.pqtl, nrow)) > 0]

# Check all pQTL peaks within a eQTL
rqtl <- do.call(rbind,lapply(qtl, function(x) data.frame(chr=unique(x$chr),
        min=min(x$pos), max=max(x$pos))))

PinE <- lapply(1:nrow(sig), function(x) inrange(chr=sig[x, "chr.snp"], start=sig[x, "min.pos"],
        end=sig[x, "max.pos"], map=rqtl, single=NULL, range=c(start="min", end="max")))
names(PinE) <- 1:nrow(sig)
unlist(lapply(PinE, nrow))[unlist(lapply(PinE, nrow)) > 0]

PinE <- do.call(rbind, lapply(names(PinE), function(x)
        data.frame(pqtl=rownames(PinE[[x]]), PinE[[x]], eqtl=rep(x, nrow(PinE[[x]])))))
rownames(PinE) <- NULL

PinE <- lapply(names(qtl), function(x) sig[as.character(PinE$eqtl[grep(x, PinE$pqtl)]),])
names(PinE) <- names(qtl)
unlist(lapply(PinE, nrow))[unlist(lapply(PinE, nrow)) > 0]

# pQTL peaks overlaping left side of eQTL
left <- lapply(qtl, function(x) inrange(chr=unique(x$chr), start=min(x$pos),
        end=max(x$pos), map=sig, single="min.pos", range=NULL))
unlist(lapply(left, nrow))[unlist(lapply(left, nrow)) > 0]

# pQTL peaks overlaping right side of eQTL
right <- lapply(qtl, function(x) inrange(chr=unique(x$chr), start=min(x$pos),
        end=max(x$pos), map=sig, single="max.pos", range=NULL))
unlist(lapply(right, nrow))[unlist(lapply(right, nrow)) > 0]

#' Merge all pQTL co-localizing with eQTL
coloc <- lapply(names(eqtl.pqtl), function(x) rbind(eqtl.pqtl[[x]],
        PinE[[x]][!rownames(PinE[[x]]) %in% rownames(eqtl.pqtl[[x]]),]))
names(coloc) <- names(eqtl.pqtl)

coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]],
        left[[x]][!rownames(left[[x]]) %in% rownames(coloc[[x]]),]))
names(coloc) <- names(eqtl.pqtl)

coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]],
        right[[x]][!rownames(right[[x]]) %in% rownames(coloc[[x]]),]))
names(coloc) <- names(eqtl.pqtl)

# Matrix of pQTL and colocalized eQTL
qtlM <- do.call(rbind, lapply(names(qtl), function(x)
        data.frame(chr=unique(qtl[[x]]$chr), lowb=min(qtl[[x]]$pos), upperb=max(qtl[[x]]$pos),
                eQTL=nrow(coloc[[x]]))))
rownames(qtlM) <- names(qtl)

# Number of pQTL with co-localized eQTL
sum(qtlM$eQTL > 0)

# Summary colocalized QTL
coloc <- coloc[unlist(lapply(coloc, nrow)) > 0]
unlist(lapply(coloc, nrow))


#' Create table with pQTL information
sum.qtl <- do.call(rbind, lapply(qtl, function(x) data.frame(chr=unique(x$chr), start=min(x$pos), end=max(x$pos),
        SNP=rownames(x)[min(x$pvalue) == x$pvalue][1], x[min(x$pvalue) == x$pvalue,][1,
        c("pos", "std.eff", "pvalue", "qvalue")], nSNP=nrow(x))))
sum.qtl <- sum.qtl[order(sum.qtl$chr),]

# Add heritability of each pQTL and number of co-localized eQTL
idx <- unlist(lapply(rownames(sum.qtl), function(x) strsplit(x, "[.]")[[1]][1]))
sum.qtl <- data.frame(sum.qtl, h2=pqtl$gblup[idx,"h2"], pval.h2=pqtl$gblup[idx,"pvalue"],
        qval.h2=pqtl$gblup[idx,"qvalue"], eQTL=unlist(lapply(rownames(sum.qtl),
        function(x) ifelse(x %in% names(coloc), nrow(coloc[[x]]), 0))))
dim(sum.qtl)

#' ---
#' Order colocalized eQTL per chromosome
coloc <- coloc[rownames(sum.qtl)]
coloc <- coloc[!is.na(names(coloc))]

# Number of eQTL co-localized with pQTL
length(unlist(lapply(coloc,rownames)))

# Number of unique eQTL co-localized with pQTL
length(unique(unlist(lapply(coloc,rownames))))

#' Add colocalized eQTL-pQTL infomation to eQTL table (regul R object)
lst <- do.call(rbind,lapply(names(coloc), function(x)
        data.frame(eqtl=rownames(coloc[[x]]), pqtl=strsplit(x,"[.]")[[1]][1])))

regul <- cbind(regul,
        # snp.effect=unlist(lapply(1:nrow(regul),
        #         function(x) eqtl$gwa[regul$SNP[x], regul$miRNA[x]])),
        qval.snp=unlist(lapply(1:nrow(regul),
               function(x) rst.gwa[(rst.gwa$miRNA==regul$miRNA[x] & rst.gwa$SNPid==regul$SNP[x]),"gwa.qval"])),
        pval.snp=unlist(lapply(1:nrow(regul),
               function(x) rst.gwa[(rst.gwa$miRNA==regul$miRNA[x] & rst.gwa$SNPid==regul$SNP[x]),"gwa.pval"])),
        summary.rst.gblup[regul$miRNA,c("h2","lrtpvalue","qvalue")],
        colocalized.pqtl=unlist(lapply(rownames(regul),
                function(x) paste(as.character(lst[lst$eqtl == x, "pqtl"]), collapse=', '))))

#' Reorder columns in regul data frame
regul <- regul[,c("chr.snp", "SNP", "pos.snp", "pval.snp", "qval.snp", "min.pos", "max.pos", "range.peak", "num.snp",
        "miRNA", "chr.miR", "start.miR", "end.miR", "range.miR", "h2", "lrtpvalue", "qvalue", "regulator", "colocalized.pqtl")]

head(regul)

# Change megabase positions to base positions
# idx <- c("min.pos", "max.pos", "range.peak", "pos.snp", "start.miR", "end.miR", "range.miR")
# regul[,idx] <- regul[,idx]

#' Order regul matrix by peak chromosome and position
tmp <- lapply(1:18, function(x) regul[regul$chr.snp == x,])
regul <- do.call(rbind, lapply(tmp, function(x) x[order(x$pos.snp),]))

head(regul)

#' ---
#' Estimate the absolute position of markers and gene expressions
#+ warning=FALSE
posann <- (total.mature.annot2$end+total.mature.annot2$start)/2000000
names(posann) <- as.character(rownames(total.mature.annot2))
mapZ <- map
mapZ$pos <- mapZ$pos/1000000
mapt <- data.frame(chr=total.mature.annot2$chr0,pos=posann)
mapt <- mapt[!is.na(as.numeric(as.character(mapt$chr))), ]
map <- rbind(mapZ, mapt)
map$chr <- as.numeric(map$chr)
map <- map[order(map$chr), ]

#' Absolute postions of genes and markers
abspos <- absmap(map)


#' Summary statistics for each phenotype and selected 174 animals
pheno <- MSUPRP_meat$pheno[rownames(MSUPRP_miRNA$pheno[,,1]),,1]
sumstat <- data.frame(N=apply(pheno,2, function(x) length(x) - sum(is.na(x))),
        Mean=apply(pheno,2, mean, na.rm=T), SD=apply(pheno,2, sd, na.rm=T))

#' ## Visualize
#' 
#' #### Manhattan Plots: miRNA eQTL colocalized with pQTL
idx <- lapply(coloc, function(x) as.character(x$miRNA))

#' Make sure the rst.gwa df is ordered by miRNA then by SNPid, so that SNPid can be used for the rownames

dim(rst.gwa)
dim(rst.gwa[order(rst.gwa[,c("miRNA","SNPid")]),])

# eQTL qvalues
qvalE <- lapply(idx, function(x) data.frame(rst.gwa[(rst.gwa$miRNA%in%x), c("miRNA","SNPid","gwa.qval")]))
str(qvalE)
traits<-names(qvalE)


qvalE<-lapply(names(qvalE), function(x) as.data.frame(split(qvalE[[x]]$"gwa.qval", qvalE[[x]]$"miRNA"),col.names=idx[[x]]))
names(qvalE)<-traits

str(qvalE)

# pQTL qvalues
nms <- unique(unlist(lapply(names(coloc), function(x) strsplit(x, "[.]")[[1]][1])))
qvalP <- qval[,nms]

#' Growth Phenotypes
idx <- c("bf10_13wk","bf10_16wk",
        "lrf_19wk")

#+ growth, fig.align='center', fig.width=8, fig.height=16, echo=FALSE
par(oma=c(2,2,2,2), mfrow=c(1,3))
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 6))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        points(abspos[rownames(qvalP)], -log10(qvalP[,i]), pch=20,
                col=c("deepskyblue","blue")[(mps%%2)+1])
        abline(h=-log10(0.05),col="red")
}

#' Carcass composition phenotypes
idx <- c("ph_24h","num_ribs","car_length")
#+ carcass, fig.align='center', fig.width=8, fig.height=16, echo=FALSE
par(oma=c(2,2,2,2), mfrow=c(1,3))
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 10))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        points(abspos[rownames(qvalP)], -log10(qvalP[,i]), pch=20,
                col=c("deepskyblue","blue")[(mps%%2)+1])
        abline(h=-log10(0.05),col="red")
}

#' Meat quality phenotypes
idx <- c("WBS","tenderness","overtend","juiciness","cook_yield","driploss","protein")
#+ quality, fig.align='center', fig.width=8, fig.height=16, echo=FALSE
par(oma=c(2,2,2,2), mfrow=c(3,3))
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 25))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        points(abspos[rownames(qvalP)], -log10(qvalP[,i]), pch=20,
                col=c("deepskyblue","blue")[(mps%%2)+1])
        abline(h=-log10(0.05),col="red")
}

#' ## Save data
#' 
#' Save the colocalized miRNA eQTL pQTL peaks:
save(regul, file="../1_miRNA_eQTL_pQTL_colocalized_peaks.Rdata")
#' Write table of colocalized miRNA eQTL pQTL peaks:
write.table(regul, quote=F, col.names=T, row.names=F, sep="\t", file="../2_miRNA_eQTL_pQTL_colocalized_peaks.txt")
#' Save pQTL table to text file
write.table(sum.qtl, quote=F, row.names=T, col.names=T, sep="\t",file="../3_pQTL_summary.txt")
#' Save summary statistics for each phenotype and selected 174 animals
write.table(sumstat, quote=F, col.names=T, row.names=T, sep="\t", file="../4_summary_phenotypes_174.txt")