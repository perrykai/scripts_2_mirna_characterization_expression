#' **Script:** `5_miRNA_eQTL_examine_peaks.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`
#' 
#' **Date:**  `7/19-27/17`
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 
#' **Input File(s):** 
#' 
#' 1. `2_gwa_results_summary.Rdata`, `3_eqtl_summary_tables_maps.Rdata`
#' 
#' 2. `2_mature_mirna_annotation.Rdata`, `3_msuprp_mirna_gpdata_pheno_counts.Rdata`, `4_normalized_dge_object_and_voom_output.Rdata`, `5_Z_G_miRNA_gblup_gwas.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`
#' 
#' **Output File(s):** 
#' 
#' 1. `7_exam_eqtl_peaks_fix_snps.Rdata`
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
#' For each of the significant miRNA eQTL peaks identifyed in the miRNA eQTL scan this code will:
#' 
#' 1. Compute the variance explained by the markers significanly associated to the miRNA expression in question
#' 
#' 2. Test the significance of the variance accounted for by the markers associated to the miRNA expression
#' 
#' 3. Fit the top significant marker per miRNA eQTL peak as a covariate in the model
#' 
#' 4. Observe resulting peaks (if any) after fixing the top significant SNP per eQTL peak
#'
#' ## Install libraries
rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")

library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(corrplot)


#' ## Load data
#' 
#' Load DV's eqtl functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

#' Load required data:
load("../2_gwa_results_summary.Rdata")
load("../3_eqtl_summary_tables_maps.Rdata")
load("../../3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")
load("../../3_build_dge_object_for_eqtl/5_Z_G_miRNA_gblup_gwas.Rdata")

#' ## Analysis
#' 
#' Create covariates for analysis:
X <- data.frame(sex=as.factor(MSUPRP_miRNA$covar$sex),
	 selcrit=as.factor(MSUPRP_miRNA$covar$growth_group))
rownames(X) <- MSUPRP_miRNA$covar$id

#' Change miRNAs and markers to characters in correct format:
sum.eqtl$miRNA<-gsub("-",".", as.character(sum.eqtl$miRNA))
sum.eqtl$SNP<-as.character(sum.eqtl$SNP)

fullsum.eqtl$miRNA<-gsub("-",".", as.character(fullsum.eqtl$miRNA))
fullsum.eqtl$SNP<-as.character(fullsum.eqtl$SNP)

rownames(map.full)<-gsub("-",".", rownames(map.full))
tail(map.full)


#' Matrix of miRNAs expressions and covariates
data <- data.frame(t(v$E), X, Z[,unique(sum.eqtl$SNP)])

colnames(wtcen)<-gsub("-",".", colnames(wtcen))
head(colnames(wtcen))

#' ### Compute variance accounted associated markers
#' 
#' Vector of genes to test
rsp <- gsub("-",".", unique(sum.eqtl$miRNA))
length(rsp)

#' Design for GBLUP
design <- c(~sex + selcrit)

#' GBLUP
#' 
#+ results='hide'
system.time({
    rst.gblup <- lapply(rsp, function(x) gblup(rsp=x, data=data,
        design=design, G=G, vdata=NULL, wt=wtcen, pos=c(T,T)))
    names(rst.gblup) <- rsp
})

#' I first received this error due to a typo in the gblup function of gwaR (notice, rsp is called resp here):

# Error in paste("name of response variable", resp, "should also be present in wt matrix") :
#   object 'resp' not found

#' Once I fixed that typo, I discovered the true error -- the names of the response variables differ between the data mx and the wt mx: 

# Error in gblup.default(rsp = x, data = data, design = design, G = G, vdata = NULL,  :
  # name of response variable ssc.let.7d.5p should also be present in wt matrix

#' Upon editing the colnames of the wtcen object, the function ran smoothly.
#' 
#' ---
#' 
#' Test Peak
system.time({
    lrt.peak <- lapply(rsp, function(x) tryCatch(test.peak(gb=rst.gblup[[x]], x=t(Z),
        peak_pos=sum.eqtl[sum.eqtl$miRNA == x, "SNP"]), error=function(e) NULL))
    names(lrt.peak) <- rsp
})

# Eliminate NULL results
length(lrt.peak)
lrt.peak <- lrt.peak[unlist(lapply(lrt.peak, length)) == 3]
length(lrt.peak)

#' Compute qvalues
pval <- unlist(lapply(lrt.peak, function(x) x$pvalue))
qval <- qvalue(pval, lambda=0)$qvalue

#' Merge the results of the LRT for each eQTL peak
varpeak <- do.call (rbind, lapply(names(lrt.peak), function(x)
		cbind(t(lrt.peak[[x]]$vars[1]), h2=lrt.peak[[x]]$vars[1, 3],
		lrt.peak[[x]]$llik, pvalue=lrt.peak[[x]]$pvalue)))

rownames(varpeak) <- names(qval)

varpeak <- cbind(varpeak, qvalue=qval)

#' ---
#' 
#' Summary variance explaned by peak
summary(varpeak$h2)

varpeak


#' ### GWA fixed top SNP per eQTL peak
#' 
#' Genome Wide Association fixing top significant SNP per peak for each miRNA expression containing one or more eQTL

#+ results='hide'
Z <- t(Z)
system.time({
	rst.gwa <- lapply(1:nrow(sum.eqtl), function(x) run.gwa(rsp=sum.eqtl$miRNA[x], data=data,
		design=as.formula(paste("~sex + selcrit +",sum.eqtl$SNP[x])), G=G, vdata=NULL,
		wt=wtcen, x=Z[!rownames(Z) %in% sum.eqtl$SNP[x],], LRT=F, threshold = 0.05, returnz = T,
        pos=c(T,T)))
names(rst.gwa)<-1:length(rst.gwa)
})

#' Calculate pvalues from GWA Zscores
gwa.pv <- lapply(rst.gwa, getpvalue, log.p=F, is.z=T)

#' Merge the results of the GWA into a matrix
rst.gwa <- do.call(cbind, rst.gwa)

#' Gene-wise Multiple Test Correction (FDR) for GWA pvalues (compute qvalues)
system.time({
	gwa.qv <- do.call(cbind, lapply(gwa.pv, function(x) qvalue(x)$qvalues))

})

#' Merge the pvalues of the GWA into a matrix
gwa.pv <- do.call(cbind, gwa.pv)

#' Merge results of eQTL analysis into one R object
nms.gwa <- sum.eqtl$miRNA
names(nms.gwa) <- colnames(gwa.pv) 
exam.eqtl <- list(varpeak=varpeak, gwa=rst.gwa, gwa.pval=gwa.pv, gwa.qval=gwa.qv, nms.gwa=nms.gwa)

#' Number of genes still containing a significant peak
threshold <- 0.05
qval <- exam.eqtl$gwa.qval
sig <- qval < threshold
gw <- colSums(sig)[colSums(sig)!=0]
length(gw)
summary(gw)

#' Remove NAs (ssc-miR-429 (column 13) and ssc-miR-7135-3p (column 16))
gw<-gw[!is.na(gw)]

#' Reduce qval matrix, retain only genes with a significant association
qval <- qval[,names(gw)]

nmt <- exam.eqtl$nms.gwa[colnames(qval)]
dim(qval)

#' Separate nmt into groups based on y-axis appropriate for Manhattan plots:
nmt
nmt8<-nmt[1:5]

nmt4<-nmt[6]

nmt10<-nmt[7]

#' ---
#' Investigate why I'm getting NAs in the rst.gwa:
SNP13<-c(sum.eqtl[13,"SNP"], rownames(rst.gwa)[is.na(rst.gwa[,13])])
corZmiR13<-cor(t(Z[SNP13,]))

SNP16<-c(rownames(rst.gwa)[is.na(rst.gwa[,16])])
corZmiR16<-cor(t(Z[SNP16,]))

#' It's because of the LD! Notice the extremely strong correlation between the SNPs in each case.
corZmiR13

corZmiR16

#' Identify the correlation of the SNPs in the two complete eQTL peaks yielding NAs:
#' 
#' miRNA-429 eQTL peak contains 93 SNPs:
snp.429<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.429","SNP"]
#' One other significant SNP exists on SSC8; this needs to be removed.
snp.429<-snp.429[snp.429!="ALGA0046283"]
length(snp.429)
#' Calculate correlation:
corsnp.429<-cor(t(Z[snp.429,]))

#' miRNA-7135-3p eQTL peak contains 15 SNP:
snp.7135<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.7135.3p","SNP"]
length(snp.7135)

#' Calculate correlation:
corsnp.7135<-cor(t(Z[snp.7135,]))

#' These correlations will be visualized using corrplot package later on.

#' ---
#' 
#' Create tables to summarize results of fixing peak eQTL:
rownames(total.mature.annot2)<-gsub("-",".", rownames(total.mature.annot2))
total.mature.annot2$Name<-gsub("-",".", total.mature.annot2$Name)

stb.nmiR<-function(nm,qval,map,annotation,Z,threshold=0.01,gene.name="genes",pergene=TRUE){
    idx <- qval < threshold
    snp <- names(qval)[idx]
    snp.effect <- ifelse(Z[idx] < 0, x<-"-", x<-"+")
    gene <- rep(nm, sum(idx))

    rst <- data.frame(miR=gene, chr.miR=annotation[gene,"chr0"],
        start.miR=annotation[gene,"start"]/1e6, end.miR=annotation[gene,"end"]/1e6,
        strand=annotation[gene,"strand"], SNP=snp, chr.snp=map[snp,"chr"], pos.snp=map[snp,"pos"], 
        snp.effect=snp.effect,
        qvalue=qval[idx], row.names=NULL)

    if(pergene){
        id<-unique(rst[,"miR"])
        x<-list()

        for (i in id){
            a<-rst[rst[,"miR"]==i,]

            if (length(unique(a[,"chr.snp"]))==1){
                a<-a[order(a[,"qvalue"])[1],]

                } else {

                b<-a[order(a[,"qvalue"]),]
                a<-list()

                    for (j in unique(b$chr.snp)){
                        a[[j]]<-b[b[,"chr.snp"]==j,][1,]
                    }

                a<-do.call(rbind,a)
                }

        x[[i]]<-a

        }

    rst<-do.call(rbind,x)
    }
    rownames(rst)<-NULL
    return(rst)
}


#' Summary Table for all gene-marker associations
sumtb.exam <- lapply(names(nmt), function(x) stb.nmiR(nm=nmt[x],
    qval=qval[,x], map=map.full, annotation=total.mature.annot2, Z=exam.eqtl$gwa[,x],
    threshold=0.05, pergene=F))
names(sumtb.exam) <- names(nmt)
length(sumtb.exam)

#' Summary table for all eQTL peaks (pergene=T)
rsumtb.exam <- lapply(names(nmt), function(x) stb.nmiR(nm=nmt[x],
    qval=qval[,x], map=map.full, annotation=total.mature.annot2, Z=exam.eqtl$gwa[,x],
    threshold=0.05, pergene=T))
names(rsumtb.exam) <- names(nmt)
length(rsumtb.exam)

#' Create data.frame of eQTL peaks for diferentiationg between local and distant regulators:

peakrngmir<-function(nmtb, sumtb){

# Positions of SNPs in eQTL peak
    nmtb <- nmtb
    map.CI <- sumtb[sumtb$miR==nmtb,c("SNP","chr.snp","pos.snp","qvalue")]

# Number of associated markers by chromosome
    chr <- table(map.CI$chr.snp)
    cat("Number of markers per chromosomal peak for",nmtb,":")
    print(chr)
    idx <- as.numeric(names(chr))

    min.pos <- unlist(lapply(idx, function(x) min(map.CI[map.CI$chr.snp == x,"pos.snp"])))
    max.pos <- unlist(lapply(idx, function(x) max(map.CI[map.CI$chr.snp == x,"pos.snp"])))

    start.miR <- rep(sumtb[sumtb$miR == nmtb,"start.miR"][1], length(idx))
    end.miR <- rep(sumtb[sumtb$miR == nmtb,"end.miR"][1], length(idx))

# Identify the position of marker extreams for each peak
    peaks <- data.frame(miRNA=rep(nmtb,length(idx)),
                chr.miR=rep(sumtb[sumtb$miR == nmtb,"chr.miR"][1], length(idx)),
                start.miR=start.miR, end.miR=end.miR, range.miR=end.miR-start.miR,
                chr.snp=idx, range.peak=max.pos-min.pos,
                min.pos=min.pos, max.pos=max.pos, num.snp=as.vector(chr))

    return(peaks)
}


#+ results='hide'
peaks.exam <- lapply(names(nmt), function(x) data.frame(peakrngmir(nmtb=nmt[x],sumtb=sumtb.exam[[x]]), 
	rsumtb.exam[[x]][,c("SNP","pos.snp")]))
names(peaks.exam) <- names(nmt)

mirpeaks$miRNA<-gsub("-",".", mirpeaks$miRNA)

#' Compare the previous peaks with those obtained after fixing the top snp
for (i in names(nmt)){
cat(nmt[i], '\n')
cat('\n','All the peaks for this miRNA:','\n')
print(mirpeaks[mirpeaks$miRNA == nmt[i],])
cat('\n','Fixed this top snp:','\n')
print(mirpeaks[names(nmt),][i,])
cat('\n','Peak after fixing top SNP:','\n')
print(peaks.exam[[i]])
cat('\n','\n','\n')
}


#' ## Visualize
#' 
#' ### Correlation of SNPs giving NAs in the gwa after fixing the peak SNP:
#' 
#' miR-429
#+ cor_plot_miR429, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(corsnp.429, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.7, tl.col="black", title="miR-429 peak SNP ASGA0094554")
#' miR-7135-3p
#+ cor_plot_miR7135, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(corsnp.7135, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", tl.srt=45, title="miR-7135-3p peak SNP MARC0056802")

#' ### Manhattan Plots
par(oma=c(2,2,2,2))
#+ man_plot_qval_y4, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
x <- lapply(names(nmt4), function(x) manhpt(nm=nmt[x], abspos=absposmap, rst=qval[,x],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,4)))
#+ man_plot_qval_y8, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
y <- lapply(names(nmt8), function(x) manhpt(nm=nmt[x], abspos=absposmap, rst=qval[,x],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,8)))
#+ man_plot_qval_y10, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
z <- lapply(names(nmt10), function(x) manhpt(nm=nmt[x], abspos=absposmap, rst=qval[,x],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,10)))

#' Checking the other miR-874 peak, where fixing the top SNP eliminated a very large peak:
#+ man_plot_miR874, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
manhpt(nm="ssc.miR.874", abspos=absposmap, rst=exam.eqtl$gwa.qval[,17], map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,10))

#' ## Save data
#' 
#' Save eQTL results summary:
save(exam.eqtl, sumtb.exam, rsumtb.exam, peaks.exam, file="../7_exam_eqtl_peaks_fix_snps.Rdata")