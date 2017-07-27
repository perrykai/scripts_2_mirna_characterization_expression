**Script:** `5_miRNA_eQTL_examine_peaks.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`

**Date:**  `7/19-27/17`

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`

2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`

**Input File(s):** 

1. `2_gwa_results_summary.Rdata`, `3_eqtl_summary_tables_maps.Rdata`

2. `2_mature_mirna_annotation.Rdata`, `3_msuprp_mirna_gpdata_pheno_counts.Rdata`, `4_normalized_dge_object_and_voom_output.Rdata`, `5_Z_G_miRNA_gblup_gwas.Rdata`

**Output File Directory:** 

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`

**Output File(s):** 

1. `7_exam_eqtl_peaks_fix_snps.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

For each of the significant miRNA eQTL peaks identifyed in the miRNA eQTL scan this code will:

1. Compute the variance explained by the markers significanly associated to the miRNA expression in question

2. Test the significance of the variance accounted for by the markers associated to the miRNA expression

3. Fit the top significant marker per miRNA eQTL peak as a covariate in the model

4. Observe resulting peaks (if any) after fixing the top significant SNP per eQTL peak

## Install libraries


```r
rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")

library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(corrplot)
```

## Load data

Load DV's eqtl functions:


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
```

Load required data:


```r
load("../2_gwa_results_summary.Rdata")
load("../3_eqtl_summary_tables_maps.Rdata")
load("../../3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")
load("../../3_build_dge_object_for_eqtl/5_Z_G_miRNA_gblup_gwas.Rdata")
```

## Analysis

Create covariates for analysis:


```r
X <- data.frame(sex=as.factor(MSUPRP_miRNA$covar$sex),
	 selcrit=as.factor(MSUPRP_miRNA$covar$growth_group))
rownames(X) <- MSUPRP_miRNA$covar$id
```

Change miRNAs and markers to characters in correct format:


```r
sum.eqtl$miRNA<-gsub("-",".", as.character(sum.eqtl$miRNA))
sum.eqtl$SNP<-as.character(sum.eqtl$SNP)

fullsum.eqtl$miRNA<-gsub("-",".", as.character(fullsum.eqtl$miRNA))
fullsum.eqtl$SNP<-as.character(fullsum.eqtl$SNP)

rownames(map.full)<-gsub("-",".", rownames(map.full))
tail(map.full)
```

```
##                 chr       pos
## ssc.miR.7135.3p   3  29052314
## ssc.miR.874       2 145381119
## ssc.miR.95        8   4277296
## ssc.miR.9785.5p   3  22233468
## ssc.miR.9810.3p   4  90746085
## ssc.miR.9843.3p   8 122371924
```

Matrix of miRNAs expressions and covariates


```r
data <- data.frame(t(v$E), X, Z[,unique(sum.eqtl$SNP)])

colnames(wtcen)<-gsub("-",".", colnames(wtcen))
head(colnames(wtcen))
```

```
## [1] "ssc.let.7a"    "ssc.let.7c"    "ssc.let.7d.3p" "ssc.let.7d.5p"
## [5] "ssc.let.7e"    "ssc.let.7f"
```

### Compute variance accounted associated markers

Vector of genes to test


```r
rsp <- gsub("-",".", unique(sum.eqtl$miRNA))
length(rsp)
```

```
## [1] 17
```

Design for GBLUP


```r
design <- c(~sex + selcrit)
```

GBLUP



```r
system.time({
    rst.gblup <- lapply(rsp, function(x) gblup(rsp=x, data=data,
        design=design, G=G, vdata=NULL, wt=wtcen, pos=c(T,T)))
    names(rst.gblup) <- rsp
})
```

I first received this error due to a typo in the gblup function of gwaR (notice, rsp is called resp here):


```r
# Error in paste("name of response variable", resp, "should also be present in wt matrix") :
#   object 'resp' not found
```

Once I fixed that typo, I discovered the true error -- the names of the response variables differ between the data mx and the wt mx: 


```r
# Error in gblup.default(rsp = x, data = data, design = design, G = G, vdata = NULL,  :
  # name of response variable ssc.let.7d.5p should also be present in wt matrix
```

Upon editing the colnames of the wtcen object, the function ran smoothly.

---

Test Peak


```r
system.time({
    lrt.peak <- lapply(rsp, function(x) tryCatch(test.peak(gb=rst.gblup[[x]], x=t(Z),
        peak_pos=sum.eqtl[sum.eqtl$miRNA == x, "SNP"]), error=function(e) NULL))
    names(lrt.peak) <- rsp
})
```

```
##    user  system elapsed 
##   6.234   0.330   6.571
```

```r
# Eliminate NULL results
length(lrt.peak)
```

```
## [1] 17
```

```r
lrt.peak <- lrt.peak[unlist(lapply(lrt.peak, length)) == 3]
length(lrt.peak)
```

```
## [1] 5
```

Compute qvalues


```r
pval <- unlist(lapply(lrt.peak, function(x) x$pvalue))
qval <- qvalue(pval, lambda=0)$qvalue
```

Merge the results of the LRT for each eQTL peak


```r
varpeak <- do.call (rbind, lapply(names(lrt.peak), function(x)
		cbind(t(lrt.peak[[x]]$vars[1]), h2=lrt.peak[[x]]$vars[1, 3],
		lrt.peak[[x]]$llik, pvalue=lrt.peak[[x]]$pvalue)))

rownames(varpeak) <- names(qval)

varpeak <- cbind(varpeak, qvalue=qval)
```

---

Summary variance explaned by peak


```r
summary(varpeak$h2)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1703  0.3339  0.3562  0.3380  0.3967  0.4331
```

```r
varpeak
```

```
##                          G         wt        G_bkg        h2       full
## ssc.miR.140.5p  0.03590145 0.05970515  0.005189245 0.3561799  143.71419
## ssc.miR.184     0.44932968 0.47236823  0.115886518 0.4330536  -41.80699
## ssc.miR.429     0.63260791 1.14439184  0.117710861 0.3338810 -106.66825
## ssc.miR.874     0.33950143 0.48053215  0.035735951 0.3967206  -31.56483
## ssc.miR.9785.5p 0.26737500 1.32484255 -0.022080329 0.1702877 -109.75448
##                        red      dif       pvalue       qvalue
## ssc.miR.140.5p   115.91356 27.80063 4.438318e-14 7.397196e-14
## ssc.miR.184      -78.91870 37.11170 3.487921e-18 8.719802e-18
## ssc.miR.429     -131.76438 25.09613 6.969841e-13 8.712301e-13
## ssc.miR.874      -68.91103 37.34620 2.750374e-18 8.719802e-18
## ssc.miR.9785.5p -124.03640 14.28192 4.532959e-08 4.532959e-08
```

### GWA fixed top SNP per eQTL peak

Genome Wide Association fixing top significant SNP per peak for each miRNA expression containing one or more eQTL


```r
Z <- t(Z)
system.time({
	rst.gwa <- lapply(1:nrow(sum.eqtl), function(x) run.gwa(rsp=sum.eqtl$miRNA[x], data=data,
		design=as.formula(paste("~sex + selcrit +",sum.eqtl$SNP[x])), G=G, vdata=NULL,
		wt=wtcen, x=Z[!rownames(Z) %in% sum.eqtl$SNP[x],], LRT=F, threshold = 0.05, returnz = T,
        pos=c(T,T)))
names(rst.gwa)<-1:length(rst.gwa)
})
```

```
## Warning in sqrt(gw[, 2]): NaNs produced

## Warning in sqrt(gw[, 2]): NaNs produced
```

Calculate pvalues from GWA Zscores


```r
gwa.pv <- lapply(rst.gwa, getpvalue, log.p=F, is.z=T)
```

Merge the results of the GWA into a matrix


```r
rst.gwa <- do.call(cbind, rst.gwa)
```

Gene-wise Multiple Test Correction (FDR) for GWA pvalues (compute qvalues)


```r
system.time({
	gwa.qv <- do.call(cbind, lapply(gwa.pv, function(x) qvalue(x)$qvalues))

})
```

```
##    user  system elapsed 
##   6.989   0.018   7.003
```

Merge the pvalues of the GWA into a matrix


```r
gwa.pv <- do.call(cbind, gwa.pv)
```

Merge results of eQTL analysis into one R object


```r
nms.gwa <- sum.eqtl$miRNA
names(nms.gwa) <- colnames(gwa.pv) 
exam.eqtl <- list(varpeak=varpeak, gwa=rst.gwa, gwa.pval=gwa.pv, gwa.qval=gwa.qv, nms.gwa=nms.gwa)
```

Number of genes still containing a significant peak


```r
threshold <- 0.05
qval <- exam.eqtl$gwa.qval
sig <- qval < threshold
gw <- colSums(sig)[colSums(sig)!=0]
length(gw)
```

```
## [1] 9
```

```r
summary(gw)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    1.00   18.50   27.00   50.71   81.00  128.00       2
```

Remove NAs (ssc-miR-429 (column 13) and ssc-miR-7135-3p (column 16))


```r
gw<-gw[!is.na(gw)]
```

Reduce qval matrix, retain only genes with a significant association


```r
qval <- qval[,names(gw)]

nmt <- exam.eqtl$nms.gwa[colnames(qval)]
dim(qval)
```

```
## [1] 38165     7
```

Separate nmt into groups based on y-axis appropriate for Manhattan plots:


```r
nmt
```

```
##                 5                 6                 8                 9 
##  "ssc.miR.140.5p"  "ssc.miR.140.5p"     "ssc.miR.184"     "ssc.miR.184" 
##                14                15                18 
##     "ssc.miR.429" "ssc.miR.6782.3p"     "ssc.miR.874"
```

```r
nmt8<-nmt[1:5]

nmt4<-nmt[6]

nmt10<-nmt[7]
```

---
Investigate why I'm getting NAs in the rst.gwa:


```r
SNP13<-c(sum.eqtl[13,"SNP"], rownames(rst.gwa)[is.na(rst.gwa[,13])])
corZmiR13<-cor(t(Z[SNP13,]))

SNP16<-c(rownames(rst.gwa)[is.na(rst.gwa[,16])])
corZmiR16<-cor(t(Z[SNP16,]))
```

It's because of the LD! Notice the extremely strong correlation between the SNPs in each case.


```r
corZmiR13
```

```
##             ASGA0094554 ASGA0028271 MARC0001121 ALGA0117083 ALGA0118145
## ASGA0094554   1.0000000  -0.8167067   0.9393866  -0.7896967   1.0000000
## ASGA0028271  -0.8167067   1.0000000  -0.7754324   0.9620977  -0.8167067
## MARC0001121   0.9393866  -0.7754324   1.0000000  -0.7552581   0.9393866
## ALGA0117083  -0.7896967   0.9620977  -0.7552581   1.0000000  -0.7896967
## ALGA0118145   1.0000000  -0.8167067   0.9393866  -0.7896967   1.0000000
## MARC0018157  -1.0000000   0.8167067  -0.9393866   0.7896967  -1.0000000
## M1GA0024787  -1.0000000   0.8167067  -0.9393866   0.7896967  -1.0000000
## ASGA0030400   0.8307748  -0.9929296   0.7914451  -0.9682450   0.8307748
## M1GA0009140   1.0000000  -0.8167067   0.9393866  -0.7896967   1.0000000
## H3GA0019376  -0.8191336   0.9787790  -0.7806128   0.9540794  -0.8191336
## MARC0020138   1.0000000  -0.8167067   0.9393866  -0.7896967   1.0000000
##             MARC0018157 M1GA0024787 ASGA0030400 M1GA0009140 H3GA0019376
## ASGA0094554  -1.0000000  -1.0000000   0.8307748   1.0000000  -0.8191336
## ASGA0028271   0.8167067   0.8167067  -0.9929296  -0.8167067   0.9787790
## MARC0001121  -0.9393866  -0.9393866   0.7914451   0.9393866  -0.7806128
## ALGA0117083   0.7896967   0.7896967  -0.9682450  -0.7896967   0.9540794
## ALGA0118145  -1.0000000  -1.0000000   0.8307748   1.0000000  -0.8191336
## MARC0018157   1.0000000   1.0000000  -0.8307748  -1.0000000   0.8191336
## M1GA0024787   1.0000000   1.0000000  -0.8307748  -1.0000000   0.8191336
## ASGA0030400  -0.8307748  -0.8307748   1.0000000   0.8307748  -0.9859219
## M1GA0009140  -1.0000000  -1.0000000   0.8307748   1.0000000  -0.8191336
## H3GA0019376   0.8191336   0.8191336  -0.9859219  -0.8191336   1.0000000
## MARC0020138  -1.0000000  -1.0000000   0.8307748   1.0000000  -0.8191336
##             MARC0020138
## ASGA0094554   1.0000000
## ASGA0028271  -0.8167067
## MARC0001121   0.9393866
## ALGA0117083  -0.7896967
## ALGA0118145   1.0000000
## MARC0018157  -1.0000000
## M1GA0024787  -1.0000000
## ASGA0030400   0.8307748
## M1GA0009140   1.0000000
## H3GA0019376  -0.8191336
## MARC0020138   1.0000000
```

```r
corZmiR16
```

```
##             MARC0056802 H3GA0009141 ASGA0014022 ASGA0014024 ASGA0014023
## MARC0056802   1.0000000  -0.8943923   1.0000000  -0.8449980   0.7548917
## H3GA0009141  -0.8943923   1.0000000  -0.8943923   0.7124859  -0.8372249
## ASGA0014022   1.0000000  -0.8943923   1.0000000  -0.8449980   0.7548917
## ASGA0014024  -0.8449980   0.7124859  -0.8449980   1.0000000  -0.8830306
## ASGA0014023   0.7548917  -0.8372249   0.7548917  -0.8830306   1.0000000
## ALGA0112651   1.0000000  -0.8943923   1.0000000  -0.8449980   0.7548917
## ALGA0124243   1.0000000  -0.8943923   1.0000000  -0.8449980   0.7548917
## ASGA0097769   1.0000000  -0.8943923   1.0000000  -0.8449980   0.7548917
## ASGA0090944  -1.0000000   0.8943923  -1.0000000   0.8449980  -0.7548917
## ALGA0018214  -1.0000000   0.8943923  -1.0000000   0.8449980  -0.7548917
## ALGA0106253   1.0000000  -0.8943923   1.0000000  -0.8449980   0.7548917
##             ALGA0112651 ALGA0124243 ASGA0097769 ASGA0090944 ALGA0018214
## MARC0056802   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## H3GA0009141  -0.8943923  -0.8943923  -0.8943923   0.8943923   0.8943923
## ASGA0014022   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ASGA0014024  -0.8449980  -0.8449980  -0.8449980   0.8449980   0.8449980
## ASGA0014023   0.7548917   0.7548917   0.7548917  -0.7548917  -0.7548917
## ALGA0112651   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ALGA0124243   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ASGA0097769   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ASGA0090944  -1.0000000  -1.0000000  -1.0000000   1.0000000   1.0000000
## ALGA0018214  -1.0000000  -1.0000000  -1.0000000   1.0000000   1.0000000
## ALGA0106253   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
##             ALGA0106253
## MARC0056802   1.0000000
## H3GA0009141  -0.8943923
## ASGA0014022   1.0000000
## ASGA0014024  -0.8449980
## ASGA0014023   0.7548917
## ALGA0112651   1.0000000
## ALGA0124243   1.0000000
## ASGA0097769   1.0000000
## ASGA0090944  -1.0000000
## ALGA0018214  -1.0000000
## ALGA0106253   1.0000000
```

Identify the correlation of the SNPs in the two complete eQTL peaks yielding NAs:

miRNA-429 eQTL peak contains 93 SNPs:


```r
snp.429<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.429","SNP"]
```

One other significant SNP exists on SSC8; this needs to be removed.


```r
snp.429<-snp.429[snp.429!="ALGA0046283"]
length(snp.429)
```

```
## [1] 93
```

Calculate correlation:


```r
corsnp.429<-cor(t(Z[snp.429,]))
```

miRNA-7135-3p eQTL peak contains 15 SNP:


```r
snp.7135<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.7135.3p","SNP"]
length(snp.7135)
```

```
## [1] 15
```

Calculate correlation:


```r
corsnp.7135<-cor(t(Z[snp.7135,]))
```

These correlations will be visualized using corrplot package later on.
---

Create tables to summarize results of fixing peak eQTL:


```r
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
```

Summary Table for all gene-marker associations


```r
sumtb.exam <- lapply(names(nmt), function(x) stb.nmiR(nm=nmt[x],
    qval=qval[,x], map=map.full, annotation=total.mature.annot2, Z=exam.eqtl$gwa[,x],
    threshold=0.05, pergene=F))
names(sumtb.exam) <- names(nmt)
length(sumtb.exam)
```

```
## [1] 7
```

Summary table for all eQTL peaks (pergene=T)


```r
rsumtb.exam <- lapply(names(nmt), function(x) stb.nmiR(nm=nmt[x],
    qval=qval[,x], map=map.full, annotation=total.mature.annot2, Z=exam.eqtl$gwa[,x],
    threshold=0.05, pergene=T))
names(rsumtb.exam) <- names(nmt)
length(rsumtb.exam)
```

```
## [1] 7
```

Create data.frame of eQTL peaks for diferentiationg between local and distant regulators:


```r
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
```

```r
peaks.exam <- lapply(names(nmt), function(x) data.frame(peakrngmir(nmtb=nmt[x],sumtb=sumtb.exam[[x]]), 
	rsumtb.exam[[x]][,c("SNP","pos.snp")]))
names(peaks.exam) <- names(nmt)

mirpeaks$miRNA<-gsub("-",".", mirpeaks$miRNA)
```

Compare the previous peaks with those obtained after fixing the top snp


```r
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
```

```
## ssc.miR.140.5p 
## 
##  All the peaks for this miRNA: 
##            miRNA chr.miR start.miR end.miR range.miR miRBase.ID chr.snp
## 5 ssc.miR.140.5p      NA        NA      NA        NA       <NA>       4
## 6 ssc.miR.140.5p      NA        NA      NA        NA       <NA>       6
##   range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 5          0  7002305  7002305       1 ASGA0017748  7002305
## 6     386244 15626776 16013020       3 H3GA0055369 15626776
## 
##  Fixed this top snp: 
##            miRNA chr.miR start.miR end.miR range.miR miRBase.ID chr.snp
## 5 ssc.miR.140.5p      NA        NA      NA        NA       <NA>       4
##   range.peak min.pos max.pos num.snp         SNP pos.snp
## 5          0 7002305 7002305       1 ASGA0017748 7002305
## 
##  Peak after fixing top SNP: 
##            miRNA chr.miR start.miR end.miR range.miR chr.snp range.peak
## 5 ssc.miR.140.5p    <NA>        NA      NA        NA       6    6077010
##   min.pos  max.pos num.snp         SNP  pos.snp
## 5 9908019 15985029      19 ASGA0095962 15557072
## 
##  
##  
## ssc.miR.140.5p 
## 
##  All the peaks for this miRNA: 
##            miRNA chr.miR start.miR end.miR range.miR miRBase.ID chr.snp
## 5 ssc.miR.140.5p      NA        NA      NA        NA       <NA>       4
## 6 ssc.miR.140.5p      NA        NA      NA        NA       <NA>       6
##   range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 5          0  7002305  7002305       1 ASGA0017748  7002305
## 6     386244 15626776 16013020       3 H3GA0055369 15626776
## 
##  Fixed this top snp: 
##            miRNA chr.miR start.miR end.miR range.miR miRBase.ID chr.snp
## 6 ssc.miR.140.5p      NA        NA      NA        NA       <NA>       6
##   range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 6     386244 15626776 16013020       3 H3GA0055369 15626776
## 
##  Peak after fixing top SNP: 
##            miRNA chr.miR start.miR end.miR range.miR chr.snp range.peak
## 6 ssc.miR.140.5p    <NA>        NA      NA        NA       4          0
##   min.pos max.pos num.snp         SNP pos.snp
## 6 7002305 7002305       1 ASGA0017748 7002305
## 
##  
##  
## ssc.miR.184 
## 
##  All the peaks for this miRNA: 
##          miRNA chr.miR start.miR  end.miR range.miR   miRBase.ID chr.snp
## 8  ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       3
## 9  ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       6
## 10 ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       7
##    range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 8   125903917   9121844 135025761       2 DBWU0000430   9121844
## 9           0 157003423 157003423       1 M1GA0026172 157003423
## 10   84868032  49815607 134683639      50 ALGA0041952  55936003
## 
##  Fixed this top snp: 
##         miRNA chr.miR start.miR  end.miR range.miR   miRBase.ID chr.snp
## 8 ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       3
##   range.peak min.pos   max.pos num.snp         SNP pos.snp
## 8  125903917 9121844 135025761       2 DBWU0000430 9121844
## 
##  Peak after fixing top SNP: 
##         miRNA chr.miR start.miR end.miR range.miR chr.snp range.peak
## 1 ssc.miR.184    chr7  53.88368 53.8837   2.1e-05       6          0
## 2 ssc.miR.184    chr7  53.88368 53.8837   2.1e-05       7   84845881
##     min.pos   max.pos num.snp         SNP   pos.snp
## 1 156993286 156993286       1 M1GA0024350 156993286
## 2  49785872 134631753      67 H3GA0021739  55905779
## 
##  
##  
## ssc.miR.184 
## 
##  All the peaks for this miRNA: 
##          miRNA chr.miR start.miR  end.miR range.miR   miRBase.ID chr.snp
## 8  ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       3
## 9  ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       6
## 10 ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       7
##    range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 8   125903917   9121844 135025761       2 DBWU0000430   9121844
## 9           0 157003423 157003423       1 M1GA0026172 157003423
## 10   84868032  49815607 134683639      50 ALGA0041952  55936003
## 
##  Fixed this top snp: 
##         miRNA chr.miR start.miR  end.miR range.miR   miRBase.ID chr.snp
## 9 ssc.miR.184       7  53883677 53883698        21 MIMAT0002127       6
##   range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 9          0 157003423 157003423       1 M1GA0026172 157003423
## 
##  Peak after fixing top SNP: 
##         miRNA chr.miR start.miR end.miR range.miR chr.snp range.peak
## 1 ssc.miR.184    chr7  53.88368 53.8837   2.1e-05       3          0
## 2 ssc.miR.184    chr7  53.88368 53.8837   2.1e-05       7   81510294
##     min.pos   max.pos num.snp         SNP   pos.snp
## 1 135025761 135025761       1 ASGA0016793 135025761
## 2  52803473 134313767      26 H3GA0021739  55905779
## 
##  
##  
## ssc.miR.429 
## 
##  All the peaks for this miRNA: 
##          miRNA chr.miR start.miR  end.miR range.miR   miRBase.ID chr.snp
## 13 ssc.miR.429       6  58044184 58044205        21 MIMAT0020591       6
## 14 ssc.miR.429       6  58044184 58044205        21 MIMAT0020591       8
##    range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 13   22928010 38548262 61476272      93 ASGA0094554 41809035
## 14          0  6551222  6551222       1 ALGA0046283  6551222
## 
##  Fixed this top snp: 
##          miRNA chr.miR start.miR  end.miR range.miR   miRBase.ID chr.snp
## 14 ssc.miR.429       6  58044184 58044205        21 MIMAT0020591       8
##    range.peak min.pos max.pos num.snp         SNP pos.snp
## 14          0 6551222 6551222       1 ALGA0046283 6551222
## 
##  Peak after fixing top SNP: 
##          miRNA chr.miR start.miR end.miR range.miR chr.snp range.peak
## 14 ssc.miR.429    chr6  58.04418 58.0442   2.1e-05       6   28371338
##     min.pos  max.pos num.snp         SNP  pos.snp
## 14 38548262 66919600      94 ASGA0094554 41809035
## 
##  
##  
## ssc.miR.6782.3p 
## 
##  All the peaks for this miRNA: 
##              miRNA chr.miR start.miR end.miR range.miR   miRBase.ID
## 15 ssc.miR.6782.3p       6    956806  956827        21 MIMAT0037064
##    chr.snp range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 15      10   11498076 20485599 31983675       4 ASGA0094215 29227834
## 
##  Fixed this top snp: 
##              miRNA chr.miR start.miR end.miR range.miR   miRBase.ID
## 15 ssc.miR.6782.3p       6    956806  956827        21 MIMAT0037064
##    chr.snp range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 15      10   11498076 20485599 31983675       4 ASGA0094215 29227834
## 
##  Peak after fixing top SNP: 
##             miRNA chr.miR start.miR  end.miR range.miR chr.snp range.peak
## 1 ssc.miR.6782.3p    chr6  0.956806 0.956827   2.1e-05       5     862846
## 2 ssc.miR.6782.3p    chr6  0.956806 0.956827   2.1e-05      10    3489062
##    min.pos  max.pos num.snp         SNP  pos.snp
## 1 61751304 62614150      11 MARC0014603 61751304
## 2 28853263 32342325       7 DIAS0000707 29229917
## 
##  
##  
## ssc.miR.874 
## 
##  All the peaks for this miRNA: 
##          miRNA chr.miR start.miR   end.miR range.miR   miRBase.ID chr.snp
## 17 ssc.miR.874       2 145381108 145381130        22 MIMAT0025384       2
## 18 ssc.miR.874       2 145381108 145381130        22 MIMAT0025384       3
##    range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 17   23140264 132997331 156137595     116 ALGA0016550 145462524
## 18          0  64022095  64022095       1 ALGA0122273  64022095
## 
##  Fixed this top snp: 
##          miRNA chr.miR start.miR   end.miR range.miR   miRBase.ID chr.snp
## 18 ssc.miR.874       2 145381108 145381130        22 MIMAT0025384       3
##    range.peak  min.pos  max.pos num.snp         SNP  pos.snp
## 18          0 64022095 64022095       1 ALGA0122273 64022095
## 
##  Peak after fixing top SNP: 
##          miRNA chr.miR start.miR  end.miR range.miR chr.snp range.peak
## 18 ssc.miR.874    chr2  145.3811 145.3811   2.2e-05       2   23140264
##      min.pos   max.pos num.snp         SNP   pos.snp
## 18 132997331 156137595     128 ALGA0016550 145462524
## 
##  
## 
```

## Visualize

### Correlation of SNPs giving NAs in the gwa after fixing the peak SNP:

miR-429


```r
corrplot(corsnp.429, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.7, tl.col="black", title="miR-429 peak SNP ASGA0094554")
```

![plot of chunk cor_plot_miR429](figure/cor_plot_miR429-1.tiff)

miR-7135-3p


```r
corrplot(corsnp.7135, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", tl.srt=45, title="miR-7135-3p peak SNP MARC0056802")
```

![plot of chunk cor_plot_miR7135](figure/cor_plot_miR7135-1.tiff)

### Manhattan Plots


```r
par(oma=c(2,2,2,2))
```

```r
x <- lapply(names(nmt4), function(x) manhpt(nm=nmt[x], abspos=absposmap, rst=qval[,x],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,4)))
```

![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-1.tiff)

```r
y <- lapply(names(nmt8), function(x) manhpt(nm=nmt[x], abspos=absposmap, rst=qval[,x],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,8)))
```

![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-1.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-2.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-3.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-4.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-5.tiff)

```r
z <- lapply(names(nmt10), function(x) manhpt(nm=nmt[x], abspos=absposmap, rst=qval[,x],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,10)))
```

![plot of chunk man_plot_qval_y10](figure/man_plot_qval_y10-1.tiff)

Checking the other miR-874 peak, where fixing the top SNP eliminated a very large peak:


```r
manhpt(nm="ssc.miR.874", abspos=absposmap, rst=exam.eqtl$gwa.qval[,17], map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,10))
```

![plot of chunk man_plot_miR874](figure/man_plot_miR874-1.tiff)

```
## $transcript
## [1] "ssc.miR.874"
## 
## $max
## [1] 0.0313
## 
## $num.snp
## [1] 0
```

## Save data

Save eQTL results summary:


```r
save(exam.eqtl, sumtb.exam, rsumtb.exam, peaks.exam, file="../7_exam_eqtl_peaks_fix_snps.Rdata")
```

