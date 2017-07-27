**Script:** `4_miRNA_eQTL_hotspot_characterization.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`

**Date:**  07/12/17

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`

2. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/`

**Input File(s):** 

1. `1_precursor_mirna_annotation.Rdata`

1. `2_mature_mirna_annotation.Rdata`

1. `3_msuprp_mirna_gpdata_pheno_counts.Rdata`

2. `4_normalized_dge_object_and_voom_output.Rdata`

3. `5_Z_G_miRNA_gblup_gwas.Rdata`

4. `gpData_PRKAG3.Rdata`

**Output File Directory:** 

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`

**Output File(s):** 

1. `6_hotspot_miRNA_correlation_results.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to investigate the putative miRNA eQTL hotspots to determine if they are truly hotspots, or spurious associations due to high correlation of miRNA expression or genotype data.

To do this, miRNA expression will be correlated between hotspot miRNAs utilizing the log-cpm (v$E) . 
Pearson correlation will be used. 

I will also investigate the genomic origins of these miRNAs and SNPs, to see if the miRNAs are coming from similar genomic regions, if they have similar seed sequences, etc.

## Install libraries


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")

rm(list=ls())

library(regress)
library(gwaR)
library(limma)
library(edgeR)
library(parallel)
library(qvalue)
library(corrplot)
```

## Load data

Load microRNA expression data


```r
load("../../3_build_dge_object_for_eqtl/1_precursor_mirna_annotation.Rdata")
load("../../3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")
load("../../3_build_dge_object_for_eqtl/5_Z_G_miRNA_gblup_gwas.Rdata")
```

Load the MSUPRP_meat gpdata object (for allele freq calculation):


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/gpData_PRKAG3.Rdata")

rm(PRKAG3)

ls()
```

```
##  [1] "dge"                  "G"                    "MSUPRP_meat"         
##  [4] "MSUPRP_miRNA"         "precursor.mirannot"   "summary_MSUPRP_miRNA"
##  [7] "total.mature.annot2"  "v"                    "wtcen"               
## [10] "Z"
```

## Analysis

Extract the names of the miRNA in the hotspots:

Three hotspots (H3GA0052416, MARC0027291, MARC0047188) were associated with the same four miRNAs:


```r
htspt.mirna4 <- as.character(c("ssc-let-7d-5p","ssc-miR-345-3p","ssc-miR-95","ssc-miR-9843-3p"))
htspt.mirna4
```

```
## [1] "ssc-let-7d-5p"   "ssc-miR-345-3p"  "ssc-miR-95"      "ssc-miR-9843-3p"
```

One hotspot (MARC0093624) was associated with five miRNAs, overlapping to some degree with the other list:


```r
htspt.mirna5 <- as.character(c("ssc-let-7d-5p","ssc-let-7g","ssc-miR-1468","ssc-miR-95","ssc-miR-9843-3p"))
htspt.mirna5
```

```
## [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-1468"    "ssc-miR-95"     
## [5] "ssc-miR-9843-3p"
```

Correlate the miRNA expression profiles using the log-cpm counts of miRNA expression:


```r
data4<-t(v$E[htspt.mirna4,])
head(data4)
```

```
##      ssc-let-7d-5p ssc-miR-345-3p ssc-miR-95 ssc-miR-9843-3p
## 1034      11.63362       6.601181   11.44954        7.893082
## 1036      11.10306       6.476953   10.90727        7.694085
## 1041      11.28869       6.288982   10.94674        7.390966
## 1049      11.33382       6.763551   10.89172        7.804735
## 1058      11.16527       6.340149   11.20300        8.222093
## 1060      11.86027       5.896578   11.26415        6.215755
```

```r
cormx4<-cor(data4)

data5<-t(v$E[htspt.mirna5,])
head(data5)
```

```
##      ssc-let-7d-5p ssc-let-7g ssc-miR-1468 ssc-miR-95 ssc-miR-9843-3p
## 1034      11.63362   13.39516     8.860415   11.44954        7.893082
## 1036      11.10306   13.15001     7.920425   10.90727        7.694085
## 1041      11.28869   13.24917     8.282074   10.94674        7.390966
## 1049      11.33382   13.15911     8.066415   10.89172        7.804735
## 1058      11.16527   12.88285     8.728817   11.20300        8.222093
## 1060      11.86027   13.64016     8.511007   11.26415        6.215755
```

```r
cormx5<-cor(data5)
```

Function to perform correlation analysis of miRNA to miRNA espression and export a summary of correlation analysis


```r
cor.exp <- function(var1, var2, data, ...){
	x <- cor.test(data[,as.character(var1)], data[,as.character(var2)], ...)
	rst <- data.frame(var1, var2, cor=x$estimate, 
		t=x$statistic, 
		pvalue=x$p.value)
	return(rst)
}
```

Correlation test of 4 miRNAs: 


```r
mir4.htspt<-list()
for(i in 1:length(colnames(data4))){
mir4.htspt[[i]]<-data.frame(var1=colnames(data4)[i], var2=colnames(data4)[-i])
}

vars4<-do.call(rbind, mir4.htspt)

mir.cor4<- do.call(rbind, mapply(cor.exp, var1=vars4[,1], var2=vars4[,2], MoreArgs=c(list(data=data4), alternative="two.sided", method="p"), SIMPLIFY=FALSE))
```

Perform multiple test correction (FDR) on correlation analysis:


```r
mir.cor4$qval<-qvalue(mir.cor4$pvalue)$qvalue
mir.cor4$pi0<-qvalue(mir.cor4$pvalue)$pi0
mir.cor4
```

```
##                  var1            var2          cor           t
## cor     ssc-let-7d-5p  ssc-miR-345-3p  0.081568516  1.07333770
## cor1    ssc-let-7d-5p      ssc-miR-95  0.283173491  3.87228309
## cor2    ssc-let-7d-5p ssc-miR-9843-3p -0.436356740 -6.36022669
## cor3   ssc-miR-345-3p   ssc-let-7d-5p  0.081568516  1.07333770
## cor4   ssc-miR-345-3p      ssc-miR-95 -0.045211686 -0.59355266
## cor5   ssc-miR-345-3p ssc-miR-9843-3p  0.453703788  6.67704916
## cor6       ssc-miR-95   ssc-let-7d-5p  0.283173491  3.87228309
## cor7       ssc-miR-95  ssc-miR-345-3p -0.045211686 -0.59355266
## cor8       ssc-miR-95 ssc-miR-9843-3p  0.001440106  0.01888683
## cor9  ssc-miR-9843-3p   ssc-let-7d-5p -0.436356740 -6.36022669
## cor10 ssc-miR-9843-3p  ssc-miR-345-3p  0.453703788  6.67704916
## cor11 ssc-miR-9843-3p      ssc-miR-95  0.001440106  0.01888683
##             pvalue         qval pi0
## cor   2.846232e-01 4.269348e-01   1
## cor1  1.529773e-04 3.059546e-04   1
## cor2  1.756288e-09 5.268864e-09   1
## cor3  2.846232e-01 4.269348e-01   1
## cor4  5.535911e-01 6.643093e-01   1
## cor5  3.231275e-10 1.938765e-09   1
## cor6  1.529773e-04 3.059546e-04   1
## cor7  5.535911e-01 6.643093e-01   1
## cor8  9.849533e-01 9.849533e-01   1
## cor9  1.756288e-09 5.268864e-09   1
## cor10 3.231275e-10 1.938765e-09   1
## cor11 9.849533e-01 9.849533e-01   1
```

Significantly correlated miRNA pairs:


```r
thres<-0.05
mir.cor4[mir.cor4$qval<thres,]
```

```
##                  var1            var2        cor         t       pvalue
## cor1    ssc-let-7d-5p      ssc-miR-95  0.2831735  3.872283 1.529773e-04
## cor2    ssc-let-7d-5p ssc-miR-9843-3p -0.4363567 -6.360227 1.756288e-09
## cor5   ssc-miR-345-3p ssc-miR-9843-3p  0.4537038  6.677049 3.231275e-10
## cor6       ssc-miR-95   ssc-let-7d-5p  0.2831735  3.872283 1.529773e-04
## cor9  ssc-miR-9843-3p   ssc-let-7d-5p -0.4363567 -6.360227 1.756288e-09
## cor10 ssc-miR-9843-3p  ssc-miR-345-3p  0.4537038  6.677049 3.231275e-10
##               qval pi0
## cor1  3.059546e-04   1
## cor2  5.268864e-09   1
## cor5  1.938765e-09   1
## cor6  3.059546e-04   1
## cor9  5.268864e-09   1
## cor10 1.938765e-09   1
```

---

Correlation test of 5 miRNAs:


```r
mir5.htspt<-list()
for(i in 1:length(colnames(data5))){
mir5.htspt[[i]]<-data.frame(var1=colnames(data5)[i], var2=colnames(data5)[-i])
}

vars5<-do.call(rbind, mir5.htspt)

mir.cor5<- do.call(rbind, mapply(cor.exp, var1=vars5[,1], var2=vars5[,2], MoreArgs=c(list(data=data5), alternative="two.sided", method="p"), SIMPLIFY=FALSE))
```

Perform multiple test correction (FDR) on correlation analysis:


```r
mir.cor5$qval<-qvalue(mir.cor5$pvalue)$qvalue
mir.cor5$pi0<-qvalue(mir.cor5$pvalue)$pi0
mir.cor5
```

```
##                  var1            var2          cor           t
## cor     ssc-let-7d-5p      ssc-let-7g  0.853618417 21.49189018
## cor1    ssc-let-7d-5p    ssc-miR-1468  0.742557670 14.53988103
## cor2    ssc-let-7d-5p      ssc-miR-95  0.283173491  3.87228309
## cor3    ssc-let-7d-5p ssc-miR-9843-3p -0.436356740 -6.36022669
## cor4       ssc-let-7g   ssc-let-7d-5p  0.853618417 21.49189018
## cor5       ssc-let-7g    ssc-miR-1468  0.624759973 10.49369013
## cor6       ssc-let-7g      ssc-miR-95  0.301946906  4.15387996
## cor7       ssc-let-7g ssc-miR-9843-3p -0.353237197 -4.95189209
## cor8     ssc-miR-1468   ssc-let-7d-5p  0.742557670 14.53988103
## cor9     ssc-miR-1468      ssc-let-7g  0.624759973 10.49369013
## cor10    ssc-miR-1468      ssc-miR-95  0.424386351  6.14675903
## cor11    ssc-miR-1468 ssc-miR-9843-3p -0.291866596 -4.00204752
## cor12      ssc-miR-95   ssc-let-7d-5p  0.283173491  3.87228309
## cor13      ssc-miR-95      ssc-let-7g  0.301946906  4.15387996
## cor14      ssc-miR-95    ssc-miR-1468  0.424386351  6.14675903
## cor15      ssc-miR-95 ssc-miR-9843-3p  0.001440106  0.01888683
## cor16 ssc-miR-9843-3p   ssc-let-7d-5p -0.436356740 -6.36022669
## cor17 ssc-miR-9843-3p      ssc-let-7g -0.353237197 -4.95189209
## cor18 ssc-miR-9843-3p    ssc-miR-1468 -0.291866596 -4.00204752
## cor19 ssc-miR-9843-3p      ssc-miR-95  0.001440106  0.01888683
##             pvalue         qval pi0
## cor   0.000000e+00 0.000000e+00   1
## cor1  0.000000e+00 0.000000e+00   1
## cor2  1.529773e-04 1.699748e-04   1
## cor3  1.756288e-09 4.390720e-09   1
## cor4  0.000000e+00 0.000000e+00   1
## cor5  0.000000e+00 0.000000e+00   1
## cor6  5.140065e-05 7.342951e-05   1
## cor7  1.744514e-06 2.907524e-06   1
## cor8  0.000000e+00 0.000000e+00   1
## cor9  0.000000e+00 0.000000e+00   1
## cor10 5.349823e-09 1.069965e-08   1
## cor11 9.319294e-05 1.164912e-04   1
## cor12 1.529773e-04 1.699748e-04   1
## cor13 5.140065e-05 7.342951e-05   1
## cor14 5.349823e-09 1.069965e-08   1
## cor15 9.849533e-01 9.849533e-01   1
## cor16 1.756288e-09 4.390720e-09   1
## cor17 1.744514e-06 2.907524e-06   1
## cor18 9.319294e-05 1.164912e-04   1
## cor19 9.849533e-01 9.849533e-01   1
```

Significantly correlated miRNA pairs:  


```r
mir.cor5[mir.cor5$qval<thres,]
```

```
##                  var1            var2        cor         t       pvalue
## cor     ssc-let-7d-5p      ssc-let-7g  0.8536184 21.491890 0.000000e+00
## cor1    ssc-let-7d-5p    ssc-miR-1468  0.7425577 14.539881 0.000000e+00
## cor2    ssc-let-7d-5p      ssc-miR-95  0.2831735  3.872283 1.529773e-04
## cor3    ssc-let-7d-5p ssc-miR-9843-3p -0.4363567 -6.360227 1.756288e-09
## cor4       ssc-let-7g   ssc-let-7d-5p  0.8536184 21.491890 0.000000e+00
## cor5       ssc-let-7g    ssc-miR-1468  0.6247600 10.493690 0.000000e+00
## cor6       ssc-let-7g      ssc-miR-95  0.3019469  4.153880 5.140065e-05
## cor7       ssc-let-7g ssc-miR-9843-3p -0.3532372 -4.951892 1.744514e-06
## cor8     ssc-miR-1468   ssc-let-7d-5p  0.7425577 14.539881 0.000000e+00
## cor9     ssc-miR-1468      ssc-let-7g  0.6247600 10.493690 0.000000e+00
## cor10    ssc-miR-1468      ssc-miR-95  0.4243864  6.146759 5.349823e-09
## cor11    ssc-miR-1468 ssc-miR-9843-3p -0.2918666 -4.002048 9.319294e-05
## cor12      ssc-miR-95   ssc-let-7d-5p  0.2831735  3.872283 1.529773e-04
## cor13      ssc-miR-95      ssc-let-7g  0.3019469  4.153880 5.140065e-05
## cor14      ssc-miR-95    ssc-miR-1468  0.4243864  6.146759 5.349823e-09
## cor16 ssc-miR-9843-3p   ssc-let-7d-5p -0.4363567 -6.360227 1.756288e-09
## cor17 ssc-miR-9843-3p      ssc-let-7g -0.3532372 -4.951892 1.744514e-06
## cor18 ssc-miR-9843-3p    ssc-miR-1468 -0.2918666 -4.002048 9.319294e-05
##               qval pi0
## cor   0.000000e+00   1
## cor1  0.000000e+00   1
## cor2  1.699748e-04   1
## cor3  4.390720e-09   1
## cor4  0.000000e+00   1
## cor5  0.000000e+00   1
## cor6  7.342951e-05   1
## cor7  2.907524e-06   1
## cor8  0.000000e+00   1
## cor9  0.000000e+00   1
## cor10 1.069965e-08   1
## cor11 1.164912e-04   1
## cor12 1.699748e-04   1
## cor13 7.342951e-05   1
## cor14 1.069965e-08   1
## cor16 4.390720e-09   1
## cor17 2.907524e-06   1
## cor18 1.164912e-04   1
```

---

Investigate correlation of genotypes between the 4 putative eQTL hotspots:


```r
snp.htspt<- c("H3GA0052416", "MARC0027291", "MARC0047188", "MARC0093624")
```

Examine allele frequencies of the hotspot SNPs

Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:


```r
geno_f2<-MSUPRP_meat$geno[!((as.numeric(rownames(MSUPRP_meat$geno))<=1000) | (as.numeric(rownames(MSUPRP_meat$geno))>=6000)),]
dim(geno_f2)
```

```
## [1]   940 45331
```

Subset the hotspot SNPs:


```r
geno.htspt<-geno_f2[,snp.htspt]
dim(geno.htspt)
```

```
## [1] 940   4
```

Subset the 174 animals:


```r
geno.htspt<-geno.htspt[rownames(data4),]
dim(geno.htspt)
```

```
## [1] 174   4
```

```r
head(geno.htspt)
```

```
##      H3GA0052416 MARC0027291 MARC0047188 MARC0093624
## 1034           0           2           0           2
## 1036           0           2           0           2
## 1041           0           2           0           2
## 1049           0           2           0           2
## 1058           0           2           0           2
## 1060           1           1           1           1
```

```r
if(sum(rownames(geno.htspt)!=rownames(data4)) != 0) stop ("Animal IDs not correct for genotype object")

geno.cor<-cor(geno.htspt)
```


Calculate allele frequency for all the F2 animals:


```r
allele_freq<-colMeans(geno.htspt,na.rm=T)/2
allele_freq
```

```
## H3GA0052416 MARC0027291 MARC0047188 MARC0093624 
##   0.1293103   0.8735632   0.1293103   0.8965517
```

Minor allele frequency:


```r
maf<-ifelse(allele_freq>0.5, (1- allele_freq), allele_freq)
maf
```

```
## H3GA0052416 MARC0027291 MARC0047188 MARC0093624 
##   0.1293103   0.1264368   0.1293103   0.1034483
```

---

Identify the associated miRNAs in the annotation file, to query miRBase and the literature for information on their effects 


```r
head(total.mature.annot2)
```

```
##                        Name  chr0     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
## ssc-let-7f       ssc-let-7f  chr3  44864810  44864831    22      + miRNA
##                      Alias          Precursors
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085
## ssc-let-7c    MIMAT0002151           MI0002445
## ssc-let-7d-3p MIMAT0025357           MI0022120
## ssc-let-7d-5p MIMAT0025356           MI0022120
## ssc-let-7e    MIMAT0013866           MI0013086
## ssc-let-7f    MIMAT0002152 MI0022121,MI0002446
```

```r
mat.annot<-total.mature.annot2[c("ssc-let-7d-5p", "ssc-let-7g","ssc-miR-1468", "ssc-miR-345-3p", "ssc-miR-95", "ssc-miR-9843-3p"),]
mat.annot$Precursors
```

```
## [1] "MI0022120" "MI0013087" "MI0022160" "MI0013117" "MI0002436" "MI0031612"
```

```r
head(precursor.mirannot)
```

```
##            Name chr0    start      end width strand
## 1  ssc-mir-9815 chr1  1882209  1882304    96      -
## 3  ssc-mir-9821 chr1  2600232  2600312    81      -
## 5  ssc-mir-9824 chr1  2816327  2816402    76      -
## 7  ssc-mir-9831 chr1  4108414  4108490    77      +
## 9  ssc-mir-9857 chr1 12666392 12666461    70      -
## 11 ssc-mir-9787 chr1 12884779 12884869    91      -
##                        type     Alias Derives_from
## 1  miRNA_primary_transcript MI0031582         <NA>
## 3  miRNA_primary_transcript MI0031588         <NA>
## 5  miRNA_primary_transcript MI0031591         <NA>
## 7  miRNA_primary_transcript MI0031600         <NA>
## 9  miRNA_primary_transcript MI0031635         <NA>
## 11 miRNA_primary_transcript MI0031548         <NA>
```

```r
prec.annot<-precursor.mirannot[precursor.mirannot$Alias %in% mat.annot$Precursors,]

if(sum(!as.character(prec.annot$Alias) %in% mat.annot$Precursors) != 0) stop ("Precursor IDs do not match between datasets")

prec.annot
```

```
##              Name  chr0     start       end width strand
## 191    ssc-let-7g chr13  37599004  37599083    80      +
## 437    ssc-let-7d  chr3  44867267  44867362    96      +
## 610 ssc-mir-345-1  chr7 128658300 128658379    80      +
## 625    ssc-mir-95  chr8   4277238   4277318    81      +
## 638  ssc-mir-9843  chr8 122371904 122371984    81      -
## 710  ssc-mir-1468  chrX  56757045  56757127    83      -
##                         type     Alias Derives_from
## 191 miRNA_primary_transcript MI0013087         <NA>
## 437 miRNA_primary_transcript MI0022120         <NA>
## 610 miRNA_primary_transcript MI0013117         <NA>
## 625 miRNA_primary_transcript MI0002436         <NA>
## 638 miRNA_primary_transcript MI0031612         <NA>
## 710 miRNA_primary_transcript MI0022160         <NA>
```

Information extracted from target prediction input:

miRNA           | Seed Sequence (Human)
--------------- | -------------
let-7-5p/98-5p  | GAGGUAG
miR-1468-5p	   | UCCGUUU
miR-345-3p	   | CCCUGAA
miR-95-3p	   | UCAACGG

Also see `miRNA_targets_common-1.tiff` to compare the targets in common expressed in this dataset between miRNAs.

## Visualize

Correlation plots of 4 and 5 hotspot-associated miRNAs:

4 miRNAs:


```r
corrplot.mixed(cormx4, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Four hotspot-associated miRNAs")
```

![plot of chunk miRNA4_hotspots](figure/miRNA4_hotspots-1.tiff)

5 miRNAs:


```r
corrplot.mixed(cormx5, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Five hotspot-associated miRNAs")
```

![plot of chunk miRNA5_hotspots](figure/miRNA5_hotspots-1.tiff)

```r
corrplot.mixed(geno.cor, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Putative Hotspot SNP Genotype Correlation")
```

![plot of chunk geno_cor_hotspots](figure/geno_cor_hotspots-1.tiff)

## Save data


```r
save(mir.cor4, mir.cor5, geno.cor, maf, file="../6_hotspot_miRNA_correlation_results.Rdata")
```

