**Script:** `2_miRNA_eqtl_summary.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`

**Date:**  2/25/16 UPDATED 6/29/16

**Input File Directory:**
1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`

2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`

**Input File(s):**

1. `1_gblup_results_summary.Rdata`

2. `2_gwa_results_summary.Rdata`

3. `3_msuprp_mirna_gpdata_pheno_counts.Rdata`

4. `4_normalized_dge_object_and_voom_output.Rdata`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`

**Output File(s):** `2_miRNA_eqtl_summary.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to summarize the number of eQTL peaks per miRNA output from the first eQTL analysis of the 174 F2 MSUPRP pig miRNA expression profiles.

## Install libraries



```r
rm(list=ls())
```

Load eqtl function Rdata containing stb function, which will summarize the eQTL results


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
ls()
```

```
##  [1] "absmap"     "add_legend" "AddPosGene" "distance"   "inrange"   
##  [6] "manhpt"     "peakrng"    "plot.GMA"   "sigpval"    "stb"       
## [11] "stb.nm"     "tbpos"      "zstandard"
```

```r
library(limma)
library(edgeR)
library(gwaR)
library(plyr)

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")
```

## Load data

Load the gblup and eQTL output:


```r
load("../1_gblup_results_summary.Rdata")

load("../2_gwa_results_summary.Rdata")
```

Load the dge object to obtain the mature miRNA annotation:


```r
load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")
```

Load the MSUPRP_miRNA gpdata object for the mapping information:


```r
load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
ls()
```

```
##  [1] "absmap"               "add_legend"           "AddPosGene"          
##  [4] "dge"                  "distance"             "inrange"             
##  [7] "manhpt"               "MSUPRP_miRNA"         "peakrng"             
## [10] "plot.GMA"             "rst.gwa"              "sigpval"             
## [13] "stb"                  "stb.nm"               "summary_MSUPRP_miRNA"
## [16] "summary.rst.gblup"    "tbpos"                "v"                   
## [19] "wtcen"                "zstandard"
```

## Analysis

### Summarize the heritability of the miRNAs, output from GBLUP:

The average heritability of all the miRNA expression profiles:


```r
mean(summary.rst.gblup$h2)
```

```
## [1] 0.1201534
```

```r
summary(summary.rst.gblup$h2)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0000000 0.0000731 0.0858800 0.1202000 0.1863000 0.6324000
```

The average heritability of the miRNAs found to be significantly heritable at FDR < 0.05:

How many significantly heritable miRNAs found at FDR < 0.05? (Should match from eQTL scan output md file)


```r
sum(summary.rst.gblup$qvalue<0.05)
```

```
## [1] 50
```

Extract the significantly heritable miRNAs from the summary.rst.gblup dataset and calculate mean h2:


```r
sigh2<-summary.rst.gblup[summary.rst.gblup$qvalue<0.05,]
dim(sigh2)
```

```
## [1] 50 10
```

```r
mean(sigh2$h2)
```

```
## [1] 0.3290721
```

```r
summary(sigh2$h2)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1752  0.2371  0.3234  0.3291  0.3650  0.6324
```

Define the minimum p-value that is not significant (based on q-value < 0.05) as the threshold for plotting significant points


```r
summary(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.008773 0.089420 0.247300 0.271400 0.500000 0.500000
```

```r
sigthresh<-min(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
sigthresh
```

```
## [1] 0.00877268
```

Plot h2 vs -log10(p-value) like before to determine trend in significance and h2:


```r
plot(summary.rst.gblup$h2, -log10(summary.rst.gblup$lrtpvalue),
        xlab = expression("Heritability"~(h^{2})),
    ylab = "-log10(p-value)",
    main = "Significance vs Heritability")
points(summary.rst.gblup$h2[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       -log10(summary.rst.gblup$lrtpvalue)[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       pch=19,col="red")
abline(a = -log10(sigthresh), b = 0, lty = 5)
```

![plot of chunk vol_sig_h2](figure/vol_sig_h2-1.tiff)

---

### The Summary Table of GWAS Results:


Assign the correct names to the different objects:


```r
map <- MSUPRP_miRNA$map
colnames(map)
```

```
## [1] "chr" "pos"
```

```r
annotation <- dge$genes
head(annotation)
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
colnames(annotation)
```

```
## [1] "Name"       "chr0"       "start"      "end"        "width"     
## [6] "strand"     "type"       "Alias"      "Precursors"
```

Annotation must contain columns "chr", "start", "end", and "strand"


```r
colnames(annotation)<-c("Name","chr","start","end","width","strand","type","Alias","Precursors")
annotation$mid.mir<-round(rowMeans(annotation[,c("start", "end")], na.rm=TRUE))
head(annotation)
```

```
##                        Name   chr     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
## ssc-let-7f       ssc-let-7f  chr3  44864810  44864831    22      + miRNA
##                      Alias          Precursors   mid.mir
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085  44864454
## ssc-let-7c    MIMAT0002151           MI0002445 191559362
## ssc-let-7d-3p MIMAT0025357           MI0022120  44867342
## ssc-let-7d-5p MIMAT0025356           MI0022120  44867288
## ssc-let-7e    MIMAT0013866           MI0013086  51858396
## ssc-let-7f    MIMAT0002152 MI0022121,MI0002446  44864820
```

Add the mid.mir to the map object for later use in manhattan plots (arrow where midpoint of transcript lies)

First, substitute the chromosome number for the "chr" column of the annotation file and to add to the chr column of the map:


```r
annotation$chr<-gsub("X", "19", annotation$chr)
annotation$chr<-as.numeric(gsub("chr", "", annotation$chr))
table(annotation$chr)
```

```
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 
## 25 25 13  9 10 19 17  4 12  7  8 23 14 10  3  4  8 14 27
```

```r
head(annotation)
```

```
##                        Name chr     start       end width strand  type
## ssc-let-7a       ssc-let-7a   3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c  13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p   3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p   3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e   6  51858385  51858406    22      + miRNA
## ssc-let-7f       ssc-let-7f   3  44864810  44864831    22      + miRNA
##                      Alias          Precursors   mid.mir
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085  44864454
## ssc-let-7c    MIMAT0002151           MI0002445 191559362
## ssc-let-7d-3p MIMAT0025357           MI0022120  44867342
## ssc-let-7d-5p MIMAT0025356           MI0022120  44867288
## ssc-let-7e    MIMAT0013866           MI0013086  51858396
## ssc-let-7f    MIMAT0002152 MI0022121,MI0002446  44864820
```

Extract the chromosome and the position of the miRNA and build a data.frame to add to the map object later:


```r
mirpos<-data.frame(chr=annotation$chr,
	pos=annotation$mid.mir, row.names=annotation$Name)
str(mirpos)
```

```
## 'data.frame':	295 obs. of  2 variables:
##  $ chr: num  3 13 3 3 6 3 13 5 17 9 ...
##  $ pos: num  4.49e+07 1.92e+08 4.49e+07 4.49e+07 5.19e+07 ...
```

```r
head(mirpos)
```

```
##               chr       pos
## ssc-let-7a      3  44864454
## ssc-let-7c     13 191559362
## ssc-let-7d-3p   3  44867342
## ssc-let-7d-5p   3  44867288
## ssc-let-7e      6  51858396
## ssc-let-7f      3  44864820
```

```r
dim(mirpos)
```

```
## [1] 295   2
```

```r
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
```

Create the summary table, using pergene = TRUE to get gene-wise eQTL peaks:


```r
sum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=TRUE)
dim(sum.eqtl)
```

```
## [1] 24 13
```

```r
head(sum.eqtl)
```

```
##             miRNA chr.miR start.miR  end.miR  mid.miR strand   miRBase.ID
## 1   ssc-let-7d-5p       3  44867277 44867298 44867288      + MIMAT0025356
## 2      ssc-let-7g      13  37599007 37599028 37599018      + MIMAT0013867
## 3     ssc-miR-128      13  22859239 22859259 22859249      + MIMAT0002157
## 4 ssc-miR-1306-3p      14  55111005 55111025 55111015      + MIMAT0013938
## 5  ssc-miR-140-5p      NA        NA       NA      NaN   <NA>         <NA>
## 6  ssc-miR-140-5p      NA        NA       NA      NaN   <NA>         <NA>
##            precursors         SNP chr.snp   pos.snp snp.sign       qvalue
## 1           MI0022120 MARC0093624      15 135538479        - 0.0117673699
## 2           MI0013087 MARC0093624      15 135538479        - 0.0117673699
## 3 MI0013094,MI0002451 ALGA0023517       4  15959123        + 0.0455542195
## 4           MI0013148 H3GA0034702      12  54736413        + 0.0065172015
## 5                <NA> ASGA0017748       4   7002305        - 0.0109425100
## 6                <NA> H3GA0055369       6  15626776        + 0.0003022932
```

How many unique miRNAs have eQTL?


```r
length(unique(sum.eqtl$miRNA))
```

```
## [1] 17
```

miR-eQTL peaks at FDR<0.05:


```r
sum.eqtl
```

```
##              miRNA chr.miR start.miR   end.miR   mid.miR strand
## 1    ssc-let-7d-5p       3  44867277  44867298  44867288      +
## 2       ssc-let-7g      13  37599007  37599028  37599018      +
## 3      ssc-miR-128      13  22859239  22859259  22859249      +
## 4  ssc-miR-1306-3p      14  55111005  55111025  55111015      +
## 5   ssc-miR-140-5p      NA        NA        NA       NaN   <NA>
## 6   ssc-miR-140-5p      NA        NA        NA       NaN   <NA>
## 7     ssc-miR-1468      19  56757096  56757117  56757106      -
## 8      ssc-miR-184       7  53883677  53883698  53883688      +
## 9      ssc-miR-184       7  53883677  53883698  53883688      +
## 10     ssc-miR-184       7  53883677  53883698  53883688      +
## 11    ssc-miR-190b       4 104457820 104457841 104457830      +
## 12  ssc-miR-345-3p       7 128658348 128658368 128658358      +
## 13     ssc-miR-429       6  58044184  58044205  58044194      -
## 14     ssc-miR-429       6  58044184  58044205  58044194      -
## 15 ssc-miR-6782-3p       6    956806    956827    956816      +
## 16 ssc-miR-7135-3p       3  29052303  29052324  29052314      +
## 17     ssc-miR-874       2 145381108 145381130 145381119      -
## 18     ssc-miR-874       2 145381108 145381130 145381119      -
## 19      ssc-miR-95       8   4277286   4277307   4277296      +
## 20 ssc-miR-9785-5p       3  22233456  22233479  22233468      -
## 21 ssc-miR-9785-5p       3  22233456  22233479  22233468      -
## 22 ssc-miR-9785-5p       3  22233456  22233479  22233468      -
## 23 ssc-miR-9810-3p       4  90746075  90746095  90746085      +
## 24 ssc-miR-9843-3p       8 122371914 122371935 122371924      -
##      miRBase.ID          precursors         SNP chr.snp   pos.snp snp.sign
## 1  MIMAT0025356           MI0022120 MARC0093624      15 135538479        -
## 2  MIMAT0013867           MI0013087 MARC0093624      15 135538479        -
## 3  MIMAT0002157 MI0013094,MI0002451 ALGA0023517       4  15959123        +
## 4  MIMAT0013938           MI0013148 H3GA0034702      12  54736413        +
## 5          <NA>                <NA> ASGA0017748       4   7002305        -
## 6          <NA>                <NA> H3GA0055369       6  15626776        +
## 7  MIMAT0025386           MI0022160 MARC0093624      15 135538479        -
## 8  MIMAT0002127           MI0002421 DBWU0000430       3   9121844        -
## 9  MIMAT0002127           MI0002421 M1GA0026172       6 157003423        +
## 10 MIMAT0002127           MI0002421 ALGA0041952       7  55936003        +
## 11 MIMAT0020588           MI0017988 ALGA0026452       4  94724175        -
## 12 MIMAT0013900           MI0013117 MARC0027291      15 135171935        -
## 13 MIMAT0020591           MI0017991 ASGA0094554       6  41809035        -
## 14 MIMAT0020591           MI0017991 ALGA0046283       8   6551222        +
## 15 MIMAT0037064           MI0031620 ASGA0094215      10  29227834        +
## 16 MIMAT0028146           MI0023568 MARC0056802       3  28401592        -
## 17 MIMAT0025384           MI0022157 ALGA0016550       2 145462524        +
## 18 MIMAT0025384           MI0022157 ALGA0122273       3  64022095        -
## 19 MIMAT0002142           MI0002436 MARC0093624      15 135538479        -
## 20 MIMAT0037002           MI0031545 DRGA0003812       3  21432908        -
## 21 MIMAT0037002           MI0031545 MARC0081878      10    223442        -
## 22 MIMAT0037002           MI0031545 ALGA0121561      17   3681417        +
## 23 MIMAT0037028           MI0031577 MARC0021620       5  16689918        +
## 24 MIMAT0037061           MI0031612 MARC0093624      15 135538479        -
##          qvalue
## 1  1.176737e-02
## 2  1.176737e-02
## 3  4.555422e-02
## 4  6.517202e-03
## 5  1.094251e-02
## 6  3.022932e-04
## 7  1.176737e-02
## 8  3.806704e-02
## 9  1.461839e-02
## 10 1.645128e-07
## 11 1.898190e-02
## 12 2.188055e-02
## 13 5.415319e-06
## 14 3.016255e-02
## 15 1.870050e-04
## 16 1.150891e-02
## 17 2.900787e-09
## 18 2.116294e-02
## 19 1.176737e-02
## 20 3.557354e-02
## 21 3.557354e-02
## 22 3.557354e-02
## 23 3.627222e-02
## 24 1.176737e-02
```

#### Summary of GWAS results at FDR < 0.05

Number of eQTL peaks per chromosome:


```r
table(sum.eqtl$chr.snp)
```

```
## 
##  2  3  4  5  6  7  8 10 12 15 17 
##  1  4  3  1  3  1  1  2  1  6  1
```

Names of associated miRNAs:


```r
table(sum.eqtl$miRNA)
```

```
## 
##   ssc-let-7d-5p      ssc-let-7g     ssc-miR-128 ssc-miR-1306-3p 
##               1               1               1               1 
##  ssc-miR-140-5p    ssc-miR-1468     ssc-miR-184    ssc-miR-190b 
##               2               1               3               1 
##  ssc-miR-345-3p     ssc-miR-429 ssc-miR-6782-3p ssc-miR-7135-3p 
##               1               2               1               1 
##     ssc-miR-874      ssc-miR-95 ssc-miR-9785-5p ssc-miR-9810-3p 
##               2               1               3               1 
## ssc-miR-9843-3p 
##               1
```

Chromosomes of associated miRNAs:


```r
table(sum.eqtl$chr.miR)
```

```
## 
##  2  3  4  6  7  8 13 14 19 
##  2  5  2  3  4  2  2  1  1
```

Names of associated peak markers:


```r
table(as.character(sum.eqtl$SNP))
```

```
## 
## ALGA0016550 ALGA0023517 ALGA0026452 ALGA0041952 ALGA0046283 ALGA0121561 
##           1           1           1           1           1           1 
## ALGA0122273 ASGA0017748 ASGA0094215 ASGA0094554 DBWU0000430 DRGA0003812 
##           1           1           1           1           1           1 
## H3GA0034702 H3GA0055369 M1GA0026172 MARC0021620 MARC0027291 MARC0056802 
##           1           1           1           1           1           1 
## MARC0081878 MARC0093624 
##           1           5
```

---

### Determining the ranges of associated SNPs per eQTL peak on SSC15 (for ISAG abstract):

First, create the summary table at FDR 5% again, this time with pergene=F to identify all markers associated with each eQTL peak:


```r
fullsum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=FALSE)
dim(fullsum.eqtl)
```

```
## [1] 338  13
```

Summarize the number of SNPs associated with each miRNA eQTL (some have multiple peaks)


```r
numsnps<-by(fullsum.eqtl, as.character(fullsum.eqtl$miRNA), nrow)
numsnps<-ldply(numsnps, fun=NULL, id=names(numsnps))
colnames(numsnps) <- c("miRNA", "numsnps")
numsnps
```

```
##              miRNA numsnps
## 1    ssc-let-7d-5p       4
## 2       ssc-let-7g       1
## 3      ssc-miR-128       1
## 4  ssc-miR-1306-3p       1
## 5   ssc-miR-140-5p       4
## 6     ssc-miR-1468       1
## 7      ssc-miR-184      53
## 8     ssc-miR-190b      11
## 9   ssc-miR-345-3p       4
## 10     ssc-miR-429      94
## 11 ssc-miR-6782-3p       4
## 12 ssc-miR-7135-3p      15
## 13     ssc-miR-874     117
## 14      ssc-miR-95       4
## 15 ssc-miR-9785-5p      18
## 16 ssc-miR-9810-3p       2
## 17 ssc-miR-9843-3p       4
```

```r
sum(numsnps$numsnps)
```

```
## [1] 338
```

---

### Extract peak range data from all miRNA eQTL peaks

I can obtain information on the range of each peak based on miRNA (adapted from Deborah's function "peakrng"):


```r
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
```

nmt = name transcript (in this case, miRNA); so, make a list of the miRNAs and loop through to get the peak range information for each miRNA with significant eQTL peaks
sumtb = output from the summary table function, with pergene = FALSE


```r
sigmirnames <- unique(as.character(sum.eqtl$miRNA))
sigmirnames
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-140-5p"  "ssc-miR-1468"   
##  [7] "ssc-miR-184"     "ssc-miR-190b"    "ssc-miR-345-3p" 
## [10] "ssc-miR-429"     "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [13] "ssc-miR-874"     "ssc-miR-95"      "ssc-miR-9785-5p"
## [16] "ssc-miR-9810-3p" "ssc-miR-9843-3p"
```

```r
mirpeaks<-data.frame(do.call(rbind, lapply(sigmirnames, peakrngmir, fullsum.eqtl)), sum.eqtl[,c("SNP", "pos.snp")])
```

```
## Number of markers per chromosomal peak for ssc-let-7d-5p :
## 15 
##  4 
## Number of markers per chromosomal peak for ssc-let-7g :
## 15 
##  1 
## Number of markers per chromosomal peak for ssc-miR-128 :
## 4 
## 1 
## Number of markers per chromosomal peak for ssc-miR-1306-3p :
## 12 
##  1 
## Number of markers per chromosomal peak for ssc-miR-140-5p :
## 4 6 
## 1 3 
## Number of markers per chromosomal peak for ssc-miR-1468 :
## 15 
##  1 
## Number of markers per chromosomal peak for ssc-miR-184 :
##  3  6  7 
##  2  1 50 
## Number of markers per chromosomal peak for ssc-miR-190b :
##  4 
## 11 
## Number of markers per chromosomal peak for ssc-miR-345-3p :
## 15 
##  4 
## Number of markers per chromosomal peak for ssc-miR-429 :
##  6  8 
## 93  1 
## Number of markers per chromosomal peak for ssc-miR-6782-3p :
## 10 
##  4 
## Number of markers per chromosomal peak for ssc-miR-7135-3p :
##  3 
## 15 
## Number of markers per chromosomal peak for ssc-miR-874 :
##   2   3 
## 116   1 
## Number of markers per chromosomal peak for ssc-miR-95 :
## 15 
##  4 
## Number of markers per chromosomal peak for ssc-miR-9785-5p :
##  3 10 17 
## 16  1  1 
## Number of markers per chromosomal peak for ssc-miR-9810-3p :
## 5 
## 2 
## Number of markers per chromosomal peak for ssc-miR-9843-3p :
## 15 
##  4
```

```r
mirpeaks
```

```
##              miRNA chr.miR start.miR   end.miR range.miR   miRBase.ID
## 1    ssc-let-7d-5p       3  44867277  44867298        21 MIMAT0025356
## 2       ssc-let-7g      13  37599007  37599028        21 MIMAT0013867
## 3      ssc-miR-128      13  22859239  22859259        20 MIMAT0002157
## 4  ssc-miR-1306-3p      14  55111005  55111025        20 MIMAT0013938
## 5   ssc-miR-140-5p      NA        NA        NA        NA         <NA>
## 6   ssc-miR-140-5p      NA        NA        NA        NA         <NA>
## 7     ssc-miR-1468      19  56757096  56757117        21 MIMAT0025386
## 8      ssc-miR-184       7  53883677  53883698        21 MIMAT0002127
## 9      ssc-miR-184       7  53883677  53883698        21 MIMAT0002127
## 10     ssc-miR-184       7  53883677  53883698        21 MIMAT0002127
## 11    ssc-miR-190b       4 104457820 104457841        21 MIMAT0020588
## 12  ssc-miR-345-3p       7 128658348 128658368        20 MIMAT0013900
## 13     ssc-miR-429       6  58044184  58044205        21 MIMAT0020591
## 14     ssc-miR-429       6  58044184  58044205        21 MIMAT0020591
## 15 ssc-miR-6782-3p       6    956806    956827        21 MIMAT0037064
## 16 ssc-miR-7135-3p       3  29052303  29052324        21 MIMAT0028146
## 17     ssc-miR-874       2 145381108 145381130        22 MIMAT0025384
## 18     ssc-miR-874       2 145381108 145381130        22 MIMAT0025384
## 19      ssc-miR-95       8   4277286   4277307        21 MIMAT0002142
## 20 ssc-miR-9785-5p       3  22233456  22233479        23 MIMAT0037002
## 21 ssc-miR-9785-5p       3  22233456  22233479        23 MIMAT0037002
## 22 ssc-miR-9785-5p       3  22233456  22233479        23 MIMAT0037002
## 23 ssc-miR-9810-3p       4  90746075  90746095        20 MIMAT0037028
## 24 ssc-miR-9843-3p       8 122371914 122371935        21 MIMAT0037061
##    chr.snp range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 1       15     366544 135171935 135538479       4 MARC0093624 135538479
## 2       15          0 135538479 135538479       1 MARC0093624 135538479
## 3        4          0  15959123  15959123       1 ALGA0023517  15959123
## 4       12          0  54736413  54736413       1 H3GA0034702  54736413
## 5        4          0   7002305   7002305       1 ASGA0017748   7002305
## 6        6     386244  15626776  16013020       3 H3GA0055369  15626776
## 7       15          0 135538479 135538479       1 MARC0093624 135538479
## 8        3  125903917   9121844 135025761       2 DBWU0000430   9121844
## 9        6          0 157003423 157003423       1 M1GA0026172 157003423
## 10       7   84868032  49815607 134683639      50 ALGA0041952  55936003
## 11       4   16724835  94724175 111449010      11 ALGA0026452  94724175
## 12      15    1285898 133948641 135234539       4 MARC0027291 135171935
## 13       6   22928010  38548262  61476272      93 ASGA0094554  41809035
## 14       8          0   6551222   6551222       1 ALGA0046283   6551222
## 15      10   11498076  20485599  31983675       4 ASGA0094215  29227834
## 16       3     723613  28401592  29125205      15 MARC0056802  28401592
## 17       2   23140264 132997331 156137595     116 ALGA0016550 145462524
## 18       3          0  64022095  64022095       1 ALGA0122273  64022095
## 19      15     366544 135171935 135538479       4 MARC0093624 135538479
## 20       3   14394187  20983778  35377965      16 DRGA0003812  21432908
## 21      10          0    223442    223442       1 MARC0081878    223442
## 22      17          0   3681417   3681417       1 ALGA0121561   3681417
## 23       5      86125  16689918  16776043       2 MARC0021620  16689918
## 24      15     366544 135171935 135538479       4 MARC0093624 135538479
```

---

### Creating Manhattan plots of the six miRNA with the highest numbers of associated SNP markers (for ISAG poster)

First, convert the map positions to absolute map positions using Deborah's function "absmap"

Notice how the map object's positions are relative to chromosome:


```r
head(map)
```

```
##             chr    pos
## MARC0044150   1 286933
## ASGA0000014   1 342481
## H3GA0000026   1 389876
## ASGA0000021   1 489855
## ALGA0000009   1 538161
## ALGA0000014   1 565627
```

```r
dim(map)
```

```
## [1] 38166     2
```

Add the map positions of the miRNA with eQTL


```r
map.full<-rbind(map, mirpos[sigmirnames,])
dim(map.full)
```

```
## [1] 38183     2
```

```r
if(nrow(map.full) - nrow(map) != length(sigmirnames)) stop ("miRNA map positions not added correctly")
```

Use the absmap function to convert the chromosomal map positions to absolute map positions:


```r
absposmap<-absmap(map.full)
head(absposmap)
```

```
## MARC0044150 ASGA0000014 H3GA0000026 ASGA0000021 ALGA0000009 ALGA0000014 
##      286933      342481      389876      489855      538161      565627
```

```r
tail(absposmap)
```

```
##  ALGA0098927  ALGA0098928  ASGA0080447  ASGA0080449  M1GA0023446 
##   2445228408   2445284532   2445355730   2445542182   2445560578 
## ssc-miR-1468 
##   2502317684
```

Notice it didn't include the miRNA with no map position (ssc-miR-140-5p)


```r
absposmap[which(names(absposmap) %in% sigmirnames)]
```

```
##     ssc-miR-874 ssc-miR-9785-5p ssc-miR-7135-3p   ssc-let-7d-5p 
##       460595619       499746054       506564900       522379874 
## ssc-miR-9810-3p    ssc-miR-190b ssc-miR-6782-3p     ssc-miR-429 
##       712423896       726135641       877456244       934543622 
##     ssc-miR-184  ssc-miR-345-3p      ssc-miR-95 ssc-miR-9843-3p 
##      1088117455      1162892125      1173194702      1291289330 
##     ssc-miR-128      ssc-let-7g ssc-miR-1306-3p    ssc-miR-1468 
##      1722239460      1736979229      1972702496      2502317684
```

Divide by 1e6 to get absolute position in Mb (nicer x-axis for plots)


```r
head(absposmap/1e6)
```

```
## MARC0044150 ASGA0000014 H3GA0000026 ASGA0000021 ALGA0000009 ALGA0000014 
##    0.286933    0.342481    0.389876    0.489855    0.538161    0.565627
```

Use sigpval function (from DV) to calculate the significant p-value cutoff for plotting manhattan plots for each miRNA

Then, provide the vector of miRNA names and the absposmap object to the manhpt function (also from Deborah's func_eqtl.Rdata) and loop through the vector of miRNA names to create the Manhattan plots:
## Visualize



```r
sigmirnames
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-140-5p"  "ssc-miR-1468"   
##  [7] "ssc-miR-184"     "ssc-miR-190b"    "ssc-miR-345-3p" 
## [10] "ssc-miR-429"     "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [13] "ssc-miR-874"     "ssc-miR-95"      "ssc-miR-9785-5p"
## [16] "ssc-miR-9810-3p" "ssc-miR-9843-3p"
```

```r
for(i in sigmirnames){
	pvalcutoff<-sigpval(i, pvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.pval"], qvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.qval"], fdr=0.05)
	pvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.pval"]
	names(pvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=pvals,map=map.full,annotation=NULL,pvalues=TRUE,cutoff=pvalcutoff, arrow=TRUE)

}
```

![plot of chunk man_plot_pval](figure/man_plot_pval-1.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-2.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-3.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-4.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-5.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-6.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-7.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-8.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-9.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-10.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-11.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-12.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-13.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-14.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-15.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-16.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-17.tiff)

The Manhattan Plots of q-values should be put on equal y-axes for comparison on poster.

Three of the peaks are strong enough signals to be on a y-axis of 0-10 (-log10(qval))

The remaining 14 peaks are less strong, on a y-axis of 0-4



```r
sigmirnames4<-sigmirnames[c(1:4,6,8,9,11,12,14:17)]
sigmirnames4
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-1468"    "ssc-miR-190b"   
##  [7] "ssc-miR-345-3p"  "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [10] "ssc-miR-95"      "ssc-miR-9785-5p" "ssc-miR-9810-3p"
## [13] "ssc-miR-9843-3p"
```

```r
sigmirnames8<-sigmirnames[c(5,7,10)]
sigmirnames8
```

```
## [1] "ssc-miR-140-5p" "ssc-miR-184"    "ssc-miR-429"
```

```r
sigmirnames10<-sigmirnames[13]
sigmirnames10
```

```
## [1] "ssc-miR-874"
```

```r
for(i in sigmirnames4){
	qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
	names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,4))

}
```

![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-1.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-2.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-3.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-4.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-5.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-6.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-7.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-8.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-9.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-10.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-11.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-12.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-13.tiff)

```r
for(i in sigmirnames8){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,8))

}
```

![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-1.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-2.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-3.tiff)

```r
for(i in sigmirnames10){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,10))

}
```

![plot of chunk man_plot_qval_y10](figure/man_plot_qval_y10-1.tiff)

## Save data


```r
save(sum.eqtl, fullsum.eqtl, absposmap, map.full, mirpeaks, file = "../3_eqtl_summary_tables_maps.Rdata")
```

