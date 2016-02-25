**Script:** `2_miRNA_eqtl_summary.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts`

**Date:**  2/25/16

**Input File Directory:**  
1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/`

2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`

**Input File(s):** 
 
1. `1_gblup_gwa_results_summary_full.Rdata`

2. `5_normalized_dge_object_and_voom_output.Rata`

3. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`

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

library(limma)
```

```
## Loading required package: methods
```

```r
library(edgeR)

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/4_miRNA_gblup_eqtl_analysis/scripts")
```

## Load data

Load the eQTL output:


```r
load("../1_gblup_gwa_results_summary_full.Rdata")
```

Load the dge object to obtain the mature miRNA annotation:


```r
load("../../3_build_dge_object_for_eqtl/5_normalized_dge_object_and_voom_output.Rata")
```

Load the final_MSUPRP_mirna gpdata object for the mapping information:


```r
load("../../3_build_dge_object_for_eqtl/4_msuprp_mirna_gpdata_pheno_counts.Rdata")
ls()
```

```
##  [1] "absmap"                     "add_legend"                
##  [3] "AddPosGene"                 "dge"                       
##  [5] "distance"                   "eqtl"                      
##  [7] "final_MSUPRP_miRNA"         "manhpt"                    
##  [9] "peakrng"                    "plot.GMA"                  
## [11] "stb"                        "summary_final_MSUPRP_miRNA"
## [13] "tbpos"                      "v"                         
## [15] "wtcen"                      "zstandard"
```

## Analysis

### Summarize the heritability of the miRNAs, output from GBLUP:

The average heritability of all the miRNA expression profiles:


```r
mean(eqtl$gblup$h2)
```

```
## [1] 0.118076
```

```r
summary(eqtl$gblup$h2)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0000000 0.0000722 0.0853400 0.1181000 0.1826000 0.6328000
```

The average heritability of the miRNAs found to be significantly heritable at FDR < 0.05:

How many significantly heritable miRNAs found at FDR < 0.05? (Should match from eQTL scan output md file)


```r
sum(eqtl$gblup$qvalue<0.05)
```

```
## [1] 47
```

Extract the significantly heritable miRNAs from the eqtl$gblup dataset and calculate mean h2:


```r
sigh2<-eqtl$gblup[eqtl$gblup$qvalue<0.05,]
dim(sigh2)
```

```
## [1] 47 10
```

```r
mean(sigh2$h2)
```

```
## [1] 0.3291513
```

```r
summary(sigh2$h2)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1738  0.2431  0.3243  0.3292  0.3666  0.6328
```

---

### The Summary Table of GWAS Results:


Assign the correct names to the different objects:


```r
map <- final_MSUPRP_miRNA$map
colnames(map)
```

```
## [1] "chr" "pos"
```

```r
annotation <- dge$genes
colnames(annotation)
```

```
## [1] "Name"       "chr0"       "start"      "end"        "width"     
## [6] "strand"     "type"       "Alias"      "Precursors"
```

Annotation must contain columns "chr", "start", "end", and "strand"


```r
colnames(annotation)<-c("Name","chr","start","end","width","strand","type","Alias","Precursors")
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
##                      Alias          Precursors
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085
## ssc-let-7c    MIMAT0002151           MI0002445
## ssc-let-7d-3p MIMAT0025357           MI0022120
## ssc-let-7d-5p MIMAT0025356           MI0022120
## ssc-let-7e    MIMAT0013866           MI0013086
## ssc-let-7f    MIMAT0002152 MI0022121,MI0002446
```

This function (stb) returns a summary table of the eQTL peaks per chromosome/per gene at FDR < 0.01:


```r
rsumtb1 <- stb(qval=eqtl$gwa.qval, map=map, annotation=annotation, Z=eqtl$gwa, threshold=0.01, gene.name="Precursors", pergene=T)
dim(rsumtb1)
```

```
## [1]  8 11
```

miR-eQTL peaks at FDR<0.01:


```r
rsumtb1
```

```
##              Gene chr.gene start.gene   end.gene strand gene.name
## 1      ssc-let-7g    chr13  37.599007  37.599028      + MI0013087
## 2 ssc-miR-1306-3p    chr14  55.111005  55.111025      + MI0013148
## 3  ssc-miR-140-5p     <NA>         NA         NA   <NA>      <NA>
## 4     ssc-miR-184     chr7  53.883677  53.883698      + MI0002421
## 5     ssc-miR-429     chr6  58.044184  58.044205      - MI0017991
## 6 ssc-miR-6782-3p     chr6   0.956806   0.956827      + MI0031620
## 7     ssc-miR-874     chr2 145.381108 145.381130      - MI0022157
## 8 ssc-miR-9843-3p     chr8 122.371914 122.371935      - MI0031612
##           SNP chr.snp   pos.snp snp.effect       qvalue
## 1 MARC0093624      15 135538479          - 5.937313e-03
## 2 H3GA0034702      12  54736413          + 6.197482e-03
## 3 H3GA0055369       6  15626776          + 2.867644e-04
## 4 ALGA0041952       7  55936003          + 1.543732e-07
## 5 ASGA0094554       6  41809035          - 5.046088e-06
## 6 ASGA0094215      10  29227834          + 1.828363e-04
## 7 ALGA0016550       2 145462524          + 2.897624e-09
## 8 MARC0093624      15 135538479          + 2.945101e-03
```

This function (stb) returns a summary table of the eQTL peaks per chromosome/per gene at FDR < 0.05:


```r
rsumtb5 <- stb(qval=eqtl$gwa.qval, map=map, annotation=annotation, Z=eqtl$gwa, threshold=0.05, gene.name="Precursors", pergene=T)
dim(rsumtb5)
```

```
## [1] 26 11
```

miR-eQTL peaks at FDR<0.05:


```r
rsumtb5
```

```
##               Gene chr.gene start.gene   end.gene strand
## 1    ssc-let-7d-5p     chr3  44.867277  44.867298      +
## 2       ssc-let-7g    chr13  37.599007  37.599028      +
## 3      ssc-miR-128    chr13  22.859239  22.859259      +
## 4  ssc-miR-1306-3p    chr14  55.111005  55.111025      +
## 5   ssc-miR-140-5p     <NA>         NA         NA   <NA>
## 6   ssc-miR-140-5p     <NA>         NA         NA   <NA>
## 7     ssc-miR-1468     chrX  56.757096  56.757117      -
## 8      ssc-miR-184     chr7  53.883677  53.883698      +
## 9      ssc-miR-184     chr7  53.883677  53.883698      +
## 10     ssc-miR-184     chr7  53.883677  53.883698      +
## 11    ssc-miR-190b     chr4 104.457820 104.457841      +
## 12    ssc-miR-190b     chr4 104.457820 104.457841      +
## 13   ssc-miR-22-3p     <NA>         NA         NA   <NA>
## 14  ssc-miR-345-3p     chr7 128.658348 128.658368      +
## 15     ssc-miR-429     chr6  58.044184  58.044205      -
## 16     ssc-miR-429     chr6  58.044184  58.044205      -
## 17 ssc-miR-6782-3p     chr6   0.956806   0.956827      +
## 18 ssc-miR-7135-3p     chr3  29.052303  29.052324      +
## 19     ssc-miR-874     chr2 145.381108 145.381130      -
## 20     ssc-miR-874     chr2 145.381108 145.381130      -
## 21      ssc-miR-95     chr8   4.277286   4.277307      +
## 22 ssc-miR-9785-5p     chr3  22.233456  22.233479      -
## 23 ssc-miR-9785-5p     chr3  22.233456  22.233479      -
## 24 ssc-miR-9785-5p     chr3  22.233456  22.233479      -
## 25 ssc-miR-9810-3p     chr4  90.746075  90.746095      +
## 26 ssc-miR-9843-3p     chr8 122.371914 122.371935      -
##              gene.name         SNP chr.snp   pos.snp snp.effect
## 1            MI0022120 MARC0093624      15 135538479          -
## 2            MI0013087 MARC0093624      15 135538479          -
## 3  MI0013094,MI0002451 ALGA0023517       4  15959123          +
## 4            MI0013148 H3GA0034702      12  54736413          +
## 5                 <NA> ASGA0017748       4   7002305          -
## 6                 <NA> H3GA0055369       6  15626776          +
## 7            MI0022160 MARC0093624      15 135538479          -
## 8            MI0002421 DBWU0000430       3   9121844          -
## 9            MI0002421 M1GA0026172       6 157003423          +
## 10           MI0002421 ALGA0041952       7  55936003          +
## 11           MI0017988 ALGA0102568       3  76844137          +
## 12           MI0017988 ALGA0026452       4  94724175          -
## 13                <NA> MARC0009333      15 134397712          -
## 14           MI0013117 MARC0027291      15 135171935          +
## 15           MI0017991 ASGA0094554       6  41809035          -
## 16           MI0017991 ALGA0046283       8   6551222          +
## 17           MI0031620 ASGA0094215      10  29227834          +
## 18           MI0023568 MARC0056802       3  28401592          -
## 19           MI0022157 ALGA0016550       2 145462524          +
## 20           MI0022157 ALGA0122273       3  64022095          -
## 21           MI0002436 MARC0027291      15 135171935          -
## 22           MI0031545 ASGA0013843       3  21495074          +
## 23           MI0031545 MARC0081878      10    223442          -
## 24           MI0031545 ALGA0121561      17   3681417          +
## 25           MI0031577 MARC0021620       5  16689918          +
## 26           MI0031612 MARC0093624      15 135538479          +
##          qvalue
## 1  1.150336e-02
## 2  5.937313e-03
## 3  4.474533e-02
## 4  6.197482e-03
## 5  1.089555e-02
## 6  2.867644e-04
## 7  2.129363e-02
## 8  3.774307e-02
## 9  1.355902e-02
## 10 1.543732e-07
## 11 4.890582e-02
## 12 1.999064e-02
## 13 4.976423e-02
## 14 2.213197e-02
## 15 5.046088e-06
## 16 2.998539e-02
## 17 1.828363e-04
## 18 1.155505e-02
## 19 2.897624e-09
## 20 2.052799e-02
## 21 3.714473e-02
## 22 3.702342e-02
## 23 3.702342e-02
## 24 3.702342e-02
## 25 3.498560e-02
## 26 2.945101e-03
```

#### Summary of GWAS results at FDR < 0.01

Number of eQTL peaks per chromosome:


```r
table(rsumtb1$chr.snp)
```

```
## 
##  2  6  7 10 12 15 
##  1  2  1  1  1  2
```

Names of associated miRNAs:


```r
table(rsumtb1$Gene)
```

```
## 
##      ssc-let-7g ssc-miR-1306-3p  ssc-miR-140-5p     ssc-miR-184 
##               1               1               1               1 
##     ssc-miR-429 ssc-miR-6782-3p     ssc-miR-874 ssc-miR-9843-3p 
##               1               1               1               1
```

Names of associated markers:


```r
table(as.character(rsumtb1$SNP))
```

```
## 
## ALGA0016550 ALGA0041952 ASGA0094215 ASGA0094554 H3GA0034702 H3GA0055369 
##           1           1           1           1           1           1 
## MARC0093624 
##           2
```

#### Summary of GWAS results at FDR < 0.05

Number of eQTL peaks per chromosome:


```r
table(rsumtb5$chr.snp)
```

```
## 
##  2  3  4  5  6  7  8 10 12 15 17 
##  1  5  3  1  3  1  1  2  1  7  1
```

Names of associated miRNAs:


```r
table(rsumtb5$Gene)
```

```
## 
##   ssc-let-7d-5p      ssc-let-7g     ssc-miR-128 ssc-miR-1306-3p 
##               1               1               1               1 
##  ssc-miR-140-5p    ssc-miR-1468     ssc-miR-184    ssc-miR-190b 
##               2               1               3               2 
##   ssc-miR-22-3p  ssc-miR-345-3p     ssc-miR-429 ssc-miR-6782-3p 
##               1               1               2               1 
## ssc-miR-7135-3p     ssc-miR-874      ssc-miR-95 ssc-miR-9785-5p 
##               1               2               1               3 
## ssc-miR-9810-3p ssc-miR-9843-3p 
##               1               1
```

Names of associated markers:


```r
table(as.character(rsumtb5$SNP))
```

```
## 
## ALGA0016550 ALGA0023517 ALGA0026452 ALGA0041952 ALGA0046283 ALGA0102568 
##           1           1           1           1           1           1 
## ALGA0121561 ALGA0122273 ASGA0013843 ASGA0017748 ASGA0094215 ASGA0094554 
##           1           1           1           1           1           1 
## DBWU0000430 H3GA0034702 H3GA0055369 M1GA0026172 MARC0009333 MARC0021620 
##           1           1           1           1           1           1 
## MARC0027291 MARC0056802 MARC0081878 MARC0093624 
##           2           1           1           4
```

## Save data


```r
save(rsumtb1, file = "../2_eqtl_summary_table_fdr1.Rdata")
save(rsumtb5, file = "../3_eqtl_summary_table_fdr5.Rdata")
```

