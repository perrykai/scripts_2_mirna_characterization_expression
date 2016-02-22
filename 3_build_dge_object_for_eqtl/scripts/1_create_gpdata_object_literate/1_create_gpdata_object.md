**Script:** `1_create_gpdata_object.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`

**Date:**  2/19/16

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
2. `/mnt/research/pigsnp/MSUPRP/carcass_quality/data/`
3. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/`

**Input File(s):** 

1. `1_mean_mature_mirna_expression_unfiltered.Rdata`
2. `MSUPRP_meat.RData`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`

**Output File(s):** 

1. `3_msuprp_mirna_gpdata.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives
The objective of this script is to create the gpdata object needed for the miRNA eQTL analysis.
This gpdata object will be filtered from the MSUPRP_meat gpdata object in the following ways:

1. The covariate data will be reduced to include only animals in this analysis (n=174), and will have the additional column of growth_group as a factor.  

2. The geno data will be reduced to include only animals in this analysis (n=174), and filtered for removal of fixed SNPs, SNPs with maf < 0.10, and SNPs on sex chromosomes. 

## This analysis conducted using R/3.2.0, not R/3.1.0

## Install libraries


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts/")

library(synbreed)
library(regress)
library (limma)
```

```
## Loading required package: methods
```

```r
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)
```

```
## Warning: replacing previous import by 'grid::arrow' when loading 'qvalue'
```

```
## Warning: replacing previous import by 'grid::unit' when loading 'qvalue'
```

## Load data


```r
rm(list=ls())
```

Load the miRNA expression data:


```r
load("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")
```

Load the MSUPRP_meat gpdata object:


```r
load("/mnt/research/pigsnp/MSUPRP/carcass_quality/data/MSUPRP_meat.RData")
```

## Analysis
### 1. A matrix of gene expression, with dimensions genes x samples (335 x 174)


```r
dim(no.zero.dfmeanrcround)
```

```
## [1] 335 174
```

```r
colnames(no.zero.dfmeanrcround)
```

```
##   [1] "1034" "1036" "1041" "1049" "1058" "1060" "1080" "1082" "1085" "1091"
##  [11] "1096" "1100" "1107" "1110" "1111" "1113" "1116" "1123" "1134" "1136"
##  [21] "1145" "1147" "1152" "1154" "1158" "1170" "1177" "1179" "1192" "1194"
##  [31] "1197" "1199" "1205" "1207" "1237" "1239" "1240" "1242" "1265" "1267"
##  [41] "1278" "1282" "1291" "1295" "1300" "1304" "1321" "1323" "1423" "1424"
##  [51] "1425" "1426" "1431" "1434" "1435" "1444" "1445" "1449" "1456" "1458"
##  [61] "1482" "1484" "1491" "1493" "1502" "1504" "1510" "1512" "1517" "1523"
##  [71] "1529" "1532" "1533" "1534" "1537" "1543" "1578" "1580" "1589" "1592"
##  [81] "1593" "1594" "1625" "1627" "1638" "1640" "1644" "1646" "1652" "1662"
##  [91] "1669" "1677" "1685" "1687" "1695" "1697" "1746" "1758" "1760" "1776"
## [101] "1778" "1782" "1784" "1785" "1789" "1793" "1798" "1800" "1818" "1819"
## [111] "1820" "1833" "1836" "1839" "1843" "1844" "1879" "1881" "1884" "1886"
## [121] "1889" "1891" "1903" "1904" "1907" "1910" "1914" "1916" "1928" "1930"
## [131] "1965" "1971" "1976" "1980" "1989" "1991" "1999" "2003" "2018" "2020"
## [141] "2022" "2024" "2026" "2027" "2029" "2030" "2064" "2071" "2073" "2076"
## [151] "2094" "2100" "2118" "2119" "2120" "2123" "2131" "2135" "2141" "2143"
## [161] "2152" "2154" "2164" "2168" "2195" "2197" "2229" "2231" "2261" "2263"
## [171] "2297" "2303" "2311" "2317"
```

```r
head(rownames(no.zero.dfmeanrcround))
```

```
## [1] "ssc-let-7a"    "ssc-let-7c"    "ssc-let-7d-3p" "ssc-let-7d-5p"
## [5] "ssc-let-7e"    "ssc-let-7f"
```

Create vector of animal IDs from column names of expression matrix:


```r
pigid <- colnames(no.zero.dfmeanrcround)
length(pigid)
```

```
## [1] 174
```

```r
pigid
```

```
##   [1] "1034" "1036" "1041" "1049" "1058" "1060" "1080" "1082" "1085" "1091"
##  [11] "1096" "1100" "1107" "1110" "1111" "1113" "1116" "1123" "1134" "1136"
##  [21] "1145" "1147" "1152" "1154" "1158" "1170" "1177" "1179" "1192" "1194"
##  [31] "1197" "1199" "1205" "1207" "1237" "1239" "1240" "1242" "1265" "1267"
##  [41] "1278" "1282" "1291" "1295" "1300" "1304" "1321" "1323" "1423" "1424"
##  [51] "1425" "1426" "1431" "1434" "1435" "1444" "1445" "1449" "1456" "1458"
##  [61] "1482" "1484" "1491" "1493" "1502" "1504" "1510" "1512" "1517" "1523"
##  [71] "1529" "1532" "1533" "1534" "1537" "1543" "1578" "1580" "1589" "1592"
##  [81] "1593" "1594" "1625" "1627" "1638" "1640" "1644" "1646" "1652" "1662"
##  [91] "1669" "1677" "1685" "1687" "1695" "1697" "1746" "1758" "1760" "1776"
## [101] "1778" "1782" "1784" "1785" "1789" "1793" "1798" "1800" "1818" "1819"
## [111] "1820" "1833" "1836" "1839" "1843" "1844" "1879" "1881" "1884" "1886"
## [121] "1889" "1891" "1903" "1904" "1907" "1910" "1914" "1916" "1928" "1930"
## [131] "1965" "1971" "1976" "1980" "1989" "1991" "1999" "2003" "2018" "2020"
## [141] "2022" "2024" "2026" "2027" "2029" "2030" "2064" "2071" "2073" "2076"
## [151] "2094" "2100" "2118" "2119" "2120" "2123" "2131" "2135" "2141" "2143"
## [161] "2152" "2154" "2164" "2168" "2195" "2197" "2229" "2231" "2261" "2263"
## [171] "2297" "2303" "2311" "2317"
```

### 2. A data frame of metadata (trait data) with dimensions samples x categories (174 x categories), including growth_group as a factor:


Discard the animals from the gpdata geno object not present in the miRNA expression matrix:


```r
todisc <- rownames(MSUPRP_meat$geno)[!rownames(MSUPRP_meat$geno) %in% pigid]

redMSU <- discard.individuals(MSUPRP_meat, todisc)
```

Remove covariate data from animals not in this analysis:


```r
todisc <- redMSU$covar$id[!redMSU$covar$id %in% rownames(redMSU$geno)]

redMSU <- discard.individuals(redMSU, todisc)
```

Create the growth_group column as a factor by combining selcrit and Status


```r
redMSU$covar <- data.frame(redMSU$covar,
        growth_group=paste(redMSU$covar$selcrit, redMSU$covar$Status, sep="-"))
```

Change the levels of the factor growth_group to have lma-L the first level


```r
redMSU$covar$growth_group<-relevel(redMSU$covar$growth_group, ref = "lma-L")

is.factor(redMSU$covar$growth_group)
```

```
## [1] TRUE
```

```r
head(redMSU$covar)
```

```
##      id phenotyped genotyped family sex litter wt_10wk bf10_22wk lma_22wk
## 26 1034       TRUE      TRUE     NA   F      5   30.84     16.51    47.55
## 28 1036       TRUE      TRUE     NA   F      5   29.48     18.80    31.16
## 33 1041       TRUE      TRUE     NA   M      5   32.66     21.08    41.94
## 41 1049       TRUE      TRUE     NA   M      5   29.03     21.59    31.94
## 50 1058       TRUE      TRUE     NA   F     10   24.95     10.16    39.10
## 52 1060       TRUE      TRUE     NA   F     10   20.87     17.02    31.03
##    slgdt_cd age_slg car_wt Status selcrit microarray_Dye microarray_file
## 26        2     160  85.49      H     lma            cy3     slide33.gpr
## 28        2     160  71.66      L     lma            cy5     slide33.gpr
## 33        2     160  80.27      H     lma            cy5     slide34.gpr
## 41        2     160  81.86      L     lma            cy3     slide34.gpr
## 50        2     159  72.11      H     lma            cy5     slide52.gpr
## 52        2     159  70.75      L     lma            cy3     slide52.gpr
##    Color Hairden earset Spots Underbelly face_color line_origin perc_duroc
## 26  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.5051327
## 28  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.5471133
## 33  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.4721328
## 41  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.4565887
## 50  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.4336773
## 52  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.5354381
##    growth_group
## 26        lma-H
## 28        lma-L
## 33        lma-H
## 41        lma-L
## 50        lma-H
## 52        lma-L
```

```r
if (sum(redMSU$covar$id!=pigid)!=0) stop ("pigids of redMSU$covar not correct")
```

The dataframe of covariates is complete.

### 3. A matrix of genotype data, filtered for removing fixed SNPs, SNPs with low maf, and SNPs on sex chromosomes. 
Extract the genotypes from the GPData object:


```r
genomat <- redMSU$geno

head(genomat[,1:10])
```

```
##      MARC0044150 ASGA0000014 H3GA0000026 ASGA0000021 ALGA0000009
## 1034           0           2           2           2           2
## 1036           0           2           2           2           2
## 1041           1           2           1           1           1
## 1049           1           2           1           1           1
## 1058           1           2           1           1           1
## 1060           2           2           0           0           0
##      ALGA0000014 H3GA0000032 ASGA0000005 M1GA0000060 ASGA0000047
## 1034           2           0           0           2           2
## 1036           2           0           0           2           2
## 1041           1           1           1           1           2
## 1049           1           1           1           1           2
## 1058           1           1           1           1           1
## 1060           0           2           2           0           0
```

```r
dim(genomat)
```

```
## [1]   174 45329
```

```r
if (sum(rownames(genomat)!=redMSU$covar$id) != 0) stop ("rows of geno not the same as ids of covar")
if (sum(colnames(genomat)!=colnames(MSUPRP_meat$geno))!=0) stop ("columns of genomat not the same as MSUPRP_meat$geno")
```

Filter the matrix of genotypes: 

Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):


```r
sdv <- apply(genomat, 2, sd)
length(sdv)
```

```
## [1] 45329
```

```r
sum(sdv == 0)
```

```
## [1] 34
```

Remove fixed SNPs from genotype matrix:


```r
genomat <- genomat[,sdv>0]

dim(genomat)
```

```
## [1]   174 45295
```

Filter SNPs for minor allele frequency (maf > 0.10 retained):

First, get the allelic frequency:


```r
af <- colMeans(genomat)/2
summary(af)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.002874 0.287400 0.511500 0.508000 0.729900 0.997100
```

How many SNPs have a af < 0.10 and >0.90?

We filter this way since we don't 'know' the minor allele in each SNP based on allele frequency:


```r
sum(af<0.10)
```

```
## [1] 3184
```

```r
sum(af>0.90)
```

```
## [1] 3486
```

```r
sum(c(af<0.10,af>0.90))
```

```
## [1] 6670
```

Retain all SNPs with af > 0.10 and <0.90 (maf < 0.10 discarded):


```r
genomat <- genomat[ ,af>0.10 & af<0.90]
```

Dimensions of remaining SNP marker matrix:


```r
dim(genomat)
```

```
## [1]   174 38625
```

```r
if (sum((colMeans(genomat)/2)<=0.10) !=0) stop ("maf filtering did not work correctly")
if (sum((colMeans(genomat)/2)>=0.90) !=0) stop ("maf filtering did not work correctly")
if (sum(rownames(genomat)!=redMSU$covar$id)!=0) stop ("rownames of marker matrix not equal to trait data")
```

Eliminate markers on sex chromosomes:


```r
table(redMSU$map$chr)
```

```
## 
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
## 5383 2967 2418 3109 2034 2843 2827 2370 2814 1588 1641 1436 3470 3388 2434 
##   16   17   18   19 
## 1669 1470 1050  418
```

```r
sexchr <- rownames(redMSU$map)[redMSU$map$chr == 19]
length(sexchr)
```

```
## [1] 418
```

```r
sum((colnames(genomat) %in% sexchr))
```

```
## [1] 306
```

```r
genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]
```

This object contains the markers that are not fixed, have a "maf" > 0.10, and are not mapped to the sex chromosomes:


```r
dim(genomatfil)
```

```
## [1]   174 38319
```

```r
if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(redMSU$covar$id != rownames(genomatfil)) != 0) stop ("rownames of trait data and genotype matrix not the same")
if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")
```

Create the todel vector, containing the names of the SNPs in the gpdata geno object NOT found in the filtered genotype matrix:


```r
todel <- colnames(redMSU$geno)[!colnames(redMSU$geno) %in% colnames(genomatfil)]
length(todel)
```

```
## [1] 7010
```

That means that, in total, we are removing 7010 SNPs from the gpdata geno object. 

Using discard.markers allows us to remove both the markers we don't want, and the map information we don't want, all in one step.


```r
final_MSUPRP_miRNA<-discard.markers(redMSU, todel)
```

The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 


```r
summary_final_MSUPRP_miRNA<-summary(final_MSUPRP_miRNA)
summary_final_MSUPRP_miRNA$geno
```

```
## $nMarkers
## [1] 38319
## 
## $genotypes
##     0     1     2 
## 0.281 0.427 0.292 
## 
## $nNA
## [1] 0
## 
## $markerChr
## 
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
## 4662 2530 1933 2703 1710 2403 2367 2086 2307 1366 1377 1278 3049 2825 2134 
##   16   17   18 
## 1496 1265  828 
## 
## $mappedMarkers
## [1] 38319
```

```r
summary_final_MSUPRP_miRNA$pedigree
```

```
## Number of 
## 	 individuals  174 
## 	 males :  86 , females :  88 
## 	 Par 1 (sire)  0 
## 	 Par 2 (dam)   0 
## 	 generations  1
```


## Save data

Save the gpdata object with filtered genotypes, map information, and covariate information:


```r
save(final_MSUPRP_miRNA, summary_final_MSUPRP_miRNA, file="../3_msuprp_mirna_gpdata.Rdata")
```

