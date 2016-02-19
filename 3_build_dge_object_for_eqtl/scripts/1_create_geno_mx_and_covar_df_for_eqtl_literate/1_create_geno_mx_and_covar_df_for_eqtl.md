**Script:** `1_create_metadata_object.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/scripts`

**Date:**  2/18/16

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
2. `/mnt/research/pigsnp/MSUPRP/carcass_quality/data/`
3. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/`

**Input File(s):** 

1. `1_mean_mature_mirna_expression_unfiltered.Rdata`
2. `MSUPRP_meat.RData`
3. `ssc.gff3`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`

**Output File(s):** 

1. `3_covardata_for_eQTL_analysis.Rdata`
2. `4_filtered_genotypes_for_G_matrix.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives
The objective of this script is to create the other components of data needed for the miRNA eQTL analysis.
These data objects are needed for this analysis:

1. A matrix of gene expression, with dimensions genes x samples (335 x 174)

2. A data frame of metadata (trait data) with dimensions samples x categories (174 x categories)

3. A matrix of genotype data, filtered for removing fixed SNPs, those SNPs with low maf, and those SNPs on sex chromosomes. 

## This analysis conducted using R/3.2.0, not R/3.1.0

## Install libraries


```r
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

Load the gp data object:


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

### 2. A data frame of metadata (trait data) with dimensions samples x categories (174 x categories)

Extract data frame of covariates (trait data) from GPData object:


```r
covaMSU<-MSUPRP_meat$covar
rownames(covaMSU)<-covaMSU$id
colnames(covaMSU)
```

```
##  [1] "id"              "phenotyped"      "genotyped"      
##  [4] "family"          "sex"             "litter"         
##  [7] "wt_10wk"         "bf10_22wk"       "lma_22wk"       
## [10] "slgdt_cd"        "age_slg"         "car_wt"         
## [13] "Status"          "selcrit"         "microarray_Dye" 
## [16] "microarray_file" "Color"           "Hairden"        
## [19] "earset"          "Spots"           "Underbelly"     
## [22] "face_color"      "line_origin"     "perc_duroc"
```

```r
head(covaMSU)
```

```
##        id phenotyped genotyped family  sex litter wt_10wk bf10_22wk
## 6070 6070      FALSE      TRUE     NA <NA>   <NA>      NA        NA
## 6071 6071      FALSE      TRUE     NA <NA>   <NA>      NA        NA
## 6086 6086      FALSE      TRUE     NA <NA>   <NA>      NA        NA
## 6088 6088      FALSE      TRUE     NA <NA>   <NA>      NA        NA
## 6092 6092      FALSE      TRUE     NA <NA>   <NA>      NA        NA
## 6323 6323      FALSE      TRUE     NA <NA>   <NA>      NA        NA
##      lma_22wk slgdt_cd age_slg car_wt Status selcrit microarray_Dye
## 6070       NA     <NA>      NA     NA   <NA>    <NA>           <NA>
## 6071       NA     <NA>      NA     NA   <NA>    <NA>           <NA>
## 6086       NA     <NA>      NA     NA   <NA>    <NA>           <NA>
## 6088       NA     <NA>      NA     NA   <NA>    <NA>           <NA>
## 6092       NA     <NA>      NA     NA   <NA>    <NA>           <NA>
## 6323       NA     <NA>      NA     NA   <NA>    <NA>           <NA>
##      microarray_file Color Hairden earset Spots Underbelly face_color
## 6070            <NA>  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 6071            <NA>  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 6086            <NA>  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 6088            <NA>  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 6092            <NA>  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 6323            <NA>  <NA>    <NA>   <NA>  <NA>       <NA>         NA
##      line_origin perc_duroc
## 6070           2          0
## 6071           2          0
## 6086           2          0
## 6088           2          0
## 6092           2          0
## 6323           1          1
```

```r
dim(covaMSU)
```

```
## [1] 1037   24
```

```r
if (sum(rownames(covaMSU) != MSUPRP_meat$covar$id) != 0) stop ("rownames not the same as MSUPRP_meat id column")
```

Subset the covaMSU dataframe, retaining only those animals in pigid:


```r
redcovaMSU<-covaMSU[(rownames(covaMSU) %in% pigid),]
dim(redcovaMSU)
```

```
## [1] 174  24
```

```r
head(redcovaMSU)
```

```
##        id phenotyped genotyped family sex litter wt_10wk bf10_22wk
## 1034 1034       TRUE      TRUE     NA   F      5   30.84     16.51
## 1036 1036       TRUE      TRUE     NA   F      5   29.48     18.80
## 1041 1041       TRUE      TRUE     NA   M      5   32.66     21.08
## 1049 1049       TRUE      TRUE     NA   M      5   29.03     21.59
## 1058 1058       TRUE      TRUE     NA   F     10   24.95     10.16
## 1060 1060       TRUE      TRUE     NA   F     10   20.87     17.02
##      lma_22wk slgdt_cd age_slg car_wt Status selcrit microarray_Dye
## 1034    47.55        2     160  85.49      H     lma            cy3
## 1036    31.16        2     160  71.66      L     lma            cy5
## 1041    41.94        2     160  80.27      H     lma            cy5
## 1049    31.94        2     160  81.86      L     lma            cy3
## 1058    39.10        2     159  72.11      H     lma            cy5
## 1060    31.03        2     159  70.75      L     lma            cy3
##      microarray_file Color Hairden earset Spots Underbelly face_color
## 1034     slide33.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1036     slide33.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1041     slide34.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1049     slide34.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1058     slide52.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1060     slide52.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
##      line_origin perc_duroc
## 1034          NA  0.5051327
## 1036          NA  0.5471133
## 1041          NA  0.4721328
## 1049          NA  0.4565887
## 1058          NA  0.4336773
## 1060          NA  0.5354381
```

```r
if (sum(rownames(redcovaMSU) != pigid) != 0) stop ("rows are not correct")
if (sum(rownames(redcovaMSU) != redcovaMSU$id) != 0) stop ("rownames do not match id column")
```

Create the growth_group column as a factor by combining Status and selcrit


```r
growth_group<-paste(redcovaMSU$selcrit, redcovaMSU$Status, sep = "-")
length(growth_group)
```

```
## [1] 174
```

```r
head(growth_group)
```

```
## [1] "lma-H" "lma-L" "lma-H" "lma-L" "lma-H" "lma-L"
```

```r
redcovaMSU$growth_group<-as.factor(growth_group)
```

Change the levels of the factor growth_group to have lma-L the first level


```r
redcovaMSU$growth_group<-relevel(redcovaMSU$growth_group, ref = "lma-L")

is.factor(redcovaMSU$growth_group)
```

```
## [1] TRUE
```

```r
head(redcovaMSU)
```

```
##        id phenotyped genotyped family sex litter wt_10wk bf10_22wk
## 1034 1034       TRUE      TRUE     NA   F      5   30.84     16.51
## 1036 1036       TRUE      TRUE     NA   F      5   29.48     18.80
## 1041 1041       TRUE      TRUE     NA   M      5   32.66     21.08
## 1049 1049       TRUE      TRUE     NA   M      5   29.03     21.59
## 1058 1058       TRUE      TRUE     NA   F     10   24.95     10.16
## 1060 1060       TRUE      TRUE     NA   F     10   20.87     17.02
##      lma_22wk slgdt_cd age_slg car_wt Status selcrit microarray_Dye
## 1034    47.55        2     160  85.49      H     lma            cy3
## 1036    31.16        2     160  71.66      L     lma            cy5
## 1041    41.94        2     160  80.27      H     lma            cy5
## 1049    31.94        2     160  81.86      L     lma            cy3
## 1058    39.10        2     159  72.11      H     lma            cy5
## 1060    31.03        2     159  70.75      L     lma            cy3
##      microarray_file Color Hairden earset Spots Underbelly face_color
## 1034     slide33.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1036     slide33.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1041     slide34.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1049     slide34.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1058     slide52.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
## 1060     slide52.gpr  <NA>    <NA>   <NA>  <NA>       <NA>         NA
##      line_origin perc_duroc growth_group
## 1034          NA  0.5051327        lma-H
## 1036          NA  0.5471133        lma-L
## 1041          NA  0.4721328        lma-H
## 1049          NA  0.4565887        lma-L
## 1058          NA  0.4336773        lma-H
## 1060          NA  0.5354381        lma-L
```

```r
tail(redcovaMSU)
```

```
##        id phenotyped genotyped family sex litter wt_10wk bf10_22wk
## 2261 2261       TRUE      TRUE     NA   M    134   24.95     21.08
## 2263 2263       TRUE      TRUE     NA   M    134   26.76     36.32
## 2297 2297       TRUE      TRUE     NA   M    139   30.84     39.62
## 2303 2303       TRUE      TRUE     NA   M    139   33.57     25.40
## 2311 2311       TRUE      TRUE     NA   M    140   24.49     31.24
## 2317 2317       TRUE      TRUE     NA   M    140   24.95     19.81
##      lma_22wk slgdt_cd age_slg car_wt Status selcrit microarray_Dye
## 2261    31.81       32     167  69.84      L      bf            cy5
## 2263    41.42       32     167  79.59      H      bf            cy3
## 2297    39.35       32     166  90.02      H      bf            cy5
## 2303    41.55       32     166  96.60      L      bf            cy3
## 2311    43.81       32     165  84.13      H      bf            cy3
## 2317    41.61       32     165  78.68      L      bf            cy5
##      microarray_file     Color      Hairden       earset  Spots
## 2261    slide173.gpr       red intermediate intermediate small 
## 2263    slide173.gpr      fawn intermediate         down small 
## 2297     slide56.gpr       red intermediate intermediate   none
## 2303     slide56.gpr light red intermediate intermediate small 
## 2311     slide53.gpr       red intermediate intermediate   none
## 2317     slide53.gpr       red intermediate        erect small 
##        Underbelly face_color line_origin perc_duroc growth_group
## 2261         none          1          NA  0.4974465         bf-L
## 2263         none          0          NA  0.3115709         bf-H
## 2297 intermediate          0          NA  0.5569193         bf-H
## 2303 intermediate          0          NA  0.4879955         bf-L
## 2311 intermediate          0          NA  0.5136357         bf-H
## 2317         none          1          NA  0.5229859         bf-L
```

```r
dim(redcovaMSU)
```

```
## [1] 174  25
```

The matrix of covariates is complete. This will be utilized in the GBLUP when combined with the matrix of normalized read counts.

### 3. A matrix of genotype data, filtered for removing fixed SNPs, SNPs with low maf, and SNPs on sex chromosomes. 

Discard the animals from the GPData Geno object not present in the miRNA expression matrix:


```r
todisc <- rownames(MSUPRP_meat$geno)[!rownames(MSUPRP_meat$geno) %in% pigid]

redgenoMSU <- discard.individuals(MSUPRP_meat, todisc)
```

Extract the genotypes from the GPData object:


```r
genomat <- redgenoMSU$geno

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
if (sum(rownames(genomat)!=rownames(redcovaMSU)) != 0) stop ("rows of geno not the same as rows of pheno")
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

Filter SNPs for minor allele frequency (maf > 0.1 retained):


```r
maf <- colMeans(genomat)/2

head(maf)
```

```
## MARC0044150 ASGA0000014 H3GA0000026 ASGA0000021 ALGA0000009 ALGA0000014 
##   0.3793103   0.8017241   0.6293103   0.6206897   0.6293103   0.6293103
```

How many SNPs have a maf <= 0.1?


```r
sum(maf<=0.1)
```

```
## [1] 3184
```

Retain all SNPs with maf > 0.1:


```r
genomat <- genomat[,maf>0.1]
```

Dimensions of remaining SNP marker matrix:


```r
dim(genomat)
```

```
## [1]   174 42111
```

```r
if (sum((colMeans(genomat)/2)<=0.1) !=0) stop ("maf filtering did not work correctly")

if (sum(rownames(genomat)!=rownames(redcovaMSU))!=0) stop ("rownames of marker matrix not equal to trait data")
```

Eliminate markers on sex chromosomes:


```r
sexchr <- rownames(redgenoMSU$map)[redgenoMSU$map$chr == 19]
length(sexchr)
```

```
## [1] 418
```

```r
sum((colnames(genomat) %in% sexchr))
```

```
## [1] 367
```

```r
genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]

dim(genomatfil)
```

```
## [1]   174 41744
```

```r
if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(rownames(redcovaMSU) != rownames(genomatfil)) != 0) stop ("rownames of trait data and genotype matrix not the same")
if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")
```

The filtered matrix of SNPs is complete. This will be used in the GBLUP and GWAS as the G matrix (genomic relationship matrix) after standardizing the SNPs. 

## Save data

Save the trait data (covariates):


```r
save(redcovaMSU,file = "../3_covardata_for_eQTL_analysis.Rdata")
```

Save the filtered genotype matrix:


```r
save(genomatfil,file = "../4_filtered_genotypes_for_G_matrix.Rdata")
```

