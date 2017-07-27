**Script:** `5_target_mRNA-pQTL-colocalization.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/scripts`

**Date:**  6/29/17

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/`

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/`

1. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/`

**Input File(s):** 

1. `data.Rdata`, `pQTL_ALL.Rdata`, `inrange_function.Rdata`

1. `10_mRNA_miRNA_correlation_output.Rdata`

1. `gpData_PRKAG3.Rdata`

**Output File Directory:** 

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/`

**Output File(s):** ``

1. `12_target_mrna_colocalized_pqtl.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to co-localize the significant target mRNAs with pQTL identified in the MSUPRP dataset

## Install libraries



```r
rm(list=ls())
library(methods)
library(limma)
library(edgeR)
```

Session Information


```r
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: CentOS release 6.8 (Final)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] edgeR_3.12.1 limma_3.26.9 knitr_1.14  
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5  tools_3.2.0   stringi_1.1.1 stringr_1.1.0 evaluate_0.9
```

## Load data

Load required R objects 


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/data.Rdata")
load("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/10_mRNA_miRNA_correlation_output.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/pQTL_ALL.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/gpData_PRKAG3.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/inrange_function.Rdata")
```

## Analysis

### Annotation of corrleated target genes


```r
genes <- dge$genes
colnames(genes)[1] <- "chr"
annot <- do.call(rbind, lapply(names(sig.mrnaR), function(x) 
	data.frame(genes[sig.mrnaR[[x]],], miRNA=rep(x, length(sig.mrnaR[[x]])), 
		rst.corR[[x]][sig.mrnaR[[x]],c("cor", "pvalue", "qvalue")])))
```

### Phenotypic QTL MSUPRP
q-values pQTL GWAS


```r
qval <- pqtl$gwa.qval
dim(qval)
```

```
## [1] 44911    67
```

p-values pQTL GWAS


```r
pval <- pqtl$gwa.pval
dim(pval)
```

```
## [1] 44911    67
```

Standardized SNP effects pQTL GWA


```r
sdEff <- pqtl$gwa
dim(sdEff)
```

```
## [1] 44911    67
```

Marker Map


```r
map <- data.frame(chr=MSUPRP_meat$map[[1]], pos=MSUPRP_meat$map[[2]])
rownames(map) <- colnames(MSUPRP_meat$geno)
dim(map)
```

```
## [1] 45331     2
```

Retain marker information for all pQTL


```r
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(map[names(sig[[x]]),], 
	std.eff=sdEff[names(sig[[x]]),x], pvalue=pval[names(sig[[x]]),x], qvalue=sig[[x]]))
names(sig) <- colnames(qval)
sig <- sig[unlist(lapply(sig, nrow)) > 0]
length(sig)
```

```
## [1] 27
```

```r
names(sig)
```

```
##  [1] "bf10_10wk"  "lrf_10wk"   "bf10_13wk"  "lrf_13wk"   "bf10_16wk" 
##  [6] "lma_16wk"   "lrf_16wk"   "bf10_19wk"  "lrf_19wk"   "bf10_22wk" 
## [11] "lrf_22wk"   "dress_ptg"  "cook_yield" "WBS"        "juiciness" 
## [16] "tenderness" "overtend"   "driploss"   "ph_24h"     "car_length"
## [21] "num_ribs"   "last_lum"   "car_bf10"   "car_lma"    "loin"      
## [26] "belly"      "protein"
```

If a pQTL contains more than one peak split each peak


```r
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
```

pQTL genomic regions


```r
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
```

```
## [1] 58  4
```

### Co-localization mRNA with pQTL 

Check all mRNA within a eQTL


```r
win <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single=NULL, range=c(start="start", end="end")))
names(win) <- rownames(qtlP)
```

mRNA overlaping left side of pQTL


```r
left <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="start", range=NULL))
names(left) <- rownames(qtlP)
```

mRNA overlaping right side of pQTL


```r
right <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="end", range=NULL))
names(right) <- rownames(qtlP)
```

Merge all mRNA co-localizing with pQTL  


```r
# Merge within and left side
coloc <- lapply(names(win), function(x) rbind(win[[x]], 
	left[[x]][!as.character(left[[x]]$geneID) %in% as.character(win[[x]]$geneID) | 
	!as.character(left[[x]]$miRNA) == as.character(win[[x]]$miRNA),]))
names(coloc) <- names(win)
```

When we create the merged coloc & right side object, 
we receive a warning. 

Check why we're getting the warning:


```r
right[[36]][!rownames(right[[36]]) %in% rownames(coloc[[36]]) | !right[[36]]$miRNA == coloc[[36]]$miRNA,]
```

```
## Warning in is.na(e1) | is.na(e2): longer object length is not a multiple of
## shorter object length
```

```
## Warning in `==.default`(right[[36]]$miRNA, coloc[[36]]$miRNA): longer
## object length is not a multiple of shorter object length
```

```
##             chr     start       end strand      geneID genes       miRNA
## XLOC_022328   7  64023230  64071596      - XLOC_022328   PML ssc-miR-874
## XLOC_022593   7 105973962 105978050      - XLOC_022593 VASH1 ssc-miR-874
##                   cor      pvalue      qvalue
## XLOC_022328  0.143337 0.006120952 0.011768918
## XLOC_022593 -0.165973 0.001502863 0.006561806
```

For the num_ribs phenotype (36), one gene overlaps the pQTL peak to the right that is not within the coloc object.
So, when building the merged object, the two objects will have different lengths, triggering the warning 
and adding the last line of the longer list to the merged object (hence the duplicate, below).


```r
# Merge coloc and right side
coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]], 
	right[[x]][!as.character(right[[x]]$geneID) %in% as.character(coloc[[x]]$geneID) |
	!as.character(right[[x]]$miRNA) == as.character(coloc[[x]]$miRNA),]))
```

```
## Warning in as.character(right[[x]]$miRNA) == as.character(coloc[[x]]
## $miRNA): longer object length is not a multiple of shorter object length
```

```r
names(coloc) <- names(win)
```

Final list of mRNA targets significantly correlated with miRNA and co-localizing with a pQTL


```r
coloc <- do.call(rbind, lapply(names(coloc), function(x) data.frame(coloc[[x]], 
	pheno=rep(strsplit(x, "[.]")[[1]][1], nrow(coloc[[x]])))))
rownames(coloc) <- NULL
dim(coloc)
```

```
## [1] 181  11
```

```r
head(coloc)
```

```
##   chr     start       end strand      geneID   genes           miRNA
## 1   1 305501875 305633144      - XLOC_002454 RAPGEF1 ssc-miR-6782-3p
## 2   1 305264778 305272661      + XLOC_001188  FAM78A     ssc-miR-874
## 3   1 307984863 308176494      + XLOC_001212  COL5A1     ssc-miR-874
## 4   1 305269495 305301268      - XLOC_002450  FAM78A     ssc-miR-874
## 5   1 307276016 307449234      - XLOC_002471    VAV2     ssc-miR-874
## 6   2   5146894   5149713      - XLOC_012991   CD248  ssc-miR-140-5p
##          cor      pvalue      qvalue    pheno
## 1  0.1107704 0.034140009 0.048159832 car_bf10
## 2 -0.1323841 0.011349059 0.016629714 car_bf10
## 3  0.1534137 0.003346881 0.009463316 car_bf10
## 4 -0.1599854 0.002216171 0.007922780 car_bf10
## 5  0.1262505 0.015758490 0.018814749 car_bf10
## 6  0.1228916 0.018762077 0.048013127 overtend
```

Investigate the occurences of multiple rows:


```r
sort(table(as.character(coloc$genes)))
```

```
## 
##   ANGEL1     BAG2     BOLL     BYSL   COL5A1     DIO2   HECTD1     HES6 
##        1        1        1        1        1        1        1        1 
##   HS6ST1   INPP5B     JPH4   METTL8   NCKAP1     NOP9    NR4A2      PML 
##        1        1        1        1        1        1        1        1 
##    POMT2    PROX2  PRPF40A    PTPN6  RAPGEF1    RAPH1    RRP36     SCLY 
##        1        1        1        1        1        1        1        1 
##    SF3B1  SLC29A1   SLC7A1  SMPDL3B    STAT1     TFEB  TINAGL1 TMEM200B 
##        1        1        1        1        1        1        1        1 
##    USP12     VAV2    VSNL1    WDR33   ZNF451   B3GAT3     CAP1    EFHD1 
##        1        1        1        1        1        2        2        2 
##     ETV7   FAM78A     GRM4     HEYL    KHNYN   KLHL30   KLHL31    LEMD2 
##        2        2        2        2        2        2        2        2 
##     MEN1    PFDN6   STOML1  TRMT112  ZSCAN20    ATG2A   CARNS1  CCDC85B 
##        2        2        2        2        2        3        3        3 
##      CCS    CD248 CDC42EP2  CDK2AP2    POLA2   PPP1CA     RELA    SYT12 
##        3        3        3        3        3        3        3        3 
##    VASH1     RGMA   RNF125    SP110 KIAA1328    OVOL1    POLD4    PTPRN 
##        3        4        4        4        5        6        6        8 
##    TPGS2      CTH   IGFBP5 
##       10       11       18
```

This is the only gene-miRNA-phenotype that is duplicated in the table; 
all other occurrences of multiple gene names differ either by phenotype, miRNA-associated, or XLOC gene ID.


```r
coloc[coloc$genes=="VASH1",]
```

```
##     chr     start       end strand      geneID genes       miRNA
## 123   7 105959003 105978240      + XLOC_021682 VASH1 ssc-miR-874
## 131   7 105973962 105978050      - XLOC_022593 VASH1 ssc-miR-874
## 133   7 105973962 105978050      - XLOC_022593 VASH1 ssc-miR-874
##            cor      pvalue      qvalue    pheno
## 123 -0.1662651 0.001474195 0.006561806 num_ribs
## 131 -0.1659730 0.001502863 0.006561806 num_ribs
## 133 -0.1659730 0.001502863 0.006561806 num_ribs
```

Extract the pertinent columns for the significant negatively-associated mRNAs overlapping pQTLs:


```r
negcoloc<-coloc[coloc$cor < 0,c("chr", "geneID", "genes", "miRNA", "cor", "pheno")]
negcoloc
```

```
##     chr      geneID   genes           miRNA         cor      pheno
## 2     1 XLOC_001188  FAM78A     ssc-miR-874 -0.13238408   car_bf10
## 4     1 XLOC_002450  FAM78A     ssc-miR-874 -0.15998540   car_bf10
## 10    2 XLOC_011808  PPP1CA     ssc-miR-874 -0.11106243   overtend
## 11    2 XLOC_011809   POLD4     ssc-miR-874 -0.10960204   overtend
## 13    2 XLOC_011840   OVOL1     ssc-miR-874 -0.16962395   overtend
## 15    2 XLOC_011874 TRMT112     ssc-miR-874 -0.17137641   overtend
## 18    2 XLOC_012974   POLD4     ssc-miR-874 -0.11953268   overtend
## 19    2 XLOC_013003 CCDC85B     ssc-miR-874 -0.09631252   overtend
## 20    2 XLOC_013007   OVOL1     ssc-miR-874 -0.10010953   overtend
## 27    2 XLOC_011808  PPP1CA     ssc-miR-874 -0.11106243        WBS
## 28    2 XLOC_011809   POLD4     ssc-miR-874 -0.10960204        WBS
## 30    2 XLOC_011840   OVOL1     ssc-miR-874 -0.16962395        WBS
## 33    2 XLOC_012974   POLD4     ssc-miR-874 -0.11953268        WBS
## 34    2 XLOC_013003 CCDC85B     ssc-miR-874 -0.09631252        WBS
## 35    2 XLOC_013007   OVOL1     ssc-miR-874 -0.10010953        WBS
## 43    2 XLOC_011808  PPP1CA     ssc-miR-874 -0.11106243 tenderness
## 44    2 XLOC_011809   POLD4     ssc-miR-874 -0.10960204 tenderness
## 46    2 XLOC_011840   OVOL1     ssc-miR-874 -0.16962395 tenderness
## 48    2 XLOC_011874 TRMT112     ssc-miR-874 -0.17137641 tenderness
## 51    2 XLOC_012974   POLD4     ssc-miR-874 -0.11953268 tenderness
## 52    2 XLOC_013003 CCDC85B     ssc-miR-874 -0.09631252 tenderness
## 53    2 XLOC_013007   OVOL1     ssc-miR-874 -0.10010953 tenderness
## 57    3 XLOC_015732   VSNL1     ssc-miR-874 -0.12960935  bf10_13wk
## 61    6 XLOC_019472 SMPDL3B     ssc-miR-874 -0.16816356   car_bf10
## 63    6 XLOC_019526 ZSCAN20     ssc-miR-874 -0.14888645   car_bf10
## 66    6 XLOC_020607 ZSCAN20     ssc-miR-874 -0.25140562   car_bf10
## 99    7 XLOC_021297   RRP36  ssc-miR-140-5p -0.15166119  dress_ptg
## 101   7 XLOC_021210  KLHL31     ssc-miR-874 -0.10536692  dress_ptg
## 102   7 XLOC_021219   PFDN6     ssc-miR-874 -0.11646586  dress_ptg
## 103   7 XLOC_021224   LEMD2     ssc-miR-874 -0.15355969  dress_ptg
## 104   7 XLOC_021280    BYSL     ssc-miR-874 -0.09251552  dress_ptg
## 106   7 XLOC_022087  KLHL31     ssc-miR-874 -0.12858708  dress_ptg
## 109   7 XLOC_022131    ETV7     ssc-miR-874 -0.13764147  dress_ptg
## 111   7 XLOC_021219   PFDN6     ssc-miR-874 -0.11646586    car_lma
## 112   7 XLOC_021224   LEMD2     ssc-miR-874 -0.15355969    car_lma
## 114   7 XLOC_022131    ETV7     ssc-miR-874 -0.13764147    car_lma
## 117   7 XLOC_022329  STOML1 ssc-miR-6782-3p -0.15253742   num_ribs
## 118   7 XLOC_022380   KHNYN ssc-miR-6782-3p -0.12420591   num_ribs
## 119   7 XLOC_022509    RGMA ssc-miR-6782-3p -0.12595838   num_ribs
## 120   7 XLOC_021474  HECTD1     ssc-miR-874 -0.09777291   num_ribs
## 123   7 XLOC_021682   VASH1     ssc-miR-874 -0.16626506   num_ribs
## 126   7 XLOC_022329  STOML1     ssc-miR-874 -0.12070099   num_ribs
## 127   7 XLOC_022380   KHNYN     ssc-miR-874 -0.17093830   num_ribs
## 128   7 XLOC_022386    NOP9     ssc-miR-874 -0.17152245   num_ribs
## 129   7 XLOC_022509    RGMA     ssc-miR-874 -0.11544359   num_ribs
## 130   7 XLOC_022579   PROX2     ssc-miR-874 -0.09908726   num_ribs
## 131   7 XLOC_022593   VASH1     ssc-miR-874 -0.16597298   num_ribs
## 133   7 XLOC_022593   VASH1     ssc-miR-874 -0.16597298   num_ribs
## 134  11 XLOC_003611   USP12     ssc-miR-874 -0.11442132  dress_ptg
## 142  15 XLOC_008772   WDR33     ssc-miR-874 -0.10230011    protein
## 143  15 XLOC_009022  IGFBP5     ssc-miR-874 -0.12508215    protein
## 144  15 XLOC_009357 PRPF40A     ssc-miR-874 -0.16334429    protein
## 147  15 XLOC_009584  IGFBP5     ssc-miR-874 -0.10098576    protein
## 148  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982    protein
## 149  15 XLOC_009656   SP110     ssc-miR-874 -0.14757211    protein
## 152  15 XLOC_009022  IGFBP5     ssc-miR-874 -0.12508215     ph_24h
## 153  15 XLOC_009584  IGFBP5     ssc-miR-874 -0.10098576     ph_24h
## 154  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982     ph_24h
## 155  15 XLOC_009656   SP110     ssc-miR-874 -0.14757211     ph_24h
## 158  15 XLOC_009022  IGFBP5     ssc-miR-874 -0.12508215 cook_yield
## 159  15 XLOC_009584  IGFBP5     ssc-miR-874 -0.10098576 cook_yield
## 160  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982 cook_yield
## 161  15 XLOC_009656   SP110     ssc-miR-874 -0.14757211 cook_yield
## 164  15 XLOC_009022  IGFBP5     ssc-miR-874 -0.12508215   driploss
## 167  15 XLOC_009584  IGFBP5     ssc-miR-874 -0.10098576   driploss
## 168  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982   driploss
## 169  15 XLOC_009656   SP110     ssc-miR-874 -0.14757211   driploss
## 170  15 XLOC_009701  KLHL30     ssc-miR-874 -0.12902519   driploss
## 173  15 XLOC_009022  IGFBP5     ssc-miR-874 -0.12508215   overtend
## 174  15 XLOC_009584  IGFBP5     ssc-miR-874 -0.10098576   overtend
## 175  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982   overtend
## 177  15 XLOC_009022  IGFBP5     ssc-miR-874 -0.12508215 tenderness
## 178  15 XLOC_009584  IGFBP5     ssc-miR-874 -0.10098576 tenderness
## 179  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982 tenderness
## 180  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982  juiciness
## 181  15 XLOC_009614   PTPRN     ssc-miR-874 -0.11631982        WBS
```

```r
as.character(unique(negcoloc$genes))
```

```
##  [1] "FAM78A"  "PPP1CA"  "POLD4"   "OVOL1"   "TRMT112" "CCDC85B" "VSNL1"  
##  [8] "SMPDL3B" "ZSCAN20" "RRP36"   "KLHL31"  "PFDN6"   "LEMD2"   "BYSL"   
## [15] "ETV7"    "STOML1"  "KHNYN"   "RGMA"    "HECTD1"  "VASH1"   "NOP9"   
## [22] "PROX2"   "USP12"   "WDR33"   "IGFBP5"  "PRPF40A" "PTPRN"   "SP110"  
## [29] "KLHL30"
```

## Save data


```r
save(coloc, negcoloc, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/12_target_mrna_colocalized_pqtl.Rdata")
```

