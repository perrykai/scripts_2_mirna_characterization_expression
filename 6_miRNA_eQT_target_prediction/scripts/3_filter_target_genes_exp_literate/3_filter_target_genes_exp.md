**Script:** `3_filter_target_genes_exp.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/scripts`

**Date:**  06/26/17

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/`

2. `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/`

**Input File(s):** 

1. `4_filtered_target_prediction_results.Rdata`

2. `data.Rdata`

**Output File Directory:** 

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/`

**Output File(s):** 

`5_filtered_targets_expression_results.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to filter the target prediction results based on genes expressed in the 168 LD samples DV obtained in her analysis.

The filtered list will dictate which genes will be included in the correlation analysis and pQTL co-localization.

## Install libraries


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/scripts")

rm(list=ls())

library(limma)
library(edgeR)
library(methods)
library(biomaRt)
```

## Load data

Load the list of filtered target prediction results:


```r
load("../4_filtered_target_prediction_results.Rdata")
```

Load DV's dge object to obtain her annotation file:


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/data.Rdata")
```

## Analysis



```r
genes<-dge$genes
dim(genes)
```

```
## [1] 18839     6
```

```r
head(genes)
```

```
##             chr0   start     end strand      geneID              genes
## XLOC_000001    1  280458  290229      + XLOC_000001               DLL1
## XLOC_000003    1  643982  645242      + XLOC_000003               <NA>
## XLOC_000005    1  935012  955845      + XLOC_000005 ENSSSCG00000030155
## XLOC_000007    1  987441 1175801      + XLOC_000007 ENSSSCG00000004008
## XLOC_000011    1 2447572 2457750      + XLOC_000011              DACT2
## XLOC_000012    1 2657897 2672308      + XLOC_000012              FRMD1
```

How many genes have no name?


```r
sum(is.na(genes$genes))
```

```
## [1] 3854
```

How many genes have ensembl IDs, but no gene symbol?


```r
length(grep("ENSSSCG", as.character(genes$genes)))
```

```
## [1] 3610
```

Filter out the NA genes and the ensembl IDs:


```r
genes<-genes[!is.na(genes$genes),]
dim(genes)
```

```
## [1] 14985     6
```

```r
genes<-genes[!grepl("ENSSSCG", as.character(genes$genes)),]
```

Dimensions of retained genes object:


```r
dim(genes)
```

```
## [1] 11375     6
```

```r
head(genes)
```

```
##             chr0   start     end strand      geneID   genes
## XLOC_000001    1  280458  290229      + XLOC_000001    DLL1
## XLOC_000011    1 2447572 2457750      + XLOC_000011   DACT2
## XLOC_000012    1 2657897 2672308      + XLOC_000012   FRMD1
## XLOC_000015    1 3405126 3662244      + XLOC_000015 RPS6KA2
## XLOC_000017    1 3716401 3734041      + XLOC_000017  SFT2D1
## XLOC_000021    1 4635847 4742908      + XLOC_000021  PDE10A
```

How many unique genes:


```r
length(unique(genes$genes))
```

```
## [1] 10392
```

```r
tst<-targets.unique$"miR-128-3p"
head(tst)
```

```
##           a_Gene_ID miRNA_family_ID Site_type ensid_version
## 25  ENST00000002165      miR-128-3p   8mer-1a             6
## 62  ENST00000003302      miR-128-3p   7mer-m8             4
## 75  ENST00000004921      miR-128-3p   8mer-1a             3
## 164 ENST00000007722      miR-128-3p   7mer-m8             7
## 211 ENST00000012049      miR-128-3p   8mer-1a             5
## 230 ENST00000012443      miR-128-3p   7mer-m8             4
##     ensembl_gene_id ensembl_transc_id external_gene_name hgnc_symbol
## 25  ENSG00000001036   ENST00000002165              FUCA2       FUCA2
## 62  ENSG00000048028   ENST00000003302              USP28       USP28
## 75             <NA>              <NA>               <NA>        <NA>
## 164 ENSG00000005884   ENST00000007722              ITGA3       ITGA3
## 211 ENSG00000011478   ENST00000012049              QPCTL       QPCTL
## 230 ENSG00000011485   ENST00000012443              PPP5C       PPP5C
##     status
## 25   KNOWN
## 62   KNOWN
## 75    <NA>
## 164  KNOWN
## 211  KNOWN
## 230  KNOWN
```

```r
sum(tst$external_gene_name %in% genes$genes)
```

```
## [1] 1558
```

```r
dim(tst[tst$external_gene_name %in% genes$genes,])
```

```
## [1] 1558    9
```

Create a list of miRNA names to loop through:


```r
miRnm<-as.character(names(targets.unique))
miRnm
```

```
##  [1] "miR-200-3p/429"   "let-7-5p/98-5p"   "miR-128-3p"      
##  [4] "miR-140-5p"       "miR-6821-3p"      "miR-6888-3p"     
##  [7] "miR-874-3p"       "miR-345-3p"       "miR-6072/6891-3p"
## [10] "miR-1306-3p"      "miR-184"          "miR-190-5p"      
## [13] "miR-1468-5p"      "miR-95-3p"
```

Filter the genes targets based on those expressed in the dataset:


```r
targets.exp<-list()
targets.exp.sum<-list()
for(i in miRnm){
targets.exp[[i]]<-targets.unique[[i]][targets.unique[[i]]$external_gene_name %in% genes$genes,]
# Create summary file:
targets.exp.sum[[i]]$summary<-data.frame(
	# How many targets were input into biomart? 
	gene.input=nrow(targets.unique[[i]]),
	# How many targets were expressed?
	gene.output=nrow(targets.exp[[i]]),
	# What fraction of the targets input were expressed/retained?
	gene.prop=nrow(targets.exp[[i]])/nrow(targets.unique[[i]]))
}

str(targets.exp)
```

```
## List of 14
##  $ miR-200-3p/429  :'data.frame':	1297 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1297] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ miRNA_family_ID   : chr [1:1297] "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" ...
##   ..$ Site_type         : chr [1:1297] "7mer-m8" "8mer-1a" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1297] "3" "5" "8" "4" ...
##   ..$ ensembl_gene_id   : chr [1:1297] "ENSG00000003056" "ENSG00000002587" "ENSG00000001630" "ENSG00000075790" ...
##   ..$ ensembl_transc_id : chr [1:1297] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ external_gene_name: chr [1:1297] "M6PR" "HS3ST1" "CYP51A1" "BCAP29" ...
##   ..$ hgnc_symbol       : chr [1:1297] "M6PR" "HS3ST1" "CYP51A1" "BCAP29" ...
##   ..$ status            : chr [1:1297] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ let-7-5p/98-5p  :'data.frame':	768 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:768] "ENST00000001008" "ENST00000002829" "ENST00000005386" "ENST00000007699" ...
##   ..$ miRNA_family_ID   : chr [1:768] "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" ...
##   ..$ Site_type         : chr [1:768] "7mer-m8" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:768] "4" "3" "3" "5" ...
##   ..$ ensembl_gene_id   : chr [1:768] "ENSG00000004478" "ENSG00000001617" "ENSG00000005175" "ENSG00000006047" ...
##   ..$ ensembl_transc_id : chr [1:768] "ENST00000001008" "ENST00000002829" "ENST00000005386" "ENST00000007699" ...
##   ..$ external_gene_name: chr [1:768] "FKBP4" "SEMA3F" "RPAP3" "YBX2" ...
##   ..$ hgnc_symbol       : chr [1:768] "FKBP4" "SEMA3F" "RPAP3" "YBX2" ...
##   ..$ status            : chr [1:768] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-128-3p      :'data.frame':	1558 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1558] "ENST00000002165" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ miRNA_family_ID   : chr [1:1558] "miR-128-3p" "miR-128-3p" "miR-128-3p" "miR-128-3p" ...
##   ..$ Site_type         : chr [1:1558] "8mer-1a" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1558] "6" "4" "7" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1558] "ENSG00000001036" "ENSG00000048028" "ENSG00000005884" "ENSG00000011478" ...
##   ..$ ensembl_transc_id : chr [1:1558] "ENST00000002165" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ external_gene_name: chr [1:1558] "FUCA2" "USP28" "ITGA3" "QPCTL" ...
##   ..$ hgnc_symbol       : chr [1:1558] "FUCA2" "USP28" "ITGA3" "QPCTL" ...
##   ..$ status            : chr [1:1558] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-140-5p      :'data.frame':	812 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:812] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000025301" ...
##   ..$ miRNA_family_ID   : chr [1:812] "miR-140-5p" "miR-140-5p" "miR-140-5p" "miR-140-5p" ...
##   ..$ Site_type         : chr [1:812] "7mer-m8" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:812] "5" "3" "3" "2" ...
##   ..$ ensembl_gene_id   : chr [1:812] "ENSG00000002587" "ENSG00000001461" "ENSG00000005175" "ENSG00000023516" ...
##   ..$ ensembl_transc_id : chr [1:812] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000025301" ...
##   ..$ external_gene_name: chr [1:812] "HS3ST1" "NIPAL3" "RPAP3" "AKAP11" ...
##   ..$ hgnc_symbol       : chr [1:812] "HS3ST1" "NIPAL3" "RPAP3" "AKAP11" ...
##   ..$ status            : chr [1:812] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6821-3p     :'data.frame':	593 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:593] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000024061" ...
##   ..$ miRNA_family_ID   : chr [1:593] "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" ...
##   ..$ Site_type         : chr [1:593] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:593] "5" "4" "7" "3" ...
##   ..$ ensembl_gene_id   : chr [1:593] "ENSG00000002587" "ENSG00000075790" "ENSG00000010270" "ENSG00000022567" ...
##   ..$ ensembl_transc_id : chr [1:593] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000024061" ...
##   ..$ external_gene_name: chr [1:593] "HS3ST1" "BCAP29" "STARD3NL" "SLC45A4" ...
##   ..$ hgnc_symbol       : chr [1:593] "HS3ST1" "BCAP29" "STARD3NL" "SLC45A4" ...
##   ..$ status            : chr [1:593] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6888-3p     :'data.frame':	1481 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1481] "ENST00000002596" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ miRNA_family_ID   : chr [1:1481] "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" ...
##   ..$ Site_type         : chr [1:1481] "7mer-m8" "8mer-1a" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1481] "5" "4" "7" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1481] "ENSG00000002587" "ENSG00000048028" "ENSG00000005884" "ENSG00000011478" ...
##   ..$ ensembl_transc_id : chr [1:1481] "ENST00000002596" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ external_gene_name: chr [1:1481] "HS3ST1" "USP28" "ITGA3" "QPCTL" ...
##   ..$ hgnc_symbol       : chr [1:1481] "HS3ST1" "USP28" "ITGA3" "QPCTL" ...
##   ..$ status            : chr [1:1481] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-874-3p      :'data.frame':	1119 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1119] "ENST00000002596" "ENST00000029410" "ENST00000046640" "ENST00000062104" ...
##   ..$ miRNA_family_ID   : chr [1:1119] "miR-874-3p" "miR-874-3p" "miR-874-3p" "miR-874-3p" ...
##   ..$ Site_type         : chr [1:1119] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1119] "5" "5" "3" "2" ...
##   ..$ ensembl_gene_id   : chr [1:1119] "ENSG00000002587" "ENSG00000027847" "ENSG00000040531" "ENSG00000053438" ...
##   ..$ ensembl_transc_id : chr [1:1119] "ENST00000002596" "ENST00000029410" "ENST00000046640" "ENST00000062104" ...
##   ..$ external_gene_name: chr [1:1119] "HS3ST1" "B4GALT7" "CTNS" "NNAT" ...
##   ..$ hgnc_symbol       : chr [1:1119] "HS3ST1" "B4GALT7" "CTNS" "NNAT" ...
##   ..$ status            : chr [1:1119] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-345-3p      :'data.frame':	1052 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1052] "ENST00000002829" "ENST00000005259" "ENST00000011653" "ENST00000025301" ...
##   ..$ miRNA_family_ID   : chr [1:1052] "miR-345-3p" "miR-345-3p" "miR-345-3p" "miR-345-3p" ...
##   ..$ Site_type         : chr [1:1052] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1052] "3" "4" "4" "2" ...
##   ..$ ensembl_gene_id   : chr [1:1052] "ENSG00000001617" "ENSG00000075790" "ENSG00000010610" "ENSG00000023516" ...
##   ..$ ensembl_transc_id : chr [1:1052] "ENST00000002829" "ENST00000005259" "ENST00000011653" "ENST00000025301" ...
##   ..$ external_gene_name: chr [1:1052] "SEMA3F" "BCAP29" "CD4" "AKAP11" ...
##   ..$ hgnc_symbol       : chr [1:1052] "SEMA3F" "BCAP29" "CD4" "AKAP11" ...
##   ..$ status            : chr [1:1052] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6072/6891-3p:'data.frame':	901 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:901] "ENST00000003912" "ENST00000007722" "ENST00000012443" "ENST00000014914" ...
##   ..$ miRNA_family_ID   : chr [1:901] "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" ...
##   ..$ Site_type         : chr [1:901] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:901] "3" "7" "4" "5" ...
##   ..$ ensembl_gene_id   : chr [1:901] "ENSG00000001461" "ENSG00000005884" "ENSG00000011485" "ENSG00000013588" ...
##   ..$ ensembl_transc_id : chr [1:901] "ENST00000003912" "ENST00000007722" "ENST00000012443" "ENST00000014914" ...
##   ..$ external_gene_name: chr [1:901] "NIPAL3" "ITGA3" "PPP5C" "GPRC5A" ...
##   ..$ hgnc_symbol       : chr [1:901] "NIPAL3" "ITGA3" "PPP5C" "GPRC5A" ...
##   ..$ status            : chr [1:901] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-1306-3p     :'data.frame':	141 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:141] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000222249" ...
##   ..$ miRNA_family_ID   : chr [1:141] "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" ...
##   ..$ Site_type         : chr [1:141] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:141] "8" "4" "3" "9" ...
##   ..$ ensembl_gene_id   : chr [1:141] "ENSG00000006453" "ENSG00000088256" "ENSG00000100344" "ENSG00000105642" ...
##   ..$ ensembl_transc_id : chr [1:141] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000222249" ...
##   ..$ external_gene_name: chr [1:141] "BAIAP2L1" "GNA11" "PNPLA3" "KCNN1" ...
##   ..$ hgnc_symbol       : chr [1:141] "BAIAP2L1" "GNA11" "PNPLA3" "KCNN1" ...
##   ..$ status            : chr [1:141] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-184         :'data.frame':	179 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:179] "ENST00000200181" "ENST00000214869" "ENST00000220772" "ENST00000222812" ...
##   ..$ miRNA_family_ID   : chr [1:179] "miR-184" "miR-184" "miR-184" "miR-184" ...
##   ..$ Site_type         : chr [1:179] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:179] "3" "2" "3" "3" ...
##   ..$ ensembl_gene_id   : chr [1:179] "ENSG00000132470" "ENSG00000099203" "ENSG00000104332" "ENSG00000106089" ...
##   ..$ ensembl_transc_id : chr [1:179] "ENST00000200181" "ENST00000214869" "ENST00000220772" "ENST00000222812" ...
##   ..$ external_gene_name: chr [1:179] "ITGB4" "TMED1" "SFRP1" "STX1A" ...
##   ..$ hgnc_symbol       : chr [1:179] "ITGB4" "TMED1" "SFRP1" "STX1A" ...
##   ..$ status            : chr [1:179] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-190-5p      :'data.frame':	520 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:520] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000192788" ...
##   ..$ miRNA_family_ID   : chr [1:520] "miR-190-5p" "miR-190-5p" "miR-190-5p" "miR-190-5p" ...
##   ..$ Site_type         : chr [1:520] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:520] "5" "4" "3" "5" ...
##   ..$ ensembl_gene_id   : chr [1:520] "ENSG00000014919" "ENSG00000054965" "ENSG00000005001" "ENSG00000065060" ...
##   ..$ ensembl_transc_id : chr [1:520] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000192788" ...
##   ..$ external_gene_name: chr [1:520] "COX15" "FAM168A" "PRSS22" "UHRF1BP1" ...
##   ..$ hgnc_symbol       : chr [1:520] "COX15" "FAM168A" "PRSS22" "UHRF1BP1" ...
##   ..$ status            : chr [1:520] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-1468-5p     :'data.frame':	198 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:198] "ENST00000173229" "ENST00000215838" "ENST00000219097" "ENST00000219252" ...
##   ..$ miRNA_family_ID   : chr [1:198] "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" ...
##   ..$ Site_type         : chr [1:198] "8mer-1a" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:198] "2" "3" "2" "5" ...
##   ..$ ensembl_gene_id   : chr [1:198] "ENSG00000065320" "ENSG00000185339" "ENSG00000091651" "ENSG00000102978" ...
##   ..$ ensembl_transc_id : chr [1:198] "ENST00000173229" "ENST00000215838" "ENST00000219097" "ENST00000219252" ...
##   ..$ external_gene_name: chr [1:198] "NTN1" "TCN2" "ORC6" "POLR2C" ...
##   ..$ hgnc_symbol       : chr [1:198] "NTN1" "TCN2" "ORC6" "POLR2C" ...
##   ..$ status            : chr [1:198] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-95-3p       :'data.frame':	105 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:105] "ENST00000216797" "ENST00000247461" "ENST00000254667" "ENST00000258538" ...
##   ..$ miRNA_family_ID   : chr [1:105] "miR-95-3p" "miR-95-3p" "miR-95-3p" "miR-95-3p" ...
##   ..$ Site_type         : chr [1:105] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:105] "5" "4" "3" "3" ...
##   ..$ ensembl_gene_id   : chr [1:105] "ENSG00000100906" "ENSG00000127022" "ENSG00000132334" "ENSG00000136052" ...
##   ..$ ensembl_transc_id : chr [1:105] "ENST00000216797" "ENST00000247461" "ENST00000254667" "ENST00000258538" ...
##   ..$ external_gene_name: chr [1:105] "NFKBIA" "CANX" "PTPRE" "SLC41A2" ...
##   ..$ hgnc_symbol       : chr [1:105] "NFKBIA" "CANX" "PTPRE" "SLC41A2" ...
##   ..$ status            : chr [1:105] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
```

```r
targets.exp.sum
```

```
## $`miR-200-3p/429`
## $`miR-200-3p/429`$summary
##   gene.input gene.output gene.prop
## 1       2301        1297  0.563668
## 
## 
## $`let-7-5p/98-5p`
## $`let-7-5p/98-5p`$summary
##   gene.input gene.output gene.prop
## 1       1398         768 0.5493562
## 
## 
## $`miR-128-3p`
## $`miR-128-3p`$summary
##   gene.input gene.output gene.prop
## 1       2982        1558 0.5224681
## 
## 
## $`miR-140-5p`
## $`miR-140-5p`$summary
##   gene.input gene.output gene.prop
## 1       1518         812 0.5349144
## 
## 
## $`miR-6821-3p`
## $`miR-6821-3p`$summary
##   gene.input gene.output gene.prop
## 1       1118         593 0.5304114
## 
## 
## $`miR-6888-3p`
## $`miR-6888-3p`$summary
##   gene.input gene.output gene.prop
## 1       2816        1481 0.5259233
## 
## 
## $`miR-874-3p`
## $`miR-874-3p`$summary
##   gene.input gene.output gene.prop
## 1       2175        1119 0.5144828
## 
## 
## $`miR-345-3p`
## $`miR-345-3p`$summary
##   gene.input gene.output gene.prop
## 1       1950        1052 0.5394872
## 
## 
## $`miR-6072/6891-3p`
## $`miR-6072/6891-3p`$summary
##   gene.input gene.output gene.prop
## 1       1680         901 0.5363095
## 
## 
## $`miR-1306-3p`
## $`miR-1306-3p`$summary
##   gene.input gene.output gene.prop
## 1        264         141 0.5340909
## 
## 
## $`miR-184`
## $`miR-184`$summary
##   gene.input gene.output gene.prop
## 1        365         179  0.490411
## 
## 
## $`miR-190-5p`
## $`miR-190-5p`$summary
##   gene.input gene.output gene.prop
## 1        979         520 0.5311542
## 
## 
## $`miR-1468-5p`
## $`miR-1468-5p`$summary
##   gene.input gene.output gene.prop
## 1        394         198 0.5025381
## 
## 
## $`miR-95-3p`
## $`miR-95-3p`$summary
##   gene.input gene.output gene.prop
## 1        196         105 0.5357143
```

---

Extract the gene symbols from the full expressed dataset to convert them to their ensembl gene IDs for use in DAVID:



```r
length(unique(genes$genes))
```

```
## [1] 10392
```

## Visualize
## Save data


```r
save(targets.exp, targets.exp.sum, file="../5_filtered_targets_expression_results.Rdata")
write.table(unique(genes$genes), file="../6_DAVID_background_gene_names_expressed.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

