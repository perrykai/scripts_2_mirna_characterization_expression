**Script:** `6_summarize_susscrofa_ensembl_blast_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Date:**  `6/14/16 UPDATED 7/5/16`

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `1_susscrofa_ensembl_blastn_output_e5.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** 

1. `1a_susscrofa_filtered_ensembl_blastn_results.Rdata`
2. `1b_susscrofa_filtered_uniqseq_ensembl_blastn_results.csv`
3. `1c_susscrofa_filtered_totseq_ensembl_blastn_results.csv`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
This script filters and summarizes the blast results from the Sus scrofa Ensembl database query.

First, if there is only one hit for a sequence, return that hit.
If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
Then, retain matches with > 96% sequence identity.
Finally, the first filter returns sequences with the maximum bitscore.

The second filter returns the sequence with greatest percent match to blast hit.
Finally, if there are still remaining sequences with multiple hits, retain only the first hit. 

## Install libraries



```r
library(biomaRt)
```

```
## Loading required package: methods
```

```r
library(parallel)
library(stringr)
rm(list=ls())
```

## Load data



```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts")

system.time({
blastresults<-read.table("../1_susscrofa_ensembl_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   5.687   0.071   5.808
```

```r
dim(blastresults)
```

```
## [1] 316347     12
```

```r
head(blastresults)
```

```
##               query_id             dbseq_id perc_identical length mismatch
## 1 seq_81539680_x686919 ENSSSCT00000024283.1            100     22        0
## 2 seq_82226599_x674294 ENSSSCT00000021744.1            100     21        0
## 3 seq_90218330_x552719 ENSSSCT00000021744.1            100     24        0
## 4 seq_92386826_x527823 ENSSSCT00000021744.1            100     21        0
## 5 seq_93438878_x519473 ENSSSCT00000021744.1            100     22        0
## 6 seq_93958351_x515658 ENSSSCT00000021744.1            100     23        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22           1        22  4e-07     44.1
## 2       0           1        21          63        83  1e-06     42.1
## 3       0           1        24          60        83  3e-08     48.1
## 4       0           1        21          61        81  1e-06     42.1
## 5       0           1        22          60        81  4e-07     44.1
## 6       0           1        23          61        83  1e-07     46.1
```

```r
str(blastresults)
```

```
## 'data.frame':	316347 obs. of  12 variables:
##  $ query_id      : Factor w/ 132139 levels "seq_104712401_x365907",..: 132134 132135 132136 132137 132138 132139 1 2 3 4 ...
##  $ dbseq_id      : Factor w/ 2044 levels "ENSSSCT00000001325.2",..: 1345 1030 1030 1030 1030 1030 1030 25 1030 1593 ...
##  $ perc_identical: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ length        : int  22 21 24 21 22 23 20 25 22 22 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  22 21 24 21 22 23 20 25 22 22 ...
##  $ dbseq_start   : int  1 63 60 61 60 61 62 1 62 51 ...
##  $ dbseq_end     : int  22 83 83 81 81 83 81 25 83 72 ...
##  $ evalue        : num  4e-07 1e-06 3e-08 1e-06 4e-07 1e-07 5e-06 7e-09 4e-07 4e-07 ...
##  $ bitscore      : num  44.1 42.1 48.1 42.1 44.1 46.1 40.1 50.1 44.1 44.1 ...
```

## Analysis



```r
blastresults$dbseq_id<-as.character(blastresults$dbseq_id)
blastresults$query_id<-as.character(blastresults$query_id)
```

Extract the ensembl id to another column:

Move the version numbers (".1", ".2" or ".3") after the ensembl ids to another column for each sequence to make it compatible with biomaRt

See website http://useast.ensembl.org/Help/View?id=181 for details on the stability of ensembl IDs


```r
blastresults$ensid<-as.character(blastresults$dbseq_id)
blastresults$ensid<-gsub("\\..*","",blastresults$ensid)
blastresults$ensid_version <- as.character(lapply(strsplit(as.character(blastresults$dbseq_id), '.', fixed = TRUE), "[", 2))

head(blastresults)
```

```
##               query_id             dbseq_id perc_identical length mismatch
## 1 seq_81539680_x686919 ENSSSCT00000024283.1            100     22        0
## 2 seq_82226599_x674294 ENSSSCT00000021744.1            100     21        0
## 3 seq_90218330_x552719 ENSSSCT00000021744.1            100     24        0
## 4 seq_92386826_x527823 ENSSSCT00000021744.1            100     21        0
## 5 seq_93438878_x519473 ENSSSCT00000021744.1            100     22        0
## 6 seq_93958351_x515658 ENSSSCT00000021744.1            100     23        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22           1        22  4e-07     44.1
## 2       0           1        21          63        83  1e-06     42.1
## 3       0           1        24          60        83  3e-08     48.1
## 4       0           1        21          61        81  1e-06     42.1
## 5       0           1        22          60        81  4e-07     44.1
## 6       0           1        23          61        83  1e-07     46.1
##                ensid ensid_version
## 1 ENSSSCT00000024283             1
## 2 ENSSSCT00000021744             1
## 3 ENSSSCT00000021744             1
## 4 ENSSSCT00000021744             1
## 5 ENSSSCT00000021744             1
## 6 ENSSSCT00000021744             1
```

```r
head(table(blastresults$ensid))
```

```
## 
## ENSSSCT00000001325 ENSSSCT00000001326 ENSSSCT00000001333 
##                 19                 17                 14 
## ENSSSCT00000001335 ENSSSCT00000001337 ENSSSCT00000001522 
##                 13                119                  3
```

Check how many unique small RNA sequences had ensembl hits:


```r
length(unique(blastresults$query_id))
```

```
## [1] 132139
```

Check how many unique ensembl database sequences had hits:


```r
length(unique(blastresults$dbseq_id))
```

```
## [1] 2044
```

Use biomaRt package (v/2.20.0) to obtain annotation information for Ensembl dataset BLAST hits:


```r
mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="dec2015.archive.ensembl.org", dataset="sscrofa_gene_ensembl")

system.time({
rst <- getBM(
attributes = c("ensembl_transcript_id","chromosome_name", "ensembl_gene_id", "rfam", "rfam_transcript_name", "version", "transcript_version", "transcript_source", "status", "transcript_status", "gene_biotype"),
filters = "ensembl_transcript_id", 
values=blastresults$ensid, 
mart = mart)
})
```

```
##    user  system elapsed 
##   0.230   0.009   7.875
```

```r
dim(rst)
```

```
## [1] 2044   11
```

```r
head(rst)
```

```
##   ensembl_transcript_id chromosome_name    ensembl_gene_id rfam
## 1    ENSSSCT00000001325               7 ENSSSCG00000001227     
## 2    ENSSSCT00000001326               7 ENSSSCG00000001227     
## 3    ENSSSCT00000001333               7 ENSSSCG00000001228     
## 4    ENSSSCT00000001335               7 ENSSSCG00000001397     
## 5    ENSSSCT00000001337               7 ENSSSCG00000001396     
## 6    ENSSSCT00000001522               7 ENSSSCG00000001396     
##   rfam_transcript_name version transcript_version transcript_source status
## 1                            3                  2           ensembl  KNOWN
## 2                            3                  2           ensembl  KNOWN
## 3                            3                  2           ensembl  KNOWN
## 4                            3                  2           ensembl  KNOWN
## 5                            3                  2           ensembl  KNOWN
## 6                            3                  2           ensembl  KNOWN
##   transcript_status         gene_biotype
## 1             KNOWN processed_transcript
## 2             KNOWN processed_transcript
## 3             KNOWN processed_transcript
## 4             KNOWN processed_transcript
## 5             KNOWN processed_transcript
## 6             KNOWN processed_transcript
```

```r
tail(rst)
```

```
##      ensembl_transcript_id chromosome_name    ensembl_gene_id rfam
## 2039    ENSSSCT00000035929               7 ENSSSCG00000001727     
## 2040    ENSSSCT00000035935               X ENSSSCG00000024809     
## 2041    ENSSSCT00000036045               1 ENSSSCG00000004943     
## 2042    ENSSSCT00000036444              17 ENSSSCG00000007280     
## 2043    ENSSSCT00000036555               1 ENSSSCG00000004334     
## 2044    ENSSSCT00000036566              13 ENSSSCG00000031014     
##      rfam_transcript_name version transcript_version transcript_source
## 2039                            3                  1            havana
## 2040                            2                  1            havana
## 2041                            3                  1            havana
## 2042                            3                  1            havana
## 2043                            3                  1            havana
## 2044                            1                  1            havana
##      status transcript_status         gene_biotype
## 2039  KNOWN             KNOWN processed_transcript
## 2040  KNOWN             KNOWN processed_transcript
## 2041  KNOWN             KNOWN processed_transcript
## 2042  KNOWN             KNOWN processed_transcript
## 2043  KNOWN             KNOWN processed_transcript
## 2044  KNOWN             KNOWN              lincRNA
```

```r
table(rst$rfam)
```

```
## 
##         RF00001 RF00002 RF00003 RF00004 RF00006 RF00007 RF00009 RF00012 
##     493     150       2      53      31       2       1       4      15 
## RF00015 RF00016 RF00017 RF00019 RF00020 RF00024 RF00026 RF00030 RF00045 
##      11       9      17       9       9       1     684       1       2 
## RF00049 RF00054 RF00055 RF00056 RF00067 RF00068 RF00069 RF00070 RF00071 
##       2       1       1       3       2       1       1       2       2 
## RF00072 RF00085 RF00086 RF00087 RF00088 RF00089 RF00090 RF00091 RF00092 
##       2       1       1       1       1       1       3       4       4 
## RF00093 RF00096 RF00099 RF00100 RF00108 RF00133 RF00134 RF00136 RF00137 
##       3       4       7     110       5       4       1       1       2 
## RF00138 RF00139 RF00147 RF00150 RF00151 RF00152 RF00153 RF00154 RF00155 
##       3       3       2       2       3       1       2       1       1 
## RF00156 RF00157 RF00158 RF00181 RF00186 RF00187 RF00188 RF00190 RF00191 
##      20       2       1       9       1       1       3       2       1 
## RF00211 RF00213 RF00217 RF00218 RF00221 RF00231 RF00263 RF00264 RF00265 
##       3       3       2       1       1       1       1       1       1 
## RF00266 RF00271 RF00272 RF00273 RF00274 RF00275 RF00276 RF00277 RF00278 
##       3       2       3       2       2       2       1       1       1 
## RF00279 RF00281 RF00283 RF00284 RF00286 RF00288 RF00289 RF00302 RF00319 
##       3       1       1       1       1       1       1       2       2 
## RF00322 RF00324 RF00325 RF00334 RF00340 RF00341 RF00342 RF00377 RF00392 
##      23       1       3       4       2       1       1       2       4 
## RF00393 RF00394 RF00396 RF00397 RF00398 RF00399 RF00401 RF00402 RF00403 
##       1       1       1       2       2       1       2       2       1 
## RF00404 RF00405 RF00406 RF00407 RF00408 RF00409 RF00410 RF00411 RF00412 
##       1       1       4       2       1       1       4       2       7 
## RF00413 RF00414 RF00415 RF00416 RF00418 RF00419 RF00420 RF00421 RF00422 
##      14       3       2       1       5       1       3       1       1 
## RF00424 RF00425 RF00426 RF00428 RF00429 RF00430 RF00431 RF00432 RF00438 
##       1       7       2       1       3       2       1       2       3 
## RF00439 RF00440 RF00443 RF00476 RF00478 RF00494 RF00548 RF00553 RF00554 
##       1       1       4       2       2       1       3       1       4 
## RF00560 RF00561 RF00563 RF00564 RF00567 RF00568 RF00569 RF00570 RF00571 
##       2       5       1       1       1       2       1       1       1 
## RF00572 RF00573 RF00574 RF00575 RF00576 RF00577 RF00578 RF00579 RF00581 
##       1       1       1       3       1       1       1       2       3 
## RF00582 RF00584 RF00586 RF00588 RF00591 RF00592 RF00593 RF00594 RF00598 
##       1       1       1       2       2       1       1       2       1 
## RF00599 RF00600 RF00601 RF00602 RF00603 RF00604 RF00606 RF00608 RF00609 
##       1       3       1       1       1       4       1       2       1 
## RF00610 RF00611 RF00613 RF00614 RF00618 RF00619 RF01156 RF01161 RF01164 
##       2       2       1       9       1      11       2       1       1 
## RF01169 RF01170 RF01173 RF01183 RF01186 RF01188 RF01191 RF01192 RF01200 
##       3       1       1       1       1       1       2       2       1 
## RF01210 RF01211 RF01229 RF01233 RF01241 RF01268 RF01277 RF01290 RF01291 
##      14       2       2       1       1       2       1       2       1 
## RF01294 RF01299 
##       1       3
```

```r
table(rst$gene_biotype)
```

```
## 
##              lincRNA                miRNA             misc_RNA 
##                   12                  404                  144 
##              Mt_rRNA              Mt_tRNA processed_transcript 
##                    2                   22                   53 
##                 rRNA               snoRNA                snRNA 
##                  152                  451                  804
```

```r
sum(rst$rfam == "")
```

```
## [1] 493
```

```r
sum(rst$gene_biotype == "")
```

```
## [1] 0
```

Use match function to obtain the gene_biotype and gene status information for each ensembl match and add it as a column to blastresults:


```r
blastresults$gene_biotype<-as.character(rst[match(blastresults$ensid, rst$ensembl_transcript_id), "gene_biotype"])
blastresults$status<-rst[match(blastresults$ensid, rst$ensembl_transcript_id), "status"]

dim(blastresults)
```

```
## [1] 316347     16
```

```r
head(blastresults)
```

```
##               query_id             dbseq_id perc_identical length mismatch
## 1 seq_81539680_x686919 ENSSSCT00000024283.1            100     22        0
## 2 seq_82226599_x674294 ENSSSCT00000021744.1            100     21        0
## 3 seq_90218330_x552719 ENSSSCT00000021744.1            100     24        0
## 4 seq_92386826_x527823 ENSSSCT00000021744.1            100     21        0
## 5 seq_93438878_x519473 ENSSSCT00000021744.1            100     22        0
## 6 seq_93958351_x515658 ENSSSCT00000021744.1            100     23        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22           1        22  4e-07     44.1
## 2       0           1        21          63        83  1e-06     42.1
## 3       0           1        24          60        83  3e-08     48.1
## 4       0           1        21          61        81  1e-06     42.1
## 5       0           1        22          60        81  4e-07     44.1
## 6       0           1        23          61        83  1e-07     46.1
##                ensid ensid_version gene_biotype status
## 1 ENSSSCT00000024283             1        miRNA  NOVEL
## 2 ENSSSCT00000021744             1         rRNA  KNOWN
## 3 ENSSSCT00000021744             1         rRNA  KNOWN
## 4 ENSSSCT00000021744             1         rRNA  KNOWN
## 5 ENSSSCT00000021744             1         rRNA  KNOWN
## 6 ENSSSCT00000021744             1         rRNA  KNOWN
```

```r
tail(blastresults)
```

```
##                query_id             dbseq_id perc_identical length
## 316342 seq_227425488_x2 ENSSSCT00000020202.1            100     22
## 316343 seq_227425488_x2 ENSSSCT00000020550.1            100     22
## 316344 seq_227425510_x2 ENSSSCT00000019658.3            100     21
## 316345 seq_227425644_x2 ENSSSCT00000019690.3            100     21
## 316346 seq_227425734_x2 ENSSSCT00000019669.1            100     20
## 316347 seq_227425774_x2 ENSSSCT00000019690.3            100     21
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316342        0       0           1        22           1        22  5e-07
## 316343        0       0           1        22           1        22  5e-07
## 316344        0       0           1        21        1140      1160  1e-06
## 316345        0       0           1        21          39        59  2e-06
## 316346        0       0           1        20          22         3  7e-06
## 316347        0       0           1        21          42        22  2e-06
##        bitscore              ensid ensid_version gene_biotype status
## 316342     44.1 ENSSSCT00000020202             1       snoRNA  NOVEL
## 316343     44.1 ENSSSCT00000020550             1       snoRNA  NOVEL
## 316344     42.1 ENSSSCT00000019658             3      Mt_rRNA  NOVEL
## 316345     42.1 ENSSSCT00000019690             3      Mt_tRNA  NOVEL
## 316346     40.1 ENSSSCT00000019669             1      Mt_tRNA  NOVEL
## 316347     42.1 ENSSSCT00000019690             3      Mt_tRNA  NOVEL
```

```r
table(blastresults$gene_biotype)
```

```
## 
##              lincRNA                miRNA             misc_RNA 
##                 1297                14397                43199 
##              Mt_rRNA              Mt_tRNA processed_transcript 
##                31467                42533                 2626 
##                 rRNA               snoRNA                snRNA 
##                40920                30997               108911
```

```r
######################################################################
```

Obtained the following code from Deborah Velez, how she filtered her BLAST results for gene annotation (/mnt/research/ernstc_lab/fat_eqtl/BLAST/code/annotation_pigoligoarray.R)


```r
######################################################################
```

Format the ensembl blast output data.frame to a list to filter the sequences with multiple ensembl hits


```r
# Number of unique sequences
n <- length(unique(blastresults$query_id))
n
```

```
## [1] 132139
```

```r
# Create the sequence list to filter
system.time({
idx1 <- mclapply(as.character(unique(blastresults$query_id))[1:round(n/5)], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx2 <- mclapply(as.character(unique(blastresults$query_id))[(round(n/5) + 1):(2*round(n/5))], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx3 <- mclapply(as.character(unique(blastresults$query_id))[((2*round(n/5)) + 1):(3*round(n/5))], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx4 <- mclapply(as.character(unique(blastresults$query_id))[((3*round(n/5)) + 1):(4*round(n/5))], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx5 <- mclapply(as.character(unique(blastresults$query_id))[((4*round(n/5)) + 1):n], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx <- c(idx1,idx2,idx3,idx4,idx5)
})
```

```
##      user    system   elapsed 
## 18042.649    82.773  2416.182
```

```r
length(idx)
```

```
## [1] 132139
```

Function to filter multiple ensembl blast hits per sequence


```r
filter <- function(seqblast){
# Sequence has only one ensembl blast hit
        if (nrow(seqblast) == 1){
                return(seqblast)
        }
# Sequence with 100% match to ensembl blast hit
        ident <- seqblast[sapply(seqblast$perc_identical, function(x) x == 100),]
                if (nrow(ident) == 1){
                        return(ident)
        }
# Sequence with greater than 96% match to ensembl blast hit
        ident2 <- seqblast[sapply(seqblast$perc_identical, function(x) x > 96),]
                if (nrow(ident2) == 1){
                        return(ident2)
        }
# Select ensembl blast hit with the greatest bitscore for a sequence
        bit <- seqblast[seqblast$bitscore == max(seqblast$bitscore),]
                        if (nrow(bit) == 1){
                        return(bit)
        }

        return(bit)
}
```

Apply filter to the sequence list


```r
system.time({
p1 <- mclapply(idx[1:round(n/5)], filter, mc.cores=10)
p2 <- mclapply(idx[(round(n/5) + 1):(2*round(n/5))], filter, mc.cores=10)
p3 <- mclapply(idx[((2*round(n/5)) + 1):(3*round(n/5))], filter, mc.cores=10)
p4 <- mclapply(idx[((3*round(n/5)) + 1):(4*round(n/5))], filter, mc.cores=10)
p5 <- mclapply(idx[((4*round(n/5)) + 1):n], filter, mc.cores=10)
stp1 <- c(p1,p2,p3,p4,p5)
})
```

```
##    user  system elapsed 
##  55.625  11.861  14.480
```

```r
length(stp1)
```

```
## [1] 132139
```

```r
head(stp1)
```

```
## [[1]]
##               query_id             dbseq_id perc_identical length mismatch
## 1 seq_81539680_x686919 ENSSSCT00000024283.1            100     22        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22           1        22  4e-07     44.1
##                ensid ensid_version gene_biotype status
## 1 ENSSSCT00000024283             1        miRNA  NOVEL
## 
## [[2]]
##               query_id             dbseq_id perc_identical length mismatch
## 2 seq_82226599_x674294 ENSSSCT00000021744.1            100     21        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 2       0           1        21          63        83  1e-06     42.1
##                ensid ensid_version gene_biotype status
## 2 ENSSSCT00000021744             1         rRNA  KNOWN
## 
## [[3]]
##               query_id             dbseq_id perc_identical length mismatch
## 3 seq_90218330_x552719 ENSSSCT00000021744.1            100     24        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 3       0           1        24          60        83  3e-08     48.1
##                ensid ensid_version gene_biotype status
## 3 ENSSSCT00000021744             1         rRNA  KNOWN
## 
## [[4]]
##               query_id             dbseq_id perc_identical length mismatch
## 4 seq_92386826_x527823 ENSSSCT00000021744.1            100     21        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 4       0           1        21          61        81  1e-06     42.1
##                ensid ensid_version gene_biotype status
## 4 ENSSSCT00000021744             1         rRNA  KNOWN
## 
## [[5]]
##               query_id             dbseq_id perc_identical length mismatch
## 5 seq_93438878_x519473 ENSSSCT00000021744.1            100     22        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 5       0           1        22          60        81  4e-07     44.1
##                ensid ensid_version gene_biotype status
## 5 ENSSSCT00000021744             1         rRNA  KNOWN
## 
## [[6]]
##               query_id             dbseq_id perc_identical length mismatch
## 6 seq_93958351_x515658 ENSSSCT00000021744.1            100     23        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 6       0           1        23          61        83  1e-07     46.1
##                ensid ensid_version gene_biotype status
## 6 ENSSSCT00000021744             1         rRNA  KNOWN
```

```r
tail(stp1)
```

```
## [[1]]
##                query_id             dbseq_id perc_identical length
## 316340 seq_227425466_x2 ENSSSCT00000022687.1            100     22
## 316341 seq_227425466_x2 ENSSSCT00000021153.1            100     22
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316340        0       0           1        22          43        64  4e-07
## 316341        0       0           1        22          43        64  4e-07
##        bitscore              ensid ensid_version gene_biotype status
## 316340     44.1 ENSSSCT00000022687             1        miRNA  NOVEL
## 316341     44.1 ENSSSCT00000021153             1        miRNA  NOVEL
## 
## [[2]]
##                query_id             dbseq_id perc_identical length
## 316342 seq_227425488_x2 ENSSSCT00000020202.1            100     22
## 316343 seq_227425488_x2 ENSSSCT00000020550.1            100     22
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316342        0       0           1        22           1        22  5e-07
## 316343        0       0           1        22           1        22  5e-07
##        bitscore              ensid ensid_version gene_biotype status
## 316342     44.1 ENSSSCT00000020202             1       snoRNA  NOVEL
## 316343     44.1 ENSSSCT00000020550             1       snoRNA  NOVEL
## 
## [[3]]
##                query_id             dbseq_id perc_identical length
## 316344 seq_227425510_x2 ENSSSCT00000019658.3            100     21
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316344        0       0           1        21        1140      1160  1e-06
##        bitscore              ensid ensid_version gene_biotype status
## 316344     42.1 ENSSSCT00000019658             3      Mt_rRNA  NOVEL
## 
## [[4]]
##                query_id             dbseq_id perc_identical length
## 316345 seq_227425644_x2 ENSSSCT00000019690.3            100     21
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316345        0       0           1        21          39        59  2e-06
##        bitscore              ensid ensid_version gene_biotype status
## 316345     42.1 ENSSSCT00000019690             3      Mt_tRNA  NOVEL
## 
## [[5]]
##                query_id             dbseq_id perc_identical length
## 316346 seq_227425734_x2 ENSSSCT00000019669.1            100     20
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316346        0       0           1        20          22         3  7e-06
##        bitscore              ensid ensid_version gene_biotype status
## 316346     40.1 ENSSSCT00000019669             1      Mt_tRNA  NOVEL
## 
## [[6]]
##                query_id             dbseq_id perc_identical length
## 316347 seq_227425774_x2 ENSSSCT00000019690.3            100     21
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 316347        0       0           1        21          42        22  2e-06
##        bitscore              ensid ensid_version gene_biotype status
## 316347     42.1 ENSSSCT00000019690             3      Mt_tRNA  NOVEL
```

How many sequences have more than one hit after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 15617
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 8271 2104 1293  668  627  293  322  144   81   62   51   57   25  443   79 
##   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
##  191  112   19   25   25   38   44  113   57   43    7    8    3    5    2 
##   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46 
##   24   17   17   10   16    5    6   11   11    4    9    1   10    8   10 
##   47   48   49   50   51   52   53   54   55   56   57   58   59   60   61 
##   13    4    9   10   12    1   13   11    7   10   10   15   10   14    7 
##   62   63   64   65   66   67   68   69   71   72   73   74   79   88   93 
##    5   14    4   11   13    5   10    3    3    1    3    2    2    1    1 
##   96  102  111  116  118  123  127  129  132  136  137  140  145  147  154 
##    1    1    1    1    1    1    2    1    1    1    1    1    1    1    1 
##  155  164  165  175 
##    1    2    1    2
```

```r
filter2 <- filter <- function(seqblast){
# Sequence has only one ensembl blast hit
        if (nrow(seqblast) == 1){
                        return(seqblast)
        }
# Select sequence with greatest percent match to blast hit
        ident <- seqblast[seqblast$perc_identical == max(seqblast$perc_identical),]
                if (nrow(ident) == 1){
                        return(ident)
        }
# Select the first blast hit if there are still multiple hits for one sequence
	return(seqblast[1,])
}
```

Apply a second filter to the filtered sequence list


```r
system.time({
q1 <- mclapply(stp1[1:round(n/5)], filter2, mc.cores=10)
q2 <- mclapply(stp1[(round(n/5) + 1):(2*round(n/5))], filter2, mc.cores=10)
q3 <- mclapply(stp1[((2*round(n/5)) + 1):(3*round(n/5))], filter2, mc.cores=10)
q4 <- mclapply(stp1[((3*round(n/5)) + 1):(4*round(n/5))], filter2, mc.cores=10)
q5 <- mclapply(stp1[((4*round(n/5)) + 1):n], filter2, mc.cores=10)
stp3 <- c(q1,q2,q3,q4,q5)
})
```

```
##    user  system elapsed 
##  24.282  10.538   8.755
```

This command sums the number of sequence IDs appearing more than once in the dataset.


```r
length(stp3[unlist(lapply(stp3,nrow)) > 1])
```

```
## [1] 0
```

```r
stp4 <- stp3[unlist(lapply(stp3,nrow)) > 1]
length(stp4)
```

```
## [1] 0
```

```r
table(unlist(lapply(stp4,nrow)))
```

```
## < table of extent 0 >
```

Summary file of small RNA sequence ensembl blast results


```r
sumblast2 <- do.call(rbind,stp3)
head(sumblast2)
```

```
##               query_id             dbseq_id perc_identical length mismatch
## 1 seq_81539680_x686919 ENSSSCT00000024283.1            100     22        0
## 2 seq_82226599_x674294 ENSSSCT00000021744.1            100     21        0
## 3 seq_90218330_x552719 ENSSSCT00000021744.1            100     24        0
## 4 seq_92386826_x527823 ENSSSCT00000021744.1            100     21        0
## 5 seq_93438878_x519473 ENSSSCT00000021744.1            100     22        0
## 6 seq_93958351_x515658 ENSSSCT00000021744.1            100     23        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22           1        22  4e-07     44.1
## 2       0           1        21          63        83  1e-06     42.1
## 3       0           1        24          60        83  3e-08     48.1
## 4       0           1        21          61        81  1e-06     42.1
## 5       0           1        22          60        81  4e-07     44.1
## 6       0           1        23          61        83  1e-07     46.1
##                ensid ensid_version gene_biotype status
## 1 ENSSSCT00000024283             1        miRNA  NOVEL
## 2 ENSSSCT00000021744             1         rRNA  KNOWN
## 3 ENSSSCT00000021744             1         rRNA  KNOWN
## 4 ENSSSCT00000021744             1         rRNA  KNOWN
## 5 ENSSSCT00000021744             1         rRNA  KNOWN
## 6 ENSSSCT00000021744             1         rRNA  KNOWN
```

```r
dim(sumblast2)
```

```
## [1] 132139     16
```

```r
length(unique(sumblast2$query_id))
```

```
## [1] 132139
```

```r
table(sumblast2$gene_biotype)
```

```
## 
##              lincRNA                miRNA             misc_RNA 
##                  432                10582                12935 
##              Mt_rRNA              Mt_tRNA processed_transcript 
##                31467                42533                 1461 
##                 rRNA               snoRNA                snRNA 
##                 8450                19427                 4852
```

```r
uniqsumblast2<-as.data.frame(table(sumblast2$gene_biotype))
colnames(uniqsumblast2)<-c("Gene_Biotype", "Freq")
uniqsumblast2
```

```
##           Gene_Biotype  Freq
## 1              lincRNA   432
## 2                miRNA 10582
## 3             misc_RNA 12935
## 4              Mt_rRNA 31467
## 5              Mt_tRNA 42533
## 6 processed_transcript  1461
## 7                 rRNA  8450
## 8               snoRNA 19427
## 9                snRNA  4852
```

Add the column of sequence count to the sumblast2 data frame


```r
sumblast2$seq_count<-as.numeric(str_split_fixed(sumblast2$query_id, "_x", 2)[,2])
head(sumblast2)
```

```
##               query_id             dbseq_id perc_identical length mismatch
## 1 seq_81539680_x686919 ENSSSCT00000024283.1            100     22        0
## 2 seq_82226599_x674294 ENSSSCT00000021744.1            100     21        0
## 3 seq_90218330_x552719 ENSSSCT00000021744.1            100     24        0
## 4 seq_92386826_x527823 ENSSSCT00000021744.1            100     21        0
## 5 seq_93438878_x519473 ENSSSCT00000021744.1            100     22        0
## 6 seq_93958351_x515658 ENSSSCT00000021744.1            100     23        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22           1        22  4e-07     44.1
## 2       0           1        21          63        83  1e-06     42.1
## 3       0           1        24          60        83  3e-08     48.1
## 4       0           1        21          61        81  1e-06     42.1
## 5       0           1        22          60        81  4e-07     44.1
## 6       0           1        23          61        83  1e-07     46.1
##                ensid ensid_version gene_biotype status seq_count
## 1 ENSSSCT00000024283             1        miRNA  NOVEL    686919
## 2 ENSSSCT00000021744             1         rRNA  KNOWN    674294
## 3 ENSSSCT00000021744             1         rRNA  KNOWN    552719
## 4 ENSSSCT00000021744             1         rRNA  KNOWN    527823
## 5 ENSSSCT00000021744             1         rRNA  KNOWN    519473
## 6 ENSSSCT00000021744             1         rRNA  KNOWN    515658
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype2<-as.matrix(by(sumblast2$seq_count, sumblast2$gene_biotype, sum))
totalsumbiotype2
```

```
##                         [,1]
## lincRNA                 1305
## miRNA                3183953
## misc_RNA             1573721
## Mt_rRNA              3080094
## Mt_tRNA              2660413
## processed_transcript  246214
## rRNA                 5576806
## snoRNA               1322467
## snRNA                 234975
```

```r
if (sum(rownames(totalsumbiotype2) != uniqsumblast2$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the miRNAs and the processed_transcript fields:


```r
sum(sumblast2$seq_count[sumblast2$gene_biotype == "miRNA"])
```

```
## [1] 3183953
```

```r
sum(sumblast2$seq_count[sumblast2$gene_biotype == "processed_transcript"])
```

```
## [1] 246214
```

```r
if (sum(sumblast2$seq_count[sumblast2$gene_biotype == "miRNA"]) != totalsumbiotype2["miRNA",]) stop ("miRNA counts not equal")
if (sum(sumblast2$seq_count[sumblast2$gene_biotype == "processed_transcript"]) != totalsumbiotype2["processed_transcript",]) stop ("processed_transcript counts not equal")
```

## Save data


```r
save(sumblast2, uniqsumblast2, totalsumbiotype2, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1a_susscrofa_filtered_ensembl_blastn_results.Rdata"))
write.csv(uniqsumblast2, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1b_susscrofa_filtered_uniqseq_ensembl_blastn_results.csv")
write.csv(totalsumbiotype2, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1c_susscrofa_filtered_totseq_ensembl_blastn_results.csv", row.names=TRUE)
```

