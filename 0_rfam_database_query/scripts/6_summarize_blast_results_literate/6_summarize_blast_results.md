**Script:** `6_summarize_ensembl_blast_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Date:**  `6/14/16 UPDATED 7/5/16`

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `1_susscrofa_ensembl_blastn_output_e6.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** `5_susscrofa_filtered_ensembl_blastn_results.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the first of the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
This script filters and summarizes the blast results from the Sus scrofa Ensembl database query.

First, if there is only one hit for a sequence, return that hit.
If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
Then, retain matches with > 96% sequence identity.
Finally, the first filter returns sequences with the maximum bitscore.

The second filter returns one BLAST hit if the gene_biotype is the same between the redundant hits.
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
blastresults<-read.table("../1_susscrofa_ensembl_blastn_output_e6.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   4.207   0.045   4.280
```

```r
dim(blastresults)
```

```
## [1] 254756     12
```

```r
head(blastresults)
```

```
##                query_id             dbseq_id perc_identical length
## 1        seq_0_x9179927 ENSSSCT00000031380.1            100     22
## 2        seq_0_x9179927 ENSSSCT00000020599.1            100     22
## 3 seq_14912851_x5116126 ENSSSCT00000020827.2            100     22
## 4 seq_20028977_x4202989 ENSSSCT00000031380.1            100     23
## 5 seq_20028977_x4202989 ENSSSCT00000020599.1            100     23
## 6 seq_31132448_x3292947 ENSSSCT00000030019.1            100     22
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 1        0       0           1        22          53        74  4e-07
## 2        0       0           1        22          64        85  4e-07
## 3        0       0           1        22          50        71  4e-07
## 4        0       0           1        23          53        75  1e-07
## 5        0       0           1        23          64        86  1e-07
## 6        0       0           1        22          44        65  4e-07
##   bitscore
## 1     44.1
## 2     44.1
## 3     44.1
## 4     46.1
## 5     46.1
## 6     44.1
```

```r
str(blastresults)
```

```
## 'data.frame':	254756 obs. of  12 variables:
##  $ query_id      : Factor w/ 127885 levels "seq_0_x9179927",..: 1 1 159 1462 1462 127845 127845 127846 127846 127847 ...
##  $ dbseq_id      : Factor w/ 1863 levels "ENSSSCT00000001325.2",..: 1729 460 566 1729 460 1644 317 825 138 889 ...
##  $ perc_identical: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ length        : int  22 22 22 23 23 22 22 22 22 22 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  22 22 22 23 23 22 22 22 22 22 ...
##  $ dbseq_start   : int  53 64 50 53 64 44 49 8 11 18 ...
##  $ dbseq_end     : int  74 85 71 75 86 65 70 29 32 39 ...
##  $ evalue        : num  4e-07 4e-07 4e-07 1e-07 1e-07 4e-07 4e-07 4e-07 4e-07 4e-07 ...
##  $ bitscore      : num  44.1 44.1 44.1 46.1 46.1 44.1 44.1 44.1 44.1 44.1 ...
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
##                query_id             dbseq_id perc_identical length
## 1        seq_0_x9179927 ENSSSCT00000031380.1            100     22
## 2        seq_0_x9179927 ENSSSCT00000020599.1            100     22
## 3 seq_14912851_x5116126 ENSSSCT00000020827.2            100     22
## 4 seq_20028977_x4202989 ENSSSCT00000031380.1            100     23
## 5 seq_20028977_x4202989 ENSSSCT00000020599.1            100     23
## 6 seq_31132448_x3292947 ENSSSCT00000030019.1            100     22
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 1        0       0           1        22          53        74  4e-07
## 2        0       0           1        22          64        85  4e-07
## 3        0       0           1        22          50        71  4e-07
## 4        0       0           1        23          53        75  1e-07
## 5        0       0           1        23          64        86  1e-07
## 6        0       0           1        22          44        65  4e-07
##   bitscore              ensid ensid_version
## 1     44.1 ENSSSCT00000031380             1
## 2     44.1 ENSSSCT00000020599             1
## 3     44.1 ENSSSCT00000020827             2
## 4     46.1 ENSSSCT00000031380             1
## 5     46.1 ENSSSCT00000020599             1
## 6     44.1 ENSSSCT00000030019             1
```

```r
head(table(blastresults$ensid))
```

```
## 
## ENSSSCT00000001325 ENSSSCT00000001326 ENSSSCT00000001333 
##                 11                  9                  9 
## ENSSSCT00000001335 ENSSSCT00000001337 ENSSSCT00000001522 
##                  9                 52                  2
```

Check how many unique small RNA sequences had ensembl hits:


```r
length(unique(blastresults$query_id))
```

```
## [1] 127885
```

Check how many unique ensembl database sequences had hits:


```r
length(unique(blastresults$dbseq_id))
```

```
## [1] 1863
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
##   0.188   0.007   6.970
```

```r
dim(rst)
```

```
## [1] 1863   11
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
## 1858    ENSSSCT00000035867              17 ENSSSCG00000007146     
## 1859    ENSSSCT00000035935               X ENSSSCG00000024809     
## 1860    ENSSSCT00000036045               1 ENSSSCG00000004943     
## 1861    ENSSSCT00000036444              17 ENSSSCG00000007280     
## 1862    ENSSSCT00000036555               1 ENSSSCG00000004334     
## 1863    ENSSSCT00000036566              13 ENSSSCG00000031014     
##      rfam_transcript_name version transcript_version transcript_source
## 1858                            3                  1            havana
## 1859                            2                  1            havana
## 1860                            3                  1            havana
## 1861                            3                  1            havana
## 1862                            3                  1            havana
## 1863                            1                  1            havana
##      status transcript_status         gene_biotype
## 1858  KNOWN             KNOWN processed_transcript
## 1859  KNOWN             KNOWN processed_transcript
## 1860  KNOWN             KNOWN processed_transcript
## 1861  KNOWN             KNOWN processed_transcript
## 1862  KNOWN             KNOWN processed_transcript
## 1863  KNOWN             KNOWN              lincRNA
```

```r
table(rst$rfam)
```

```
## 
##         RF00001 RF00002 RF00003 RF00004 RF00006 RF00007 RF00009 RF00012 
##     511     131       2      49      29       1       1       4      14 
## RF00015 RF00016 RF00017 RF00019 RF00020 RF00024 RF00026 RF00030 RF00045 
##      10       9      16       9       7       1     565       1       2 
## RF00049 RF00054 RF00055 RF00056 RF00067 RF00068 RF00069 RF00070 RF00071 
##       2       1       1       3       2       1       1       2       2 
## RF00072 RF00085 RF00086 RF00087 RF00088 RF00089 RF00090 RF00091 RF00092 
##       2       1       1       1       1       1       3       4       4 
## RF00093 RF00096 RF00099 RF00100 RF00108 RF00133 RF00134 RF00136 RF00137 
##       3       4       7     100       4       4       1       1       2 
## RF00138 RF00139 RF00147 RF00150 RF00151 RF00152 RF00153 RF00154 RF00155 
##       3       3       2       2       3       1       2       1       1 
## RF00156 RF00157 RF00158 RF00181 RF00186 RF00187 RF00188 RF00190 RF00191 
##      13       2       1       8       1       1       3       1       1 
## RF00211 RF00213 RF00217 RF00218 RF00221 RF00231 RF00263 RF00264 RF00265 
##       3       3       2       1       1       1       1       1       1 
## RF00266 RF00271 RF00272 RF00273 RF00274 RF00275 RF00276 RF00277 RF00278 
##       3       2       3       2       2       2       1       1       1 
## RF00279 RF00281 RF00283 RF00284 RF00286 RF00288 RF00289 RF00302 RF00319 
##       3       1       1       1       1       1       1       2       2 
## RF00322 RF00324 RF00325 RF00334 RF00340 RF00341 RF00342 RF00377 RF00392 
##      20       1       3       4       2       1       1       2       4 
## RF00393 RF00394 RF00396 RF00397 RF00398 RF00399 RF00401 RF00402 RF00403 
##       1       1       1       2       2       1       2       2       1 
## RF00404 RF00405 RF00406 RF00407 RF00408 RF00409 RF00410 RF00411 RF00412 
##       1       1       4       2       1       1       4       2       7 
## RF00413 RF00414 RF00415 RF00416 RF00418 RF00419 RF00420 RF00421 RF00422 
##      13       3       2       1       5       1       3       1       1 
## RF00424 RF00425 RF00426 RF00428 RF00429 RF00430 RF00431 RF00432 RF00438 
##       1       6       2       1       3       2       1       2       2 
## RF00439 RF00440 RF00443 RF00476 RF00478 RF00494 RF00548 RF00553 RF00554 
##       1       1       2       2       2       1       1       1       4 
## RF00560 RF00561 RF00563 RF00564 RF00567 RF00568 RF00569 RF00570 RF00571 
##       2       5       1       1       1       1       1       1       1 
## RF00572 RF00573 RF00574 RF00575 RF00576 RF00577 RF00578 RF00579 RF00581 
##       1       1       1       3       1       1       1       2       3 
## RF00582 RF00584 RF00586 RF00588 RF00591 RF00592 RF00593 RF00598 RF00599 
##       1       1       1       2       2       1       1       1       1 
## RF00600 RF00601 RF00602 RF00603 RF00604 RF00606 RF00608 RF00609 RF00610 
##       1       1       1       1       4       1       2       1       2 
## RF00611 RF00613 RF00614 RF00618 RF00619 RF01156 RF01161 RF01164 RF01169 
##       2       1       6       1       7       2       1       1       3 
## RF01170 RF01173 RF01183 RF01186 RF01188 RF01191 RF01192 RF01200 RF01210 
##       1       1       1       1       1       2       2       1       8 
## RF01211 RF01229 RF01233 RF01241 RF01268 RF01277 RF01290 RF01291 RF01294 
##       2       2       1       1       2       1       2       1       1 
## RF01299 
##       2
```

```r
table(rst$gene_biotype)
```

```
## 
##              lincRNA                miRNA             misc_RNA 
##                    9                  431                  132 
##              Mt_rRNA              Mt_tRNA processed_transcript 
##                    2                   22                   47 
##                 rRNA               snoRNA                snRNA 
##                  133                  417                  670
```

```r
sum(rst$rfam == "")
```

```
## [1] 511
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
## [1] 254756     16
```

```r
head(blastresults)
```

```
##                query_id             dbseq_id perc_identical length
## 1        seq_0_x9179927 ENSSSCT00000031380.1            100     22
## 2        seq_0_x9179927 ENSSSCT00000020599.1            100     22
## 3 seq_14912851_x5116126 ENSSSCT00000020827.2            100     22
## 4 seq_20028977_x4202989 ENSSSCT00000031380.1            100     23
## 5 seq_20028977_x4202989 ENSSSCT00000020599.1            100     23
## 6 seq_31132448_x3292947 ENSSSCT00000030019.1            100     22
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 1        0       0           1        22          53        74  4e-07
## 2        0       0           1        22          64        85  4e-07
## 3        0       0           1        22          50        71  4e-07
## 4        0       0           1        23          53        75  1e-07
## 5        0       0           1        23          64        86  1e-07
## 6        0       0           1        22          44        65  4e-07
##   bitscore              ensid ensid_version gene_biotype status
## 1     44.1 ENSSSCT00000031380             1        miRNA  NOVEL
## 2     44.1 ENSSSCT00000020599             1        miRNA  KNOWN
## 3     44.1 ENSSSCT00000020827             2        miRNA  KNOWN
## 4     46.1 ENSSSCT00000031380             1        miRNA  NOVEL
## 5     46.1 ENSSSCT00000020599             1        miRNA  KNOWN
## 6     44.1 ENSSSCT00000030019             1        miRNA  KNOWN
```

```r
tail(blastresults)
```

```
##                query_id             dbseq_id perc_identical length
## 254751 seq_227425466_x2 ENSSSCT00000021153.1            100     22
## 254752 seq_227425488_x2 ENSSSCT00000020202.1            100     22
## 254753 seq_227425488_x2 ENSSSCT00000020550.1            100     22
## 254754 seq_227425572_x2 ENSSSCT00000019802.2            100     25
## 254755 seq_227425572_x2 ENSSSCT00000024765.1            100     24
## 254756 seq_227425572_x2 ENSSSCT00000021115.1            100     24
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254751        0       0           1        22          43        64  4e-07
## 254752        0       0           1        22           1        22  5e-07
## 254753        0       0           1        22           1        22  5e-07
## 254754        0       0           4        28           9        33  8e-09
## 254755        0       0           5        28          10        33  3e-08
## 254756        0       0           5        28          10        33  3e-08
##        bitscore              ensid ensid_version gene_biotype status
## 254751     44.1 ENSSSCT00000021153             1        miRNA  NOVEL
## 254752     44.1 ENSSSCT00000020202             1       snoRNA  NOVEL
## 254753     44.1 ENSSSCT00000020550             1       snoRNA  NOVEL
## 254754     50.1 ENSSSCT00000019802             2        miRNA  KNOWN
## 254755     48.1 ENSSSCT00000024765             1        miRNA  KNOWN
## 254756     48.1 ENSSSCT00000021115             1        miRNA  KNOWN
```

```r
table(blastresults$gene_biotype)
```

```
## 
##              lincRNA                miRNA             misc_RNA 
##                  615                59382                29385 
##              Mt_rRNA              Mt_tRNA processed_transcript 
##                23187                30282                 1302 
##                 rRNA               snoRNA                snRNA 
##                19835                22542                68226
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
## [1] 127885
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
## 14236.281    48.939  2071.037
```

```r
length(idx)
```

```
## [1] 127885
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
##  83.309  13.185  18.507
```

```r
length(stp1)
```

```
## [1] 127885
```

```r
head(stp1)
```

```
## [[1]]
##         query_id             dbseq_id perc_identical length mismatch
## 1 seq_0_x9179927 ENSSSCT00000031380.1            100     22        0
## 2 seq_0_x9179927 ENSSSCT00000020599.1            100     22        0
##   gapopen query_start query_end dbseq_start dbseq_end evalue bitscore
## 1       0           1        22          53        74  4e-07     44.1
## 2       0           1        22          64        85  4e-07     44.1
##                ensid ensid_version gene_biotype status
## 1 ENSSSCT00000031380             1        miRNA  NOVEL
## 2 ENSSSCT00000020599             1        miRNA  KNOWN
## 
## [[2]]
##                query_id             dbseq_id perc_identical length
## 3 seq_14912851_x5116126 ENSSSCT00000020827.2            100     22
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 3        0       0           1        22          50        71  4e-07
##   bitscore              ensid ensid_version gene_biotype status
## 3     44.1 ENSSSCT00000020827             2        miRNA  KNOWN
## 
## [[3]]
##                query_id             dbseq_id perc_identical length
## 4 seq_20028977_x4202989 ENSSSCT00000031380.1            100     23
## 5 seq_20028977_x4202989 ENSSSCT00000020599.1            100     23
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 4        0       0           1        23          53        75  1e-07
## 5        0       0           1        23          64        86  1e-07
##   bitscore              ensid ensid_version gene_biotype status
## 4     46.1 ENSSSCT00000031380             1        miRNA  NOVEL
## 5     46.1 ENSSSCT00000020599             1        miRNA  KNOWN
## 
## [[4]]
##                query_id             dbseq_id perc_identical length
## 6 seq_31132448_x3292947 ENSSSCT00000030019.1            100     22
## 7 seq_31132448_x3292947 ENSSSCT00000020296.2            100     22
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 6        0       0           1        22          44        65  4e-07
## 7        0       0           1        22          49        70  4e-07
##   bitscore              ensid ensid_version gene_biotype status
## 6     44.1 ENSSSCT00000030019             1        miRNA  KNOWN
## 7     44.1 ENSSSCT00000020296             2        miRNA  KNOWN
## 
## [[5]]
##                query_id             dbseq_id perc_identical length
## 8 seq_37418191_x2862885 ENSSSCT00000021393.2            100     22
## 9 seq_37418191_x2862885 ENSSSCT00000019914.2            100     22
##   mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 8        0       0           1        22           8        29  4e-07
## 9        0       0           1        22          11        32  4e-07
##   bitscore              ensid ensid_version gene_biotype status
## 8     44.1 ENSSSCT00000021393             2        miRNA  KNOWN
## 9     44.1 ENSSSCT00000019914             2        miRNA  KNOWN
## 
## [[6]]
##                 query_id             dbseq_id perc_identical length
## 10 seq_40281076_x2762402 ENSSSCT00000021537.1            100     22
##    mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 10        0       0           1        22          18        39  4e-07
##    bitscore              ensid ensid_version gene_biotype status
## 10     44.1 ENSSSCT00000021537             1        miRNA  NOVEL
```

```r
tail(stp1)
```

```
## [[1]]
##                query_id             dbseq_id perc_identical length
## 254569 seq_227425356_x2 ENSSSCT00000021744.1          96.15     26
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254569        1       0           1        26         128       153  5e-07
##        bitscore              ensid ensid_version gene_biotype status
## 254569     44.1 ENSSSCT00000021744             1         rRNA  KNOWN
## 
## [[2]]
##                query_id             dbseq_id perc_identical length
## 254570 seq_227425360_x2 ENSSSCT00000020085.1            100     24
## 254571 seq_227425360_x2 ENSSSCT00000020721.1            100     24
## 254572 seq_227425360_x2 ENSSSCT00000030290.1            100     24
## 254573 seq_227425360_x2 ENSSSCT00000020727.1            100     24
## 254574 seq_227425360_x2 ENSSSCT00000028304.1            100     24
## 254575 seq_227425360_x2 ENSSSCT00000026881.1            100     24
## 254576 seq_227425360_x2 ENSSSCT00000022019.1            100     24
## 254577 seq_227425360_x2 ENSSSCT00000023875.1            100     24
## 254578 seq_227425360_x2 ENSSSCT00000020140.1            100     24
## 254579 seq_227425360_x2 ENSSSCT00000029572.1            100     24
## 254580 seq_227425360_x2 ENSSSCT00000022989.1            100     24
## 254581 seq_227425360_x2 ENSSSCT00000024709.1            100     24
## 254582 seq_227425360_x2 ENSSSCT00000020928.1            100     24
## 254583 seq_227425360_x2 ENSSSCT00000027304.1            100     24
## 254584 seq_227425360_x2 ENSSSCT00000020634.1            100     24
## 254585 seq_227425360_x2 ENSSSCT00000032431.1            100     24
## 254586 seq_227425360_x2 ENSSSCT00000020027.1            100     24
## 254587 seq_227425360_x2 ENSSSCT00000032035.1            100     24
## 254588 seq_227425360_x2 ENSSSCT00000020569.1            100     24
## 254589 seq_227425360_x2 ENSSSCT00000029527.1            100     24
## 254590 seq_227425360_x2 ENSSSCT00000025500.1            100     24
## 254591 seq_227425360_x2 ENSSSCT00000020492.2            100     24
## 254592 seq_227425360_x2 ENSSSCT00000023907.1            100     24
## 254593 seq_227425360_x2 ENSSSCT00000022234.1            100     24
## 254594 seq_227425360_x2 ENSSSCT00000024431.1            100     24
## 254595 seq_227425360_x2 ENSSSCT00000020802.1            100     24
## 254596 seq_227425360_x2 ENSSSCT00000028162.1            100     24
## 254597 seq_227425360_x2 ENSSSCT00000021249.1            100     24
## 254598 seq_227425360_x2 ENSSSCT00000028426.1            100     24
## 254599 seq_227425360_x2 ENSSSCT00000020580.1            100     24
## 254600 seq_227425360_x2 ENSSSCT00000029976.1            100     24
## 254601 seq_227425360_x2 ENSSSCT00000030339.1            100     24
## 254602 seq_227425360_x2 ENSSSCT00000031547.1            100     24
## 254603 seq_227425360_x2 ENSSSCT00000024693.1            100     24
## 254604 seq_227425360_x2 ENSSSCT00000029517.1            100     24
## 254605 seq_227425360_x2 ENSSSCT00000020903.1            100     24
## 254606 seq_227425360_x2 ENSSSCT00000030525.1            100     24
## 254607 seq_227425360_x2 ENSSSCT00000020443.1            100     24
## 254608 seq_227425360_x2 ENSSSCT00000030465.1            100     24
## 254609 seq_227425360_x2 ENSSSCT00000032399.1            100     24
## 254610 seq_227425360_x2 ENSSSCT00000023871.1            100     24
## 254611 seq_227425360_x2 ENSSSCT00000030427.1            100     24
## 254612 seq_227425360_x2 ENSSSCT00000025238.1            100     24
## 254613 seq_227425360_x2 ENSSSCT00000020053.1            100     24
## 254614 seq_227425360_x2 ENSSSCT00000019782.1            100     24
## 254615 seq_227425360_x2 ENSSSCT00000029229.1            100     24
## 254616 seq_227425360_x2 ENSSSCT00000020788.1            100     24
## 254617 seq_227425360_x2 ENSSSCT00000023783.1            100     24
## 254618 seq_227425360_x2 ENSSSCT00000026194.1            100     24
## 254619 seq_227425360_x2 ENSSSCT00000023839.1            100     24
## 254620 seq_227425360_x2 ENSSSCT00000020564.1            100     24
## 254621 seq_227425360_x2 ENSSSCT00000025321.1            100     24
## 254622 seq_227425360_x2 ENSSSCT00000020332.1            100     24
## 254623 seq_227425360_x2 ENSSSCT00000024846.1            100     24
## 254624 seq_227425360_x2 ENSSSCT00000020779.2            100     24
## 254625 seq_227425360_x2 ENSSSCT00000031226.1            100     24
## 254626 seq_227425360_x2 ENSSSCT00000023636.1            100     24
## 254627 seq_227425360_x2 ENSSSCT00000028064.1            100     24
## 254628 seq_227425360_x2 ENSSSCT00000031808.1            100     24
## 254629 seq_227425360_x2 ENSSSCT00000019899.1            100     24
## 254630 seq_227425360_x2 ENSSSCT00000025902.1            100     24
## 254631 seq_227425360_x2 ENSSSCT00000032067.1            100     24
## 254632 seq_227425360_x2 ENSSSCT00000025313.1            100     24
## 254633 seq_227425360_x2 ENSSSCT00000022665.1            100     24
## 254634 seq_227425360_x2 ENSSSCT00000032438.1            100     24
## 254635 seq_227425360_x2 ENSSSCT00000027479.1            100     24
## 254636 seq_227425360_x2 ENSSSCT00000020915.1            100     24
## 254637 seq_227425360_x2 ENSSSCT00000032481.1            100     24
## 254638 seq_227425360_x2 ENSSSCT00000027633.1            100     24
## 254639 seq_227425360_x2 ENSSSCT00000019836.1            100     24
## 254640 seq_227425360_x2 ENSSSCT00000030530.1            100     24
## 254641 seq_227425360_x2 ENSSSCT00000028645.1            100     24
## 254642 seq_227425360_x2 ENSSSCT00000029138.1            100     24
## 254643 seq_227425360_x2 ENSSSCT00000024649.1            100     24
## 254644 seq_227425360_x2 ENSSSCT00000021567.1            100     24
## 254645 seq_227425360_x2 ENSSSCT00000019890.1            100     24
## 254646 seq_227425360_x2 ENSSSCT00000028493.1            100     24
## 254647 seq_227425360_x2 ENSSSCT00000020147.2            100     24
## 254648 seq_227425360_x2 ENSSSCT00000027754.1            100     24
## 254649 seq_227425360_x2 ENSSSCT00000030471.1            100     24
## 254650 seq_227425360_x2 ENSSSCT00000023693.1            100     24
## 254651 seq_227425360_x2 ENSSSCT00000031524.1            100     24
## 254652 seq_227425360_x2 ENSSSCT00000029851.1            100     24
## 254653 seq_227425360_x2 ENSSSCT00000032346.1            100     24
## 254654 seq_227425360_x2 ENSSSCT00000019982.1            100     24
## 254655 seq_227425360_x2 ENSSSCT00000021433.1            100     24
## 254656 seq_227425360_x2 ENSSSCT00000025838.1            100     24
## 254657 seq_227425360_x2 ENSSSCT00000020465.1            100     24
## 254658 seq_227425360_x2 ENSSSCT00000027505.1            100     24
## 254659 seq_227425360_x2 ENSSSCT00000020920.1            100     24
## 254660 seq_227425360_x2 ENSSSCT00000021281.1            100     24
## 254661 seq_227425360_x2 ENSSSCT00000020176.1            100     24
## 254662 seq_227425360_x2 ENSSSCT00000020720.2            100     24
## 254663 seq_227425360_x2 ENSSSCT00000025951.1            100     24
## 254664 seq_227425360_x2 ENSSSCT00000022999.1            100     24
## 254665 seq_227425360_x2 ENSSSCT00000028930.1            100     24
## 254666 seq_227425360_x2 ENSSSCT00000020650.1            100     24
## 254667 seq_227425360_x2 ENSSSCT00000032218.1            100     24
## 254668 seq_227425360_x2 ENSSSCT00000021093.1            100     24
## 254669 seq_227425360_x2 ENSSSCT00000028966.1            100     24
## 254670 seq_227425360_x2 ENSSSCT00000020292.1            100     24
## 254671 seq_227425360_x2 ENSSSCT00000032191.1            100     24
## 254672 seq_227425360_x2 ENSSSCT00000025974.1            100     24
## 254673 seq_227425360_x2 ENSSSCT00000032436.1            100     24
## 254674 seq_227425360_x2 ENSSSCT00000025798.1            100     24
## 254675 seq_227425360_x2 ENSSSCT00000031329.1            100     24
## 254676 seq_227425360_x2 ENSSSCT00000028404.1            100     24
## 254677 seq_227425360_x2 ENSSSCT00000020804.1            100     24
## 254678 seq_227425360_x2 ENSSSCT00000024611.1            100     24
## 254679 seq_227425360_x2 ENSSSCT00000020542.2            100     24
## 254680 seq_227425360_x2 ENSSSCT00000024518.1            100     24
## 254681 seq_227425360_x2 ENSSSCT00000022431.1            100     24
## 254682 seq_227425360_x2 ENSSSCT00000020100.2            100     24
## 254683 seq_227425360_x2 ENSSSCT00000021020.1            100     24
## 254684 seq_227425360_x2 ENSSSCT00000020397.1            100     24
## 254685 seq_227425360_x2 ENSSSCT00000023499.1            100     24
## 254686 seq_227425360_x2 ENSSSCT00000021613.1            100     24
## 254687 seq_227425360_x2 ENSSSCT00000025999.1            100     24
## 254688 seq_227425360_x2 ENSSSCT00000020536.1            100     24
## 254689 seq_227425360_x2 ENSSSCT00000020035.1            100     24
## 254690 seq_227425360_x2 ENSSSCT00000024813.1            100     24
## 254691 seq_227425360_x2 ENSSSCT00000027173.1            100     24
## 254692 seq_227425360_x2 ENSSSCT00000020767.1            100     24
## 254693 seq_227425360_x2 ENSSSCT00000020598.1            100     24
## 254694 seq_227425360_x2 ENSSSCT00000021297.2            100     24
## 254695 seq_227425360_x2 ENSSSCT00000019882.2            100     24
## 254696 seq_227425360_x2 ENSSSCT00000021259.1            100     24
## 254697 seq_227425360_x2 ENSSSCT00000026280.1            100     24
## 254698 seq_227425360_x2 ENSSSCT00000024484.1            100     24
## 254699 seq_227425360_x2 ENSSSCT00000021166.1            100     24
## 254700 seq_227425360_x2 ENSSSCT00000022664.1            100     24
## 254701 seq_227425360_x2 ENSSSCT00000028085.1            100     24
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254570        0       0           1        24          39        62  3e-08
## 254571        0       0           1        24          39        62  3e-08
## 254572        0       0           1        24          39        62  3e-08
## 254573        0       0           1        24          39        62  3e-08
## 254574        0       0           1        24          39        62  3e-08
## 254575        0       0           1        24          39        62  3e-08
## 254576        0       0           1        24          31        54  3e-08
## 254577        0       0           1        24          39        62  3e-08
## 254578        0       0           1        24          39        62  3e-08
## 254579        0       0           1        24          39        62  3e-08
## 254580        0       0           1        24          39        62  3e-08
## 254581        0       0           1        24          39        62  3e-08
## 254582        0       0           1        24          38        61  3e-08
## 254583        0       0           1        24          39        62  3e-08
## 254584        0       0           1        24          43        66  3e-08
## 254585        0       0           1        24          20        43  3e-08
## 254586        0       0           1        24          38        61  3e-08
## 254587        0       0           1        24          39        62  3e-08
## 254588        0       0           1        24          14        37  3e-08
## 254589        0       0           1        24          39        62  3e-08
## 254590        0       0           1        24          39        62  3e-08
## 254591        0       0           1        24          39        62  3e-08
## 254592        0       0           1        24          35        58  3e-08
## 254593        0       0           1        24          31        54  3e-08
## 254594        0       0           1        24          39        62  3e-08
## 254595        0       0           1        24          40        63  3e-08
## 254596        0       0           1        24          39        62  3e-08
## 254597        0       0           1        24          36        59  3e-08
## 254598        0       0           1        24          39        62  3e-08
## 254599        0       0           1        24          39        62  3e-08
## 254600        0       0           1        24          39        62  3e-08
## 254601        0       0           1        24          39        62  3e-08
## 254602        0       0           1        24          39        62  3e-08
## 254603        0       0           1        24          39        62  3e-08
## 254604        0       0           1        24          39        62  3e-08
## 254605        0       0           1        24          37        60  3e-08
## 254606        0       0           1        24          39        62  3e-08
## 254607        0       0           1        24          41        64  3e-08
## 254608        0       0           1        24          39        62  3e-08
## 254609        0       0           1        24          39        62  3e-08
## 254610        0       0           1        24          39        62  3e-08
## 254611        0       0           1        24          39        62  3e-08
## 254612        0       0           1        24          39        62  3e-08
## 254613        0       0           1        24          31        54  3e-08
## 254614        0       0           1        24          39        62  3e-08
## 254615        0       0           1        24          39        62  3e-08
## 254616        0       0           1        24          35        58  3e-08
## 254617        0       0           1        24          39        62  3e-08
## 254618        0       0           1        24          39        62  3e-08
## 254619        0       0           1        24          39        62  3e-08
## 254620        0       0           1        24          40        63  3e-08
## 254621        0       0           1        24          39        62  3e-08
## 254622        0       0           1        24          35        58  3e-08
## 254623        0       0           1        24          39        62  3e-08
## 254624        0       0           1        24          39        62  3e-08
## 254625        0       0           1        24          39        62  3e-08
## 254626        0       0           1        24          39        62  3e-08
## 254627        0       0           1        24          39        62  3e-08
## 254628        0       0           1        24          39        62  3e-08
## 254629        0       0           1        24          35        58  3e-08
## 254630        0       0           1        24          39        62  3e-08
## 254631        0       0           1        24          39        62  3e-08
## 254632        0       0           1        24          39        62  3e-08
## 254633        0       0           1        24          39        62  3e-08
## 254634        0       0           1        24          39        62  3e-08
## 254635        0       0           1        24          39        62  3e-08
## 254636        0       0           1        24          39        62  3e-08
## 254637        0       0           1        24          39        62  3e-08
## 254638        0       0           1        24          39        62  3e-08
## 254639        0       0           1        24          39        62  3e-08
## 254640        0       0           1        24          43        66  3e-08
## 254641        0       0           1        24          39        62  3e-08
## 254642        0       0           1        24          39        62  3e-08
## 254643        0       0           1        24          39        62  3e-08
## 254644        0       0           1        24          39        62  3e-08
## 254645        0       0           1        24          39        62  3e-08
## 254646        0       0           1        24          39        62  3e-08
## 254647        0       0           1        24          39        62  3e-08
## 254648        0       0           1        24          39        62  3e-08
## 254649        0       0           1        24          39        62  3e-08
## 254650        0       0           1        24          39        62  3e-08
## 254651        0       0           1        24          39        62  3e-08
## 254652        0       0           1        24          39        62  3e-08
## 254653        0       0           1        24          39        62  3e-08
## 254654        0       0           1        24          14        37  3e-08
## 254655        0       0           1        24          41        64  3e-08
## 254656        0       0           1        24          39        62  3e-08
## 254657        0       0           1        24          38        61  3e-08
## 254658        0       0           1        24          39        62  3e-08
## 254659        0       0           1        24          39        62  3e-08
## 254660        0       0           1        24          39        62  3e-08
## 254661        0       0           1        24          39        62  3e-08
## 254662        0       0           1        24          39        62  3e-08
## 254663        0       0           1        24          39        62  3e-08
## 254664        0       0           1        24          39        62  3e-08
## 254665        0       0           1        24          39        62  3e-08
## 254666        0       0           1        24          35        58  3e-08
## 254667        0       0           1        24          39        62  3e-08
## 254668        0       0           1        24          35        58  3e-08
## 254669        0       0           1        24          39        62  3e-08
## 254670        0       0           1        24          39        62  3e-08
## 254671        0       0           1        24          39        62  3e-08
## 254672        0       0           1        24          39        62  3e-08
## 254673        0       0           1        24          39        62  3e-08
## 254674        0       0           1        24          39        62  3e-08
## 254675        0       0           1        24          39        62  3e-08
## 254676        0       0           1        24          39        62  3e-08
## 254677        0       0           1        24          40        63  3e-08
## 254678        0       0           1        24          39        62  3e-08
## 254679        0       0           1        24          39        62  3e-08
## 254680        0       0           1        24          39        62  3e-08
## 254681        0       0           1        24          39        62  3e-08
## 254682        0       0           1        24          39        62  3e-08
## 254683        0       0           1        24          39        62  3e-08
## 254684        0       0           1        24          39        62  3e-08
## 254685        0       0           1        24          39        62  3e-08
## 254686        0       0           1        24          17        40  3e-08
## 254687        0       0           1        24          39        62  3e-08
## 254688        0       0           1        24          39        62  3e-08
## 254689        0       0           1        24          14        37  3e-08
## 254690        0       0           1        24          39        62  3e-08
## 254691        0       0           1        24          39        62  3e-08
## 254692        0       0           1        24          37        60  3e-08
## 254693        0       0           1        24          38        61  3e-08
## 254694        0       0           1        24          39        62  3e-08
## 254695        0       0           1        24          40        63  3e-08
## 254696        0       0           1        24          38        61  3e-08
## 254697        0       0           1        24          39        62  3e-08
## 254698        0       0           1        24          39        62  3e-08
## 254699        0       0           1        24          38        61  3e-08
## 254700        0       0           1        24          39        62  3e-08
## 254701        0       0           1        24          39        62  3e-08
##        bitscore              ensid ensid_version gene_biotype status
## 254570     48.1 ENSSSCT00000020085             1        snRNA  NOVEL
## 254571     48.1 ENSSSCT00000020721             1        snRNA  NOVEL
## 254572     48.1 ENSSSCT00000030290             1        snRNA  NOVEL
## 254573     48.1 ENSSSCT00000020727             1        snRNA  NOVEL
## 254574     48.1 ENSSSCT00000028304             1        snRNA  NOVEL
## 254575     48.1 ENSSSCT00000026881             1        snRNA  NOVEL
## 254576     48.1 ENSSSCT00000022019             1        snRNA  NOVEL
## 254577     48.1 ENSSSCT00000023875             1        snRNA  NOVEL
## 254578     48.1 ENSSSCT00000020140             1        snRNA  NOVEL
## 254579     48.1 ENSSSCT00000029572             1        snRNA  NOVEL
## 254580     48.1 ENSSSCT00000022989             1        snRNA  NOVEL
## 254581     48.1 ENSSSCT00000024709             1        snRNA  NOVEL
## 254582     48.1 ENSSSCT00000020928             1        snRNA  NOVEL
## 254583     48.1 ENSSSCT00000027304             1        snRNA  NOVEL
## 254584     48.1 ENSSSCT00000020634             1        snRNA  NOVEL
## 254585     48.1 ENSSSCT00000032431             1        snRNA  NOVEL
## 254586     48.1 ENSSSCT00000020027             1        snRNA  NOVEL
## 254587     48.1 ENSSSCT00000032035             1        snRNA  NOVEL
## 254588     48.1 ENSSSCT00000020569             1        snRNA  NOVEL
## 254589     48.1 ENSSSCT00000029527             1        snRNA  NOVEL
## 254590     48.1 ENSSSCT00000025500             1        snRNA  NOVEL
## 254591     48.1 ENSSSCT00000020492             2        snRNA  NOVEL
## 254592     48.1 ENSSSCT00000023907             1        snRNA  NOVEL
## 254593     48.1 ENSSSCT00000022234             1        snRNA  NOVEL
## 254594     48.1 ENSSSCT00000024431             1        snRNA  NOVEL
## 254595     48.1 ENSSSCT00000020802             1        snRNA  NOVEL
## 254596     48.1 ENSSSCT00000028162             1        snRNA  NOVEL
## 254597     48.1 ENSSSCT00000021249             1        snRNA  NOVEL
## 254598     48.1 ENSSSCT00000028426             1        snRNA  NOVEL
## 254599     48.1 ENSSSCT00000020580             1        snRNA  NOVEL
## 254600     48.1 ENSSSCT00000029976             1        snRNA  NOVEL
## 254601     48.1 ENSSSCT00000030339             1        snRNA  NOVEL
## 254602     48.1 ENSSSCT00000031547             1        snRNA  NOVEL
## 254603     48.1 ENSSSCT00000024693             1        snRNA  NOVEL
## 254604     48.1 ENSSSCT00000029517             1        snRNA  NOVEL
## 254605     48.1 ENSSSCT00000020903             1        snRNA  NOVEL
## 254606     48.1 ENSSSCT00000030525             1        snRNA  NOVEL
## 254607     48.1 ENSSSCT00000020443             1        snRNA  NOVEL
## 254608     48.1 ENSSSCT00000030465             1        snRNA  NOVEL
## 254609     48.1 ENSSSCT00000032399             1        snRNA  NOVEL
## 254610     48.1 ENSSSCT00000023871             1        snRNA  NOVEL
## 254611     48.1 ENSSSCT00000030427             1        snRNA  NOVEL
## 254612     48.1 ENSSSCT00000025238             1        snRNA  NOVEL
## 254613     48.1 ENSSSCT00000020053             1        snRNA  NOVEL
## 254614     48.1 ENSSSCT00000019782             1        snRNA  NOVEL
## 254615     48.1 ENSSSCT00000029229             1        snRNA  NOVEL
## 254616     48.1 ENSSSCT00000020788             1        snRNA  NOVEL
## 254617     48.1 ENSSSCT00000023783             1        snRNA  NOVEL
## 254618     48.1 ENSSSCT00000026194             1        snRNA  NOVEL
## 254619     48.1 ENSSSCT00000023839             1        snRNA  NOVEL
## 254620     48.1 ENSSSCT00000020564             1        snRNA  NOVEL
## 254621     48.1 ENSSSCT00000025321             1        snRNA  NOVEL
## 254622     48.1 ENSSSCT00000020332             1        snRNA  NOVEL
## 254623     48.1 ENSSSCT00000024846             1        snRNA  NOVEL
## 254624     48.1 ENSSSCT00000020779             2        snRNA  NOVEL
## 254625     48.1 ENSSSCT00000031226             1        snRNA  NOVEL
## 254626     48.1 ENSSSCT00000023636             1        snRNA  NOVEL
## 254627     48.1 ENSSSCT00000028064             1        snRNA  NOVEL
## 254628     48.1 ENSSSCT00000031808             1        snRNA  NOVEL
## 254629     48.1 ENSSSCT00000019899             1        snRNA  NOVEL
## 254630     48.1 ENSSSCT00000025902             1        snRNA  NOVEL
## 254631     48.1 ENSSSCT00000032067             1        snRNA  NOVEL
## 254632     48.1 ENSSSCT00000025313             1        snRNA  NOVEL
## 254633     48.1 ENSSSCT00000022665             1        snRNA  NOVEL
## 254634     48.1 ENSSSCT00000032438             1        snRNA  NOVEL
## 254635     48.1 ENSSSCT00000027479             1        snRNA  NOVEL
## 254636     48.1 ENSSSCT00000020915             1        snRNA  NOVEL
## 254637     48.1 ENSSSCT00000032481             1        snRNA  NOVEL
## 254638     48.1 ENSSSCT00000027633             1        snRNA  NOVEL
## 254639     48.1 ENSSSCT00000019836             1        snRNA  NOVEL
## 254640     48.1 ENSSSCT00000030530             1        snRNA  NOVEL
## 254641     48.1 ENSSSCT00000028645             1        snRNA  NOVEL
## 254642     48.1 ENSSSCT00000029138             1        snRNA  NOVEL
## 254643     48.1 ENSSSCT00000024649             1        snRNA  NOVEL
## 254644     48.1 ENSSSCT00000021567             1        snRNA  NOVEL
## 254645     48.1 ENSSSCT00000019890             1        snRNA  NOVEL
## 254646     48.1 ENSSSCT00000028493             1        snRNA  NOVEL
## 254647     48.1 ENSSSCT00000020147             2        snRNA  NOVEL
## 254648     48.1 ENSSSCT00000027754             1        snRNA  NOVEL
## 254649     48.1 ENSSSCT00000030471             1        snRNA  NOVEL
## 254650     48.1 ENSSSCT00000023693             1        snRNA  NOVEL
## 254651     48.1 ENSSSCT00000031524             1        snRNA  NOVEL
## 254652     48.1 ENSSSCT00000029851             1        snRNA  NOVEL
## 254653     48.1 ENSSSCT00000032346             1        snRNA  NOVEL
## 254654     48.1 ENSSSCT00000019982             1        snRNA  NOVEL
## 254655     48.1 ENSSSCT00000021433             1        snRNA  NOVEL
## 254656     48.1 ENSSSCT00000025838             1        snRNA  NOVEL
## 254657     48.1 ENSSSCT00000020465             1        snRNA  NOVEL
## 254658     48.1 ENSSSCT00000027505             1        snRNA  NOVEL
## 254659     48.1 ENSSSCT00000020920             1        snRNA  NOVEL
## 254660     48.1 ENSSSCT00000021281             1        snRNA  NOVEL
## 254661     48.1 ENSSSCT00000020176             1        snRNA  NOVEL
## 254662     48.1 ENSSSCT00000020720             2        snRNA  NOVEL
## 254663     48.1 ENSSSCT00000025951             1        snRNA  NOVEL
## 254664     48.1 ENSSSCT00000022999             1        snRNA  NOVEL
## 254665     48.1 ENSSSCT00000028930             1        snRNA  NOVEL
## 254666     48.1 ENSSSCT00000020650             1        snRNA  NOVEL
## 254667     48.1 ENSSSCT00000032218             1        snRNA  NOVEL
## 254668     48.1 ENSSSCT00000021093             1        snRNA  NOVEL
## 254669     48.1 ENSSSCT00000028966             1        snRNA  NOVEL
## 254670     48.1 ENSSSCT00000020292             1        snRNA  NOVEL
## 254671     48.1 ENSSSCT00000032191             1        snRNA  NOVEL
## 254672     48.1 ENSSSCT00000025974             1        snRNA  NOVEL
## 254673     48.1 ENSSSCT00000032436             1        snRNA  NOVEL
## 254674     48.1 ENSSSCT00000025798             1        snRNA  NOVEL
## 254675     48.1 ENSSSCT00000031329             1        snRNA  NOVEL
## 254676     48.1 ENSSSCT00000028404             1        snRNA  NOVEL
## 254677     48.1 ENSSSCT00000020804             1        snRNA  NOVEL
## 254678     48.1 ENSSSCT00000024611             1        snRNA  NOVEL
## 254679     48.1 ENSSSCT00000020542             2        snRNA  NOVEL
## 254680     48.1 ENSSSCT00000024518             1        snRNA  NOVEL
## 254681     48.1 ENSSSCT00000022431             1        snRNA  NOVEL
## 254682     48.1 ENSSSCT00000020100             2        snRNA  NOVEL
## 254683     48.1 ENSSSCT00000021020             1        snRNA  NOVEL
## 254684     48.1 ENSSSCT00000020397             1        snRNA  NOVEL
## 254685     48.1 ENSSSCT00000023499             1        snRNA  NOVEL
## 254686     48.1 ENSSSCT00000021613             1        snRNA  NOVEL
## 254687     48.1 ENSSSCT00000025999             1        snRNA  NOVEL
## 254688     48.1 ENSSSCT00000020536             1        snRNA  NOVEL
## 254689     48.1 ENSSSCT00000020035             1        snRNA  NOVEL
## 254690     48.1 ENSSSCT00000024813             1        snRNA  NOVEL
## 254691     48.1 ENSSSCT00000027173             1        snRNA  NOVEL
## 254692     48.1 ENSSSCT00000020767             1        snRNA  NOVEL
## 254693     48.1 ENSSSCT00000020598             1        snRNA  NOVEL
## 254694     48.1 ENSSSCT00000021297             2        snRNA  NOVEL
## 254695     48.1 ENSSSCT00000019882             2        snRNA  NOVEL
## 254696     48.1 ENSSSCT00000021259             1        snRNA  NOVEL
## 254697     48.1 ENSSSCT00000026280             1        snRNA  NOVEL
## 254698     48.1 ENSSSCT00000024484             1        snRNA  NOVEL
## 254699     48.1 ENSSSCT00000021166             1        snRNA  NOVEL
## 254700     48.1 ENSSSCT00000022664             1        snRNA  NOVEL
## 254701     48.1 ENSSSCT00000028085             1        snRNA  NOVEL
## 
## [[3]]
##                query_id             dbseq_id perc_identical length
## 254748 seq_227425458_x2 ENSSSCT00000030800.1            100     24
## 254749 seq_227425458_x2 ENSSSCT00000023836.1            100     24
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254748        0       0           1        24          52        75  3e-08
## 254749        0       0           1        24          52        75  3e-08
##        bitscore              ensid ensid_version gene_biotype status
## 254748     48.1 ENSSSCT00000030800             1        snRNA  NOVEL
## 254749     48.1 ENSSSCT00000023836             1        snRNA  NOVEL
## 
## [[4]]
##                query_id             dbseq_id perc_identical length
## 254750 seq_227425466_x2 ENSSSCT00000022687.1            100     22
## 254751 seq_227425466_x2 ENSSSCT00000021153.1            100     22
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254750        0       0           1        22          43        64  4e-07
## 254751        0       0           1        22          43        64  4e-07
##        bitscore              ensid ensid_version gene_biotype status
## 254750     44.1 ENSSSCT00000022687             1        miRNA  NOVEL
## 254751     44.1 ENSSSCT00000021153             1        miRNA  NOVEL
## 
## [[5]]
##                query_id             dbseq_id perc_identical length
## 254752 seq_227425488_x2 ENSSSCT00000020202.1            100     22
## 254753 seq_227425488_x2 ENSSSCT00000020550.1            100     22
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254752        0       0           1        22           1        22  5e-07
## 254753        0       0           1        22           1        22  5e-07
##        bitscore              ensid ensid_version gene_biotype status
## 254752     44.1 ENSSSCT00000020202             1       snoRNA  NOVEL
## 254753     44.1 ENSSSCT00000020550             1       snoRNA  NOVEL
## 
## [[6]]
##                query_id             dbseq_id perc_identical length
## 254754 seq_227425572_x2 ENSSSCT00000019802.2            100     25
##        mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 254754        0       0           4        28           9        33  8e-09
##        bitscore              ensid ensid_version gene_biotype status
## 254754     50.1 ENSSSCT00000019802             2        miRNA  KNOWN
```

How many sequences have more than one hit after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 26554
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##     2     3     4     5     6     7     8     9    10    11    12    13 
## 20556  2294   924   416   462   181   204    97    56    39    34    34 
##    14    15    16    17    18    19    20    21    22    23    24    25 
##    16   357    68   150    85    16    19    10    29    32    91    35 
##    26    27    28    29    30    31    32    33    34    35    36    37 
##    19     2     5     3     6     2    24    12    13    11     9     5 
##    38    39    40    41    42    43    44    45    46    47    48    49 
##     6     8     4     4     9     1    10     8     9     9     4     7 
##    50    51    52    53    54    55    56    57    58    59    60    61 
##     7     9     1    11     4     5     8     7    15    10    14     7 
##    62    63    64    65    66    67    68    79    88    93    96   102 
##     3    13     2    11    10     2     1     2     1     1     1     1 
##   111   116   118   123   127   129   132   136   137   145   147   154 
##     1     1     1     1     1     1     1     1     1     1     1     1 
##   155 
##     1
```

```r
filter2 <- filter <- function(seqblast){
# Sequence has only one ensembl blast hit
        if (nrow(seqblast) == 1){
                        return(seqblast)
        }
# Select sequence with greatest percent match to blast hit
        ident <- seqblast[seqblast$pident == max(seqblast$pident),]
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
##  43.847  11.951  12.593
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
sumblast <- do.call(rbind,stp3)
head(sumblast)
```

```
##                 query_id             dbseq_id perc_identical length
## 1         seq_0_x9179927 ENSSSCT00000031380.1            100     22
## 3  seq_14912851_x5116126 ENSSSCT00000020827.2            100     22
## 4  seq_20028977_x4202989 ENSSSCT00000031380.1            100     23
## 6  seq_31132448_x3292947 ENSSSCT00000030019.1            100     22
## 8  seq_37418191_x2862885 ENSSSCT00000021393.2            100     22
## 10 seq_40281076_x2762402 ENSSSCT00000021537.1            100     22
##    mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 1         0       0           1        22          53        74  4e-07
## 3         0       0           1        22          50        71  4e-07
## 4         0       0           1        23          53        75  1e-07
## 6         0       0           1        22          44        65  4e-07
## 8         0       0           1        22           8        29  4e-07
## 10        0       0           1        22          18        39  4e-07
##    bitscore              ensid ensid_version gene_biotype status
## 1      44.1 ENSSSCT00000031380             1        miRNA  NOVEL
## 3      44.1 ENSSSCT00000020827             2        miRNA  KNOWN
## 4      46.1 ENSSSCT00000031380             1        miRNA  NOVEL
## 6      44.1 ENSSSCT00000030019             1        miRNA  KNOWN
## 8      44.1 ENSSSCT00000021393             2        miRNA  KNOWN
## 10     44.1 ENSSSCT00000021537             1        miRNA  NOVEL
```

```r
dim(sumblast)
```

```
## [1] 127885     16
```

```r
length(unique(sumblast$query_id))
```

```
## [1] 127885
```

```r
table(sumblast$gene_biotype)
```

```
## 
##              lincRNA                miRNA             misc_RNA 
##                  228                38845                10248 
##              Mt_rRNA              Mt_tRNA processed_transcript 
##                23187                30282                  807 
##                 rRNA               snoRNA                snRNA 
##                 5230                15397                 3661
```

```r
uniqsumblast<-as.data.frame(table(sumblast$gene_biotype))
colnames(uniqsumblast)<-c("Gene_Biotype", "Freq")
uniqsumblast
```

```
##           Gene_Biotype  Freq
## 1              lincRNA   228
## 2                miRNA 38845
## 3             misc_RNA 10248
## 4              Mt_rRNA 23187
## 5              Mt_tRNA 30282
## 6 processed_transcript   807
## 7                 rRNA  5230
## 8               snoRNA 15397
## 9                snRNA  3661
```

Add the column of sequence count to the sumblast data frame


```r
sumblast$seq_count<-as.numeric(str_split_fixed(sumblast$query_id, "_x", 2)[,2])
head(sumblast)
```

```
##                 query_id             dbseq_id perc_identical length
## 1         seq_0_x9179927 ENSSSCT00000031380.1            100     22
## 3  seq_14912851_x5116126 ENSSSCT00000020827.2            100     22
## 4  seq_20028977_x4202989 ENSSSCT00000031380.1            100     23
## 6  seq_31132448_x3292947 ENSSSCT00000030019.1            100     22
## 8  seq_37418191_x2862885 ENSSSCT00000021393.2            100     22
## 10 seq_40281076_x2762402 ENSSSCT00000021537.1            100     22
##    mismatch gapopen query_start query_end dbseq_start dbseq_end evalue
## 1         0       0           1        22          53        74  4e-07
## 3         0       0           1        22          50        71  4e-07
## 4         0       0           1        23          53        75  1e-07
## 6         0       0           1        22          44        65  4e-07
## 8         0       0           1        22           8        29  4e-07
## 10        0       0           1        22          18        39  4e-07
##    bitscore              ensid ensid_version gene_biotype status seq_count
## 1      44.1 ENSSSCT00000031380             1        miRNA  NOVEL   9179927
## 3      44.1 ENSSSCT00000020827             2        miRNA  KNOWN   5116126
## 4      46.1 ENSSSCT00000031380             1        miRNA  NOVEL   4202989
## 6      44.1 ENSSSCT00000030019             1        miRNA  KNOWN   3292947
## 8      44.1 ENSSSCT00000021393             2        miRNA  KNOWN   2862885
## 10     44.1 ENSSSCT00000021537             1        miRNA  NOVEL   2762402
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype<-as.matrix(by(sumblast$seq_count, sumblast$gene_biotype, sum))
totalsumbiotype
```

```
##                          [,1]
## lincRNA                   652
## miRNA                88921566
## misc_RNA              1503043
## Mt_rRNA               2337654
## Mt_tRNA               2095395
## processed_transcript   231995
## rRNA                  3225800
## snoRNA                1156318
## snRNA                  188218
```

```r
if (sum(rownames(totalsumbiotype) != uniqsumblast$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the miRNAs and the processed_transcript fields:


```r
sum(sumblast$seq_count[sumblast$gene_biotype == "miRNA"])
```

```
## [1] 88921566
```

```r
sum(sumblast$seq_count[sumblast$gene_biotype == "processed_transcript"])
```

```
## [1] 231995
```

```r
if (sum(sumblast$seq_count[sumblast$gene_biotype == "miRNA"]) != totalsumbiotype["miRNA",]) stop ("miRNA counts not equal")
if (sum(sumblast$seq_count[sumblast$gene_biotype == "processed_transcript"]) != totalsumbiotype["processed_transcript",]) stop ("processed_transcript counts not equal")
```

## Save data


```r
save(sumblast, uniqsumblast, totalsumbiotype, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/5_susscrofa_filtered_ensembl_blastn_results.Rdata"))
write.csv(uniqsumblast, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/5_1_susscrofa_filtered_uniqseq_ensembl_blastn_results.csv", col.names=TRUE)
```

```
## Warning in write.csv(uniqsumblast, file = "/mnt/research/pigeqtl/analyses/
## microRNA/2_mirna_characterization_expression/0_rfam_database_query/
## 5_1_susscrofa_filtered_uniqseq_ensembl_blastn_results.csv", : attempt to
## set 'col.names' ignored
```

```r
write.csv(totalsumbiotype, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/5_2_susscrofa_filtered_totseq_ensembl_blastn_results.csv", row.names=TRUE)
```

