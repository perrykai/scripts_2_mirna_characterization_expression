**Script:** `15_summarize_blast_rfam_homo_sapiens_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`

**Date:**  7/7/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `3_homosapiens_rfam_blastn_output_e6.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** `7_homosapiens_filtered_rfam_blastn_results.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the third of the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
This script filters and summarizes the blast results from the Homo sapiens Rfam database query.

First, if there is only one hit for a sequence, return that hit.
If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
Then, retain matches with > 96% sequence identity.
Finally, the first filter returns sequences with the maximum bitscore.

The second filter returns one BLAST hit with the maximum percent identical between redundant hits.
Finally, if there are still remaining sequences with multiple hits, retain only the first hit. 

## Install libraries


```r
library(parallel)
library(stringr)
rm(list=ls())
```

## Load data


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts")

system.time({
hsa<-read.table("../3_homosapiens_rfam_blastn_output_e6.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   0.213   0.005   0.218
```

```r
dim(hsa)
```

```
## [1] 19293    12
```

```r
head(hsa)
```

```
##                query_id                               dbseq_id
## 1 seq_111010246_x310332 RF00005;tRNA;AL356957.27/185738-185668
## 2  seq_137541419_x92361 RF00005;tRNA;AL356957.27/185738-185668
## 3  seq_142952433_x74353 RF00005;tRNA;AL356957.27/185738-185668
## 4  seq_155498075_x38230 RF00005;tRNA;AL356957.27/185738-185668
## 5  seq_157679625_x34535 RF00005;tRNA;AL356957.27/178432-178361
## 6  seq_158827396_x32828 RF00005;tRNA;AL356957.27/185738-185668
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1         100.00     26        0       0           5        30           5
## 2         100.00     24        0       0           5        28           5
## 3         100.00     23        0       0           5        27           5
## 4         100.00     25        0       0           5        29           5
## 5         100.00     22        0       0           1        22          51
## 6          96.15     26        1       0           5        30           5
##   dbseq_end evalue bitscore
## 1        30  2e-09     52.0
## 2        28  2e-08     48.1
## 3        27  9e-08     46.1
## 4        29  6e-09     50.1
## 5        72  3e-07     44.1
## 6        30  4e-07     44.1
```

```r
str(hsa)
```

```
## 'data.frame':	19293 obs. of  12 variables:
##  $ query_id      : Factor w/ 13328 levels "seq_111010246_x310332",..: 1 2 3 4 5 6 7 8 9 10 ...
##  $ dbseq_id      : Factor w/ 475 levels "RF00001;5S_rRNA;AADB02003081.1/1831236-1831345",..: 157 157 157 157 156 157 446 14 446 156 ...
##  $ perc_identical: num  100 100 100 100 100 ...
##  $ length        : int  26 24 23 25 22 26 23 21 23 24 ...
##  $ mismatch      : int  0 0 0 0 0 1 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  5 5 5 5 1 5 4 1 3 1 ...
##  $ query_end     : int  30 28 27 29 22 30 26 21 25 24 ...
##  $ dbseq_start   : int  5 5 5 5 51 5 32 99 32 49 ...
##  $ dbseq_end     : int  30 28 27 29 72 30 10 119 10 72 ...
##  $ evalue        : num  2e-09 2e-08 9e-08 6e-09 3e-07 4e-07 8e-08 9e-07 8e-08 2e-08 ...
##  $ bitscore      : num  52 48.1 46.1 50.1 44.1 44.1 46.1 42.1 46.1 48.1 ...
```

## Analysis


```r
hsa$dbseq_id<-as.character(hsa$dbseq_id)
hsa$query_id<-as.character(hsa$query_id)
```

Summarize the Rfam database hits by separating the components of the "dbseq_id" column at the semicolons.


```r
hsa$rfam_accession <- as.character(lapply(strsplit(as.character(hsa$dbseq_id), ';', fixed = TRUE), "[", 1))
hsa$rfam_biotype <- as.character(lapply(strsplit(as.character(hsa$dbseq_id), ';', fixed = TRUE), "[", 2))
hsa$rfam_seqnamestartend <- as.character(lapply(strsplit(as.character(hsa$dbseq_id), ';', fixed = TRUE), "[", 3))
head(hsa)
```

```
##                query_id                               dbseq_id
## 1 seq_111010246_x310332 RF00005;tRNA;AL356957.27/185738-185668
## 2  seq_137541419_x92361 RF00005;tRNA;AL356957.27/185738-185668
## 3  seq_142952433_x74353 RF00005;tRNA;AL356957.27/185738-185668
## 4  seq_155498075_x38230 RF00005;tRNA;AL356957.27/185738-185668
## 5  seq_157679625_x34535 RF00005;tRNA;AL356957.27/178432-178361
## 6  seq_158827396_x32828 RF00005;tRNA;AL356957.27/185738-185668
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1         100.00     26        0       0           5        30           5
## 2         100.00     24        0       0           5        28           5
## 3         100.00     23        0       0           5        27           5
## 4         100.00     25        0       0           5        29           5
## 5         100.00     22        0       0           1        22          51
## 6          96.15     26        1       0           5        30           5
##   dbseq_end evalue bitscore rfam_accession rfam_biotype
## 1        30  2e-09     52.0        RF00005         tRNA
## 2        28  2e-08     48.1        RF00005         tRNA
## 3        27  9e-08     46.1        RF00005         tRNA
## 4        29  6e-09     50.1        RF00005         tRNA
## 5        72  3e-07     44.1        RF00005         tRNA
## 6        30  4e-07     44.1        RF00005         tRNA
##        rfam_seqnamestartend
## 1 AL356957.27/185738-185668
## 2 AL356957.27/185738-185668
## 3 AL356957.27/185738-185668
## 4 AL356957.27/185738-185668
## 5 AL356957.27/178432-178361
## 6 AL356957.27/185738-185668
```

The number of unique sequences with homo sapiens Rfam hits in this dataset


```r
length(unique(hsa$query_id))
```

```
## [1] 13328
```

Create a subset data frame containing the pertinent information for filtering


```r
rfamsubset<-hsa[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
```

```
## [1] 19293     6
```

```r
head(rfamsubset)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332         100.00  2e-09     52.0        RF00005
## 2  seq_137541419_x92361         100.00  2e-08     48.1        RF00005
## 3  seq_142952433_x74353         100.00  9e-08     46.1        RF00005
## 4  seq_155498075_x38230         100.00  6e-09     50.1        RF00005
## 5  seq_157679625_x34535         100.00  3e-07     44.1        RF00005
## 6  seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype
## 1         tRNA
## 2         tRNA
## 3         tRNA
## 4         tRNA
## 5         tRNA
## 6         tRNA
```

Create an index of the sequences to filter through


```r
system.time({
idx <- mclapply(unique(rfamsubset$query_id), function(x)
        rfamsubset[as.character(rfamsubset$query_id) == x,], mc.cores=10)
})
```

```
##    user  system elapsed 
##  48.383   0.267   8.170
```

```r
length(idx)
```

```
## [1] 13328
```

Function to filter multiple Rfam blast hits per sequence


```r
filter <- function(seqblast){
# Sequence has only one Rfam blast hit
        if (nrow(seqblast) == 1){
                return(seqblast)
        }
# Sequence with 100% match to Rfam blast hit
        ident <- seqblast[sapply(seqblast$perc_identical, function(x) x == 100),]
                if (nrow(ident) == 1){
                        return(ident)
        }
# Sequence with greater than 96% match to Rfam blast hit
        ident2 <- seqblast[sapply(seqblast$perc_identical, function(x) x > 96),]
                if (nrow(ident2) == 1){
                        return(ident2)
        }
# Select Rfam blast hit with the greatest bitscore for a sequence
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
stp1 <- mclapply(idx, filter, mc.cores=10)
})
```

```
##    user  system elapsed 
##   1.909   0.151   0.491
```

```r
length(stp1)
```

```
## [1] 13328
```

```r
head(stp1)
```

```
## [[1]]
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332            100  2e-09       52        RF00005
##   rfam_biotype
## 1         tRNA
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 2 seq_137541419_x92361            100  2e-08     48.1        RF00005
##   rfam_biotype
## 2         tRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 3 seq_142952433_x74353            100  9e-08     46.1        RF00005
##   rfam_biotype
## 3         tRNA
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 4 seq_155498075_x38230            100  6e-09     50.1        RF00005
##   rfam_biotype
## 4         tRNA
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 5 seq_157679625_x34535            100  3e-07     44.1        RF00005
##   rfam_biotype
## 5         tRNA
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 6 seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype
## 6         tRNA
```

```r
tail(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 19283 seq_227422450_x2          96.15  4e-07     44.1        RF00005
##       rfam_biotype
## 19283         tRNA
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 19284 seq_227422912_x2          93.33  4e-07     44.1        RF01960
## 19285 seq_227422912_x2          93.33  4e-07     44.1        RF01960
##           rfam_biotype
## 19284 SSU_rRNA_eukarya
## 19285 SSU_rRNA_eukarya
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 19286 seq_227423044_x2            100  9e-08     46.1        RF01960
##           rfam_biotype
## 19286 SSU_rRNA_eukarya
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 19287 seq_227423092_x2            100  9e-07     42.1        RF00017
## 19288 seq_227423092_x2            100  9e-07     42.1        RF00017
## 19289 seq_227423092_x2            100  9e-07     42.1        RF00017
## 19290 seq_227423092_x2            100  9e-07     42.1        RF00017
##       rfam_biotype
## 19287  Metazoa_SRP
## 19288  Metazoa_SRP
## 19289  Metazoa_SRP
## 19290  Metazoa_SRP
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 19291 seq_227424220_x2          96.15  4e-07     44.1        RF00005
##       rfam_biotype
## 19291         tRNA
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 19292 seq_227424364_x2            100  9e-07     42.1        RF01960
## 19293 seq_227424364_x2            100  9e-07     42.1        RF01960
##           rfam_biotype
## 19292 SSU_rRNA_eukarya
## 19293 SSU_rRNA_eukarya
```

How many sequences have more than one annotation after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 2339
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 2038   58   48   38   16   53   17   17   12    8    8    5    1    2    4 
##   17   19   32   33   35   36   38   60   67 
##    5    2    1    1    1    1    1    1    1
```

```r
filter2 <- filter <- function(seqblast){
# Sequence has only one Rfam blast hit
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
stp3 <- mclapply(stp1, filter2, mc.cores=10)
})
```

```
##    user  system elapsed 
##   0.902   0.123   0.286
```

How many sequences have more than one annotation after the second filtering


```r
stp4 <- stp3[unlist(lapply(stp3,nrow)) > 1]
length(stp4)
```

```
## [1] 0
```

Summary file of small RNA sequence ensembl blast results


```r
sumblast3 <- do.call(rbind,stp3)
head(sumblast3)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332         100.00  2e-09     52.0        RF00005
## 2  seq_137541419_x92361         100.00  2e-08     48.1        RF00005
## 3  seq_142952433_x74353         100.00  9e-08     46.1        RF00005
## 4  seq_155498075_x38230         100.00  6e-09     50.1        RF00005
## 5  seq_157679625_x34535         100.00  3e-07     44.1        RF00005
## 6  seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype
## 1         tRNA
## 2         tRNA
## 3         tRNA
## 4         tRNA
## 5         tRNA
## 6         tRNA
```

```r
dim(sumblast2)
```

```
## Error in eval(expr, envir, enclos): object 'sumblast2' not found
```

```r
length(unique(sumblast3$query_id))
```

```
## [1] 13328
```

```r
table(sumblast3$rfam_biotype)
```

```
## 
##        5_8S_rRNA          5S_rRNA              7SK            H19_3 
##              107              801               15                5 
##         Histone3          HIV_PBS       IRES_c-myc      Metazoa_SRP 
##                3              283                5              165 
##          mir-132          mir-142          mir-154          mir-551 
##                8               46               43                1 
##          mir-562          mir-628          mir-649          mir-689 
##                1                8                1               27 
##         SCARNA15          SNORA13           SNORA3          SNORA42 
##                2                1                1                1 
##         SNORD113           SNORD2          SNORD43           snoU13 
##              124                5                5                1 
##          snoZ278 SSU_rRNA_eukarya             tRNA               U1 
##                2             8711             2521               85 
##              U11               U2               U3               U4 
##               26               18               15               72 
##               U5               U6           U6atac               U7 
##                5               14                2               10 
##               U8            Y_RNA 
##                3              185
```

```r
uniqsumblast3<-as.data.frame(table(sumblast3$rfam_biotype))
colnames(uniqsumblast3)<-c("Gene_Biotype", "Freq")
uniqsumblast3
```

```
##        Gene_Biotype Freq
## 1         5_8S_rRNA  107
## 2           5S_rRNA  801
## 3               7SK   15
## 4             H19_3    5
## 5          Histone3    3
## 6           HIV_PBS  283
## 7        IRES_c-myc    5
## 8       Metazoa_SRP  165
## 9           mir-132    8
## 10          mir-142   46
## 11          mir-154   43
## 12          mir-551    1
## 13          mir-562    1
## 14          mir-628    8
## 15          mir-649    1
## 16          mir-689   27
## 17         SCARNA15    2
## 18          SNORA13    1
## 19           SNORA3    1
## 20          SNORA42    1
## 21         SNORD113  124
## 22           SNORD2    5
## 23          SNORD43    5
## 24           snoU13    1
## 25          snoZ278    2
## 26 SSU_rRNA_eukarya 8711
## 27             tRNA 2521
## 28               U1   85
## 29              U11   26
## 30               U2   18
## 31               U3   15
## 32               U4   72
## 33               U5    5
## 34               U6   14
## 35           U6atac    2
## 36               U7   10
## 37               U8    3
## 38            Y_RNA  185
```

Add the column of sequence count to the sumblast data frame


```r
sumblast3$seq_count<-as.numeric(str_split_fixed(sumblast3$query_id, "_x", 2)[,2])
head(sumblast3)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332         100.00  2e-09     52.0        RF00005
## 2  seq_137541419_x92361         100.00  2e-08     48.1        RF00005
## 3  seq_142952433_x74353         100.00  9e-08     46.1        RF00005
## 4  seq_155498075_x38230         100.00  6e-09     50.1        RF00005
## 5  seq_157679625_x34535         100.00  3e-07     44.1        RF00005
## 6  seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype seq_count
## 1         tRNA    310332
## 2         tRNA     92361
## 3         tRNA     74353
## 4         tRNA     38230
## 5         tRNA     34535
## 6         tRNA     32828
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype3<-as.matrix(by(sumblast3$seq_count, sumblast3$rfam_biotype, sum))
totalsumbiotype3
```

```
##                    [,1]
## 5_8S_rRNA          9560
## 5S_rRNA          126620
## 7SK                 112
## H19_3                13
## Histone3              9
## HIV_PBS          124513
## IRES_c-myc           10
## Metazoa_SRP        1926
## mir-132             357
## mir-142           33437
## mir-154            6745
## mir-551              38
## mir-562               2
## mir-628            1879
## mir-649               4
## mir-689             688
## SCARNA15              4
## SNORA13               7
## SNORA3                6
## SNORA42               3
## SNORD113          15242
## SNORD2               38
## SNORD43             138
## snoU13                3
## snoZ278               4
## SSU_rRNA_eukarya 384870
## tRNA             929468
## U1                 2452
## U11                  74
## U2                  103
## U3                  106
## U4                 8701
## U5                   16
## U6                  127
## U6atac               22
## U7                   40
## U8                  622
## Y_RNA              7414
```

```r
if (sum(rownames(totalsumbiotype3) != uniqsumblast3$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the 5s_rRNAs and the tRNA fields:


```r
sum(sumblast3$seq_count[sumblast3$rfam_biotype == "5S_rRNA"])
```

```
## [1] 126620
```

```r
sum(sumblast3$seq_count[sumblast3$rfam_biotype == "tRNA"])
```

```
## [1] 929468
```

```r
if (sum(sumblast3$seq_count[sumblast3$rfam_biotype == "5S_rRNA"]) != totalsumbiotype3["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast3$seq_count[sumblast3$rfam_biotype == "tRNA"]) != totalsumbiotype3["tRNA",]) stop ("tRNA counts not equal")
```

## Save data


```r
save(sumblast3, uniqsumblast3, totalsumbiotype3, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/7_homosapiens_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast3, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/7_1_homosapiens_filtered_uniqseq_rfam_blastn_results.csv", col.names=TRUE)
```

```
## Warning in write.csv(uniqsumblast3, file = "/mnt/research/pigeqtl/analyses/
## microRNA/2_mirna_characterization_expression/0_rfam_database_query/
## 7_1_homosapiens_filtered_uniqseq_rfam_blastn_results.csv", : attempt to set
## 'col.names' ignored
```

```r
write.csv(totalsumbiotype3, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/7_2_homosapiens_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)
```

