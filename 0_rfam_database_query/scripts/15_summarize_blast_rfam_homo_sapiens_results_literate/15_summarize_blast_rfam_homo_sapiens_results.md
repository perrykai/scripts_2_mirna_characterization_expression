**Script:** `15_summarize_blast_rfam_homo_sapiens_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`

**Date:**  7/13/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `3_homosapiens_rfam_blastn_output_e5.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** 

1. `3a_homosapiens_filtered_rfam_blastn_results.Rdata`
2. `3b_homosapiens_filtered_uniqseq_rfam_blastn_results.csv`
3. `3c_homosapiens_filtered_totseq_rfam_blastn_results.csv`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
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
hsa<-read.table("../3_homosapiens_rfam_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   0.257   0.007   0.328
```

```r
dim(hsa)
```

```
## [1] 27226    12
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
## 6  seq_157679625_x34535   RF00005;tRNA;AL356957.27/52701-52630
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     26        0       0           5        30           5
## 2            100     24        0       0           5        28           5
## 3            100     23        0       0           5        27           5
## 4            100     25        0       0           5        29           5
## 5            100     22        0       0           1        22          51
## 6            100     21        0       0           2        22          52
##   dbseq_end evalue bitscore
## 1        30  2e-09     52.0
## 2        28  2e-08     48.1
## 3        27  9e-08     46.1
## 4        29  6e-09     50.1
## 5        72  3e-07     44.1
## 6        72  1e-06     42.1
```

```r
str(hsa)
```

```
## 'data.frame':	27226 obs. of  12 variables:
##  $ query_id      : Factor w/ 19530 levels "seq_111010246_x310332",..: 1 2 3 4 5 5 6 7 8 9 ...
##  $ dbseq_id      : Factor w/ 236 levels "RF00001;5S_rRNA;AADB02002495.1/594380-594480",..: 104 104 104 104 103 105 104 218 13 218 ...
##  $ perc_identical: num  100 100 100 100 100 ...
##  $ length        : int  26 24 23 25 22 21 26 23 21 23 ...
##  $ mismatch      : int  0 0 0 0 0 0 1 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  5 5 5 5 1 2 5 4 1 3 ...
##  $ query_end     : int  30 28 27 29 22 22 30 26 21 25 ...
##  $ dbseq_start   : int  5 5 5 5 51 52 5 32 99 32 ...
##  $ dbseq_end     : int  30 28 27 29 72 72 30 10 119 10 ...
##  $ evalue        : num  2e-09 2e-08 9e-08 6e-09 3e-07 1e-06 4e-07 8e-08 9e-07 8e-08 ...
##  $ bitscore      : num  52 48.1 46.1 50.1 44.1 42.1 44.1 46.1 42.1 46.1 ...
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
## 6  seq_157679625_x34535   RF00005;tRNA;AL356957.27/52701-52630
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     26        0       0           5        30           5
## 2            100     24        0       0           5        28           5
## 3            100     23        0       0           5        27           5
## 4            100     25        0       0           5        29           5
## 5            100     22        0       0           1        22          51
## 6            100     21        0       0           2        22          52
##   dbseq_end evalue bitscore rfam_accession rfam_biotype
## 1        30  2e-09     52.0        RF00005         tRNA
## 2        28  2e-08     48.1        RF00005         tRNA
## 3        27  9e-08     46.1        RF00005         tRNA
## 4        29  6e-09     50.1        RF00005         tRNA
## 5        72  3e-07     44.1        RF00005         tRNA
## 6        72  1e-06     42.1        RF00005         tRNA
##        rfam_seqnamestartend
## 1 AL356957.27/185738-185668
## 2 AL356957.27/185738-185668
## 3 AL356957.27/185738-185668
## 4 AL356957.27/185738-185668
## 5 AL356957.27/178432-178361
## 6   AL356957.27/52701-52630
```

The number of unique sequences with homo sapiens Rfam hits in this dataset


```r
length(unique(hsa$query_id))
```

```
## [1] 19530
```

Create a subset data frame containing the pertinent information for filtering


```r
rfamsubset<-hsa[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
```

```
## [1] 27226     6
```

```r
head(rfamsubset)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332            100  2e-09     52.0        RF00005
## 2  seq_137541419_x92361            100  2e-08     48.1        RF00005
## 3  seq_142952433_x74353            100  9e-08     46.1        RF00005
## 4  seq_155498075_x38230            100  6e-09     50.1        RF00005
## 5  seq_157679625_x34535            100  3e-07     44.1        RF00005
## 6  seq_157679625_x34535            100  1e-06     42.1        RF00005
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
## 100.955   0.333  18.039
```

```r
length(idx)
```

```
## [1] 19530
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
##   2.544   0.162   0.732
```

```r
length(stp1)
```

```
## [1] 19530
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
## 7 seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype
## 7         tRNA
```

```r
tail(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27220 seq_227423044_x2            100  9e-08     46.1        RF01960
##           rfam_biotype
## 27220 SSU_rRNA_eukarya
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27221 seq_227423100_x2            100  4e-06     40.1        RF00005
##       rfam_biotype
## 27221         tRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27222 seq_227423942_x2             96  1e-06     42.1        RF01960
##           rfam_biotype
## 27222 SSU_rRNA_eukarya
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27223 seq_227424182_x2          95.83  6e-06     40.1        RF00005
##       rfam_biotype
## 27223         tRNA
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27224 seq_227424220_x2          96.15  4e-07     44.1        RF00005
##       rfam_biotype
## 27224         tRNA
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27225 seq_227424364_x2            100  9e-07     42.1        RF01960
## 27226 seq_227424364_x2            100  9e-07     42.1        RF01960
##           rfam_biotype
## 27225 SSU_rRNA_eukarya
## 27226 SSU_rRNA_eukarya
```

How many sequences have more than one annotation after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 3116
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##    2    3    4    5    6    7    8 
## 3056   27   14    4    2   10    3
```

```r
filter2 <- filter <- function(seqblast){
# Sequence has only one Rfam blast hit
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
stp3 <- mclapply(stp1, filter2, mc.cores=10)
})
```

```
##    user  system elapsed 
##   1.078   0.222   0.295
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
sumblast4 <- do.call(rbind,stp3)
head(sumblast4)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332         100.00  2e-09     52.0        RF00005
## 2  seq_137541419_x92361         100.00  2e-08     48.1        RF00005
## 3  seq_142952433_x74353         100.00  9e-08     46.1        RF00005
## 4  seq_155498075_x38230         100.00  6e-09     50.1        RF00005
## 5  seq_157679625_x34535         100.00  3e-07     44.1        RF00005
## 7  seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype
## 1         tRNA
## 2         tRNA
## 3         tRNA
## 4         tRNA
## 5         tRNA
## 7         tRNA
```

```r
dim(sumblast4)
```

```
## [1] 19530     6
```

```r
length(unique(sumblast4$query_id))
```

```
## [1] 19530
```

```r
table(sumblast4$rfam_biotype)
```

```
## 
##        5_8S_rRNA          5S_rRNA              7SK            H19_3 
##              104              850                2                8 
##         Histone3          HIV_PBS       IRES_c-myc        IRES_HIF1 
##                3              420                5                1 
##      Metazoa_SRP          mir-154          mir-562          mir-649 
##               43               62                1                2 
##          mir-689         SCARNA15          SNORA74         SNORD113 
##               65                1                1              135 
##          SNORD24          SNORD38          snoZ278 SSU_rRNA_eukarya 
##                2                1                2            13293 
##             tRNA               U1              U11               U2 
##             4270                2               47                2 
##               U3               U4               U7            Y_RNA 
##                1              117               16               74
```

```r
uniqsumblast4<-as.data.frame(table(sumblast4$rfam_biotype))
colnames(uniqsumblast4)<-c("Gene_Biotype", "Freq")
uniqsumblast4
```

```
##        Gene_Biotype  Freq
## 1         5_8S_rRNA   104
## 2           5S_rRNA   850
## 3               7SK     2
## 4             H19_3     8
## 5          Histone3     3
## 6           HIV_PBS   420
## 7        IRES_c-myc     5
## 8         IRES_HIF1     1
## 9       Metazoa_SRP    43
## 10          mir-154    62
## 11          mir-562     1
## 12          mir-649     2
## 13          mir-689    65
## 14         SCARNA15     1
## 15          SNORA74     1
## 16         SNORD113   135
## 17          SNORD24     2
## 18          SNORD38     1
## 19          snoZ278     2
## 20 SSU_rRNA_eukarya 13293
## 21             tRNA  4270
## 22               U1     2
## 23              U11    47
## 24               U2     2
## 25               U3     1
## 26               U4   117
## 27               U7    16
## 28            Y_RNA    74
```

Add the column of sequence count to the sumblast data frame


```r
sumblast4$seq_count<-as.numeric(str_split_fixed(sumblast4$query_id, "_x", 2)[,2])
head(sumblast4)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_111010246_x310332         100.00  2e-09     52.0        RF00005
## 2  seq_137541419_x92361         100.00  2e-08     48.1        RF00005
## 3  seq_142952433_x74353         100.00  9e-08     46.1        RF00005
## 4  seq_155498075_x38230         100.00  6e-09     50.1        RF00005
## 5  seq_157679625_x34535         100.00  3e-07     44.1        RF00005
## 7  seq_158827396_x32828          96.15  4e-07     44.1        RF00005
##   rfam_biotype seq_count
## 1         tRNA    310332
## 2         tRNA     92361
## 3         tRNA     74353
## 4         tRNA     38230
## 5         tRNA     34535
## 7         tRNA     32828
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype4<-as.matrix(by(sumblast4$seq_count, sumblast4$rfam_biotype, sum))
totalsumbiotype4
```

```
##                     [,1]
## 5_8S_rRNA           7036
## 5S_rRNA           155890
## 7SK                    4
## H19_3                 23
## Histone3               9
## HIV_PBS           127344
## IRES_c-myc            10
## IRES_HIF1              2
## Metazoa_SRP          426
## mir-154             6893
## mir-562                2
## mir-649                6
## mir-689             1694
## SCARNA15               2
## SNORA74                2
## SNORD113           12803
## SNORD24                4
## SNORD38                2
## snoZ278                4
## SSU_rRNA_eukarya  625254
## tRNA             1029531
## U1                     5
## U11                  175
## U2                     5
## U3                     2
## U4                  8258
## U7                    56
## Y_RNA                611
```

```r
if (sum(rownames(totalsumbiotype4) != uniqsumblast4$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the 5s_rRNAs and the tRNA fields:


```r
sum(sumblast4$seq_count[sumblast4$rfam_biotype == "5S_rRNA"])
```

```
## [1] 155890
```

```r
sum(sumblast4$seq_count[sumblast4$rfam_biotype == "tRNA"])
```

```
## [1] 1029531
```

```r
if (sum(sumblast4$seq_count[sumblast4$rfam_biotype == "5S_rRNA"]) != totalsumbiotype4["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast4$seq_count[sumblast4$rfam_biotype == "tRNA"]) != totalsumbiotype4["tRNA",]) stop ("tRNA counts not equal")
```

## Save data


```r
save(sumblast4, uniqsumblast4, totalsumbiotype4, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/3a_homosapiens_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast4, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/3b_homosapiens_filtered_uniqseq_rfam_blastn_results.csv")
write.csv(totalsumbiotype4, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/3c_homosapiens_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)
```

