**Script:** `16_summarize_blast_rfam_mus_musculus_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`

**Date:**  7/13/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `4_musmusculus_rfam_blastn_output_e5.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** 

1. `4a_musmusculus_filtered_rfam_blastn_results.Rdata`
2. `4b_musmusculus_filtered_uniqseq_rfam_blastn_results.csv`
3. `4c_musmusculus_filtered_totseq_rfam_blastn_results.csv`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the fourth (and final) of the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
This script filters and summarizes the blast results from the Mus musculus Rfam database query.

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
mmu<-read.table("../4_musmusculus_rfam_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   0.249   0.007   0.337
```

```r
dim(mmu)
```

```
## [1] 27267    12
```

```r
head(mmu)
```

```
##               query_id                                    dbseq_id
## 1 seq_137069853_x95209    RF00783;mir-484;AEKR01206416.1/2916-2850
## 2 seq_144363225_x67330    RF00783;mir-484;AEKR01206416.1/2916-2850
## 3 seq_157575655_x34724    RF00783;mir-484;AEKR01206416.1/2916-2850
## 4 seq_162056863_x28478    RF00783;mir-484;AEKR01206416.1/2916-2850
## 5 seq_167564473_x21752 RF01942;mir-1937;CAAA01156883.1/50076-50196
## 6 seq_167564473_x21752     RF01942;mir-1937;AAHY01278511.1/601-482
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     22        0       0           1        22           5
## 2            100     22        0       0           1        22           5
## 3            100     21        0       0           1        21           5
## 4            100     23        0       0           1        23           5
## 5            100     20        0       0           3        22          27
## 6            100     20        0       0           3        22          19
##   dbseq_end evalue bitscore
## 1        26  4e-07     44.1
## 2        26  4e-07     44.1
## 3        25  1e-06     42.1
## 4        27  1e-07     46.1
## 5        46  6e-06     40.1
## 6        38  6e-06     40.1
```

```r
str(mmu)
```

```
## 'data.frame':	27267 obs. of  12 variables:
##  $ query_id      : Factor w/ 17078 levels "seq_137069853_x95209",..: 1 2 3 4 5 5 5 5 6 7 ...
##  $ dbseq_id      : Factor w/ 496 levels "RF00001;5S_rRNA;AAHY01011173.1/676-776",..: 397 397 397 397 485 422 413 435 495 495 ...
##  $ perc_identical: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ length        : int  22 22 21 23 20 20 20 20 24 21 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 3 3 3 3 1 1 ...
##  $ query_end     : int  22 22 21 23 22 22 22 22 24 21 ...
##  $ dbseq_start   : int  5 5 5 5 27 19 27 27 1 1 ...
##  $ dbseq_end     : int  26 26 25 27 46 38 46 46 24 21 ...
##  $ evalue        : num  4e-07 4e-07 1e-06 1e-07 6e-06 6e-06 6e-06 6e-06 3e-08 1e-06 ...
##  $ bitscore      : num  44.1 44.1 42.1 46.1 40.1 40.1 40.1 40.1 48.1 42.1 ...
```

## Analysis


```r
mmu$dbseq_id<-as.character(mmu$dbseq_id)
mmu$query_id<-as.character(mmu$query_id)
```

Summarize the Rfam database hits by separating the components of the "dbseq_id" column at the semicolons.


```r
mmu$rfam_accession <- as.character(lapply(strsplit(as.character(mmu$dbseq_id), ';', fixed = TRUE), "[", 1))
mmu$rfam_biotype <- as.character(lapply(strsplit(as.character(mmu$dbseq_id), ';', fixed = TRUE), "[", 2))
mmu$rfam_seqnamestartend <- as.character(lapply(strsplit(as.character(mmu$dbseq_id), ';', fixed = TRUE), "[", 3))
head(mmu)
```

```
##               query_id                                    dbseq_id
## 1 seq_137069853_x95209    RF00783;mir-484;AEKR01206416.1/2916-2850
## 2 seq_144363225_x67330    RF00783;mir-484;AEKR01206416.1/2916-2850
## 3 seq_157575655_x34724    RF00783;mir-484;AEKR01206416.1/2916-2850
## 4 seq_162056863_x28478    RF00783;mir-484;AEKR01206416.1/2916-2850
## 5 seq_167564473_x21752 RF01942;mir-1937;CAAA01156883.1/50076-50196
## 6 seq_167564473_x21752     RF01942;mir-1937;AAHY01278511.1/601-482
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     22        0       0           1        22           5
## 2            100     22        0       0           1        22           5
## 3            100     21        0       0           1        21           5
## 4            100     23        0       0           1        23           5
## 5            100     20        0       0           3        22          27
## 6            100     20        0       0           3        22          19
##   dbseq_end evalue bitscore rfam_accession rfam_biotype
## 1        26  4e-07     44.1        RF00783      mir-484
## 2        26  4e-07     44.1        RF00783      mir-484
## 3        25  1e-06     42.1        RF00783      mir-484
## 4        27  1e-07     46.1        RF00783      mir-484
## 5        46  6e-06     40.1        RF01942     mir-1937
## 6        38  6e-06     40.1        RF01942     mir-1937
##         rfam_seqnamestartend
## 1   AEKR01206416.1/2916-2850
## 2   AEKR01206416.1/2916-2850
## 3   AEKR01206416.1/2916-2850
## 4   AEKR01206416.1/2916-2850
## 5 CAAA01156883.1/50076-50196
## 6     AAHY01278511.1/601-482
```

The number of unique sequences with mus musculus Rfam hits in this dataset


```r
length(unique(mmu$query_id))
```

```
## [1] 17078
```

Create a subset data frame containing the pertinent information for filtering


```r
rfamsubset<-mmu[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
```

```
## [1] 27267     6
```

```r
head(rfamsubset)
```

```
##               query_id perc_identical evalue bitscore rfam_accession
## 1 seq_137069853_x95209            100  4e-07     44.1        RF00783
## 2 seq_144363225_x67330            100  4e-07     44.1        RF00783
## 3 seq_157575655_x34724            100  1e-06     42.1        RF00783
## 4 seq_162056863_x28478            100  1e-07     46.1        RF00783
## 5 seq_167564473_x21752            100  6e-06     40.1        RF01942
## 6 seq_167564473_x21752            100  6e-06     40.1        RF01942
##   rfam_biotype
## 1      mir-484
## 2      mir-484
## 3      mir-484
## 4      mir-484
## 5     mir-1937
## 6     mir-1937
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
##  97.354   0.382  14.276
```

```r
length(idx)
```

```
## [1] 17078
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
##   3.617   0.176   0.733
```

```r
length(stp1)
```

```
## [1] 17078
```

```r
head(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 1 seq_137069853_x95209            100  4e-07     44.1        RF00783
##   rfam_biotype
## 1      mir-484
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 2 seq_144363225_x67330            100  4e-07     44.1        RF00783
##   rfam_biotype
## 2      mir-484
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 3 seq_157575655_x34724            100  1e-06     42.1        RF00783
##   rfam_biotype
## 3      mir-484
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 4 seq_162056863_x28478            100  1e-07     46.1        RF00783
##   rfam_biotype
## 4      mir-484
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 5 seq_167564473_x21752            100  6e-06     40.1        RF01942
## 6 seq_167564473_x21752            100  6e-06     40.1        RF01942
## 7 seq_167564473_x21752            100  6e-06     40.1        RF01942
## 8 seq_167564473_x21752            100  6e-06     40.1        RF01942
##   rfam_biotype
## 5     mir-1937
## 6     mir-1937
## 7     mir-1937
## 8     mir-1937
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 9 seq_169281334_x19740            100  3e-08     48.1        RF01960
##       rfam_biotype
## 9 SSU_rRNA_eukarya
```

```r
tail(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27258 seq_227423612_x2            100  6e-06     40.1        RF01942
## 27259 seq_227423612_x2            100  6e-06     40.1        RF01942
## 27260 seq_227423612_x2            100  6e-06     40.1        RF01942
## 27261 seq_227423612_x2            100  6e-06     40.1        RF01942
##       rfam_biotype
## 27258     mir-1937
## 27259     mir-1937
## 27260     mir-1937
## 27261     mir-1937
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27262 seq_227423696_x2            100  6e-06     40.1        RF01960
##           rfam_biotype
## 27262 SSU_rRNA_eukarya
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27263 seq_227423936_x2             96  2e-06     42.1        RF01960
##           rfam_biotype
## 27263 SSU_rRNA_eukarya
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27264 seq_227424478_x2            100  7e-06     40.1        RF00783
##       rfam_biotype
## 27264      mir-484
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27265 seq_227424846_x2             96  2e-06     42.1        RF01960
## 27266 seq_227424846_x2             96  2e-06     42.1        RF01960
##           rfam_biotype
## 27265 SSU_rRNA_eukarya
## 27266 SSU_rRNA_eukarya
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27267 seq_227425450_x2            100  2e-06     42.1        RF01960
##           rfam_biotype
## 27267 SSU_rRNA_eukarya
```

How many sequences have more than one annotation after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 5327
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##    2    3    4    5    6    7    8    9   10   12   13   14   15   23 
## 4911  244   96   19    6    5   16    5    3    3    9    5    1    4
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
##   1.411   0.173   0.379
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
sumblast5 <- do.call(rbind,stp3)
head(sumblast5)
```

```
##               query_id perc_identical evalue bitscore rfam_accession
## 1 seq_137069853_x95209            100  4e-07     44.1        RF00783
## 2 seq_144363225_x67330            100  4e-07     44.1        RF00783
## 3 seq_157575655_x34724            100  1e-06     42.1        RF00783
## 4 seq_162056863_x28478            100  1e-07     46.1        RF00783
## 5 seq_167564473_x21752            100  6e-06     40.1        RF01942
## 9 seq_169281334_x19740            100  3e-08     48.1        RF01960
##       rfam_biotype
## 1          mir-484
## 2          mir-484
## 3          mir-484
## 4          mir-484
## 5         mir-1937
## 9 SSU_rRNA_eukarya
```

```r
dim(sumblast5)
```

```
## [1] 17078     6
```

```r
length(unique(sumblast5$query_id))
```

```
## [1] 17078
```

```r
table(sumblast5$rfam_biotype)
```

```
## 
##        5_8S_rRNA          5S_rRNA              7SK         Histone3 
##               51               99                9               37 
##       IRES_L-myc      Metazoa_SRP          mir-101          mir-192 
##               11                9               13                1 
##          mir-193         mir-1937          mir-299          mir-484 
##                2              393                2              554 
##          mir-598          mir-628          mir-689          SECIS_1 
##               12                1                4                1 
##          SNORA21          SNORA28          SNORA73          SNORA74 
##                1               19                9                1 
##          SNORD58          SNORD91           snoU85 SSU_rRNA_eukarya 
##                2                4               65            13397 
##             tRNA               U1               U2               U6 
##             2363                8                7                1 
##             XIST            Y_RNA 
##                1                1
```

```r
uniqsumblast5<-as.data.frame(table(sumblast5$rfam_biotype))
colnames(uniqsumblast5)<-c("Gene_Biotype", "Freq")
uniqsumblast5
```

```
##        Gene_Biotype  Freq
## 1         5_8S_rRNA    51
## 2           5S_rRNA    99
## 3               7SK     9
## 4          Histone3    37
## 5        IRES_L-myc    11
## 6       Metazoa_SRP     9
## 7           mir-101    13
## 8           mir-192     1
## 9           mir-193     2
## 10         mir-1937   393
## 11          mir-299     2
## 12          mir-484   554
## 13          mir-598    12
## 14          mir-628     1
## 15          mir-689     4
## 16          SECIS_1     1
## 17          SNORA21     1
## 18          SNORA28    19
## 19          SNORA73     9
## 20          SNORA74     1
## 21          SNORD58     2
## 22          SNORD91     4
## 23           snoU85    65
## 24 SSU_rRNA_eukarya 13397
## 25             tRNA  2363
## 26               U1     8
## 27               U2     7
## 28               U6     1
## 29             XIST     1
## 30            Y_RNA     1
```

Add the column of sequence count to the sumblast data frame


```r
sumblast5$seq_count<-as.numeric(str_split_fixed(sumblast5$query_id, "_x", 2)[,2])
head(sumblast5)
```

```
##               query_id perc_identical evalue bitscore rfam_accession
## 1 seq_137069853_x95209            100  4e-07     44.1        RF00783
## 2 seq_144363225_x67330            100  4e-07     44.1        RF00783
## 3 seq_157575655_x34724            100  1e-06     42.1        RF00783
## 4 seq_162056863_x28478            100  1e-07     46.1        RF00783
## 5 seq_167564473_x21752            100  6e-06     40.1        RF01942
## 9 seq_169281334_x19740            100  3e-08     48.1        RF01960
##       rfam_biotype seq_count
## 1          mir-484     95209
## 2          mir-484     67330
## 3          mir-484     34724
## 4          mir-484     28478
## 5         mir-1937     21752
## 9 SSU_rRNA_eukarya     19740
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype5<-as.matrix(by(sumblast5$seq_count, sumblast5$rfam_biotype, sum))
totalsumbiotype5
```

```
##                    [,1]
## 5_8S_rRNA          3066
## 5S_rRNA            2387
## 7SK                  30
## Histone3            322
## IRES_L-myc           63
## Metazoa_SRP          54
## mir-101             341
## mir-192               2
## mir-193               4
## mir-1937          51215
## mir-299               7
## mir-484          291492
## mir-598              29
## mir-628               2
## mir-689              50
## SECIS_1               3
## SNORA21               2
## SNORA28             253
## SNORA73              39
## SNORA74               2
## SNORD58               9
## SNORD91              22
## snoU85              789
## SSU_rRNA_eukarya 618355
## tRNA              95734
## U1                   21
## U2                   29
## U6                    2
## XIST                  2
## Y_RNA                 2
```

```r
if (sum(rownames(totalsumbiotype5) != uniqsumblast5$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the 5s_rRNAs and the tRNA fields:


```r
sum(sumblast5$seq_count[sumblast5$rfam_biotype == "5S_rRNA"])
```

```
## [1] 2387
```

```r
sum(sumblast5$seq_count[sumblast5$rfam_biotype == "tRNA"])
```

```
## [1] 95734
```

```r
if (sum(sumblast5$seq_count[sumblast5$rfam_biotype == "5S_rRNA"]) != totalsumbiotype5["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast5$seq_count[sumblast5$rfam_biotype == "tRNA"]) != totalsumbiotype5["tRNA",]) stop ("tRNA counts not equal")
```

## Save data


```r
save(sumblast5, uniqsumblast5, totalsumbiotype5, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/4a_musmusculus_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast5, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/4b_musmusculus_filtered_uniqseq_rfam_blastn_results.csv")
write.csv(totalsumbiotype5, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/4c_musmusculus_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)
```

