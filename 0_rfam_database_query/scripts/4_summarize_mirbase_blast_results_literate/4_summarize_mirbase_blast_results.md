**Script:** `4_summarize_mirbase_blast_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/`

**Date:**  `7/11/16`

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `0_susscrofa_mirbase_blastn_output_e5.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** 

1. `0a_susscrofa_filtered_mirbase_blastn_results.txt`
2. `0b_susscrofa_filtered_uniqseq_miRBase_blastn_results.csv`
3. `0c_susscrofa_filtered_totseq_miRBase_blastn_results.csv`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the first of the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
This script filters and summarizes the blast results from the Sus scrofa miRBase database query.

First, if there is only one hit for a sequence, return that hit.
If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
Then, retain matches with > 96% sequence identity.
Finally, the first filter returns sequences with the maximum bitscore.

The second filter returns one BLAST hit if the sequence with greatest percent match to blast hit.
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
ssc<-read.table("../0_susscrofa_mirbase_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   2.403   0.031   2.502
```

```r
dim(ssc)
```

```
## [1] 178255     12
```

```r
head(ssc)
```

```
##                query_id    dbseq_id perc_identical length mismatch gapopen
## 1        seq_0_x9179927   ssc-miR-1            100     21        0       0
## 2  seq_9179927_x5732924   ssc-miR-1            100     21        0       0
## 3 seq_14912851_x5116126 ssc-miR-206            100     22        0       0
## 4 seq_20028977_x4202989   ssc-miR-1            100     21        0       0
## 5 seq_24231966_x3505397   ssc-miR-1            100     20        0       0
## 6 seq_27737363_x3395085   ssc-miR-1            100     21        0       0
##   query_start query_end dbseq_start dbseq_end evalue bitscore
## 1           1        21           1        21  2e-08     42.1
## 2           1        21           1        21  2e-08     42.1
## 3           1        22           1        22  4e-09     44.1
## 4           1        21           1        21  2e-08     42.1
## 5           1        20           1        20  6e-08     40.1
## 6           1        21           1        21  2e-08     42.1
```

```r
str(ssc)
```

```
## 'data.frame':	178255 obs. of  12 variables:
##  $ query_id      : Factor w/ 147657 levels "seq_0_x9179927",..: 1 147643 276 2901 147592 147593 147594 147594 147595 147596 ...
##  $ dbseq_id      : Factor w/ 328 levels "ssc-let-7a","ssc-let-7c",..: 9 9 118 9 9 9 212 213 9 1 ...
##  $ perc_identical: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ length        : int  21 21 22 21 20 21 22 20 21 22 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  21 21 22 21 20 21 22 20 21 22 ...
##  $ dbseq_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ dbseq_end     : int  21 21 22 21 20 21 22 20 21 22 ...
##  $ evalue        : num  2e-08 2e-08 4e-09 2e-08 6e-08 2e-08 4e-09 7e-08 2e-08 4e-09 ...
##  $ bitscore      : num  42.1 42.1 44.1 42.1 40.1 42.1 44.1 40.1 42.1 44.1 ...
```

## Analysis


```r
ssc$dbseq_id<-as.character(ssc$dbseq_id)
ssc$query_id<-as.character(ssc$query_id)
```

Create a subset data frame containing the pertinent information for filtering


```r
mirbasesubset<-ssc[,c("query_id", "dbseq_id", "perc_identical", "evalue", "bitscore")]
dim(mirbasesubset)
```

```
## [1] 178255      5
```

```r
head(mirbasesubset)
```

```
##                query_id    dbseq_id perc_identical evalue bitscore
## 1        seq_0_x9179927   ssc-miR-1            100  2e-08     42.1
## 2  seq_9179927_x5732924   ssc-miR-1            100  2e-08     42.1
## 3 seq_14912851_x5116126 ssc-miR-206            100  4e-09     44.1
## 4 seq_20028977_x4202989   ssc-miR-1            100  2e-08     42.1
## 5 seq_24231966_x3505397   ssc-miR-1            100  6e-08     40.1
## 6 seq_27737363_x3395085   ssc-miR-1            100  2e-08     42.1
```

The number of unique sequences with susscrofa miRBase hits in this dataset


```r
length(unique(mirbasesubset$query_id))
```

```
## [1] 147657
```

The number of unique miRBase miRNAs identified by this BLAST query


```r
length(unique(mirbasesubset$dbseq_id))
```

```
## [1] 328
```

```r
######################################################################
```

Obtained the following code from Deborah Velez, how she filtered her BLAST results for gene annotation (/mnt/research/ernstc_lab/fat_eqtl/BLAST/code/annotation_pigoligoarray.R)


```r
######################################################################
```

Format the miRBase blast output data.frame to a list to filter the sequences with multiple miRBase hits


```r
# Number of unique sequences
n <- length(unique(mirbasesubset$query_id))
n
```

```
## [1] 147657
```

```r
# Create the sequence list to filter
system.time({
idx1 <- mclapply(as.character(unique(mirbasesubset$query_id))[1:round(n/5)], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx2 <- mclapply(as.character(unique(mirbasesubset$query_id))[(round(n/5) + 1):(2*round(n/5))], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx3 <- mclapply(as.character(unique(mirbasesubset$query_id))[((2*round(n/5)) + 1):(3*round(n/5))], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx4 <- mclapply(as.character(unique(mirbasesubset$query_id))[((3*round(n/5)) + 1):(4*round(n/5))], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx5 <- mclapply(as.character(unique(mirbasesubset$query_id))[((4*round(n/5)) + 1):n], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx <- c(idx1,idx2,idx3,idx4,idx5)
})
```

```
##     user   system  elapsed 
## 4916.470   22.560  666.396
```

```r
length(idx)
```

```
## [1] 147657
```

Function to filter multiple miRBase blast hits per sequence


```r
filter <- function(seqblast){
# Sequence has only one miRBase blast hit
        if (nrow(seqblast) == 1){
                return(seqblast)
        }
# Sequence with 100% match to miRBase blast hit
        ident <- seqblast[sapply(seqblast$perc_identical, function(x) x == 100),]
                if (nrow(ident) == 1){
                        return(ident)
        }
# Sequence with greater than 96% match to miRBase blast hit
        ident2 <- seqblast[sapply(seqblast$perc_identical, function(x) x > 96),]
                if (nrow(ident2) == 1){
                        return(ident2)
        }
# Select miRBase blast hit with the greatest bitscore for a sequence
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
##  48.128   8.977  10.201
```

```r
length(stp1)
```

```
## [1] 147657
```

```r
head(stp1)
```

```
## [[1]]
##         query_id  dbseq_id perc_identical evalue bitscore
## 1 seq_0_x9179927 ssc-miR-1            100  2e-08     42.1
## 
## [[2]]
##               query_id  dbseq_id perc_identical evalue bitscore
## 2 seq_9179927_x5732924 ssc-miR-1            100  2e-08     42.1
## 
## [[3]]
##                query_id    dbseq_id perc_identical evalue bitscore
## 3 seq_14912851_x5116126 ssc-miR-206            100  4e-09     44.1
## 
## [[4]]
##                query_id  dbseq_id perc_identical evalue bitscore
## 4 seq_20028977_x4202989 ssc-miR-1            100  2e-08     42.1
## 
## [[5]]
##                query_id  dbseq_id perc_identical evalue bitscore
## 5 seq_24231966_x3505397 ssc-miR-1            100  6e-08     40.1
## 
## [[6]]
##                query_id  dbseq_id perc_identical evalue bitscore
## 6 seq_27737363_x3395085 ssc-miR-1            100  2e-08     42.1
```

```r
tail(stp1)
```

```
## [[1]]
##                query_id    dbseq_id perc_identical evalue bitscore
## 178249 seq_227425548_x2 ssc-miR-184            100  8e-07     36.2
## 
## [[2]]
##                query_id    dbseq_id perc_identical evalue bitscore
## 178250 seq_227425570_x2 ssc-miR-99a          95.24  4e-06     34.2
## 
## [[3]]
##                query_id       dbseq_id perc_identical evalue bitscore
## 178251 seq_227425572_x2 ssc-miR-30c-5p            100  2e-09     46.1
## 
## [[4]]
##                query_id       dbseq_id perc_identical evalue bitscore
## 178253 seq_227425682_x2 ssc-miR-423-3p          95.45  1e-06     36.2
## 
## [[5]]
##                query_id     dbseq_id perc_identical evalue bitscore
## 178254 seq_227425702_x2 ssc-miR-181b            100  1e-06     36.2
## 
## [[6]]
##                query_id       dbseq_id perc_identical evalue bitscore
## 178255 seq_227425704_x2 ssc-miR-425-3p          95.24  4e-06     34.2
```

How many sequences have more than one hit after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 6749
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##    2    3 
## 6466  283
```

```r
filter2 <- filter <- function(seqblast){
# Sequence has only one miRBase blast hit
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
##   8.091   4.263   3.935
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

Summary file of small RNA sequence miRBase blast results


```r
sumblast <- do.call(rbind,stp3)
head(sumblast)
```

```
##                query_id    dbseq_id perc_identical evalue bitscore
## 1        seq_0_x9179927   ssc-miR-1            100  2e-08     42.1
## 2  seq_9179927_x5732924   ssc-miR-1            100  2e-08     42.1
## 3 seq_14912851_x5116126 ssc-miR-206            100  4e-09     44.1
## 4 seq_20028977_x4202989   ssc-miR-1            100  2e-08     42.1
## 5 seq_24231966_x3505397   ssc-miR-1            100  6e-08     40.1
## 6 seq_27737363_x3395085   ssc-miR-1            100  2e-08     42.1
```

```r
dim(sumblast)
```

```
## [1] 147657      5
```

```r
length(unique(sumblast$query_id))
```

```
## [1] 147657
```

```r
table(sumblast$dbseq_id)
```

```
## 
##       ssc-let-7a       ssc-let-7c    ssc-let-7d-3p    ssc-let-7d-5p 
##             2278             1340              288              696 
##       ssc-let-7e       ssc-let-7f       ssc-let-7g       ssc-let-7i 
##              661             1708             1389              567 
##        ssc-miR-1      ssc-miR-100      ssc-miR-101      ssc-miR-103 
##            19897              675             1771              345 
##    ssc-miR-105-1    ssc-miR-105-2     ssc-miR-106a      ssc-miR-107 
##                6                5               55              398 
##   ssc-miR-10a-3p   ssc-miR-10a-5p      ssc-miR-10b      ssc-miR-122 
##               63              891             3710              431 
##     ssc-miR-1224     ssc-miR-1249     ssc-miR-124a     ssc-miR-125a 
##                1              220                2              740 
##     ssc-miR-125b   ssc-miR-126-3p   ssc-miR-126-5p      ssc-miR-127 
##             1001             1898              951             1063 
##     ssc-miR-1271     ssc-miR-1277      ssc-miR-128     ssc-miR-1285 
##               40               35             1043              833 
##  ssc-miR-1296-3p  ssc-miR-1296-5p  ssc-miR-129a-3p     ssc-miR-129b 
##                2              178               27               44 
##  ssc-miR-1306-3p  ssc-miR-1306-5p     ssc-miR-1307     ssc-miR-130a 
##               95              196              323              234 
##     ssc-miR-130b      ssc-miR-132  ssc-miR-133a-3p  ssc-miR-133a-5p 
##              158               92             3604              816 
##     ssc-miR-133b     ssc-miR-1343      ssc-miR-135      ssc-miR-136 
##             3837              125               14              206 
##      ssc-miR-137      ssc-miR-138   ssc-miR-139-3p   ssc-miR-139-5p 
##               10               27               33              177 
##   ssc-miR-140-3p   ssc-miR-140-5p   ssc-miR-142-3p   ssc-miR-142-5p 
##             1125              217              442              214 
##   ssc-miR-143-3p   ssc-miR-143-5p      ssc-miR-144   ssc-miR-145-3p 
##             2916              132               60              231 
##   ssc-miR-145-5p     ssc-miR-1468  ssc-miR-146a-5p     ssc-miR-146b 
##              538              257              643              339 
##  ssc-miR-148a-3p  ssc-miR-148a-5p  ssc-miR-148b-3p  ssc-miR-148b-5p 
##              768               64              370               84 
##      ssc-miR-149      ssc-miR-150   ssc-miR-151-3p   ssc-miR-151-5p 
##              171              368             1136             1160 
##      ssc-miR-152      ssc-miR-153   ssc-miR-155-5p      ssc-miR-15a 
##              665               12              307              768 
##      ssc-miR-15b       ssc-miR-16    ssc-miR-17-3p    ssc-miR-17-5p 
##              209             1541               95              358 
##     ssc-miR-181a     ssc-miR-181b     ssc-miR-181c  ssc-miR-181d-3p 
##             1761             1331              300                3 
##  ssc-miR-181d-5p      ssc-miR-182      ssc-miR-183  ssc-miR-1839-3p 
##              266              131               28                5 
##  ssc-miR-1839-5p      ssc-miR-184      ssc-miR-185      ssc-miR-186 
##               15               65              141              763 
##      ssc-miR-187      ssc-miR-18a      ssc-miR-18b     ssc-miR-190a 
##               16               94               40              429 
##     ssc-miR-190b      ssc-miR-191      ssc-miR-192  ssc-miR-193a-3p 
##               49              471              171              158 
##  ssc-miR-193a-5p     ssc-miR-194a  ssc-miR-194b-5p      ssc-miR-195 
##              114              134                8              511 
##     ssc-miR-196a  ssc-miR-196b-3p  ssc-miR-196b-5p  ssc-miR-199a-3p 
##              454               35              527             1575 
##  ssc-miR-199a-5p  ssc-miR-199b-3p  ssc-miR-199b-5p      ssc-miR-19a 
##              467              541              349              201 
##      ssc-miR-19b   ssc-miR-202-5p      ssc-miR-204      ssc-miR-205 
##              321               20              103               44 
##      ssc-miR-206     ssc-miR-208b      ssc-miR-20a      ssc-miR-20b 
##             3687              368              391               54 
##       ssc-miR-21      ssc-miR-210      ssc-miR-212      ssc-miR-214 
##             1581              159                8              667 
##      ssc-miR-215      ssc-miR-216      ssc-miR-217      ssc-miR-218 
##                6               12               14               77 
##   ssc-miR-218-3p     ssc-miR-218b     ssc-miR-219a  ssc-miR-219b-3p 
##               14               93               15                6 
##   ssc-miR-221-3p   ssc-miR-221-5p      ssc-miR-222    ssc-miR-22-3p 
##              248               29              218             1459 
##      ssc-miR-224    ssc-miR-22-5p  ssc-miR-2320-3p  ssc-miR-2320-5p 
##                4              337               59              409 
##     ssc-miR-2366      ssc-miR-23a      ssc-miR-23b     ssc-miR-2411 
##               23              409              708               26 
##  ssc-miR-24-1-5p  ssc-miR-24-2-5p    ssc-miR-24-3p     ssc-miR-2483 
##              178              204             1752               44 
##      ssc-miR-26a      ssc-miR-27a   ssc-miR-27b-3p   ssc-miR-27b-5p 
##             2710              462             1175               91 
##    ssc-miR-28-3p    ssc-miR-28-5p   ssc-miR-296-3p   ssc-miR-296-5p 
##              465              439              158               72 
##      ssc-miR-299      ssc-miR-29a      ssc-miR-29b      ssc-miR-29c 
##              299              991              271              671 
##      ssc-miR-301   ssc-miR-30a-3p   ssc-miR-30a-5p   ssc-miR-30b-3p 
##               48              343             1607               58 
##   ssc-miR-30b-5p ssc-miR-30c-1-3p   ssc-miR-30c-3p   ssc-miR-30c-5p 
##              595               40               51             1146 
##      ssc-miR-30d   ssc-miR-30e-3p   ssc-miR-30e-5p       ssc-miR-31 
##             1675              487             1906               13 
##       ssc-miR-32      ssc-miR-320      ssc-miR-323      ssc-miR-324 
##              184             1818              145              108 
##      ssc-miR-325      ssc-miR-326      ssc-miR-328   ssc-miR-331-3p 
##                3               76              458              392 
##   ssc-miR-331-5p      ssc-miR-335      ssc-miR-338      ssc-miR-339 
##              156              219              231              327 
##   ssc-miR-339-3p   ssc-miR-339-5p      ssc-miR-340      ssc-miR-342 
##               63              200              379              283 
##   ssc-miR-345-3p   ssc-miR-345-5p      ssc-miR-34a      ssc-miR-34c 
##              165               77              337                7 
##     ssc-miR-3613   ssc-miR-361-3p   ssc-miR-361-5p      ssc-miR-362 
##              306              219              364              293 
##      ssc-miR-363   ssc-miR-365-3p   ssc-miR-365-5p      ssc-miR-369 
##               83              315               11              271 
##      ssc-miR-370  ssc-miR-374a-3p  ssc-miR-374a-5p  ssc-miR-374b-3p 
##              209              317             1107               70 
##  ssc-miR-374b-5p  ssc-miR-376a-3p  ssc-miR-376a-5p     ssc-miR-376c 
##              345              108               63               38 
##      ssc-miR-378  ssc-miR-378b-3p      ssc-miR-381      ssc-miR-382 
##             3705             1050               25              136 
##      ssc-miR-383      ssc-miR-411   ssc-miR-421-3p   ssc-miR-421-5p 
##                2               90              102               12 
##   ssc-miR-423-3p   ssc-miR-423-5p   ssc-miR-424-3p   ssc-miR-424-5p 
##              913              466              135              373 
##   ssc-miR-425-3p   ssc-miR-425-5p      ssc-miR-429   ssc-miR-432-3p 
##              285              862               15                4 
##   ssc-miR-432-5p     ssc-miR-4331     ssc-miR-4332     ssc-miR-4337 
##              207               66              135               11 
##     ssc-miR-4338     ssc-miR-450a  ssc-miR-450b-3p  ssc-miR-450b-5p 
##                2               64               68              396 
##  ssc-miR-450c-3p  ssc-miR-450c-5p      ssc-miR-451      ssc-miR-452 
##               83              311              433               31 
##   ssc-miR-455-3p   ssc-miR-455-5p      ssc-miR-486     ssc-miR-487b 
##              140              145             3538               19 
##      ssc-miR-489      ssc-miR-490   ssc-miR-490-5p      ssc-miR-491 
##              115               10                5               74 
##   ssc-miR-493-3p   ssc-miR-493-5p      ssc-miR-497   ssc-miR-499-3p 
##               72               99              269              197 
##   ssc-miR-499-5p      ssc-miR-500      ssc-miR-503      ssc-miR-504 
##              645              375              489              165 
##      ssc-miR-505   ssc-miR-532-3p   ssc-miR-532-5p   ssc-miR-542-3p 
##              195              138              555              553 
##   ssc-miR-542-5p   ssc-miR-545-3p   ssc-miR-545-5p     ssc-miR-551a 
##              213              126               86                5 
##      ssc-miR-574      ssc-miR-582      ssc-miR-615      ssc-miR-628 
##              227               34               27              109 
##      ssc-miR-652   ssc-miR-664-3p   ssc-miR-664-5p   ssc-miR-671-3p 
##               25              128              137               28 
##   ssc-miR-671-5p   ssc-miR-676-3p   ssc-miR-676-5p  ssc-miR-6782-3p 
##               97              158                2              194 
##        ssc-miR-7   ssc-miR-708-3p   ssc-miR-708-5p  ssc-miR-7134-3p 
##              236               86              241              853 
##  ssc-miR-7134-5p  ssc-miR-7135-3p  ssc-miR-7135-5p  ssc-miR-7136-3p 
##              333               31               23                3 
##  ssc-miR-7136-5p  ssc-miR-7137-3p  ssc-miR-7137-5p  ssc-miR-7138-3p 
##                8              209               46               11 
##  ssc-miR-7138-5p  ssc-miR-7139-3p  ssc-miR-7139-5p  ssc-miR-7140-3p 
##               14               10               12                7 
##  ssc-miR-7140-5p  ssc-miR-7142-3p  ssc-miR-7142-5p  ssc-miR-7143-3p 
##                1              166                9                5 
##  ssc-miR-7143-5p  ssc-miR-7144-5p      ssc-miR-744      ssc-miR-758 
##                1               12               83              168 
##   ssc-miR-769-3p   ssc-miR-769-5p  ssc-miR-7857-3p      ssc-miR-874 
##               61              225               16              121 
##   ssc-miR-885-3p   ssc-miR-885-5p        ssc-miR-9      ssc-miR-9-2 
##              218              301               39               16 
##      ssc-miR-92a   ssc-miR-92b-3p   ssc-miR-92b-5p      ssc-miR-935 
##              874              105                6              101 
##       ssc-miR-95    ssc-miR-96-5p  ssc-miR-9785-5p  ssc-miR-9788-3p 
##              438                6               26                2 
##       ssc-miR-98  ssc-miR-9810-3p  ssc-miR-9820-5p  ssc-miR-9841-3p 
##              456               31               10               90 
##  ssc-miR-9843-3p  ssc-miR-9851-3p  ssc-miR-9858-5p  ssc-miR-9860-5p 
##              122               45                5               65 
##      ssc-miR-99a      ssc-miR-99b 
##             1025              500
```

```r
uniqsumblast<-as.data.frame(table(sumblast$dbseq_id))
colnames(uniqsumblast)<-c("miRNA", "Freq")
uniqsumblast
```

```
##                miRNA  Freq
## 1         ssc-let-7a  2278
## 2         ssc-let-7c  1340
## 3      ssc-let-7d-3p   288
## 4      ssc-let-7d-5p   696
## 5         ssc-let-7e   661
## 6         ssc-let-7f  1708
## 7         ssc-let-7g  1389
## 8         ssc-let-7i   567
## 9          ssc-miR-1 19897
## 10       ssc-miR-100   675
## 11       ssc-miR-101  1771
## 12       ssc-miR-103   345
## 13     ssc-miR-105-1     6
## 14     ssc-miR-105-2     5
## 15      ssc-miR-106a    55
## 16       ssc-miR-107   398
## 17    ssc-miR-10a-3p    63
## 18    ssc-miR-10a-5p   891
## 19       ssc-miR-10b  3710
## 20       ssc-miR-122   431
## 21      ssc-miR-1224     1
## 22      ssc-miR-1249   220
## 23      ssc-miR-124a     2
## 24      ssc-miR-125a   740
## 25      ssc-miR-125b  1001
## 26    ssc-miR-126-3p  1898
## 27    ssc-miR-126-5p   951
## 28       ssc-miR-127  1063
## 29      ssc-miR-1271    40
## 30      ssc-miR-1277    35
## 31       ssc-miR-128  1043
## 32      ssc-miR-1285   833
## 33   ssc-miR-1296-3p     2
## 34   ssc-miR-1296-5p   178
## 35   ssc-miR-129a-3p    27
## 36      ssc-miR-129b    44
## 37   ssc-miR-1306-3p    95
## 38   ssc-miR-1306-5p   196
## 39      ssc-miR-1307   323
## 40      ssc-miR-130a   234
## 41      ssc-miR-130b   158
## 42       ssc-miR-132    92
## 43   ssc-miR-133a-3p  3604
## 44   ssc-miR-133a-5p   816
## 45      ssc-miR-133b  3837
## 46      ssc-miR-1343   125
## 47       ssc-miR-135    14
## 48       ssc-miR-136   206
## 49       ssc-miR-137    10
## 50       ssc-miR-138    27
## 51    ssc-miR-139-3p    33
## 52    ssc-miR-139-5p   177
## 53    ssc-miR-140-3p  1125
## 54    ssc-miR-140-5p   217
## 55    ssc-miR-142-3p   442
## 56    ssc-miR-142-5p   214
## 57    ssc-miR-143-3p  2916
## 58    ssc-miR-143-5p   132
## 59       ssc-miR-144    60
## 60    ssc-miR-145-3p   231
## 61    ssc-miR-145-5p   538
## 62      ssc-miR-1468   257
## 63   ssc-miR-146a-5p   643
## 64      ssc-miR-146b   339
## 65   ssc-miR-148a-3p   768
## 66   ssc-miR-148a-5p    64
## 67   ssc-miR-148b-3p   370
## 68   ssc-miR-148b-5p    84
## 69       ssc-miR-149   171
## 70       ssc-miR-150   368
## 71    ssc-miR-151-3p  1136
## 72    ssc-miR-151-5p  1160
## 73       ssc-miR-152   665
## 74       ssc-miR-153    12
## 75    ssc-miR-155-5p   307
## 76       ssc-miR-15a   768
## 77       ssc-miR-15b   209
## 78        ssc-miR-16  1541
## 79     ssc-miR-17-3p    95
## 80     ssc-miR-17-5p   358
## 81      ssc-miR-181a  1761
## 82      ssc-miR-181b  1331
## 83      ssc-miR-181c   300
## 84   ssc-miR-181d-3p     3
## 85   ssc-miR-181d-5p   266
## 86       ssc-miR-182   131
## 87       ssc-miR-183    28
## 88   ssc-miR-1839-3p     5
## 89   ssc-miR-1839-5p    15
## 90       ssc-miR-184    65
## 91       ssc-miR-185   141
## 92       ssc-miR-186   763
## 93       ssc-miR-187    16
## 94       ssc-miR-18a    94
## 95       ssc-miR-18b    40
## 96      ssc-miR-190a   429
## 97      ssc-miR-190b    49
## 98       ssc-miR-191   471
## 99       ssc-miR-192   171
## 100  ssc-miR-193a-3p   158
## 101  ssc-miR-193a-5p   114
## 102     ssc-miR-194a   134
## 103  ssc-miR-194b-5p     8
## 104      ssc-miR-195   511
## 105     ssc-miR-196a   454
## 106  ssc-miR-196b-3p    35
## 107  ssc-miR-196b-5p   527
## 108  ssc-miR-199a-3p  1575
## 109  ssc-miR-199a-5p   467
## 110  ssc-miR-199b-3p   541
## 111  ssc-miR-199b-5p   349
## 112      ssc-miR-19a   201
## 113      ssc-miR-19b   321
## 114   ssc-miR-202-5p    20
## 115      ssc-miR-204   103
## 116      ssc-miR-205    44
## 117      ssc-miR-206  3687
## 118     ssc-miR-208b   368
## 119      ssc-miR-20a   391
## 120      ssc-miR-20b    54
## 121       ssc-miR-21  1581
## 122      ssc-miR-210   159
## 123      ssc-miR-212     8
## 124      ssc-miR-214   667
## 125      ssc-miR-215     6
## 126      ssc-miR-216    12
## 127      ssc-miR-217    14
## 128      ssc-miR-218    77
## 129   ssc-miR-218-3p    14
## 130     ssc-miR-218b    93
## 131     ssc-miR-219a    15
## 132  ssc-miR-219b-3p     6
## 133   ssc-miR-221-3p   248
## 134   ssc-miR-221-5p    29
## 135      ssc-miR-222   218
## 136    ssc-miR-22-3p  1459
## 137      ssc-miR-224     4
## 138    ssc-miR-22-5p   337
## 139  ssc-miR-2320-3p    59
## 140  ssc-miR-2320-5p   409
## 141     ssc-miR-2366    23
## 142      ssc-miR-23a   409
## 143      ssc-miR-23b   708
## 144     ssc-miR-2411    26
## 145  ssc-miR-24-1-5p   178
## 146  ssc-miR-24-2-5p   204
## 147    ssc-miR-24-3p  1752
## 148     ssc-miR-2483    44
## 149      ssc-miR-26a  2710
## 150      ssc-miR-27a   462
## 151   ssc-miR-27b-3p  1175
## 152   ssc-miR-27b-5p    91
## 153    ssc-miR-28-3p   465
## 154    ssc-miR-28-5p   439
## 155   ssc-miR-296-3p   158
## 156   ssc-miR-296-5p    72
## 157      ssc-miR-299   299
## 158      ssc-miR-29a   991
## 159      ssc-miR-29b   271
## 160      ssc-miR-29c   671
## 161      ssc-miR-301    48
## 162   ssc-miR-30a-3p   343
## 163   ssc-miR-30a-5p  1607
## 164   ssc-miR-30b-3p    58
## 165   ssc-miR-30b-5p   595
## 166 ssc-miR-30c-1-3p    40
## 167   ssc-miR-30c-3p    51
## 168   ssc-miR-30c-5p  1146
## 169      ssc-miR-30d  1675
## 170   ssc-miR-30e-3p   487
## 171   ssc-miR-30e-5p  1906
## 172       ssc-miR-31    13
## 173       ssc-miR-32   184
## 174      ssc-miR-320  1818
## 175      ssc-miR-323   145
## 176      ssc-miR-324   108
## 177      ssc-miR-325     3
## 178      ssc-miR-326    76
## 179      ssc-miR-328   458
## 180   ssc-miR-331-3p   392
## 181   ssc-miR-331-5p   156
## 182      ssc-miR-335   219
## 183      ssc-miR-338   231
## 184      ssc-miR-339   327
## 185   ssc-miR-339-3p    63
## 186   ssc-miR-339-5p   200
## 187      ssc-miR-340   379
## 188      ssc-miR-342   283
## 189   ssc-miR-345-3p   165
## 190   ssc-miR-345-5p    77
## 191      ssc-miR-34a   337
## 192      ssc-miR-34c     7
## 193     ssc-miR-3613   306
## 194   ssc-miR-361-3p   219
## 195   ssc-miR-361-5p   364
## 196      ssc-miR-362   293
## 197      ssc-miR-363    83
## 198   ssc-miR-365-3p   315
## 199   ssc-miR-365-5p    11
## 200      ssc-miR-369   271
## 201      ssc-miR-370   209
## 202  ssc-miR-374a-3p   317
## 203  ssc-miR-374a-5p  1107
## 204  ssc-miR-374b-3p    70
## 205  ssc-miR-374b-5p   345
## 206  ssc-miR-376a-3p   108
## 207  ssc-miR-376a-5p    63
## 208     ssc-miR-376c    38
## 209      ssc-miR-378  3705
## 210  ssc-miR-378b-3p  1050
## 211      ssc-miR-381    25
## 212      ssc-miR-382   136
## 213      ssc-miR-383     2
## 214      ssc-miR-411    90
## 215   ssc-miR-421-3p   102
## 216   ssc-miR-421-5p    12
## 217   ssc-miR-423-3p   913
## 218   ssc-miR-423-5p   466
## 219   ssc-miR-424-3p   135
## 220   ssc-miR-424-5p   373
## 221   ssc-miR-425-3p   285
## 222   ssc-miR-425-5p   862
## 223      ssc-miR-429    15
## 224   ssc-miR-432-3p     4
## 225   ssc-miR-432-5p   207
## 226     ssc-miR-4331    66
## 227     ssc-miR-4332   135
## 228     ssc-miR-4337    11
## 229     ssc-miR-4338     2
## 230     ssc-miR-450a    64
## 231  ssc-miR-450b-3p    68
## 232  ssc-miR-450b-5p   396
## 233  ssc-miR-450c-3p    83
## 234  ssc-miR-450c-5p   311
## 235      ssc-miR-451   433
## 236      ssc-miR-452    31
## 237   ssc-miR-455-3p   140
## 238   ssc-miR-455-5p   145
## 239      ssc-miR-486  3538
## 240     ssc-miR-487b    19
## 241      ssc-miR-489   115
## 242      ssc-miR-490    10
## 243   ssc-miR-490-5p     5
## 244      ssc-miR-491    74
## 245   ssc-miR-493-3p    72
## 246   ssc-miR-493-5p    99
## 247      ssc-miR-497   269
## 248   ssc-miR-499-3p   197
## 249   ssc-miR-499-5p   645
## 250      ssc-miR-500   375
## 251      ssc-miR-503   489
## 252      ssc-miR-504   165
## 253      ssc-miR-505   195
## 254   ssc-miR-532-3p   138
## 255   ssc-miR-532-5p   555
## 256   ssc-miR-542-3p   553
## 257   ssc-miR-542-5p   213
## 258   ssc-miR-545-3p   126
## 259   ssc-miR-545-5p    86
## 260     ssc-miR-551a     5
## 261      ssc-miR-574   227
## 262      ssc-miR-582    34
## 263      ssc-miR-615    27
## 264      ssc-miR-628   109
## 265      ssc-miR-652    25
## 266   ssc-miR-664-3p   128
## 267   ssc-miR-664-5p   137
## 268   ssc-miR-671-3p    28
## 269   ssc-miR-671-5p    97
## 270   ssc-miR-676-3p   158
## 271   ssc-miR-676-5p     2
## 272  ssc-miR-6782-3p   194
## 273        ssc-miR-7   236
## 274   ssc-miR-708-3p    86
## 275   ssc-miR-708-5p   241
## 276  ssc-miR-7134-3p   853
## 277  ssc-miR-7134-5p   333
## 278  ssc-miR-7135-3p    31
## 279  ssc-miR-7135-5p    23
## 280  ssc-miR-7136-3p     3
## 281  ssc-miR-7136-5p     8
## 282  ssc-miR-7137-3p   209
## 283  ssc-miR-7137-5p    46
## 284  ssc-miR-7138-3p    11
## 285  ssc-miR-7138-5p    14
## 286  ssc-miR-7139-3p    10
## 287  ssc-miR-7139-5p    12
## 288  ssc-miR-7140-3p     7
## 289  ssc-miR-7140-5p     1
## 290  ssc-miR-7142-3p   166
## 291  ssc-miR-7142-5p     9
## 292  ssc-miR-7143-3p     5
## 293  ssc-miR-7143-5p     1
## 294  ssc-miR-7144-5p    12
## 295      ssc-miR-744    83
## 296      ssc-miR-758   168
## 297   ssc-miR-769-3p    61
## 298   ssc-miR-769-5p   225
## 299  ssc-miR-7857-3p    16
## 300      ssc-miR-874   121
## 301   ssc-miR-885-3p   218
## 302   ssc-miR-885-5p   301
## 303        ssc-miR-9    39
## 304      ssc-miR-9-2    16
## 305      ssc-miR-92a   874
## 306   ssc-miR-92b-3p   105
## 307   ssc-miR-92b-5p     6
## 308      ssc-miR-935   101
## 309       ssc-miR-95   438
## 310    ssc-miR-96-5p     6
## 311  ssc-miR-9785-5p    26
## 312  ssc-miR-9788-3p     2
## 313       ssc-miR-98   456
## 314  ssc-miR-9810-3p    31
## 315  ssc-miR-9820-5p    10
## 316  ssc-miR-9841-3p    90
## 317  ssc-miR-9843-3p   122
## 318  ssc-miR-9851-3p    45
## 319  ssc-miR-9858-5p     5
## 320  ssc-miR-9860-5p    65
## 321      ssc-miR-99a  1025
## 322      ssc-miR-99b   500
```

Add the column of sequence count to the sumblast data frame


```r
sumblast$seq_count<-as.numeric(str_split_fixed(sumblast$query_id, "_x", 2)[,2])
head(sumblast)
```

```
##                query_id    dbseq_id perc_identical evalue bitscore
## 1        seq_0_x9179927   ssc-miR-1            100  2e-08     42.1
## 2  seq_9179927_x5732924   ssc-miR-1            100  2e-08     42.1
## 3 seq_14912851_x5116126 ssc-miR-206            100  4e-09     44.1
## 4 seq_20028977_x4202989   ssc-miR-1            100  2e-08     42.1
## 5 seq_24231966_x3505397   ssc-miR-1            100  6e-08     40.1
## 6 seq_27737363_x3395085   ssc-miR-1            100  2e-08     42.1
##   seq_count
## 1   9179927
## 2   5732924
## 3   5116126
## 4   4202989
## 5   3505397
## 6   3395085
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype<-as.matrix(by(sumblast$seq_count, sumblast$dbseq_id, sum))
totalsumbiotype
```

```
##                      [,1]
## ssc-let-7a        5072114
## ssc-let-7c        3430928
## ssc-let-7d-3p       51978
## ssc-let-7d-5p      493235
## ssc-let-7e         260745
## ssc-let-7f        3765039
## ssc-let-7g        1843605
## ssc-let-7i         975379
## ssc-miR-1        47655223
## ssc-miR-100        479119
## ssc-miR-101       1681080
## ssc-miR-103        189086
## ssc-miR-105-1          53
## ssc-miR-105-2          21
## ssc-miR-106a         3541
## ssc-miR-107         70680
## ssc-miR-10a-3p       3805
## ssc-miR-10a-5p     560989
## ssc-miR-10b       6985290
## ssc-miR-122        115564
## ssc-miR-1224            9
## ssc-miR-1249        28259
## ssc-miR-124a            8
## ssc-miR-125a       416054
## ssc-miR-125b       994307
## ssc-miR-126-3p    2760992
## ssc-miR-126-5p    1072816
## ssc-miR-127        331255
## ssc-miR-1271         1561
## ssc-miR-1277         2958
## ssc-miR-128        596192
## ssc-miR-1285       354830
## ssc-miR-1296-3p        17
## ssc-miR-1296-5p     17455
## ssc-miR-129a-3p       743
## ssc-miR-129b         1636
## ssc-miR-1306-3p      3981
## ssc-miR-1306-5p     30819
## ssc-miR-1307        37815
## ssc-miR-130a        26447
## ssc-miR-130b        16053
## ssc-miR-132          3806
## ssc-miR-133a-3p   5364699
## ssc-miR-133a-5p    404245
## ssc-miR-133b      4448595
## ssc-miR-1343         7447
## ssc-miR-135           475
## ssc-miR-136         27363
## ssc-miR-137           114
## ssc-miR-138           441
## ssc-miR-139-3p       1111
## ssc-miR-139-5p      15439
## ssc-miR-140-3p     246322
## ssc-miR-140-5p      42975
## ssc-miR-142-3p     122218
## ssc-miR-142-5p      52184
## ssc-miR-143-3p    5235493
## ssc-miR-143-5p       7733
## ssc-miR-144         34608
## ssc-miR-145-3p      44572
## ssc-miR-145-5p     133576
## ssc-miR-1468        57169
## ssc-miR-146a-5p    121603
## ssc-miR-146b        52909
## ssc-miR-148a-3p    812453
## ssc-miR-148a-5p     10206
## ssc-miR-148b-3p    166894
## ssc-miR-148b-5p      5247
## ssc-miR-149         15784
## ssc-miR-150        118425
## ssc-miR-151-3p     499161
## ssc-miR-151-5p    1116309
## ssc-miR-152        487103
## ssc-miR-153           180
## ssc-miR-155-5p      37302
## ssc-miR-15a        386599
## ssc-miR-15b         60580
## ssc-miR-16        1400729
## ssc-miR-17-3p        5724
## ssc-miR-17-5p       67959
## ssc-miR-181a      2072932
## ssc-miR-181b       493491
## ssc-miR-181c        56248
## ssc-miR-181d-3p        18
## ssc-miR-181d-5p     30410
## ssc-miR-182         11115
## ssc-miR-183           796
## ssc-miR-1839-3p        23
## ssc-miR-1839-5p       498
## ssc-miR-184          3245
## ssc-miR-185         17706
## ssc-miR-186        383248
## ssc-miR-187           127
## ssc-miR-18a         13365
## ssc-miR-18b          2207
## ssc-miR-190a       118234
## ssc-miR-190b         3377
## ssc-miR-191        127759
## ssc-miR-192         24955
## ssc-miR-193a-3p     20893
## ssc-miR-193a-5p     10205
## ssc-miR-194a        14480
## ssc-miR-194b-5p        42
## ssc-miR-195        318509
## ssc-miR-196a       138625
## ssc-miR-196b-3p      1257
## ssc-miR-196b-5p    151788
## ssc-miR-199a-3p   1790989
## ssc-miR-199a-5p    301561
## ssc-miR-199b-3p    466941
## ssc-miR-199b-5p    204274
## ssc-miR-19a         34440
## ssc-miR-19b        150396
## ssc-miR-202-5p       1144
## ssc-miR-204         10969
## ssc-miR-205          2839
## ssc-miR-206       8194252
## ssc-miR-208b       125821
## ssc-miR-20a        135729
## ssc-miR-20b          3608
## ssc-miR-21        3380630
## ssc-miR-210          6625
## ssc-miR-212           197
## ssc-miR-214         53986
## ssc-miR-215            33
## ssc-miR-216           153
## ssc-miR-217           230
## ssc-miR-218         11552
## ssc-miR-218-3p        133
## ssc-miR-218b        11017
## ssc-miR-219a          132
## ssc-miR-219b-3p        80
## ssc-miR-221-3p      33515
## ssc-miR-221-5p        526
## ssc-miR-222         17290
## ssc-miR-22-3p     1187902
## ssc-miR-224            32
## ssc-miR-22-5p      181130
## ssc-miR-2320-3p      1817
## ssc-miR-2320-5p     25273
## ssc-miR-2366          292
## ssc-miR-23a        152832
## ssc-miR-23b        146534
## ssc-miR-2411          206
## ssc-miR-24-1-5p     14410
## ssc-miR-24-2-5p     22309
## ssc-miR-24-3p     1352514
## ssc-miR-2483         1579
## ssc-miR-26a       4396286
## ssc-miR-27a        199105
## ssc-miR-27b-3p    1287473
## ssc-miR-27b-5p       5541
## ssc-miR-28-3p      105343
## ssc-miR-28-5p      139093
## ssc-miR-296-3p       8592
## ssc-miR-296-5p       4371
## ssc-miR-299         37653
## ssc-miR-29a        901005
## ssc-miR-29b         45410
## ssc-miR-29c        433610
## ssc-miR-301          9025
## ssc-miR-30a-3p      87320
## ssc-miR-30a-5p    2457086
## ssc-miR-30b-3p       3629
## ssc-miR-30b-5p     506838
## ssc-miR-30c-1-3p     2201
## ssc-miR-30c-3p       3362
## ssc-miR-30c-5p    1323874
## ssc-miR-30d       1866212
## ssc-miR-30e-3p     188019
## ssc-miR-30e-5p    2375573
## ssc-miR-31            233
## ssc-miR-32          44798
## ssc-miR-320        528479
## ssc-miR-323         19439
## ssc-miR-324          7405
## ssc-miR-325            17
## ssc-miR-326          4338
## ssc-miR-328         65876
## ssc-miR-331-3p      68418
## ssc-miR-331-5p       5086
## ssc-miR-335         53062
## ssc-miR-338         34033
## ssc-miR-339         53376
## ssc-miR-339-3p       1392
## ssc-miR-339-5p      40027
## ssc-miR-340        312085
## ssc-miR-342         39868
## ssc-miR-345-3p      17249
## ssc-miR-345-5p      10309
## ssc-miR-34a         48999
## ssc-miR-34c           141
## ssc-miR-3613       259600
## ssc-miR-361-3p      19786
## ssc-miR-361-5p      87132
## ssc-miR-362         51737
## ssc-miR-363          5316
## ssc-miR-365-3p      73385
## ssc-miR-365-5p        105
## ssc-miR-369        169376
## ssc-miR-370          8114
## ssc-miR-374a-3p    115119
## ssc-miR-374a-5p   1213920
## ssc-miR-374b-3p      2880
## ssc-miR-374b-5p     86740
## ssc-miR-376a-3p     12777
## ssc-miR-376a-5p      2518
## ssc-miR-376c         1866
## ssc-miR-378       8927859
## ssc-miR-378b-3p    824944
## ssc-miR-381           541
## ssc-miR-382         17557
## ssc-miR-383            14
## ssc-miR-411         13524
## ssc-miR-421-3p       5381
## ssc-miR-421-5p         61
## ssc-miR-423-3p     218821
## ssc-miR-423-5p      78061
## ssc-miR-424-3p      11105
## ssc-miR-424-5p     201549
## ssc-miR-425-3p      25038
## ssc-miR-425-5p     323183
## ssc-miR-429           289
## ssc-miR-432-3p         19
## ssc-miR-432-5p      17128
## ssc-miR-4331          305
## ssc-miR-4332         3012
## ssc-miR-4337           67
## ssc-miR-4338            4
## ssc-miR-450a        16510
## ssc-miR-450b-3p     18011
## ssc-miR-450b-5p    120754
## ssc-miR-450c-3p      7292
## ssc-miR-450c-5p    101896
## ssc-miR-451        239854
## ssc-miR-452           905
## ssc-miR-455-3p      10728
## ssc-miR-455-5p      23292
## ssc-miR-486       4645422
## ssc-miR-487b          277
## ssc-miR-489          6493
## ssc-miR-490           499
## ssc-miR-490-5p         21
## ssc-miR-491          4413
## ssc-miR-493-3p       2385
## ssc-miR-493-5p       7122
## ssc-miR-497         65871
## ssc-miR-499-3p      18025
## ssc-miR-499-5p     663073
## ssc-miR-500         34434
## ssc-miR-503         50718
## ssc-miR-504         15552
## ssc-miR-505         15776
## ssc-miR-532-3p      17485
## ssc-miR-532-5p     278033
## ssc-miR-542-3p     180363
## ssc-miR-542-5p      46371
## ssc-miR-545-3p      11069
## ssc-miR-545-5p       8576
## ssc-miR-551a           83
## ssc-miR-574         48642
## ssc-miR-582          1080
## ssc-miR-615           750
## ssc-miR-628          8250
## ssc-miR-652           339
## ssc-miR-664-3p       8055
## ssc-miR-664-5p      10481
## ssc-miR-671-3p        657
## ssc-miR-671-5p       3453
## ssc-miR-676-3p      24683
## ssc-miR-676-5p         17
## ssc-miR-6782-3p     15583
## ssc-miR-7           84265
## ssc-miR-708-3p       3750
## ssc-miR-708-5p      45309
## ssc-miR-7134-3p    745698
## ssc-miR-7134-5p     34153
## ssc-miR-7135-3p       950
## ssc-miR-7135-5p       239
## ssc-miR-7136-3p         8
## ssc-miR-7136-5p        36
## ssc-miR-7137-3p     11260
## ssc-miR-7137-5p       905
## ssc-miR-7138-3p        49
## ssc-miR-7138-5p       227
## ssc-miR-7139-3p        51
## ssc-miR-7139-5p       177
## ssc-miR-7140-3p        45
## ssc-miR-7140-5p         5
## ssc-miR-7142-3p     42825
## ssc-miR-7142-5p        38
## ssc-miR-7143-3p        29
## ssc-miR-7143-5p         2
## ssc-miR-7144-5p       135
## ssc-miR-744          5636
## ssc-miR-758         13448
## ssc-miR-769-3p       2129
## ssc-miR-769-5p      29988
## ssc-miR-7857-3p       139
## ssc-miR-874          5615
## ssc-miR-885-3p      18395
## ssc-miR-885-5p      98324
## ssc-miR-9            2913
## ssc-miR-9-2          1590
## ssc-miR-92a        492183
## ssc-miR-92b-3p       3191
## ssc-miR-92b-5p         28
## ssc-miR-935          7544
## ssc-miR-95         339465
## ssc-miR-96-5p         279
## ssc-miR-9785-5p       337
## ssc-miR-9788-3p         7
## ssc-miR-98         396367
## ssc-miR-9810-3p       397
## ssc-miR-9820-5p        96
## ssc-miR-9841-3p      2232
## ssc-miR-9843-3p     33445
## ssc-miR-9851-3p       925
## ssc-miR-9858-5p        51
## ssc-miR-9860-5p       908
## ssc-miR-99a        955356
## ssc-miR-99b        254680
```

```r
if (sum(rownames(totalsumbiotype) != uniqsumblast$dbseq_id)) stop ("miRNA names not equal")
```

As a check, manually sum the ssc-miR-1 and the ssc-miR-339 counts:


```r
sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-1"])
```

```
## [1] 47655223
```

```r
sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-339"])
```

```
## [1] 53376
```

```r
if (sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-1"]) != totalsumbiotype["ssc-miR-1",]) stop ("ssc-miR-1 counts not equal")
if (sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-339"]) != totalsumbiotype["ssc-miR-339",]) stop ("ssc-miR-339 counts not equal")
```

## Save data


```r
save(sumblast, uniqsumblast, totalsumbiotype, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/0a_susscrofa_filtered_miRBase_blastn_results.Rdata"))
write.csv(uniqsumblast, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/0b_susscrofa_filtered_uniqseq_miRBase_blastn_results.csv")
write.csv(totalsumbiotype, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/0c_susscrofa_filtered_totseq_miRBase_blastn_results.csv", row.names=TRUE)
```

