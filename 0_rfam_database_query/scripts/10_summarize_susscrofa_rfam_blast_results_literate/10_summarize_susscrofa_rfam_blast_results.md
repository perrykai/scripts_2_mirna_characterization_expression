**Script:** `10_summarize_susscrofa_rfam_blast_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`

**Date:**  7/12/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `2_susscrofa_rfam_blastn_output_e5.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** 

1. `2a_susscrofa_filtered_rfam_blastn_results.Rdata`
2. `2b_susscrofa_filtered_uniqseq_rfam_blastn_results.csv`
3. `2c_susscrofa_filtered_totseq_rfam_blastn_results.csv`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
This script filters and summarizes the blast results from the Sus scrofa Rfam database query.

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
ssc<-read.table("../2_susscrofa_rfam_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   0.333   0.010   0.425
```

```r
dim(ssc)
```

```
## [1] 33097    12
```

```r
head(ssc)
```

```
##               query_id                                   dbseq_id
## 1 seq_165791065_x23751   RF00001;5S_rRNA;CU928929.5/144654-144538
## 2 seq_165791065_x23751   RF00001;5S_rRNA;AEMK01132601.1/1867-1749
## 3 seq_165791065_x23751   RF00001;5S_rRNA;CR956386.6/192019-192106
## 4 seq_165791065_x23751     RF00001;5S_rRNA;FP565628.6/17108-17222
## 5 seq_165791065_x23751 RF00001;5S_rRNA;AEMK01149424.1/14650-14750
## 6 seq_165791065_x23751  RF00001;5S_rRNA;CT009542.22/163378-163269
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     19        0       0           1        19          71
## 2            100     18        0       0           2        19          74
## 3            100     18        0       0           1        18          63
## 4            100     18        0       0           1        18          73
## 5            100     18        0       0           1        18          64
## 6            100     18        0       0           1        18          74
##   dbseq_end evalue bitscore
## 1        89  2e-06     38.2
## 2        91  8e-06     36.2
## 3        80  8e-06     36.2
## 4        90  8e-06     36.2
## 5        81  8e-06     36.2
## 6        91  8e-06     36.2
```

```r
str(ssc)
```

```
## 'data.frame':	33097 obs. of  12 variables:
##  $ query_id      : Factor w/ 6552 levels "seq_165791065_x23751",..: 1 1 1 1 1 1 2 3 3 4 ...
##  $ dbseq_id      : Factor w/ 286 levels "RF00001;5S_rRNA;AEMK01017908.1/4690-4578",..: 21 7 11 24 9 13 111 7 21 280 ...
##  $ perc_identical: num  100 100 100 100 100 ...
##  $ length        : int  19 18 18 18 18 18 19 18 18 22 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 1 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 2 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  19 19 18 18 18 18 19 18 18 22 ...
##  $ dbseq_start   : int  71 74 63 73 64 74 47 74 72 26 ...
##  $ dbseq_end     : int  89 91 80 90 81 91 65 91 89 47 ...
##  $ evalue        : num  2e-06 8e-06 8e-06 8e-06 8e-06 8e-06 2e-06 7e-06 7e-06 1e-05 ...
##  $ bitscore      : num  38.2 36.2 36.2 36.2 36.2 36.2 38.2 36.2 36.2 36.2 ...
```

## Analysis


```r
ssc$dbseq_id<-as.character(ssc$dbseq_id)
ssc$query_id<-as.character(ssc$query_id)
```

Summarize the Rfam database hits by separating the components of the "dbseq_id" column at the semicolons.


```r
ssc$rfam_accession <- as.character(lapply(strsplit(as.character(ssc$dbseq_id), ';', fixed = TRUE), "[", 1))
ssc$rfam_biotype <- as.character(lapply(strsplit(as.character(ssc$dbseq_id), ';', fixed = TRUE), "[", 2))
ssc$rfam_seqnamestartend <- as.character(lapply(strsplit(as.character(ssc$dbseq_id), ';', fixed = TRUE), "[", 3))
head(ssc)
```

```
##               query_id                                   dbseq_id
## 1 seq_165791065_x23751   RF00001;5S_rRNA;CU928929.5/144654-144538
## 2 seq_165791065_x23751   RF00001;5S_rRNA;AEMK01132601.1/1867-1749
## 3 seq_165791065_x23751   RF00001;5S_rRNA;CR956386.6/192019-192106
## 4 seq_165791065_x23751     RF00001;5S_rRNA;FP565628.6/17108-17222
## 5 seq_165791065_x23751 RF00001;5S_rRNA;AEMK01149424.1/14650-14750
## 6 seq_165791065_x23751  RF00001;5S_rRNA;CT009542.22/163378-163269
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     19        0       0           1        19          71
## 2            100     18        0       0           2        19          74
## 3            100     18        0       0           1        18          63
## 4            100     18        0       0           1        18          73
## 5            100     18        0       0           1        18          64
## 6            100     18        0       0           1        18          74
##   dbseq_end evalue bitscore rfam_accession rfam_biotype
## 1        89  2e-06     38.2        RF00001      5S_rRNA
## 2        91  8e-06     36.2        RF00001      5S_rRNA
## 3        80  8e-06     36.2        RF00001      5S_rRNA
## 4        90  8e-06     36.2        RF00001      5S_rRNA
## 5        81  8e-06     36.2        RF00001      5S_rRNA
## 6        91  8e-06     36.2        RF00001      5S_rRNA
##         rfam_seqnamestartend
## 1   CU928929.5/144654-144538
## 2   AEMK01132601.1/1867-1749
## 3   CR956386.6/192019-192106
## 4     FP565628.6/17108-17222
## 5 AEMK01149424.1/14650-14750
## 6  CT009542.22/163378-163269
```

The number of unique sequences with susscrofa Rfam hits in this dataset


```r
length(unique(ssc$query_id))
```

```
## [1] 6552
```

Create a subset data frame containing the pertinent information for filtering


```r
rfamsubset<-ssc[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
```

```
## [1] 33097     6
```

```r
head(rfamsubset)
```

```
##               query_id perc_identical evalue bitscore rfam_accession
## 1 seq_165791065_x23751            100  2e-06     38.2        RF00001
## 2 seq_165791065_x23751            100  8e-06     36.2        RF00001
## 3 seq_165791065_x23751            100  8e-06     36.2        RF00001
## 4 seq_165791065_x23751            100  8e-06     36.2        RF00001
## 5 seq_165791065_x23751            100  8e-06     36.2        RF00001
## 6 seq_165791065_x23751            100  8e-06     36.2        RF00001
##   rfam_biotype
## 1      5S_rRNA
## 2      5S_rRNA
## 3      5S_rRNA
## 4      5S_rRNA
## 5      5S_rRNA
## 6      5S_rRNA
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
##  39.509   0.334   8.181
```

```r
length(idx)
```

```
## [1] 6552
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
##   1.004   0.192   0.269
```

```r
length(stp1)
```

```
## [1] 6552
```

```r
head(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 1 seq_165791065_x23751            100  2e-06     38.2        RF00001
##   rfam_biotype
## 1      5S_rRNA
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 7 seq_171683787_x17170            100  2e-06     38.2        RF00005
##   rfam_biotype
## 7         tRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 8 seq_171735276_x17125            100  7e-06     36.2        RF00001
## 9 seq_171735276_x17125            100  7e-06     36.2        RF00001
##   rfam_biotype
## 8      5S_rRNA
## 9      5S_rRNA
## 
## [[4]]
##                query_id perc_identical evalue bitscore rfam_accession
## 10 seq_174929842_x14404          95.45  1e-05     36.2        RF01897
##    rfam_biotype
## 10      mir-188
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 11 seq_181717095_x9462            100  7e-06     36.2        RF00012
##    rfam_biotype
## 11           U3
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 12 seq_181886722_x9388          95.65  3e-06     38.2        RF00001
##    rfam_biotype
## 12      5S_rRNA
```

```r
tail(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 33087 seq_227422786_x2            100  3e-06     38.2        RF00001
## 33088 seq_227422786_x2            100  3e-06     38.2        RF00001
## 33089 seq_227422786_x2            100  3e-06     38.2        RF00001
## 33090 seq_227422786_x2            100  3e-06     38.2        RF00001
## 33091 seq_227422786_x2            100  3e-06     38.2        RF00001
## 33092 seq_227422786_x2            100  3e-06     38.2        RF00001
##       rfam_biotype
## 33087      5S_rRNA
## 33088      5S_rRNA
## 33089      5S_rRNA
## 33090      5S_rRNA
## 33091      5S_rRNA
## 33092      5S_rRNA
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 33093 seq_227423202_x2            100  5e-08     44.1        RF00005
##       rfam_biotype
## 33093         tRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 33094 seq_227423340_x2          95.45  1e-05     36.2        RF01945
##       rfam_biotype
## 33094     mir-1388
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 33095 seq_227423736_x2            100  3e-06     38.2        RF00005
##       rfam_biotype
## 33095         tRNA
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 33096 seq_227424446_x2            100  7e-06     36.2        RF00069
##       rfam_biotype
## 33096      SNORD24
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 33097 seq_227425010_x2            100  1e-06     40.1        RF00005
##       rfam_biotype
## 33097         tRNA
```

How many sequences have more than one annotation after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 1081
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
## 583 185  53  31  22  21  14   6   4   2   1  12  14  15   4   2   2   3 
##  26  27  28  29  30  31  32  34  35 
##   4  48  10   9   6  19   4   6   1
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
##   0.439   0.190   0.161
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
## 1  seq_165791065_x23751         100.00  2e-06     38.2        RF00001
## 7  seq_171683787_x17170         100.00  2e-06     38.2        RF00005
## 8  seq_171735276_x17125         100.00  7e-06     36.2        RF00001
## 10 seq_174929842_x14404          95.45  1e-05     36.2        RF01897
## 11  seq_181717095_x9462         100.00  7e-06     36.2        RF00012
## 12  seq_181886722_x9388          95.65  3e-06     38.2        RF00001
##    rfam_biotype
## 1       5S_rRNA
## 7          tRNA
## 8       5S_rRNA
## 10      mir-188
## 11           U3
## 12      5S_rRNA
```

```r
dim(sumblast3)
```

```
## [1] 6552    6
```

```r
length(unique(sumblast3$query_id))
```

```
## [1] 6552
```

```r
table(sumblast3$rfam_biotype)
```

```
## 
##         5S_rRNA             7SK           ACA64         DLEU1_2 
##             440             267              27             152 
##        Histone3        HOTTIP_3        IRES_Bip     Metazoa_SRP 
##              85               1              10             138 
##        mir-1388         mir-154         mir-188         mir-223 
##              11               2             174              84 
##         mir-302         mir-322         mir-324         mir-544 
##               5               1               7               1 
##         mir-611         mir-676         mir-692       mtDNA_ssA 
##              14               2               9             357 
##         NEAT1_2         SCARNA2         SCARNA7         SECIS_1 
##               3               1             112               5 
##          snoR38         SNORA11         SNORA17         SNORA18 
##              46               1               2               6 
##         SNORA20         SNORA21         SNORA27         SNORA31 
##               1              32               1              16 
##         SNORA33         SNORA36         SNORA38         SNORA43 
##               3               1               2               3 
##         SNORA48         SNORA58         SNORA64         SNORA66 
##               2               5               5               3 
##         SNORA70         SNORA71         SNORA75         SNORA79 
##               1               3               7               1 
##        SNORD113         SNORD12       SNORD121A         SNORD14 
##              37              34              12              17 
##         SNORD22         SNORD23         SNORD24         SNORD26 
##              12               4             172              25 
##         SNORD27         SNORD33         SNORD36         SNORD39 
##              18              32             178               9 
##         SNORD42         SNORD62         SNORD65         SNORD82 
##              29              13              44              39 
##         SNORD83         SNORD90    snosnR60_Z15          snoU13 
##              17              43              12               1 
##          snoU54          snoZ17 Telomerase-vert            tRNA 
##               8              21               8            3216 
##              U1              U2              U3              U5 
##              15              99             296               3 
##              U6           Vault            XIST     ZNFX1-AS1_1 
##              83               4               1               1
```

```r
uniqsumblast3<-as.data.frame(table(sumblast3$rfam_biotype))
colnames(uniqsumblast3)<-c("Gene_Biotype", "Freq")
uniqsumblast3
```

```
##       Gene_Biotype Freq
## 1          5S_rRNA  440
## 2              7SK  267
## 3            ACA64   27
## 4          DLEU1_2  152
## 5         Histone3   85
## 6         HOTTIP_3    1
## 7         IRES_Bip   10
## 8      Metazoa_SRP  138
## 9         mir-1388   11
## 10         mir-154    2
## 11         mir-188  174
## 12         mir-223   84
## 13         mir-302    5
## 14         mir-322    1
## 15         mir-324    7
## 16         mir-544    1
## 17         mir-611   14
## 18         mir-676    2
## 19         mir-692    9
## 20       mtDNA_ssA  357
## 21         NEAT1_2    3
## 22         SCARNA2    1
## 23         SCARNA7  112
## 24         SECIS_1    5
## 25          snoR38   46
## 26         SNORA11    1
## 27         SNORA17    2
## 28         SNORA18    6
## 29         SNORA20    1
## 30         SNORA21   32
## 31         SNORA27    1
## 32         SNORA31   16
## 33         SNORA33    3
## 34         SNORA36    1
## 35         SNORA38    2
## 36         SNORA43    3
## 37         SNORA48    2
## 38         SNORA58    5
## 39         SNORA64    5
## 40         SNORA66    3
## 41         SNORA70    1
## 42         SNORA71    3
## 43         SNORA75    7
## 44         SNORA79    1
## 45        SNORD113   37
## 46         SNORD12   34
## 47       SNORD121A   12
## 48         SNORD14   17
## 49         SNORD22   12
## 50         SNORD23    4
## 51         SNORD24  172
## 52         SNORD26   25
## 53         SNORD27   18
## 54         SNORD33   32
## 55         SNORD36  178
## 56         SNORD39    9
## 57         SNORD42   29
## 58         SNORD62   13
## 59         SNORD65   44
## 60         SNORD82   39
## 61         SNORD83   17
## 62         SNORD90   43
## 63    snosnR60_Z15   12
## 64          snoU13    1
## 65          snoU54    8
## 66          snoZ17   21
## 67 Telomerase-vert    8
## 68            tRNA 3216
## 69              U1   15
## 70              U2   99
## 71              U3  296
## 72              U5    3
## 73              U6   83
## 74           Vault    4
## 75            XIST    1
## 76     ZNFX1-AS1_1    1
```

Add the column of sequence count to the sumblast data frame


```r
sumblast3$seq_count<-as.numeric(str_split_fixed(sumblast3$query_id, "_x", 2)[,2])
head(sumblast3)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1  seq_165791065_x23751         100.00  2e-06     38.2        RF00001
## 7  seq_171683787_x17170         100.00  2e-06     38.2        RF00005
## 8  seq_171735276_x17125         100.00  7e-06     36.2        RF00001
## 10 seq_174929842_x14404          95.45  1e-05     36.2        RF01897
## 11  seq_181717095_x9462         100.00  7e-06     36.2        RF00012
## 12  seq_181886722_x9388          95.65  3e-06     38.2        RF00001
##    rfam_biotype seq_count
## 1       5S_rRNA     23751
## 7          tRNA     17170
## 8       5S_rRNA     17125
## 10      mir-188     14404
## 11           U3      9462
## 12      5S_rRNA      9388
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype3<-as.matrix(by(sumblast3$seq_count, sumblast3$rfam_biotype, sum))
totalsumbiotype3
```

```
##                   [,1]
## 5S_rRNA          85508
## 7SK               1699
## ACA64              250
## DLEU1_2            669
## Histone3           595
## HOTTIP_3             2
## IRES_Bip            24
## Metazoa_SRP       2472
## mir-1388            56
## mir-154             11
## mir-188          28598
## mir-223           1722
## mir-302             23
## mir-322              3
## mir-324             40
## mir-544              2
## mir-611            147
## mir-676              6
## mir-692             18
## mtDNA_ssA         3246
## NEAT1_2              6
## SCARNA2              2
## SCARNA7           1136
## SECIS_1             12
## snoR38             554
## SNORA11              3
## SNORA17              9
## SNORA18             32
## SNORA20              2
## SNORA21            579
## SNORA27              2
## SNORA31            151
## SNORA33             11
## SNORA36              2
## SNORA38              4
## SNORA43            106
## SNORA48              5
## SNORA58             65
## SNORA64             11
## SNORA66             10
## SNORA70              2
## SNORA71              9
## SNORA75             26
## SNORA79              6
## SNORD113           521
## SNORD12           1410
## SNORD121A          490
## SNORD14           1202
## SNORD22            324
## SNORD23             51
## SNORD24           5050
## SNORD26            287
## SNORD27            350
## SNORD33            397
## SNORD36          17784
## SNORD39            220
## SNORD42            470
## SNORD62             57
## SNORD65           1789
## SNORD82           3128
## SNORD83            972
## SNORD90           3073
## snosnR60_Z15       235
## snoU13               2
## snoU54             194
## snoZ17             459
## Telomerase-vert     69
## tRNA            157921
## U1                  82
## U2                 903
## U3               18941
## U5                  12
## U6                1215
## Vault               25
## XIST                 2
## ZNFX1-AS1_1          2
```

```r
if (sum(rownames(totalsumbiotype3) != uniqsumblast3$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the 5s_rRNAs and the tRNA fields:


```r
sum(sumblast3$seq_count[sumblast3$rfam_biotype == "5S_rRNA"])
```

```
## [1] 85508
```

```r
sum(sumblast3$seq_count[sumblast3$rfam_biotype == "tRNA"])
```

```
## [1] 157921
```

```r
if (sum(sumblast3$seq_count[sumblast3$rfam_biotype == "5S_rRNA"]) != totalsumbiotype3["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast3$seq_count[sumblast3$rfam_biotype == "tRNA"]) != totalsumbiotype3["tRNA",]) stop ("tRNA counts not equal")
```

## Save data


```r
save(sumblast3, uniqsumblast3, totalsumbiotype3, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/2a_susscrofa_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast3, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/2b_susscrofa_filtered_uniqseq_rfam_blastn_results.csv")
write.csv(totalsumbiotype3, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/2c_susscrofa_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)
```

