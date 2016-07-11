**Script:** `14_summarize_blast_rfam_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`

**Date:**  6/14/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `2_susscrofa_rfam_blastn_output_e6.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** `6_susscrofa_filtered_rfam_blastn_results.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to filter and summarize the output from the second of the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
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
ssc<-read.table("../2_susscrofa_rfam_blastn_output_e6.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   0.224   0.003   0.277
```

```r
dim(ssc)
```

```
## [1] 27908    12
```

```r
head(ssc)
```

```
##                query_id                                   dbseq_id
## 1 seq_128408805_x137170   RF00737;mir-322;CU468475.7/117556-117473
## 2  seq_160647653_x30109   RF00001;5S_rRNA;CU928929.5/144654-144538
## 3  seq_161536921_x29165   RF00001;5S_rRNA;CU928929.5/144654-144538
## 4  seq_161536921_x29165   RF00001;5S_rRNA;CR956386.6/192019-192106
## 5  seq_161536921_x29165     RF00001;5S_rRNA;FP565628.6/17108-17222
## 6  seq_161536921_x29165 RF00001;5S_rRNA;AEMK01149424.1/14650-14750
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     21        0       0           1        21          11
## 2            100     20        0       0           1        20          70
## 3            100     21        0       0           1        21          69
## 4            100     20        0       0           1        20          61
## 5            100     20        0       0           1        20          71
## 6            100     20        0       0           1        20          62
##   dbseq_end evalue bitscore
## 1        31  1e-07     42.1
## 2        89  6e-07     40.1
## 3        89  1e-07     42.1
## 4        80  6e-07     40.1
## 5        90  6e-07     40.1
## 6        81  6e-07     40.1
```

```r
str(ssc)
```

```
## 'data.frame':	27908 obs. of  12 variables:
##  $ query_id      : Factor w/ 5854 levels "seq_128408805_x137170",..: 1 2 3 3 3 3 3 4 5 6 ...
##  $ dbseq_id      : Factor w/ 257 levels "RF00001;5S_rRNA;AEMK01035741.1/836-715",..: 236 19 19 9 21 7 11 235 236 235 ...
##  $ perc_identical: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ length        : int  21 20 21 20 20 20 20 21 22 20 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  21 20 21 20 20 20 20 21 22 20 ...
##  $ dbseq_start   : int  11 70 69 61 71 62 72 18 11 18 ...
##  $ dbseq_end     : int  31 89 89 80 90 81 91 38 32 37 ...
##  $ evalue        : num  1e-07 6e-07 1e-07 6e-07 6e-07 6e-07 6e-07 1e-07 4e-08 6e-07 ...
##  $ bitscore      : num  42.1 40.1 42.1 40.1 40.1 40.1 40.1 42.1 44.1 40.1 ...
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
##                query_id                                   dbseq_id
## 1 seq_128408805_x137170   RF00737;mir-322;CU468475.7/117556-117473
## 2  seq_160647653_x30109   RF00001;5S_rRNA;CU928929.5/144654-144538
## 3  seq_161536921_x29165   RF00001;5S_rRNA;CU928929.5/144654-144538
## 4  seq_161536921_x29165   RF00001;5S_rRNA;CR956386.6/192019-192106
## 5  seq_161536921_x29165     RF00001;5S_rRNA;FP565628.6/17108-17222
## 6  seq_161536921_x29165 RF00001;5S_rRNA;AEMK01149424.1/14650-14750
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     21        0       0           1        21          11
## 2            100     20        0       0           1        20          70
## 3            100     21        0       0           1        21          69
## 4            100     20        0       0           1        20          61
## 5            100     20        0       0           1        20          71
## 6            100     20        0       0           1        20          62
##   dbseq_end evalue bitscore rfam_accession rfam_biotype
## 1        31  1e-07     42.1        RF00737      mir-322
## 2        89  6e-07     40.1        RF00001      5S_rRNA
## 3        89  1e-07     42.1        RF00001      5S_rRNA
## 4        80  6e-07     40.1        RF00001      5S_rRNA
## 5        90  6e-07     40.1        RF00001      5S_rRNA
## 6        81  6e-07     40.1        RF00001      5S_rRNA
##         rfam_seqnamestartend
## 1   CU468475.7/117556-117473
## 2   CU928929.5/144654-144538
## 3   CU928929.5/144654-144538
## 4   CR956386.6/192019-192106
## 5     FP565628.6/17108-17222
## 6 AEMK01149424.1/14650-14750
```

The number of unique sequences with susscrofa Rfam hits in this dataset


```r
length(unique(ssc$query_id))
```

```
## [1] 5854
```

Create a subset data frame containing the pertinent information for filtering


```r
rfamsubset<-ssc[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
```

```
## [1] 27908     6
```

```r
head(rfamsubset)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_128408805_x137170            100  1e-07     42.1        RF00737
## 2  seq_160647653_x30109            100  6e-07     40.1        RF00001
## 3  seq_161536921_x29165            100  1e-07     42.1        RF00001
## 4  seq_161536921_x29165            100  6e-07     40.1        RF00001
## 5  seq_161536921_x29165            100  6e-07     40.1        RF00001
## 6  seq_161536921_x29165            100  6e-07     40.1        RF00001
##   rfam_biotype
## 1      mir-322
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
##  29.834   0.236   3.906
```

```r
length(idx)
```

```
## [1] 5854
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
##   0.851   0.147   0.175
```

```r
length(stp1)
```

```
## [1] 5854
```

```r
head(stp1)
```

```
## [[1]]
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_128408805_x137170            100  1e-07     42.1        RF00737
##   rfam_biotype
## 1      mir-322
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 2 seq_160647653_x30109            100  6e-07     40.1        RF00001
##   rfam_biotype
## 2      5S_rRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 3 seq_161536921_x29165            100  1e-07     42.1        RF00001
##   rfam_biotype
## 3      5S_rRNA
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 8 seq_164820205_x24792            100  1e-07     42.1        RF00708
##   rfam_biotype
## 8      mir-450
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 9 seq_167301252_x22072            100  4e-08     44.1        RF00737
##   rfam_biotype
## 9      mir-322
## 
## [[6]]
##                query_id perc_identical evalue bitscore rfam_accession
## 10 seq_171050786_x17914            100  6e-07     40.1        RF00708
##    rfam_biotype
## 10      mir-450
```

```r
tail(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27899 seq_227422044_x2            100  6e-07     40.1        RF00966
##       rfam_biotype
## 27899      mir-676
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27900 seq_227422048_x2            100  2e-07     42.1        RF00005
##       rfam_biotype
## 27900         tRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27901 seq_227422156_x2            100  2e-07     42.1        RF00001
##       rfam_biotype
## 27901      5S_rRNA
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27906 seq_227422686_x2            100  1e-07     42.1        RF00005
##       rfam_biotype
## 27906         tRNA
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27907 seq_227423202_x2            100  5e-08     44.1        RF00005
##       rfam_biotype
## 27907         tRNA
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 27908 seq_227424962_x2            100  6e-07     40.1        RF00086
##       rfam_biotype
## 27908      SNORD27
```

How many sequences have more than one annotation after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 895
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##   2   3   4   5   6   7   8   9  10  11  13  14  15  16  26  27  28  29 
## 498 148  40  22  16   1  15   9   4   5  14  30  15   5  14  46   1   4 
##  30  31 
##   3   5
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
##   0.410   0.174   0.107
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
sumblast2 <- do.call(rbind,stp3)
head(sumblast2)
```

```
##                 query_id perc_identical evalue bitscore rfam_accession
## 1  seq_128408805_x137170            100  1e-07     42.1        RF00737
## 2   seq_160647653_x30109            100  6e-07     40.1        RF00001
## 3   seq_161536921_x29165            100  1e-07     42.1        RF00001
## 8   seq_164820205_x24792            100  1e-07     42.1        RF00708
## 9   seq_167301252_x22072            100  4e-08     44.1        RF00737
## 10  seq_171050786_x17914            100  6e-07     40.1        RF00708
##    rfam_biotype
## 1       mir-322
## 2       5S_rRNA
## 3       5S_rRNA
## 8       mir-450
## 9       mir-322
## 10      mir-450
```

```r
dim(sumblast2)
```

```
## [1] 5854    6
```

```r
length(unique(sumblast2$query_id))
```

```
## [1] 5854
```

```r
table(sumblast2$rfam_biotype)
```

```
## 
##         5S_rRNA             7SK           ACA64         DLEU1_2 
##             461             223              23             139 
##        Histone3        IRES_Bip     Metazoa_SRP        mir-1388 
##              76               6             110              16 
##         mir-188         mir-223         mir-302         mir-322 
##              86              61               9             197 
##         mir-324         mir-450         mir-500         mir-544 
##               6              75              41               1 
##         mir-652         mir-676         mir-692       mtDNA_ssA 
##               8              46               5             260 
##         NEAT1_2         SCARNA7         SECIS_1          snoR38 
##               2              86               3              58 
##         SNORA17         SNORA18         SNORA21         SNORA25 
##               4               6              32               3 
##         SNORA31         SNORA33         SNORA36         SNORA38 
##              32              12               1               3 
##         SNORA43         SNORA48         SNORA58         SNORA62 
##               7               3               4               2 
##         SNORA63         SNORA64         SNORA66         SNORA71 
##               1               3               3               6 
##         SNORA75        SNORD107        SNORD113         SNORD12 
##               5               2              43              22 
##       SNORD121A         SNORD14         SNORD22         SNORD23 
##               6               7               6               2 
##         SNORD24         SNORD26         SNORD27         SNORD33 
##             136              20              22              28 
##         SNORD36         SNORD39         SNORD42         SNORD62 
##             109               5              22               7 
##         SNORD65         SNORD82         SNORD83         SNORD90 
##              62              32              12              29 
##    snosnR60_Z15          snoU54          snoZ17 Telomerase-vert 
##              22               7              16              10 
##            tRNA              U1             U11              U2 
##            2738               7               2              82 
##              U3              U5              U6          U6atac 
##             206               1              61               2 
##              U8            XIST 
##               2               1
```

```r
uniqsumblast2<-as.data.frame(table(sumblast2$rfam_biotype))
colnames(uniqsumblast2)<-c("Gene_Biotype", "Freq")
uniqsumblast2
```

```
##       Gene_Biotype Freq
## 1          5S_rRNA  461
## 2              7SK  223
## 3            ACA64   23
## 4          DLEU1_2  139
## 5         Histone3   76
## 6         IRES_Bip    6
## 7      Metazoa_SRP  110
## 8         mir-1388   16
## 9          mir-188   86
## 10         mir-223   61
## 11         mir-302    9
## 12         mir-322  197
## 13         mir-324    6
## 14         mir-450   75
## 15         mir-500   41
## 16         mir-544    1
## 17         mir-652    8
## 18         mir-676   46
## 19         mir-692    5
## 20       mtDNA_ssA  260
## 21         NEAT1_2    2
## 22         SCARNA7   86
## 23         SECIS_1    3
## 24          snoR38   58
## 25         SNORA17    4
## 26         SNORA18    6
## 27         SNORA21   32
## 28         SNORA25    3
## 29         SNORA31   32
## 30         SNORA33   12
## 31         SNORA36    1
## 32         SNORA38    3
## 33         SNORA43    7
## 34         SNORA48    3
## 35         SNORA58    4
## 36         SNORA62    2
## 37         SNORA63    1
## 38         SNORA64    3
## 39         SNORA66    3
## 40         SNORA71    6
## 41         SNORA75    5
## 42        SNORD107    2
## 43        SNORD113   43
## 44         SNORD12   22
## 45       SNORD121A    6
## 46         SNORD14    7
## 47         SNORD22    6
## 48         SNORD23    2
## 49         SNORD24  136
## 50         SNORD26   20
## 51         SNORD27   22
## 52         SNORD33   28
## 53         SNORD36  109
## 54         SNORD39    5
## 55         SNORD42   22
## 56         SNORD62    7
## 57         SNORD65   62
## 58         SNORD82   32
## 59         SNORD83   12
## 60         SNORD90   29
## 61    snosnR60_Z15   22
## 62          snoU54    7
## 63          snoZ17   16
## 64 Telomerase-vert   10
## 65            tRNA 2738
## 66              U1    7
## 67             U11    2
## 68              U2   82
## 69              U3  206
## 70              U5    1
## 71              U6   61
## 72          U6atac    2
## 73              U8    2
## 74            XIST    1
```

Add the column of sequence count to the sumblast data frame


```r
sumblast2$seq_count<-as.numeric(str_split_fixed(sumblast2$query_id, "_x", 2)[,2])
head(sumblast2)
```

```
##                 query_id perc_identical evalue bitscore rfam_accession
## 1  seq_128408805_x137170            100  1e-07     42.1        RF00737
## 2   seq_160647653_x30109            100  6e-07     40.1        RF00001
## 3   seq_161536921_x29165            100  1e-07     42.1        RF00001
## 8   seq_164820205_x24792            100  1e-07     42.1        RF00708
## 9   seq_167301252_x22072            100  4e-08     44.1        RF00737
## 10  seq_171050786_x17914            100  6e-07     40.1        RF00708
##    rfam_biotype seq_count
## 1       mir-322    137170
## 2       5S_rRNA     30109
## 3       5S_rRNA     29165
## 8       mir-450     24792
## 9       mir-322     22072
## 10      mir-450     17914
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype2<-as.matrix(by(sumblast2$seq_count, sumblast2$rfam_biotype, sum))
totalsumbiotype2
```

```
##                   [,1]
## 5S_rRNA         111476
## 7SK               1516
## ACA64              200
## DLEU1_2            555
## Histone3           559
## IRES_Bip            14
## Metazoa_SRP       5950
## mir-1388           378
## mir-188          21483
## mir-223           6088
## mir-302            536
## mir-322         202022
## mir-324             38
## mir-450          47744
## mir-500           5949
## mir-544              4
## mir-652             50
## mir-676          23574
## mir-692             10
## mtDNA_ssA         2003
## NEAT1_2              4
## SCARNA7            879
## SECIS_1              7
## snoR38             522
## SNORA17             15
## SNORA18             26
## SNORA21            349
## SNORA25              9
## SNORA31            795
## SNORA33             60
## SNORA36              7
## SNORA38              8
## SNORA43            240
## SNORA48              6
## SNORA58             55
## SNORA62             14
## SNORA63              2
## SNORA64              6
## SNORA66             10
## SNORA71             87
## SNORA75             15
## SNORD107             5
## SNORD113           414
## SNORD12            696
## SNORD121A           28
## SNORD14           2096
## SNORD22            307
## SNORD23             30
## SNORD24           3327
## SNORD26            152
## SNORD27            775
## SNORD33           1871
## SNORD36           9234
## SNORD39            273
## SNORD42            833
## SNORD62             48
## SNORD65          12750
## SNORD82            473
## SNORD83            711
## SNORD90           2682
## snosnR60_Z15       288
## snoU54              49
## snoZ17             388
## Telomerase-vert    152
## tRNA            110263
## U1                  36
## U11                  4
## U2                 604
## U3               11948
## U5                   2
## U6                 675
## U6atac               4
## U8                   4
## XIST                 2
```

```r
if (sum(rownames(totalsumbiotype2) != uniqsumblast2$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the 5s_rRNAs and the tRNA fields:


```r
sum(sumblast2$seq_count[sumblast2$rfam_biotype == "5S_rRNA"])
```

```
## [1] 111476
```

```r
sum(sumblast2$seq_count[sumblast2$rfam_biotype == "tRNA"])
```

```
## [1] 110263
```

```r
if (sum(sumblast2$seq_count[sumblast2$rfam_biotype == "5S_rRNA"]) != totalsumbiotype2["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast2$seq_count[sumblast2$rfam_biotype == "tRNA"]) != totalsumbiotype2["tRNA",]) stop ("tRNA counts not equal")
```

## Save data


```r
save(sumblast2, uniqsumblast2, totalsumbiotype2, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/6_susscrofa_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast2, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/6_1_susscrofa_filtered_uniqseq_rfam_blastn_results.csv", col.names=TRUE)
```

```
## Warning in write.csv(uniqsumblast2, file = "/mnt/research/pigeqtl/analyses/
## microRNA/2_mirna_characterization_expression/0_rfam_database_query/
## 6_1_susscrofa_filtered_uniqseq_rfam_blastn_results.csv", : attempt to set
## 'col.names' ignored
```

```r
write.csv(totalsumbiotype2, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/6_2_susscrofa_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)
```

