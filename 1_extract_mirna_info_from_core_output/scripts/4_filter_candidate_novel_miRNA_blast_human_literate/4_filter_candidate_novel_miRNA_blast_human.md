**Script:** `4_filter_candidate_novel_miRNA_blast_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts/`

**Date:**  `10/24/16`

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/`

**Input File(s):** 

1. `1_candidate_novel_miRNA_precursor_blastn_human_output_e5.txt`
2. `1_candidate_novel_miRNA_mature_blastn_human_output_e5.txt`
3. `2_candidate_novel_miRNA_precursor_blastn_human_output_eval1.txt`
4. `2_candidate_novel_miRNA_mature_blastn_human_output_eval1.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/`

**Output File(s):** 

1. `3_filtered_candidate_novel_miRNA_precursor_abundance_e5.txt`
2. `4_filtered_candidate_novel_precursor_eval1.txt`
3. `5_filtered_candidate_novel_mature_eval1.txt`
4. `6_filtered_candidate_novel_miRNA_precursor_ids_e5.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to summarize the output from the BLAST query in order to characterize the candidate novel miRNAs present in the small RNA seq data.
Additionally, I want to filter the miRNA BLAST results with the e-value = 1.0, to remove some of the redundant hits and make the results more managable. 

So, need to load the precursor BLAST dataset and the full dataset of candidate novel miRNA and compare the names of sequences to see if the most abundant miRNA had BLAST results.

## Install libraries


```r
rm(list=ls())
```

## Load data


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts")

hsa.blast.precursor<-read.table("../8_blastn_candidate_novel_miRNA_output/1_candidate_novel_miRNA_precursor_blastn_human_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
eval1.hsa.blast.precursor<-read.table("../8_blastn_candidate_novel_miRNA_output/2_candidate_novel_miRNA_precursor_blastn_human_output_eval1.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
eval1.hsa.blast.mature<-read.table("../8_blastn_candidate_novel_miRNA_output/2_candidate_novel_miRNA_mature_blastn_human_output_eval1.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))

dim(hsa.blast.precursor)
```

```
## [1] 41 12
```

```r
head(hsa.blast.precursor)
```

```
##                                          query_id      dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-mir-26b
## 2          seq|X_39130|candidatenovelpercursormiR   hsa-mir-660
## 3           seq|3_7507|candidatenovelpercursormiR   hsa-mir-590
## 4         seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-1
## 5         seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-2
## 6 seq|GL894231.1_44092|candidatenovelpercursormiR   hsa-mir-655
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1         100.00     49        0       0           1        49          12
## 2          97.67     43        1       0           1        43          16
## 3          96.00     25        1       0           1        25          16
## 4          96.08     51        2       0           1        51          10
## 5          95.83     24        1       0           1        24          14
## 6         100.00     38        0       0           1        38          23
##   dbseq_end evalue bitscore
## 1        60  2e-23     97.6
## 2        58  2e-17     77.8
## 3        40  1e-06     42.1
## 4        60  7e-20     85.7
## 5        37  4e-06     40.1
## 6        60  7e-17     75.8
```

```r
str(hsa.blast.precursor)
```

```
## 'data.frame':	41 obs. of  12 variables:
##  $ query_id      : Factor w/ 27 levels "seq|1_1168|candidatenovelpercursormiR",..: 15 24 12 6 6 19 17 16 3 21 ...
##  $ dbseq_id      : Factor w/ 37 levels "hsa-mir-122",..: 12 33 30 10 11 32 24 22 31 21 ...
##  $ perc_identical: num  100 97.7 96 96.1 95.8 ...
##  $ length        : int  49 43 25 51 24 38 52 54 49 48 ...
##  $ mismatch      : int  0 1 1 2 1 0 0 1 1 2 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  49 43 25 51 24 38 52 54 49 48 ...
##  $ dbseq_start   : int  12 16 16 10 14 23 9 7 12 13 ...
##  $ dbseq_end     : int  60 58 40 60 37 60 60 60 60 60 ...
##  $ evalue        : num  2e-23 2e-17 1e-06 7e-20 4e-06 ...
##  $ bitscore      : num  97.6 77.8 42.1 85.7 40.1 75.8 103 99.6 89.7 79.8 ...
```

```r
hsa.blast.precursor$dbseq_id<-as.character(hsa.blast.precursor$dbseq_id)
hsa.blast.precursor$query_id<-as.character(hsa.blast.precursor$query_id)

dim(eval1.hsa.blast.precursor)
```

```
## [1] 652  12
```

```r
head(eval1.hsa.blast.precursor)
```

```
##                                          query_id      dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-mir-26b
## 2 seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-26a-1
## 3 seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-26a-2
## 4 seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-4473
## 5 seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-1287
## 6 seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-6775
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     49        0       0           1        49          12
## 2             92     25        2       0           1        25          10
## 3             95     20        1       0           1        20          14
## 4            100     13        0       0          38        50          31
## 5            100     12        0       0          24        35           3
## 6            100     11        0       0          28        38          34
##   dbseq_end  evalue bitscore
## 1        60 2.0e-23     97.6
## 2        34 2.0e-04     34.2
## 3        33 8.0e-04     32.2
## 4        19 5.2e-02     26.3
## 5        14 2.1e-01     24.3
## 6        24 8.2e-01     22.3
```

```r
str(eval1.hsa.blast.precursor)
```

```
## 'data.frame':	652 obs. of  12 variables:
##  $ query_id      : Factor w/ 121 levels "seq|1_1036|candidatenovelpercursormiR",..: 92 92 92 92 92 92 92 92 92 92 ...
##  $ dbseq_id      : Factor w/ 455 levels "hsa-let-7b","hsa-mir-105-1",..: 72 70 71 186 31 368 148 7 48 1 ...
##  $ perc_identical: num  100 92 95 100 100 ...
##  $ length        : int  49 25 20 13 12 11 15 15 11 11 ...
##  $ mismatch      : int  0 2 1 0 0 0 1 1 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 38 24 28 23 35 27 17 ...
##  $ query_end     : int  49 25 20 50 35 38 37 49 37 27 ...
##  $ dbseq_start   : int  12 10 14 31 3 34 51 9 12 14 ...
##  $ dbseq_end     : int  60 34 33 19 14 24 37 23 2 24 ...
##  $ evalue        : num  2.0e-23 2.0e-04 8.0e-04 5.2e-02 2.1e-01 ...
##  $ bitscore      : num  97.6 34.2 32.2 26.3 24.3 22.3 22.3 22.3 22.3 22.3 ...
```

```r
eval1.hsa.blast.precursor$dbseq_id<-as.character(eval1.hsa.blast.precursor$dbseq_id)
eval1.hsa.blast.precursor$query_id<-as.character(eval1.hsa.blast.precursor$query_id)

dim(eval1.hsa.blast.mature)
```

```
## [1] 235  12
```

```r
head(eval1.hsa.blast.mature)
```

```
##                                       query_id        dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-26b-5p
## 2 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-26a-5p
## 3 seq|GL892871.2_43622|candidatenovelmaturemiR    hsa-miR-1297
## 4 seq|GL892871.2_43622|candidatenovelmaturemiR hsa-miR-3140-5p
## 5 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-876-3p
## 6          seq|X_39130|candidatenovelmaturemiR  hsa-miR-660-5p
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     21        0       0           1        21           1
## 2             95     20        1       0           1        20           1
## 3            100     15        0       0           1        15           1
## 4            100     10        0       0           6        15          11
## 5            100     10        0       0           4        13          13
## 6            100     22        0       0           1        22           1
##   dbseq_end  evalue bitscore
## 1        21 9.0e-08     42.1
## 2        20 9.0e-05     32.2
## 3        15 3.0e-04     30.2
## 4         2 3.2e-01     20.3
## 5        22 3.2e-01     20.3
## 6        22 2.0e-08     44.1
```

```r
str(eval1.hsa.blast.mature)
```

```
## 'data.frame':	235 obs. of  12 variables:
##  $ query_id      : Factor w/ 92 levels "seq|1_1036|candidatenovelmaturemiR",..: 67 67 67 67 67 84 36 36 13 13 ...
##  $ dbseq_id      : Factor w/ 205 levels "hsa-miR-1208",..: 34 33 10 40 199 141 125 126 33 34 ...
##  $ perc_identical: num  100 95 100 100 100 ...
##  $ length        : int  21 20 15 10 10 22 21 15 22 20 ...
##  $ mismatch      : int  0 1 0 0 0 0 0 1 1 2 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 6 4 1 1 4 1 1 ...
##  $ query_end     : int  21 20 15 15 13 22 21 18 22 20 ...
##  $ dbseq_start   : int  1 1 1 11 13 1 1 16 1 1 ...
##  $ dbseq_end     : int  21 20 15 2 22 22 21 2 22 20 ...
##  $ evalue        : num  9.0e-08 9.0e-05 3.0e-04 3.2e-01 3.2e-01 2.0e-08 8.0e-08 7.6e-02 5.0e-06 2.1e-02 ...
##  $ bitscore      : num  42.1 32.2 30.2 20.3 20.3 44.1 42.1 22.3 36.2 24.3 ...
```

```r
eval1.hsa.blast.mature$dbseq_id<-as.character(eval1.hsa.blast.mature$dbseq_id)
eval1.hsa.blast.mature$query_id<-as.character(eval1.hsa.blast.mature$query_id)

load("../5_putative_novel_miRNA_filtered_candidates.Rdata")

dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 132  17
```

```r
rownames(novelmir10sigrandfoldmincounts)<- novelmir10sigrandfoldmincounts$provisional.id
head(novelmir10sigrandfoldmincounts)
```

```
##                    provisional.id miRDeep2.score
## GL892871.2_43622 GL892871.2_43622       572746.3
## X_39130                   X_39130       136780.7
## GL896425.1_44884 GL896425.1_44884       128796.6
## 3_7507                     3_7507        57547.6
## 13_28851                 13_28851        40621.2
## GL894231.1_44092 GL894231.1_44092        12387.7
##                  estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## GL892871.2_43622                                                         91 +/- 2%
## X_39130                                                                  91 +/- 2%
## GL896425.1_44884                                                         91 +/- 2%
## 3_7507                                                                   91 +/- 2%
## 13_28851                                                                 91 +/- 2%
## GL894231.1_44092                                                         91 +/- 2%
##                  rfam.alert total.read.count mature.read.count
## GL892871.2_43622          -          1123410           1119829
## X_39130                   -           268282            268080
## GL896425.1_44884          -           252627            249026
## 3_7507                    -           112871            104821
## 13_28851                  -            79667             79633
## GL894231.1_44092          -            24290             24219
##                  loop.read.count star.read.count
## GL892871.2_43622               0            3581
## X_39130                        0             202
## GL896425.1_44884               0            3601
## 3_7507                         0            8050
## 13_28851                       0              34
## GL894231.1_44092               0              71
##                  significant.randfold.p.value miRBase.miRNA
## GL892871.2_43622                          yes             -
## X_39130                                   yes             -
## GL896425.1_44884                          yes             -
## 3_7507                                    yes             -
## 13_28851                                  yes             -
## GL894231.1_44092                          yes             -
##                  example.miRBase.miRNA.with.the.same.seed UCSC.browser
## GL892871.2_43622                           hsa-miR-26a-5p            -
## X_39130                                    hsa-miR-660-5p            -
## GL896425.1_44884                                        -            -
## 3_7507                                     hsa-miR-590-3p            -
## 13_28851                                   hsa-miR-26a-5p            -
## GL894231.1_44092                           hsa-miR-655-3p            -
##                  NCBI.blastn consensus.mature.sequence
## GL892871.2_43622           -    uucaaguaauucaggauagguu
## X_39130                    -    uacccauugcauaucggaguug
## GL896425.1_44884           -    uggugccugacgucuuggcagu
## 3_7507                     -     uaauuuuauguauaagcuagu
## 13_28851                   -    uucaaguaacccaggauaggcu
## GL894231.1_44092           -    auaauacaugguuaaccucuuu
##                  consensus.star.sequence
## GL892871.2_43622  ccuguucuccauuacuuggcuc
## X_39130           accuccuaugugcaugguuuac
## GL896425.1_44884 agccagggcugcaggcacugaca
## 3_7507            gagcuuauucauaaaaguacag
## 13_28851          ccuauucuugguuacuugcacg
## GL894231.1_44092  agagguuauccguguuauguuc
##                                                   consensus.precursor.sequence
## GL892871.2_43622     uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## X_39130             uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## GL896425.1_44884 agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 3_7507           gagcuuauucauaaaaguacaguauaauccaguaaaccuguaauuuuauguauaagcuagu
## 13_28851         uucaaguaacccaggauaggcugugcaggucccaaggggccuauucuugguuacuugcacg
## GL894231.1_44092  agagguuauccguguuauguucgcuucauucaucaugaauaauacaugguuaaccucuuu
##                       precursor.coordinate
## GL892871.2_43622 GL892871.2:56019..56076:+
## X_39130             X:48640826..48640884:+
## GL896425.1_44884   GL896425.1:1750..1811:-
## 3_7507              3:11024999..11025060:+
## 13_28851           13:24885255..24885316:-
## GL894231.1_44092 GL894231.1:23071..23131:-
```

## Analysis

### 1. The goal is to compare the candidate novel miRNAs with BLAST results to the most abundant candidate novel miRNAs.

First, return the name of the sequences to the way they were before BLASTing at e-value = 1x10^-5.


```r
hsa.blast.precursor$seqname<-sapply(strsplit(hsa.blast.precursor$query_id, "|", fixed=TRUE),'[',2)
head(hsa.blast.precursor)
```

```
##                                          query_id      dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-mir-26b
## 2          seq|X_39130|candidatenovelpercursormiR   hsa-mir-660
## 3           seq|3_7507|candidatenovelpercursormiR   hsa-mir-590
## 4         seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-1
## 5         seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-2
## 6 seq|GL894231.1_44092|candidatenovelpercursormiR   hsa-mir-655
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1         100.00     49        0       0           1        49          12
## 2          97.67     43        1       0           1        43          16
## 3          96.00     25        1       0           1        25          16
## 4          96.08     51        2       0           1        51          10
## 5          95.83     24        1       0           1        24          14
## 6         100.00     38        0       0           1        38          23
##   dbseq_end evalue bitscore          seqname
## 1        60  2e-23     97.6 GL892871.2_43622
## 2        58  2e-17     77.8          X_39130
## 3        40  1e-06     42.1           3_7507
## 4        60  7e-20     85.7         13_28851
## 5        37  4e-06     40.1         13_28851
## 6        60  7e-17     75.8 GL894231.1_44092
```

Identify the unique number of sequences with BLAST hits at eval = 1x10^-5


```r
seqids<-unique(hsa.blast.precursor$seqname)
length(seqids)
```

```
## [1] 27
```

```r
seqids
```

```
##  [1] "GL892871.2_43622" "X_39130"          "3_7507"          
##  [4] "13_28851"         "GL894231.1_44092" "GL894231.1_44082"
##  [7] "GL894231.1_44074" "12_25369"         "GL894231.1_44108"
## [10] "1_1168"           "X_39122"          "X_39256"         
## [13] "7_18912"          "13_28053"         "6_14565"         
## [16] "JH118928.1_41756" "GL894231.1_44094" "17_36762"        
## [19] "1_3654"           "14_31131"         "X_40048"         
## [22] "1_851"            "X_39264"          "12_26433"        
## [25] "GL894231.1_44088" "1_1860"           "2_5120"
```

Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset


```r
match(seqids, rownames(novelmir10sigrandfoldmincounts))
```

```
##  [1]   1   2   4   5   6  10  12  13  15  16  17  20  22  26  32  39  46
## [18]  49  55  56  67  69  73  74  79  92 122
```

```r
rownames(novelmir10sigrandfoldmincounts)%in%seqids
```

```
##   [1]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE FALSE
##  [12]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE FALSE  TRUE
##  [23] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
##  [34] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
##  [45] FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE
##  [56]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
##  [67]  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE
##  [78] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
##  [89] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [100] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [111] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [122]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```

Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.


```r
novelmir.abundance<-novelmir10sigrandfoldmincounts[seqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]
```

Is this object ordered by miRDeep2 score or by total.read.count?


```r
sum(rownames(novelmir.abundance[order(novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(novelmir.abundance))
```

```
## [1] 0
```

```r
sum(rownames(novelmir.abundance[order(novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(novelmir.abundance))
```

```
## [1] 0
```

Combine the information; are the predicted miRNA precursors with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?

This object (candidate novel precursors with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.


```r
novelmir.abundance<-cbind(novelmir.abundance[match(hsa.blast.precursor$seqname, rownames(novelmir.abundance)),], hsa.blast.precursor)
sum(novelmir.abundance$seqname != novelmir.abundance$provisional.id)
```

```
## [1] 0
```

```r
head(novelmir.abundance)
```

```
##                    provisional.id miRDeep2.score total.read.count
## GL892871.2_43622 GL892871.2_43622       572746.3          1123410
## X_39130                   X_39130       136780.7           268282
## 3_7507                     3_7507        57547.6           112871
## 13_28851                 13_28851        40621.2            79667
## 13_28851.1               13_28851        40621.2            79667
## GL894231.1_44092 GL894231.1_44092        12387.7            24290
##                  mature.read.count star.read.count
## GL892871.2_43622           1119829            3581
## X_39130                     268080             202
## 3_7507                      104821            8050
## 13_28851                     79633              34
## 13_28851.1                   79633              34
## GL894231.1_44092             24219              71
##                  example.miRBase.miRNA.with.the.same.seed
## GL892871.2_43622                           hsa-miR-26a-5p
## X_39130                                    hsa-miR-660-5p
## 3_7507                                     hsa-miR-590-3p
## 13_28851                                   hsa-miR-26a-5p
## 13_28851.1                                 hsa-miR-26a-5p
## GL894231.1_44092                           hsa-miR-655-3p
##                       precursor.coordinate
## GL892871.2_43622 GL892871.2:56019..56076:+
## X_39130             X:48640826..48640884:+
## 3_7507              3:11024999..11025060:+
## 13_28851           13:24885255..24885316:-
## 13_28851.1         13:24885255..24885316:-
## GL894231.1_44092 GL894231.1:23071..23131:-
##                                                         query_id
## GL892871.2_43622 seq|GL892871.2_43622|candidatenovelpercursormiR
## X_39130                   seq|X_39130|candidatenovelpercursormiR
## 3_7507                     seq|3_7507|candidatenovelpercursormiR
## 13_28851                 seq|13_28851|candidatenovelpercursormiR
## 13_28851.1               seq|13_28851|candidatenovelpercursormiR
## GL894231.1_44092 seq|GL894231.1_44092|candidatenovelpercursormiR
##                       dbseq_id perc_identical length mismatch gapopen
## GL892871.2_43622   hsa-mir-26b         100.00     49        0       0
## X_39130            hsa-mir-660          97.67     43        1       0
## 3_7507             hsa-mir-590          96.00     25        1       0
## 13_28851         hsa-mir-26a-1          96.08     51        2       0
## 13_28851.1       hsa-mir-26a-2          95.83     24        1       0
## GL894231.1_44092   hsa-mir-655         100.00     38        0       0
##                  query_start query_end dbseq_start dbseq_end evalue
## GL892871.2_43622           1        49          12        60  2e-23
## X_39130                    1        43          16        58  2e-17
## 3_7507                     1        25          16        40  1e-06
## 13_28851                   1        51          10        60  7e-20
## 13_28851.1                 1        24          14        37  4e-06
## GL894231.1_44092           1        38          23        60  7e-17
##                  bitscore          seqname
## GL892871.2_43622     97.6 GL892871.2_43622
## X_39130              77.8          X_39130
## 3_7507               42.1           3_7507
## 13_28851             85.7         13_28851
## 13_28851.1           40.1         13_28851
## GL894231.1_44092     75.8 GL894231.1_44092
```

```r
novelmir.abundance<-novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(novelmir.abundance)
```

```
## [1] 41  8
```

```r
head(novelmir.abundance)
```

```
##                    provisional.id miRDeep2.score total.read.count
## GL892871.2_43622 GL892871.2_43622       572746.3          1123410
## X_39130                   X_39130       136780.7           268282
## 3_7507                     3_7507        57547.6           112871
## 13_28851                 13_28851        40621.2            79667
## 13_28851.1               13_28851        40621.2            79667
## GL894231.1_44092 GL894231.1_44092        12387.7            24290
##                  example.miRBase.miRNA.with.the.same.seed      dbseq_id
## GL892871.2_43622                           hsa-miR-26a-5p   hsa-mir-26b
## X_39130                                    hsa-miR-660-5p   hsa-mir-660
## 3_7507                                     hsa-miR-590-3p   hsa-mir-590
## 13_28851                                   hsa-miR-26a-5p hsa-mir-26a-1
## 13_28851.1                                 hsa-miR-26a-5p hsa-mir-26a-2
## GL894231.1_44092                           hsa-miR-655-3p   hsa-mir-655
##                  perc_identical evalue      precursor.coordinate
## GL892871.2_43622         100.00  2e-23 GL892871.2:56019..56076:+
## X_39130                   97.67  2e-17    X:48640826..48640884:+
## 3_7507                    96.00  1e-06    3:11024999..11025060:+
## 13_28851                  96.08  7e-20   13:24885255..24885316:-
## 13_28851.1                95.83  4e-06   13:24885255..24885316:-
## GL894231.1_44092         100.00  7e-17 GL894231.1:23071..23131:-
```

```r
novelmir.provisional.ids<-unique(novelmir.abundance$provisional.id)
```

---------------------------------------
### 2. Create a subset data frame containing the pertinent information for filtering the miRNA precursor sequences BLASTed at eval = 1.0


```r
eval1.precursor.subset<-eval1.hsa.blast.precursor[,c("query_id", "dbseq_id", "perc_identical", "evalue", "bitscore")]
dim(eval1.precursor.subset)
```

```
## [1] 652   5
```

```r
head(eval1.precursor.subset)
```

```
##                                          query_id      dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-mir-26b
## 2 seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-26a-1
## 3 seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-26a-2
## 4 seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-4473
## 5 seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-1287
## 6 seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-6775
##   perc_identical  evalue bitscore
## 1            100 2.0e-23     97.6
## 2             92 2.0e-04     34.2
## 3             95 8.0e-04     32.2
## 4            100 5.2e-02     26.3
## 5            100 2.1e-01     24.3
## 6            100 8.2e-01     22.3
```

The number of unique sequences with homo sapien miRBase hits in this dataset at eval = 1.0


```r
length(unique(eval1.precursor.subset$query_id))
```

```
## [1] 121
```

The number of unique homo sapiens miRBase miRNAs identified by this BLAST query at eval = 1.0


```r
length(unique(eval1.precursor.subset$dbseq_id))
```

```
## [1] 455
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
length(unique(eval1.precursor.subset$query_id))
```

```
## [1] 121
```

```r
# Create the sequence list to filter
idx <- lapply(as.character(unique(eval1.precursor.subset$query_id)), function(x)
        eval1.precursor.subset[as.character(eval1.precursor.subset$query_id) == x,])

length(idx)
```

```
## [1] 121
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
#        bit <- seqblast[seqblast$bitscore == max(seqblast$bitscore),]
#                        if (nrow(bit) == 1){
#                        return(bit)
#        }

        return(ident2)
}
```

Apply filter to the candidate novel precursor sequence list at eval = 1.0


```r
stp1 <- lapply(idx, filter)

length(stp1)
```

```
## [1] 121
```

```r
head(stp1)
```

```
## [[1]]
##                                           query_id     dbseq_id
## 1  seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-26b
## 4  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-4473
## 5  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-1287
## 6  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-6775
## 9  seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-191
## 10 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-let-7b
##    perc_identical  evalue bitscore
## 1             100 2.0e-23     97.6
## 4             100 5.2e-02     26.3
## 5             100 2.1e-01     24.3
## 6             100 8.2e-01     22.3
## 9             100 8.2e-01     22.3
## 10            100 8.2e-01     22.3
## 
## [[2]]
##                                  query_id     dbseq_id perc_identical
## 11 seq|X_39130|candidatenovelpercursormiR  hsa-mir-660          97.67
## 12 seq|X_39130|candidatenovelpercursormiR hsa-mir-9500         100.00
## 13 seq|X_39130|candidatenovelpercursormiR hsa-mir-5571         100.00
## 14 seq|X_39130|candidatenovelpercursormiR hsa-mir-4723         100.00
##     evalue bitscore
## 11 2.0e-17     77.8
## 12 8.3e-01     22.3
## 13 8.3e-01     22.3
## 14 8.3e-01     22.3
## 
## [[3]]
##                                           query_id     dbseq_id
## 15 seq|GL896425.1_44884|candidatenovelpercursormiR hsa-mir-6501
## 17 seq|GL896425.1_44884|candidatenovelpercursormiR hsa-mir-4269
## 18 seq|GL896425.1_44884|candidatenovelpercursormiR hsa-mir-1202
## 19 seq|GL896425.1_44884|candidatenovelpercursormiR  hsa-mir-593
##    perc_identical evalue bitscore
## 15            100   0.89     22.3
## 17            100   0.89     22.3
## 18            100   0.89     22.3
## 19            100   0.89     22.3
## 
## [[4]]
##                                 query_id     dbseq_id perc_identical
## 22 seq|3_7507|candidatenovelpercursormiR hsa-mir-301a            100
##    evalue bitscore
## 22   0.89     22.3
## 
## [[5]]
##                                   query_id      dbseq_id perc_identical
## 23 seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-1          96.08
## 25 seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-2         100.00
## 27 seq|13_28851|candidatenovelpercursormiR  hsa-mir-4749         100.00
##     evalue bitscore
## 23 7.0e-20     85.7
## 25 8.9e-01     22.3
## 27 8.9e-01     22.3
## 
## [[6]]
##                                           query_id    dbseq_id
## 28 seq|GL894231.1_44092|candidatenovelpercursormiR hsa-mir-655
## 29 seq|GL894231.1_44092|candidatenovelpercursormiR hsa-mir-369
## 30 seq|GL894231.1_44092|candidatenovelpercursormiR hsa-mir-154
##    perc_identical  evalue bitscore
## 28            100 7.0e-17     75.8
## 29            100 4.0e-03     30.2
## 30            100 1.4e-02     28.2
```

```r
tail(stp1)
```

```
## [[1]]
##                                  query_id     dbseq_id perc_identical
## 636 seq|2_4405|candidatenovelpercursormiR hsa-mir-4259            100
## 637 seq|2_4405|candidatenovelpercursormiR hsa-mir-4322            100
## 638 seq|2_4405|candidatenovelpercursormiR  hsa-mir-183            100
##     evalue bitscore
## 636   0.87     22.3
## 637   0.87     22.3
## 638   0.87     22.3
## 
## [[2]]
##                                    query_id     dbseq_id perc_identical
## 639 seq|18_38289|candidatenovelpercursormiR hsa-mir-544b            100
##     evalue bitscore
## 639   0.87     22.3
## 
## [[3]]
##                                    query_id     dbseq_id perc_identical
## 640 seq|14_31686|candidatenovelpercursormiR hsa-mir-7854            100
## 641 seq|14_31686|candidatenovelpercursormiR hsa-mir-7111            100
## 642 seq|14_31686|candidatenovelpercursormiR hsa-mir-200c            100
##     evalue bitscore
## 640   0.78     22.3
## 641   0.78     22.3
## 642   0.78     22.3
## 
## [[4]]
##                                    query_id     dbseq_id perc_identical
## 643 seq|13_28017|candidatenovelpercursormiR hsa-mir-3927            100
## 644 seq|13_28017|candidatenovelpercursormiR  hsa-mir-337            100
##     evalue bitscore
## 643   0.85     22.3
## 644   0.85     22.3
## 
## [[5]]
##                                   query_id     dbseq_id perc_identical
## 645 seq|6_16454|candidatenovelpercursormiR  hsa-mir-659            100
## 646 seq|6_16454|candidatenovelpercursormiR hsa-mir-7705            100
## 647 seq|6_16454|candidatenovelpercursormiR hsa-mir-6746            100
## 648 seq|6_16454|candidatenovelpercursormiR hsa-mir-4764            100
## 649 seq|6_16454|candidatenovelpercursormiR hsa-mir-4641            100
## 650 seq|6_16454|candidatenovelpercursormiR hsa-mir-3074            100
## 651 seq|6_16454|candidatenovelpercursormiR hsa-mir-24-1            100
##     evalue bitscore
## 645   0.22     24.3
## 646   0.87     22.3
## 647   0.87     22.3
## 648   0.87     22.3
## 649   0.87     22.3
## 650   0.87     22.3
## 651   0.87     22.3
## 
## [[6]]
##                                 query_id     dbseq_id perc_identical
## 652 seq|1_432|candidatenovelpercursormiR hsa-mir-4801            100
##     evalue bitscore
## 652   0.89     22.3
```

How many candidate novel precursor sequences have more than one hit after the filtering?


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 105
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##  2  3  4  5  6  7  8  9 10 12 13 
## 20 17 19 13 10  6  9  4  2  3  2
```

```r
head(stp2)
```

```
## [[1]]
##                                           query_id     dbseq_id
## 1  seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-26b
## 4  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-4473
## 5  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-1287
## 6  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-6775
## 9  seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-191
## 10 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-let-7b
##    perc_identical  evalue bitscore
## 1             100 2.0e-23     97.6
## 4             100 5.2e-02     26.3
## 5             100 2.1e-01     24.3
## 6             100 8.2e-01     22.3
## 9             100 8.2e-01     22.3
## 10            100 8.2e-01     22.3
## 
## [[2]]
##                                  query_id     dbseq_id perc_identical
## 11 seq|X_39130|candidatenovelpercursormiR  hsa-mir-660          97.67
## 12 seq|X_39130|candidatenovelpercursormiR hsa-mir-9500         100.00
## 13 seq|X_39130|candidatenovelpercursormiR hsa-mir-5571         100.00
## 14 seq|X_39130|candidatenovelpercursormiR hsa-mir-4723         100.00
##     evalue bitscore
## 11 2.0e-17     77.8
## 12 8.3e-01     22.3
## 13 8.3e-01     22.3
## 14 8.3e-01     22.3
## 
## [[3]]
##                                           query_id     dbseq_id
## 15 seq|GL896425.1_44884|candidatenovelpercursormiR hsa-mir-6501
## 17 seq|GL896425.1_44884|candidatenovelpercursormiR hsa-mir-4269
## 18 seq|GL896425.1_44884|candidatenovelpercursormiR hsa-mir-1202
## 19 seq|GL896425.1_44884|candidatenovelpercursormiR  hsa-mir-593
##    perc_identical evalue bitscore
## 15            100   0.89     22.3
## 17            100   0.89     22.3
## 18            100   0.89     22.3
## 19            100   0.89     22.3
## 
## [[4]]
##                                   query_id      dbseq_id perc_identical
## 23 seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-1          96.08
## 25 seq|13_28851|candidatenovelpercursormiR hsa-mir-26a-2         100.00
## 27 seq|13_28851|candidatenovelpercursormiR  hsa-mir-4749         100.00
##     evalue bitscore
## 23 7.0e-20     85.7
## 25 8.9e-01     22.3
## 27 8.9e-01     22.3
## 
## [[5]]
##                                           query_id    dbseq_id
## 28 seq|GL894231.1_44092|candidatenovelpercursormiR hsa-mir-655
## 29 seq|GL894231.1_44092|candidatenovelpercursormiR hsa-mir-369
## 30 seq|GL894231.1_44092|candidatenovelpercursormiR hsa-mir-154
##    perc_identical  evalue bitscore
## 28            100 7.0e-17     75.8
## 29            100 4.0e-03     30.2
## 30            100 1.4e-02     28.2
## 
## [[6]]
##                                   query_id       dbseq_id perc_identical
## 35 seq|13_27466|candidatenovelpercursormiR   hsa-mir-4749            100
## 36 seq|13_27466|candidatenovelpercursormiR   hsa-mir-6805            100
## 37 seq|13_27466|candidatenovelpercursormiR   hsa-mir-6727            100
## 38 seq|13_27466|candidatenovelpercursormiR   hsa-mir-6722            100
## 39 seq|13_27466|candidatenovelpercursormiR   hsa-mir-6718            100
## 40 seq|13_27466|candidatenovelpercursormiR hsa-mir-1233-2            100
## 41 seq|13_27466|candidatenovelpercursormiR   hsa-mir-1304            100
## 42 seq|13_27466|candidatenovelpercursormiR hsa-mir-1233-1            100
## 43 seq|13_27466|candidatenovelpercursormiR    hsa-mir-874            100
##    evalue bitscore
## 35  0.014     28.2
## 36  0.220     24.3
## 37  0.850     22.3
## 38  0.850     22.3
## 39  0.850     22.3
## 40  0.850     22.3
## 41  0.850     22.3
## 42  0.850     22.3
## 43  0.850     22.3
```

Summary file of candidate novel miRNA human miRBase blast results at eval = 1.0


```r
sumblast <- do.call(rbind,stp1)
head(sumblast)
```

```
##                                           query_id     dbseq_id
## 1  seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-26b
## 4  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-4473
## 5  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-1287
## 6  seq|GL892871.2_43622|candidatenovelpercursormiR hsa-mir-6775
## 9  seq|GL892871.2_43622|candidatenovelpercursormiR  hsa-mir-191
## 10 seq|GL892871.2_43622|candidatenovelpercursormiR   hsa-let-7b
##    perc_identical  evalue bitscore
## 1             100 2.0e-23     97.6
## 4             100 5.2e-02     26.3
## 5             100 2.1e-01     24.3
## 6             100 8.2e-01     22.3
## 9             100 8.2e-01     22.3
## 10            100 8.2e-01     22.3
```

```r
dim(sumblast)
```

```
## [1] 540   5
```

How many unqiue candidate novel miRNA precursors had miRBase hits at eval = 1.0


```r
length(unique(sumblast$query_id))
```

```
## [1] 121
```

How many unique human miRBase miRNAs are identified in this dataset at eval = 1.0?


```r
length(unique(sumblast$dbseq_id))
```

```
## [1] 392
```

```r
head(table(sumblast$query_id))
```

```
## 
##   seq|1_1036|candidatenovelpercursormiR 
##                                       8 
##   seq|1_1168|candidatenovelpercursormiR 
##                                       5 
## seq|11_24386|candidatenovelpercursormiR 
##                                       4 
## seq|11_24537|candidatenovelpercursormiR 
##                                       3 
## seq|11_24928|candidatenovelpercursormiR 
##                                       3 
## seq|11_25034|candidatenovelpercursormiR 
##                                       3
```

```r
head(table(sumblast$dbseq_id))
```

```
## 
##   hsa-let-7b  hsa-mir-1-2 hsa-mir-1202 hsa-mir-1203 hsa-mir-1205 
##            1            1            1            1            1 
##  hsa-mir-122 
##            1
```

---------------------------------------

### 3. Repeat the same filtering step on the candidate novel mature miRNA at eval = 1.0
Create a subset data frame containing the pertinent information for filtering the novel mature miRNA sequences BLASTed at eval = 1.0


```r
length(unique(eval1.hsa.blast.mature$query_id))
```

```
## [1] 92
```

```r
eval1.mature.subset<-eval1.hsa.blast.mature[,c("query_id", "dbseq_id", "perc_identical", "evalue", "bitscore")]
dim(eval1.mature.subset)
```

```
## [1] 235   5
```

```r
head(eval1.mature.subset)
```

```
##                                       query_id        dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-26b-5p
## 2 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-26a-5p
## 3 seq|GL892871.2_43622|candidatenovelmaturemiR    hsa-miR-1297
## 4 seq|GL892871.2_43622|candidatenovelmaturemiR hsa-miR-3140-5p
## 5 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-876-3p
## 6          seq|X_39130|candidatenovelmaturemiR  hsa-miR-660-5p
##   perc_identical  evalue bitscore
## 1            100 9.0e-08     42.1
## 2             95 9.0e-05     32.2
## 3            100 3.0e-04     30.2
## 4            100 3.2e-01     20.3
## 5            100 3.2e-01     20.3
## 6            100 2.0e-08     44.1
```

The number of unique candidate novel mature miRNA sequences with human miRBase hits in this dataset at eval = 1.0


```r
length(unique(eval1.mature.subset$query_id))
```

```
## [1] 92
```

The number of unique human miRBase miRNAs identified by this BLAST query at eval = 1.0


```r
length(unique(eval1.mature.subset$dbseq_id))
```

```
## [1] 205
```

Format the miRBase blast output data.frame to a list to filter the sequences with multiple miRBase hits


```r
# Number of unique candidate mature miRNA sequences at eval = 1.0
length(unique(eval1.mature.subset$query_id))
```

```
## [1] 92
```

```r
# Create the unique candidate mature miRNA sequence list to filter at eval = 1.0
idx2 <- lapply(as.character(unique(eval1.mature.subset$query_id)), function(x)
        eval1.mature.subset[as.character(eval1.mature.subset$query_id) == x,])

length(idx2)
```

```
## [1] 92
```

Apply filter to the unique candidate mature miRNA sequence list to filter at eval = 1.0


```r
stp1m <- lapply(idx2, filter)

length(stp1m)
```

```
## [1] 92
```

```r
head(stp1m)
```

```
## [[1]]
##                                       query_id        dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-26b-5p
## 3 seq|GL892871.2_43622|candidatenovelmaturemiR    hsa-miR-1297
## 4 seq|GL892871.2_43622|candidatenovelmaturemiR hsa-miR-3140-5p
## 5 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-876-3p
##   perc_identical  evalue bitscore
## 1            100 9.0e-08     42.1
## 3            100 3.0e-04     30.2
## 4            100 3.2e-01     20.3
## 5            100 3.2e-01     20.3
## 
## [[2]]
##                              query_id       dbseq_id perc_identical evalue
## 6 seq|X_39130|candidatenovelmaturemiR hsa-miR-660-5p            100  2e-08
##   bitscore
## 6     44.1
## 
## [[3]]
##                             query_id       dbseq_id perc_identical evalue
## 7 seq|3_7507|candidatenovelmaturemiR hsa-miR-590-3p            100  8e-08
##   bitscore
## 7     42.1
## 
## [[4]]
## [1] query_id       dbseq_id       perc_identical evalue        
## [5] bitscore      
## <0 rows> (or 0-length row.names)
## 
## [[5]]
##                                        query_id       dbseq_id
## 12 seq|GL894231.1_44092|candidatenovelmaturemiR hsa-miR-655-3p
## 13 seq|GL894231.1_44092|candidatenovelmaturemiR hsa-miR-329-3p
## 14 seq|GL894231.1_44092|candidatenovelmaturemiR hsa-miR-369-3p
##    perc_identical evalue bitscore
## 12            100  2e-08     44.1
## 13            100  1e-03     28.2
## 14            100  5e-03     26.3
## 
## [[6]]
##                                query_id        dbseq_id perc_identical
## 15 seq|13_27466|candidatenovelmaturemiR hsa-miR-6805-3p            100
## 16 seq|13_27466|candidatenovelmaturemiR hsa-miR-6728-3p            100
## 17 seq|13_27466|candidatenovelmaturemiR hsa-miR-5010-5p            100
## 18 seq|13_27466|candidatenovelmaturemiR  hsa-miR-511-5p            100
##    evalue bitscore
## 15  0.022     24.3
## 16  0.350     20.3
## 17  0.350     20.3
## 18  0.350     20.3
```

```r
tail(stp1m)
```

```
## [[1]]
##                                query_id        dbseq_id perc_identical
## 223 seq|7_19335|candidatenovelmaturemiR hsa-miR-7156-3p            100
## 224 seq|7_19335|candidatenovelmaturemiR    hsa-miR-1291            100
## 225 seq|7_19335|candidatenovelmaturemiR hsa-miR-6838-5p            100
##     evalue bitscore
## 223  0.088     22.3
## 224  0.088     22.3
## 225  0.350     20.3
## 
## [[2]]
##                               query_id     dbseq_id perc_identical  evalue
## 226 seq|2_5120|candidatenovelmaturemiR hsa-miR-3660            100 1.0e-06
## 227 seq|2_5120|candidatenovelmaturemiR hsa-miR-7974            100 8.2e-02
##     bitscore
## 226     38.2
## 227     22.3
## 
## [[3]]
##                               query_id     dbseq_id perc_identical evalue
## 228 seq|2_4405|candidatenovelmaturemiR hsa-miR-4327            100    0.3
##     bitscore
## 228     20.3
## 
## [[4]]
##                                 query_id       dbseq_id perc_identical
## 229 seq|13_28017|candidatenovelmaturemiR hsa-miR-33a-5p            100
##     evalue bitscore
## 229   0.35     20.3
## 
## [[5]]
##                                query_id        dbseq_id perc_identical
## 230 seq|6_16454|candidatenovelmaturemiR    hsa-miR-3937            100
## 231 seq|6_16454|candidatenovelmaturemiR hsa-miR-6746-3p            100
## 232 seq|6_16454|candidatenovelmaturemiR hsa-miR-6727-5p            100
## 233 seq|6_16454|candidatenovelmaturemiR    hsa-miR-1208            100
##     evalue bitscore
## 230  0.005     26.3
## 231  0.076     22.3
## 232  0.300     20.3
## 233  0.300     20.3
## 
## [[6]]
##                                 query_id        dbseq_id perc_identical
## 234 seq|15_33472|candidatenovelmaturemiR hsa-miR-6826-5p            100
## 235 seq|15_33472|candidatenovelmaturemiR hsa-miR-6755-3p            100
##     evalue bitscore
## 234   0.35     20.3
## 235   0.35     20.3
```

How many sequences have more than one hit after the filtering


```r
stp2m <- stp1m[unlist(lapply(stp1m,nrow)) > 1]
length(stp2m)
```

```
## [1] 54
```

```r
table(unlist(lapply(stp2m,nrow)))
```

```
## 
##  2  3  4  5  6 10 
## 27  9 10  6  1  1
```

Summary file of small RNA sequence miRBase blast results


```r
sumblast.mature <- do.call(rbind,stp1m)
head(sumblast.mature)
```

```
##                                       query_id        dbseq_id
## 1 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-26b-5p
## 3 seq|GL892871.2_43622|candidatenovelmaturemiR    hsa-miR-1297
## 4 seq|GL892871.2_43622|candidatenovelmaturemiR hsa-miR-3140-5p
## 5 seq|GL892871.2_43622|candidatenovelmaturemiR  hsa-miR-876-3p
## 6          seq|X_39130|candidatenovelmaturemiR  hsa-miR-660-5p
## 7           seq|3_7507|candidatenovelmaturemiR  hsa-miR-590-3p
##   perc_identical  evalue bitscore
## 1            100 9.0e-08     42.1
## 3            100 3.0e-04     30.2
## 4            100 3.2e-01     20.3
## 5            100 3.2e-01     20.3
## 6            100 2.0e-08     44.1
## 7            100 8.0e-08     42.1
```

```r
dim(sumblast.mature)
```

```
## [1] 203   5
```

How many unique homo sapiens miRNAs were identified in this BLAST query at eval = 1.0


```r
length(unique(sumblast.mature$dbseq_id))
```

```
## [1] 178
```

How many unique candidate novel mature miRNAs were identified in this BLAST query at eval = 1.0


```r
length(unique(sumblast.mature$query_id))
```

```
## [1] 90
```

Notice that there are two fewer candidate novel miRNAs than prior to filtering;

This is because there are two BLAST hits here that are less than 96% identical and were removed in the filtering process.

These include the sequence 13_28851 and 14_31131

Just to check:


```r
eval1.mature.subset[9:11,]
```

```
##                                query_id         dbseq_id perc_identical
## 9  seq|13_28851|candidatenovelmaturemiR   hsa-miR-26a-5p          95.45
## 10 seq|13_28851|candidatenovelmaturemiR   hsa-miR-26b-5p          90.00
## 11 seq|13_28851|candidatenovelmaturemiR hsa-miR-26a-1-3p          88.89
##     evalue bitscore
## 9  5.0e-06     36.2
## 10 2.1e-02     24.3
## 11 3.2e-01     20.3
```

```r
eval1.mature.subset[111:113,]
```

```
##                                 query_id         dbseq_id perc_identical
## 111 seq|14_31131|candidatenovelmaturemiR hsa-miR-6715b-5p          94.74
## 112 seq|14_31131|candidatenovelmaturemiR hsa-miR-6715a-3p          94.74
## 113 seq|14_31131|candidatenovelmaturemiR hsa-miR-6715b-3p          88.89
##     evalue bitscore
## 111  3e-04     30.2
## 112  3e-04     30.2
## 113  3e-01     20.3
```

## Save data


```r
write.table(novelmir.abundance, file="../8_blastn_candidate_novel_miRNA_output/3_candidate_novel_miRNA_precursor_abundance_e5.txt")
write.table(sumblast, file="../8_blastn_candidate_novel_miRNA_output/4_filtered_candidate_novel_precursor_eval1.txt")
write.table(sumblast.mature, file="../8_blastn_candidate_novel_miRNA_output/5_filtered_candidate_novel_mature_eval1.txt")
```

Create a list of the pertinent precursor provisional.ids to extract the correct pdfs for a supplemental figure


```r
write.table(novelmir.provisional.ids, file="../8_blastn_candidate_novel_miRNA_output/6_filtered_candidate_novel_miRNA_precursor_ids_e5.txt", row.names=FALSE, col.names=FALSE)
```

