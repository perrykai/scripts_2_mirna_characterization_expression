**Script:** `2_novel_mirna_sequence_analysis.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts`

**Date:**  10/11/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output`

**Input File(s):** `2_extracted_mirdeep2_core_predicted_novel_mirna.csv`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output`

**Output File(s):** 

1. `5_putative_novel_miRNA_filtered_candidates.Rdata`
2. `6_candidate_novel_mature_mir.fa`
3. `7_candidate_novel_precursor_mir.fa`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to filter the extracted putative novel miRNAs for the following characteristics:
 
1. Novel miRNA candidate must have > 90% probability of being a true positive (filter by miRDeep2 score and then estimated probability that the miRNA candidate is a true positive)
2. Hairpins must have significant Randfold p-values
3. Minimum read counts for putative mature and star strand sequences

Additionally, the novel candidate miRNAs will be written into fasta file format for BLAST against human and mouse miRBase databases (mature seq and precursor seq separately)

## Install libraries


```r
library("ShortRead")
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: methods
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:stats':
## 
##     xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
```

```
## Loading required package: BiocParallel
```

```
## Loading required package: Biostrings
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: XVector
```

```
## Loading required package: Rsamtools
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomicAlignments
```

## Load data



```r
novelmir<-read.table("../2_extracted_mirdeep2_core_predicted_novel_mirna.csv", sep="\t", header=TRUE)
dim(novelmir)
```

```
## [1] 1199   17
```

```r
colnames(novelmir)
```

```
##  [1] "provisional.id"                                                   
##  [2] "miRDeep2.score"                                                   
##  [3] "estimated.probability.that.the.miRNA.candidate.is.a.true.positive"
##  [4] "rfam.alert"                                                       
##  [5] "total.read.count"                                                 
##  [6] "mature.read.count"                                                
##  [7] "loop.read.count"                                                  
##  [8] "star.read.count"                                                  
##  [9] "significant.randfold.p.value"                                     
## [10] "miRBase.miRNA"                                                    
## [11] "example.miRBase.miRNA.with.the.same.seed"                         
## [12] "UCSC.browser"                                                     
## [13] "NCBI.blastn"                                                      
## [14] "consensus.mature.sequence"                                        
## [15] "consensus.star.sequence"                                          
## [16] "consensus.precursor.sequence"                                     
## [17] "precursor.coordinate"
```

```r
head(novelmir)
```

```
##     provisional.id miRDeep2.score
## 1 GL892871.2_43622       572746.3
## 2          X_39130       136780.7
## 3 GL896425.1_44884       128796.6
## 4         16_36270       113215.6
## 5           3_7507        57547.6
## 6         13_28851        40621.2
##   estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1                                                         91 +/- 2%
## 2                                                         91 +/- 2%
## 3                                                         91 +/- 2%
## 4                                                         91 +/- 2%
## 5                                                         91 +/- 2%
## 6                                                         91 +/- 2%
##   rfam.alert total.read.count mature.read.count loop.read.count
## 1          -          1123410           1119829               0
## 2          -           268282            268080               0
## 3          -           252627            249026               0
## 4          -           222068            222028               0
## 5          -           112871            104821               0
## 6          -            79667             79633               0
##   star.read.count significant.randfold.p.value miRBase.miRNA
## 1            3581                          yes             -
## 2             202                          yes             -
## 3            3601                          yes             -
## 4              40                           no             -
## 5            8050                          yes             -
## 6              34                          yes             -
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                           hsa-miR-26a-5p            -           -
## 2                           hsa-miR-660-5p            -           -
## 3                                        -            -           -
## 4                            hsa-let-7a-5p            -           -
## 5                           hsa-miR-590-3p            -           -
## 6                           hsa-miR-26a-5p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 3    uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 4      ugagguaguaggcugugugg     aauagcaacaacaacaaca
## 5     uaauuuuauguauaagcuagu  gagcuuauucauaaaaguacag
## 6    uucaaguaacccaggauaggcu  ccuauucuugguuacuugcacg
##                                    consensus.precursor.sequence
## 1     uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3 agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 4              aauagcaacaacaacaacaacaaaagaaugagguaguaggcugugugg
## 5 gagcuuauucauaaaaguacaguauaauccaguaaaccuguaauuuuauguauaagcuagu
## 6 uucaaguaacccaggauaggcugugcaggucccaaggggccuauucuugguuacuugcacg
##        precursor.coordinate
## 1 GL892871.2:56019..56076:+
## 2    X:48640826..48640884:+
## 3   GL896425.1:1750..1811:-
## 4     16:6629237..6629285:-
## 5    3:11024999..11025060:+
## 6   13:24885255..24885316:-
```

## Analysis

### 1. Filter by the miRDeep2 score and the estimated probability of the novel miRNA being true positives:


```r
sum(novelmir$miRDeep2.score >= 10)
```

```
## [1] 352
```

```r
table(novelmir$estimated.probability.that.the.miRNA.candidate.is.a.true.positive)
```

```
## 
## 20 +/- 2% 49 +/- 2% 52 +/- 3% 56 +/- 2% 80 +/- 2% 91 +/- 1% 91 +/- 2% 
##       152       114        26       325        49       121       412
```

```r
novelmir[novelmir$estimated.probability.that.the.miRNA.candidate.is.a.true.positive == "91 +/- 1%", "miRDeep2.score"]
```

```
##   [1] 5.9 5.9 5.9 5.9 5.9 5.9 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.7 5.7 5.7 5.7
##  [18] 5.7 5.7 5.7 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6
##  [35] 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.4 5.4 5.4
##  [52] 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.3 5.3 5.3 5.3 5.3 5.3 5.3
##  [69] 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.2 5.2
##  [86] 5.2 5.2 5.2 5.2 5.2 5.2 5.2 5.2 5.2 5.1 5.1 5.1 5.1 5.1 5.1 5.1 5.1
## [103] 5.1 5.1 5.1 5.1 5.1 5.1 5.1 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0
## [120] 5.0 5.0
```

First filter by miRDeep2 score, then by estimated probability:


```r
novelmir10<-novelmir[novelmir$miRDeep2.score >= 10,]
dim(novelmir10)
```

```
## [1] 352  17
```

```r
table(novelmir10$estimated.probability.that.the.miRNA.candidate.is.a.true.positive)
```

```
## 
## 20 +/- 2% 49 +/- 2% 52 +/- 3% 56 +/- 2% 80 +/- 2% 91 +/- 1% 91 +/- 2% 
##         0         0         0         0         0         0       352
```

```r
sum(novelmir10$estimated.probability.that.the.miRNA.candidate.is.a.true.positive != "91 +/- 2%")
```

```
## [1] 0
```

Filtering by miRDeep2 score removed any miRNA with estimated probability < "91 +/- 2%"

### 2. Hairpins must have a significant Randfold p-value


```r
sum(novelmir10$significant.randfold.p.value == "yes")
```

```
## [1] 271
```

271 potential miRNA candidates have significant Randfold p-value


```r
novelmir10sigrandfold<-novelmir10[novelmir10$significant.randfold.p.value == "yes",]
dim(novelmir10sigrandfold)
```

```
## [1] 271  17
```

```r
head(novelmir10sigrandfold)
```

```
##     provisional.id miRDeep2.score
## 1 GL892871.2_43622       572746.3
## 2          X_39130       136780.7
## 3 GL896425.1_44884       128796.6
## 5           3_7507        57547.6
## 6         13_28851        40621.2
## 7          5_13151        12394.0
##   estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1                                                         91 +/- 2%
## 2                                                         91 +/- 2%
## 3                                                         91 +/- 2%
## 5                                                         91 +/- 2%
## 6                                                         91 +/- 2%
## 7                                                         91 +/- 2%
##   rfam.alert total.read.count mature.read.count loop.read.count
## 1          -          1123410           1119829               0
## 2          -           268282            268080               0
## 3          -           252627            249026               0
## 5          -           112871            104821               0
## 6          -            79667             79633               0
## 7          -            24310             24307               0
##   star.read.count significant.randfold.p.value miRBase.miRNA
## 1            3581                          yes             -
## 2             202                          yes             -
## 3            3601                          yes             -
## 5            8050                          yes             -
## 6              34                          yes             -
## 7               3                          yes             -
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                           hsa-miR-26a-5p            -           -
## 2                           hsa-miR-660-5p            -           -
## 3                                        -            -           -
## 5                           hsa-miR-590-3p            -           -
## 6                           hsa-miR-26a-5p            -           -
## 7                          hsa-miR-5011-5p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 3    uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 5     uaauuuuauguauaagcuagu  gagcuuauucauaaaaguacag
## 6    uucaaguaacccaggauaggcu  ccuauucuugguuacuugcacg
## 7     uauauauauauauguucguau    augaguauauauauauauau
##                                                          consensus.precursor.sequence
## 1                           uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2                          uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3                       agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 5                       gagcuuauucauaaaaguacaguauaauccaguaaaccuguaauuuuauguauaagcuagu
## 6                       uucaaguaacccaggauaggcugugcaggucccaaggggccuauucuugguuacuugcacg
## 7 uauauauauauauguucguauauucacauacaaucuucacaauuacccaauuaaauaagugcuaugaguauauauauauauau
##        precursor.coordinate
## 1 GL892871.2:56019..56076:+
## 2    X:48640826..48640884:+
## 3   GL896425.1:1750..1811:-
## 5    3:11024999..11025060:+
## 6   13:24885255..24885316:-
## 7    5:34081933..34082016:-
```

Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences


```r
table(novelmir10sigrandfold$rfam.alert)
```

```
## 
##         -      rRNA rRNA/tRNA 
##       271         0         0
```

```r
sum(novelmir10sigrandfold$rfam.alert != "-")
```

```
## [1] 0
```

### 3. Minimum read counts for putative mature and star strand sequences: require at least 10 counts for each


```r
summary(novelmir10sigrandfold$mature.read.count)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       1      25      88    7759     388 1120000
```

```r
summary(novelmir10sigrandfold$star.read.count)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     1.0     4.0    13.0   149.3    42.5  8050.0
```

```r
sum(novelmir10sigrandfold$mature.read.count <= 10)
```

```
## [1] 16
```

```r
novelmir10sigrandfoldmincounts<-novelmir10sigrandfold[novelmir10sigrandfold$mature.read.count > 10,]
dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 255  17
```

```r
sum(novelmir10sigrandfoldmincounts$star.read.count <= 10)
```

```
## [1] 123
```

```r
novelmir10sigrandfoldmincounts<-novelmir10sigrandfoldmincounts[novelmir10sigrandfoldmincounts$star.read.count > 10,]
dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 132  17
```

```r
sum(novelmir10sigrandfoldmincounts$star.read.count <=10)
```

```
## [1] 0
```

```r
sum(novelmir10sigrandfoldmincounts$mature.read.count <=10)
```

```
## [1] 0
```

Investigate the final output:


```r
dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 132  17
```

```r
head(novelmir10sigrandfoldmincounts)
```

```
##     provisional.id miRDeep2.score
## 1 GL892871.2_43622       572746.3
## 2          X_39130       136780.7
## 3 GL896425.1_44884       128796.6
## 5           3_7507        57547.6
## 6         13_28851        40621.2
## 8 GL894231.1_44092        12387.7
##   estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1                                                         91 +/- 2%
## 2                                                         91 +/- 2%
## 3                                                         91 +/- 2%
## 5                                                         91 +/- 2%
## 6                                                         91 +/- 2%
## 8                                                         91 +/- 2%
##   rfam.alert total.read.count mature.read.count loop.read.count
## 1          -          1123410           1119829               0
## 2          -           268282            268080               0
## 3          -           252627            249026               0
## 5          -           112871            104821               0
## 6          -            79667             79633               0
## 8          -            24290             24219               0
##   star.read.count significant.randfold.p.value miRBase.miRNA
## 1            3581                          yes             -
## 2             202                          yes             -
## 3            3601                          yes             -
## 5            8050                          yes             -
## 6              34                          yes             -
## 8              71                          yes             -
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                           hsa-miR-26a-5p            -           -
## 2                           hsa-miR-660-5p            -           -
## 3                                        -            -           -
## 5                           hsa-miR-590-3p            -           -
## 6                           hsa-miR-26a-5p            -           -
## 8                           hsa-miR-655-3p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 3    uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 5     uaauuuuauguauaagcuagu  gagcuuauucauaaaaguacag
## 6    uucaaguaacccaggauaggcu  ccuauucuugguuacuugcacg
## 8    auaauacaugguuaaccucuuu  agagguuauccguguuauguuc
##                                    consensus.precursor.sequence
## 1     uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3 agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 5 gagcuuauucauaaaaguacaguauaauccaguaaaccuguaauuuuauguauaagcuagu
## 6 uucaaguaacccaggauaggcugugcaggucccaaggggccuauucuugguuacuugcacg
## 8  agagguuauccguguuauguucgcuucauucaucaugaauaauacaugguuaaccucuuu
##        precursor.coordinate
## 1 GL892871.2:56019..56076:+
## 2    X:48640826..48640884:+
## 3   GL896425.1:1750..1811:-
## 5    3:11024999..11025060:+
## 6   13:24885255..24885316:-
## 8 GL894231.1:23071..23131:-
```

How many of the novel miRNA candidates have a human miRBase miRNA with the same seed sequence


```r
sum(novelmir10sigrandfoldmincounts$example.miRBase.miRNA.with.the.same.seed != "-")
```

```
## [1] 45
```

```r
sum(novelmir10sigrandfoldmincounts$example.miRBase.miRNA.with.the.same.seed == "-")
```

```
## [1] 87
```

What is the summary of the mature and star read counts now?


```r
summary(novelmir10sigrandfoldmincounts$mature.read.count)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##      11.0      54.0     183.5   15330.0    1096.0 1120000.0
```

```r
summary(novelmir10sigrandfoldmincounts$star.read.count)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    11.0    22.0    44.0   300.7   120.8  8050.0
```

Subset the provisional.id, consensus.mature.sequence and consensus.precursor.sequence for use in BLASTN against known human and mouse miRBase sequences:


```r
novelmircandidateBLAST<-novelmir10sigrandfoldmincounts[,c("provisional.id", "consensus.mature.sequence", "consensus.precursor.sequence")]
dim(novelmircandidateBLAST)
```

```
## [1] 132   3
```

```r
head(novelmircandidateBLAST)
```

```
##     provisional.id consensus.mature.sequence
## 1 GL892871.2_43622    uucaaguaauucaggauagguu
## 2          X_39130    uacccauugcauaucggaguug
## 3 GL896425.1_44884    uggugccugacgucuuggcagu
## 5           3_7507     uaauuuuauguauaagcuagu
## 6         13_28851    uucaaguaacccaggauaggcu
## 8 GL894231.1_44092    auaauacaugguuaaccucuuu
##                                    consensus.precursor.sequence
## 1     uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3 agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 5 gagcuuauucauaaaaguacaguauaauccaguaaaccuguaauuuuauguauaagcuagu
## 6 uucaaguaacccaggauaggcugugcaggucccaaggggccuauucuugguuacuugcacg
## 8  agagguuauccguguuauguucgcuucauucaucaugaauaauacaugguuaaccucuuu
```

### 4. The novel candidate miRNAs (both precursor and mature) will be written into fasta file format for BLAST against human and mouse miRBase databases (mature seq and precursor seq separately)

#### First, prepare the sequence names for the candidate novel mature miR


```r
seqnamesmature<-paste("seq|",novelmircandidateBLAST$provisional.id, "|candidate novel mature miR")
seqnamesmature<-gsub(" ", "", seqnamesmature, fixed = TRUE)
head(seqnamesmature)
```

```
## [1] "seq|GL892871.2_43622|candidatenovelmaturemiR"
## [2] "seq|X_39130|candidatenovelmaturemiR"         
## [3] "seq|GL896425.1_44884|candidatenovelmaturemiR"
## [4] "seq|3_7507|candidatenovelmaturemiR"          
## [5] "seq|13_28851|candidatenovelmaturemiR"        
## [6] "seq|GL894231.1_44092|candidatenovelmaturemiR"
```

Create the BStringSet object with the mature sequence names


```r
matureids<-BStringSet(seqnamesmature)
head(matureids)
```

```
##   A BStringSet instance of length 6
##     width seq
## [1]    44 seq|GL892871.2_43622|candidatenovelmaturemiR
## [2]    35 seq|X_39130|candidatenovelmaturemiR
## [3]    44 seq|GL896425.1_44884|candidatenovelmaturemiR
## [4]    34 seq|3_7507|candidatenovelmaturemiR
## [5]    36 seq|13_28851|candidatenovelmaturemiR
## [6]    44 seq|GL894231.1_44092|candidatenovelmaturemiR
```

Prepare the candidate novel mature sequences for compatability with DNAStringSet: Substitute "t" for "u" as DNAStringSet only uses DNA alphabet


```r
novelmircandidateBLAST$consensus.mature.seq.tadjust <- gsub("u","t", novelmircandidateBLAST$consensus.mature.sequence)
```

Create the DNAStringSet object with the candidate novel mature sequence reads


```r
novelmatureseqreads<-DNAStringSet(novelmircandidateBLAST$consensus.mature.seq.tadjust)
head(novelmatureseqreads)
```

```
##   A DNAStringSet instance of length 6
##     width seq
## [1]    22 TTCAAGTAATTCAGGATAGGTT
## [2]    22 TACCCATTGCATATCGGAGTTG
## [3]    22 TGGTGCCTGACGTCTTGGCAGT
## [4]    21 TAATTTTATGTATAAGCTAGT
## [5]    22 TTCAAGTAACCCAGGATAGGCT
## [6]    22 ATAATACATGGTTAACCTCTTT
```

Use the ShortRead command to combine the candidate novel mature sequence IDs and the sequences


```r
candidatenovelmaturefasta<-ShortRead(sread=novelmatureseqreads ,id=matureids)
```

#### Then, prepare the sequence names for the candidate novel precursor miR:


```r
seqnamesprecursor<-paste("seq|", novelmircandidateBLAST$provisional.id, "|candidate novel percursor miR")
seqnamesprecursor<-gsub(" ", "", seqnamesprecursor, fixed = TRUE)

head(seqnamesprecursor)
```

```
## [1] "seq|GL892871.2_43622|candidatenovelpercursormiR"
## [2] "seq|X_39130|candidatenovelpercursormiR"         
## [3] "seq|GL896425.1_44884|candidatenovelpercursormiR"
## [4] "seq|3_7507|candidatenovelpercursormiR"          
## [5] "seq|13_28851|candidatenovelpercursormiR"        
## [6] "seq|GL894231.1_44092|candidatenovelpercursormiR"
```

Create the BStringSet object with the precursor sequence names


```r
precursorids<-BStringSet(seqnamesprecursor)
head(precursorids)
```

```
##   A BStringSet instance of length 6
##     width seq
## [1]    47 seq|GL892871.2_43622|candidatenovelpercursormiR
## [2]    38 seq|X_39130|candidatenovelpercursormiR
## [3]    47 seq|GL896425.1_44884|candidatenovelpercursormiR
## [4]    37 seq|3_7507|candidatenovelpercursormiR
## [5]    39 seq|13_28851|candidatenovelpercursormiR
## [6]    47 seq|GL894231.1_44092|candidatenovelpercursormiR
```

Prepare the candidate precursor sequences for compatability with DNAStringSet: Substitute "t" for "u" as DNAStringSet only uses DNA alphabet


```r
novelmircandidateBLAST$consensus.precursor.seq.tadjust<- gsub("u","t", novelmircandidateBLAST$consensus.precursor.sequence)
```

Create the DNAStringSet object with the novel candidate precursor sequence reads


```r
novelprecursorseqreads<-DNAStringSet(novelmircandidateBLAST$consensus.precursor.seq.tadjust)
head(novelprecursorseqreads)
```

```
##   A DNAStringSet instance of length 6
##     width seq
## [1]    57 TTCAAGTAATTCAGGATAGGTTGTGTGCTGTCCAGCCTGTTCTCCATTACTTGGCTC
## [2]    58 TACCCATTGCATATCGGAGTTGTGAATTCTCAAAGCACCTCCTATGTGCATGGTTTAC
## [3]    61 AGCCAGGGCTGCAGGCACTGACATTCACCCATGGTATTGTGGTGCCTGACGTCTTGGCAGT
## [4]    61 GAGCTTATTCATAAAAGTACAGTATAATCCAGTAAACCTGTAATTTTATGTATAAGCTAGT
## [5]    61 TTCAAGTAACCCAGGATAGGCTGTGCAGGTCCCAAGGGGCCTATTCTTGGTTACTTGCACG
## [6]    60 AGAGGTTATCCGTGTTATGTTCGCTTCATTCATCATGAATAATACATGGTTAACCTCTTT
```

Use the ShortRead command to combine the candidate novel precursor sequence IDs and the sequences


```r
candidatenovelprecursorfasta<-ShortRead(sread=novelprecursorseqreads ,id=precursorids)
```

## Visualize

## Save data



```r
save(novelmir10sigrandfoldmincounts, file="../5_putative_novel_miRNA_filtered_candidates.Rdata")
writeFasta(candidatenovelmaturefasta, file="../6_candidate_novel_mature_mir.fa")
writeFasta(candidatenovelprecursorfasta, file="../7_candidate_novel_precursor_mir.fa")
```

