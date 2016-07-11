**Script:** `16_summarize_blast_rfam_mus_musculus_results.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`

**Date:**  7/7/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Input File(s):** `4_musmusculus_rfam_blastn_output_e6.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`

**Output File(s):** `8_musmusculus_filtered_rfam_blastn_results.Rdata`

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
mmu<-read.table("../4_musmusculus_rfam_blastn_output_e6.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})
```

```
##    user  system elapsed 
##   0.214   0.006   0.239
```

```r
dim(mmu)
```

```
## [1] 20788    12
```

```r
head(mmu)
```

```
##                query_id                                        dbseq_id
## 1 seq_125916817_x153124        RF00683;mir-143;AEKR01113265.1/9034-8929
## 2  seq_137069853_x95209        RF00783;mir-484;AEKR01206416.1/2916-2850
## 3  seq_144363225_x67330        RF00783;mir-484;AEKR01206416.1/2916-2850
## 4  seq_162056863_x28478        RF00783;mir-484;AEKR01206416.1/2916-2850
## 5  seq_164216929_x25364 RF01960;SSU_rRNA_eukarya;CT010467.6/67702-65856
## 6  seq_165552505_x23995 RF01960;SSU_rRNA_eukarya;CT010467.6/67702-65856
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     22        0       0           1        22          61
## 2            100     22        0       0           1        22           5
## 3            100     22        0       0           1        22           5
## 4            100     23        0       0           1        23           5
## 5            100     25        0       0           1        25         844
## 6            100     24        0       0           1        24         845
##   dbseq_end evalue bitscore
## 1        82  4e-07     44.1
## 2        26  4e-07     44.1
## 3        26  4e-07     44.1
## 4        27  1e-07     46.1
## 5       868  8e-09     50.1
## 6       868  3e-08     48.1
```

```r
str(mmu)
```

```
## 'data.frame':	20788 obs. of  12 variables:
##  $ query_id      : Factor w/ 14134 levels "seq_125916817_x153124",..: 1 2 3 4 5 6 7 8 9 9 ...
##  $ dbseq_id      : Factor w/ 284 levels "RF00001;5S_rRNA;AAHY01033596.1/32099-31993",..: 256 261 261 261 284 284 284 258 59 51 ...
##  $ perc_identical: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ length        : int  22 22 22 23 25 24 24 23 25 25 ...
##  $ mismatch      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  22 22 22 23 25 24 24 23 25 25 ...
##  $ dbseq_start   : int  61 5 5 5 844 845 1 6 1 1 ...
##  $ dbseq_end     : int  82 26 26 27 868 868 24 28 25 25 ...
##  $ evalue        : num  4e-07 4e-07 4e-07 1e-07 8e-09 3e-08 3e-08 1e-07 8e-09 8e-09 ...
##  $ bitscore      : num  44.1 44.1 44.1 46.1 50.1 48.1 48.1 46.1 50.1 50.1 ...
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
##                query_id                                        dbseq_id
## 1 seq_125916817_x153124        RF00683;mir-143;AEKR01113265.1/9034-8929
## 2  seq_137069853_x95209        RF00783;mir-484;AEKR01206416.1/2916-2850
## 3  seq_144363225_x67330        RF00783;mir-484;AEKR01206416.1/2916-2850
## 4  seq_162056863_x28478        RF00783;mir-484;AEKR01206416.1/2916-2850
## 5  seq_164216929_x25364 RF01960;SSU_rRNA_eukarya;CT010467.6/67702-65856
## 6  seq_165552505_x23995 RF01960;SSU_rRNA_eukarya;CT010467.6/67702-65856
##   perc_identical length mismatch gapopen query_start query_end dbseq_start
## 1            100     22        0       0           1        22          61
## 2            100     22        0       0           1        22           5
## 3            100     22        0       0           1        22           5
## 4            100     23        0       0           1        23           5
## 5            100     25        0       0           1        25         844
## 6            100     24        0       0           1        24         845
##   dbseq_end evalue bitscore rfam_accession     rfam_biotype
## 1        82  4e-07     44.1        RF00683          mir-143
## 2        26  4e-07     44.1        RF00783          mir-484
## 3        26  4e-07     44.1        RF00783          mir-484
## 4        27  1e-07     46.1        RF00783          mir-484
## 5       868  8e-09     50.1        RF01960 SSU_rRNA_eukarya
## 6       868  3e-08     48.1        RF01960 SSU_rRNA_eukarya
##       rfam_seqnamestartend
## 1 AEKR01113265.1/9034-8929
## 2 AEKR01206416.1/2916-2850
## 3 AEKR01206416.1/2916-2850
## 4 AEKR01206416.1/2916-2850
## 5   CT010467.6/67702-65856
## 6   CT010467.6/67702-65856
```

The number of unique sequences with mus musculus Rfam hits in this dataset


```r
length(unique(mmu$query_id))
```

```
## [1] 14134
```

Create a subset data frame containing the pertinent information for filtering


```r
rfamsubset<-mmu[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
```

```
## [1] 20788     6
```

```r
head(rfamsubset)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_125916817_x153124            100  4e-07     44.1        RF00683
## 2  seq_137069853_x95209            100  4e-07     44.1        RF00783
## 3  seq_144363225_x67330            100  4e-07     44.1        RF00783
## 4  seq_162056863_x28478            100  1e-07     46.1        RF00783
## 5  seq_164216929_x25364            100  8e-09     50.1        RF01960
## 6  seq_165552505_x23995            100  3e-08     48.1        RF01960
##       rfam_biotype
## 1          mir-143
## 2          mir-484
## 3          mir-484
## 4          mir-484
## 5 SSU_rRNA_eukarya
## 6 SSU_rRNA_eukarya
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
##  55.199   0.307   9.424
```

```r
length(idx)
```

```
## [1] 14134
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
##   2.109   0.120   0.621
```

```r
length(stp1)
```

```
## [1] 14134
```

```r
head(stp1)
```

```
## [[1]]
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_125916817_x153124            100  4e-07     44.1        RF00683
##   rfam_biotype
## 1      mir-143
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 2 seq_137069853_x95209            100  4e-07     44.1        RF00783
##   rfam_biotype
## 2      mir-484
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 3 seq_144363225_x67330            100  4e-07     44.1        RF00783
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
## 5 seq_164216929_x25364            100  8e-09     50.1        RF01960
##       rfam_biotype
## 5 SSU_rRNA_eukarya
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 6 seq_165552505_x23995            100  3e-08     48.1        RF01960
##       rfam_biotype
## 6 SSU_rRNA_eukarya
```

```r
tail(stp1)
```

```
## [[1]]
##               query_id perc_identical evalue bitscore rfam_accession
## 20782 seq_227422184_x2          96.55  1e-08     50.1        RF01960
##           rfam_biotype
## 20782 SSU_rRNA_eukarya
## 
## [[2]]
##               query_id perc_identical evalue bitscore rfam_accession
## 20783 seq_227422198_x2          93.33  7e-07     44.1        RF00005
##       rfam_biotype
## 20783         tRNA
## 
## [[3]]
##               query_id perc_identical evalue bitscore rfam_accession
## 20784 seq_227422454_x2          96.15  5e-07     44.1        RF01960
##           rfam_biotype
## 20784 SSU_rRNA_eukarya
## 
## [[4]]
##               query_id perc_identical evalue bitscore rfam_accession
## 20785 seq_227423370_x2           96.3  1e-07     46.1        RF01960
##           rfam_biotype
## 20785 SSU_rRNA_eukarya
## 
## [[5]]
##               query_id perc_identical evalue bitscore rfam_accession
## 20786 seq_227423398_x2          96.43  4e-08     48.1        RF01960
##           rfam_biotype
## 20786 SSU_rRNA_eukarya
## 
## [[6]]
##               query_id perc_identical evalue bitscore rfam_accession
## 20787 seq_227424182_x2          96.43  4e-08     48.1        RF00005
##       rfam_biotype
## 20787         tRNA
```

How many sequences have more than one annotation after the filtering


```r
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
```

```
## [1] 4332
```

```r
table(unlist(lapply(stp2,nrow)))
```

```
## 
##    2    3    4    5    6    7    8    9   11   12   13   14 
## 4115  128   48    9    2    1    3    4    1   12    8    1
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
##   1.899   0.186   0.445
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
## 1 seq_125916817_x153124            100  4e-07     44.1        RF00683
## 2  seq_137069853_x95209            100  4e-07     44.1        RF00783
## 3  seq_144363225_x67330            100  4e-07     44.1        RF00783
## 4  seq_162056863_x28478            100  1e-07     46.1        RF00783
## 5  seq_164216929_x25364            100  8e-09     50.1        RF01960
## 6  seq_165552505_x23995            100  3e-08     48.1        RF01960
##       rfam_biotype
## 1          mir-143
## 2          mir-484
## 3          mir-484
## 4          mir-484
## 5 SSU_rRNA_eukarya
## 6 SSU_rRNA_eukarya
```

```r
dim(sumblast4)
```

```
## [1] 14134     6
```

```r
length(unique(sumblast4$query_id))
```

```
## [1] 14134
```

```r
table(sumblast4$rfam_biotype)
```

```
## 
##        5_8S_rRNA          5S_rRNA         Histone3       IRES_L-myc 
##               38               85                8                2 
##      Metazoa_SRP          mir-101          mir-143          mir-146 
##                3                6               49               18 
##          mir-150           mir-17          mir-188          mir-193 
##                2                4                2                1 
##         mir-1937          mir-194           mir-21          mir-299 
##              114                1               19               20 
##          mir-374          mir-484          mir-503          mir-598 
##                1              229               87               10 
##          mir-689           mir-92          SNORA28          SNORA73 
##                4               30               13                4 
##          SNORD58           snoU85 SSU_rRNA_eukarya             tRNA 
##                1               37            11557             1777 
##               U1               U2               U6            Y_RNA 
##                1                6                1                4
```

```r
uniqsumblast4<-as.data.frame(table(sumblast4$rfam_biotype))
colnames(uniqsumblast4)<-c("Gene_Biotype", "Freq")
uniqsumblast4
```

```
##        Gene_Biotype  Freq
## 1         5_8S_rRNA    38
## 2           5S_rRNA    85
## 3          Histone3     8
## 4        IRES_L-myc     2
## 5       Metazoa_SRP     3
## 6           mir-101     6
## 7           mir-143    49
## 8           mir-146    18
## 9           mir-150     2
## 10           mir-17     4
## 11          mir-188     2
## 12          mir-193     1
## 13         mir-1937   114
## 14          mir-194     1
## 15           mir-21    19
## 16          mir-299    20
## 17          mir-374     1
## 18          mir-484   229
## 19          mir-503    87
## 20          mir-598    10
## 21          mir-689     4
## 22           mir-92    30
## 23          SNORA28    13
## 24          SNORA73     4
## 25          SNORD58     1
## 26           snoU85    37
## 27 SSU_rRNA_eukarya 11557
## 28             tRNA  1777
## 29               U1     1
## 30               U2     6
## 31               U6     1
## 32            Y_RNA     4
```

Add the column of sequence count to the sumblast data frame


```r
sumblast4$seq_count<-as.numeric(str_split_fixed(sumblast4$query_id, "_x", 2)[,2])
head(sumblast4)
```

```
##                query_id perc_identical evalue bitscore rfam_accession
## 1 seq_125916817_x153124            100  4e-07     44.1        RF00683
## 2  seq_137069853_x95209            100  4e-07     44.1        RF00783
## 3  seq_144363225_x67330            100  4e-07     44.1        RF00783
## 4  seq_162056863_x28478            100  1e-07     46.1        RF00783
## 5  seq_164216929_x25364            100  8e-09     50.1        RF01960
## 6  seq_165552505_x23995            100  3e-08     48.1        RF01960
##       rfam_biotype seq_count
## 1          mir-143    153124
## 2          mir-484     95209
## 3          mir-484     67330
## 4          mir-484     28478
## 5 SSU_rRNA_eukarya     25364
## 6 SSU_rRNA_eukarya     23995
```

Use the "by" function to sum the sequence counts by their gene biotypes


```r
totalsumbiotype4<-as.matrix(by(sumblast4$seq_count, sumblast4$rfam_biotype, sum))
totalsumbiotype4
```

```
##                    [,1]
## 5_8S_rRNA          2105
## 5S_rRNA            1341
## Histone3             47
## IRES_L-myc            7
## Metazoa_SRP           8
## mir-101             465
## mir-143          168068
## mir-146             182
## mir-150              54
## mir-17               43
## mir-188              38
## mir-193               9
## mir-1937           5336
## mir-194               7
## mir-21              997
## mir-299            5607
## mir-374               3
## mir-484          220533
## mir-503           24682
## mir-598              25
## mir-689             122
## mir-92            13510
## SNORA28             188
## SNORA73              17
## SNORD58               3
## snoU85              507
## SSU_rRNA_eukarya 582009
## tRNA             124439
## U1                    4
## U2                   18
## U6                    2
## Y_RNA                 8
```

```r
if (sum(rownames(totalsumbiotype4) != uniqsumblast4$Gene_Biotype)) stop ("Gene_Biotypes not equal")
```

As a check, manually sum the 5s_rRNAs and the tRNA fields:


```r
sum(sumblast4$seq_count[sumblast4$rfam_biotype == "5S_rRNA"])
```

```
## [1] 1341
```

```r
sum(sumblast4$seq_count[sumblast4$rfam_biotype == "tRNA"])
```

```
## [1] 124439
```

```r
if (sum(sumblast4$seq_count[sumblast4$rfam_biotype == "5S_rRNA"]) != totalsumbiotype4["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast4$seq_count[sumblast4$rfam_biotype == "tRNA"]) != totalsumbiotype4["tRNA",]) stop ("tRNA counts not equal")
```

## Save data


```r
save(sumblast4, uniqsumblast4, totalsumbiotype4, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/8_musmusculus_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast4, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/8_1_musmusculus_filtered_uniqseq_rfam_blastn_results.csv", col.names=TRUE)
```

```
## Warning in write.csv(uniqsumblast4, file = "/mnt/research/pigeqtl/analyses/
## microRNA/2_mirna_characterization_expression/0_rfam_database_query/
## 8_1_musmusculus_filtered_uniqseq_rfam_blastn_results.csv", : attempt to set
## 'col.names' ignored
```

```r
write.csv(totalsumbiotype4, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/8_2_musmusculus_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)
```

