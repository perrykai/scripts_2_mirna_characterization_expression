**File:** `1_extract_mirna_info_from_core_output.R`

**Directory Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts/`

**Date:**  1/26/16

**Description:**  This code extracts the miRDeep2 core module output from the `result_21_01_2016_t_20_01_10.csv` file
                This .csv file contains the miRDeep2 score distribution for the novel miRNA prediction step, 
                the novel miRNAs predicted by miRDeep2, the mature miRBase miRNAs detected by miRDeep2, and the 
                miRBase miRNAs not detected by miRDeep2. The first three of these items are extracted from this .csv file using this script.
                The objective here is to characterize the known and novel miRNAs present in this dataset, isolate the sequences meeting miRDeep2 score of 10 or greater,
                to isolate the sequences at that score cutoff having homologous seed sequence with a human miRBase miRNA,
                and to estimate the false discovery rate of the miRDeep2 prediction step at miRDeep2 score 10. 

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output`

**Input File(s):**  `result_21_01_2016_t_20_01_10.csv`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_mirna_info_extracted_core_output/`

**Output File(s):** 

-  `1_extracted_mirdeep2_core_output_stats.csv`
-  `2_extracted_mirdeep2_core_predicted_novel_mirna.csv`
-  `3_extracted_mirdeep2_core_mature_mirna_detected.csv`

To isolate the first section of the csv containing the miRDeep2 distribution scores:


```r
sts<-read.table("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/result_21_01_2016_t_20_01_10.csv", sep = "\t", nrows=21, header = TRUE)

head(kable(sts))
```

```
## [1] "| miRDeep2.score| novel.miRNAs.reported.by.miRDeep2|novel.miRNAs..estimated.false.positives |novel.miRNAs..estimated.true.positives | known.miRNAs.in.species| known.miRNAs.in.data|known.miRNAs.detected.by.miRDeep2 | estimated.signal.to.noise| excision.gearing|"
## [2] "|--------------:|---------------------------------:|:---------------------------------------|:--------------------------------------|-----------------------:|--------------------:|:---------------------------------|-------------------------:|----------------:|"
## [3] "|             10|                               352|32 +/- 6                                |320 +/- 6 (91 +/- 2%)                  |                     411|                  316|259 (82%)                         |                      15.7|                4|"
## [4] "|              9|                               360|33 +/- 6                                |327 +/- 6 (91 +/- 2%)                  |                     411|                  316|259 (82%)                         |                      15.6|                4|"
## [5] "|              8|                               379|34 +/- 6                                |345 +/- 6 (91 +/- 2%)                  |                     411|                  316|259 (82%)                         |                      15.7|                4|"
## [6] "|              7|                               397|35 +/- 6                                |362 +/- 6 (91 +/- 2%)                  |                     411|                  316|259 (82%)                         |                      15.7|                4|"
```

```r
write.table(sts, file="../1_extracted_mirdeep2_core_output_stats.csv", sep = "\t", col.names = TRUE)
```

To isolate the second section of the csv containing the novel predicted miRNAs:


```r
cv<-read.table("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/result_21_01_2016_t_20_01_10.csv", sep = "\t", skip=26, nrow=1199,  header = TRUE, fill=TRUE)

head(kable(cv))
```

```
## [1] "|provisional.id   | miRDeep2.score|estimated.probability.that.the.miRNA.candidate.is.a.true.positive |rfam.alert | total.read.count| mature.read.count| loop.read.count| star.read.count|significant.randfold.p.value |miRBase.miRNA |example.miRBase.miRNA.with.the.same.seed |UCSC.browser |NCBI.blastn |consensus.mature.sequence |consensus.star.sequence        |consensus.precursor.sequence                                                                    |precursor.coordinate        |"
## [2] "|:----------------|--------------:|:-----------------------------------------------------------------|:----------|----------------:|-----------------:|---------------:|---------------:|:----------------------------|:-------------|:----------------------------------------|:------------|:-----------|:-------------------------|:------------------------------|:-----------------------------------------------------------------------------------------------|:---------------------------|"
## [3] "|GL892871.2_43622 |       572746.3|91 +/- 2%                                                         |-          |          1123410|           1119829|               0|            3581|yes                          |-             |hsa-miR-26a-5p                           |-            |-           |uucaaguaauucaggauagguu    |ccuguucuccauuacuuggcuc         |uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc                                       |GL892871.2:56019..56076:+   |"
## [4] "|X_39130          |       136780.7|91 +/- 2%                                                         |-          |           268282|            268080|               0|             202|yes                          |-             |hsa-miR-660-5p                           |-            |-           |uacccauugcauaucggaguug    |accuccuaugugcaugguuuac         |uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac                                      |X:48640826..48640884:+      |"
## [5] "|GL896425.1_44884 |       128796.6|91 +/- 2%                                                         |-          |           252627|            249026|               0|            3601|yes                          |-             |-                                        |-            |-           |uggugccugacgucuuggcagu    |agccagggcugcaggcacugaca        |agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu                                   |GL896425.1:1750..1811:-     |"
## [6] "|16_36270         |       113215.6|91 +/- 2%                                                         |-          |           222068|            222028|               0|              40|no                           |-             |hsa-let-7a-5p                            |-            |-           |ugagguaguaggcugugugg      |aauagcaacaacaacaaca            |aauagcaacaacaacaacaacaaaagaaugagguaguaggcugugugg                                                |16:6629237..6629285:-       |"
```

```r
write.table(cv, file="../2_extracted_mirdeep2_core_predicted_novel_mirna.csv", sep = "\t", col.names=TRUE)
```

To isolate the third section of the csv containing the miRBase miRNAs detected by miRDeep2:


```r
md<-read.table("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/result_21_01_2016_t_20_01_10.csv", sep = "\t", skip = 1230, nrow = 306, header = TRUE, fill = TRUE)

head(kable(md))
```

```
## [1] "|tag.id           | miRDeep2.score|estimated.probability.that.the.miRNA.is.a.true.positive |rfam.alert | total.read.count| mature.read.count| loop.read.count| star.read.count|significant.randfold.p.value |mature.miRBase.miRNA |example.miRBase.miRNA.with.the.same.seed |UCSC.browser |NCBI.blastn |consensus.mature.sequence |consensus.star.sequence    |consensus.precursor.sequence                                                                |precursor.coordinate        |"
## [2] "|:----------------|--------------:|:-------------------------------------------------------|:----------|----------------:|-----------------:|---------------:|---------------:|:----------------------------|:--------------------|:----------------------------------------|:------------|:-----------|:-------------------------|:--------------------------|:-------------------------------------------------------------------------------------------|:---------------------------|"
## [3] "|6_15007          |     21446337.3|91 +/- 2%                                               |-          |         42066058|          42042148|              96|           23814|yes                          |ssc-miR-1            |hsa-miR-1-3p                             |-            |-           |uggaauguaaagaaguauguau    |acauacuucuuuauguacccaua    |acauacuucuuuauguacccauaugaacauacaaugcuauggaauguaaagaaguauguau                               |6:99481945..99482006:+      |"
## [4] "|17_37750         |     21435035.7|91 +/- 2%                                               |-          |         42043890|          42032429|            1981|            9480|yes                          |ssc-miR-1            |hsa-miR-1-3p                             |-            |-           |uggaauguaaagaaguauguau    |acauacuucuuuaugugcccaua    |acauacuucuuuaugugcccauauggaccugcuaagcuauggaauguaaagaaguauguau                               |17:69285415..69285476:-     |"
## [5] "|6_15009          |      4844835.2|91 +/- 2%                                               |-          |          9502934|           9103074|             448|          399412|yes                          |ssc-miR-133a-5p      |-                                        |-            |-           |uugguccccuucaaccagcugu    |agcugguaaaauggaaccaaau     |agcugguaaaauggaaccaaaucgccucuucaauggauuugguccccuucaaccagcugu                                |6:99485215..99485275:+      |"
## [6] "|17_37738         |      4841209.8|91 +/- 2%                                               |-          |          9495823|           9103046|             282|          392495|yes                          |ssc-miR-133a-5p      |-                                        |-            |-           |uugguccccuucaaccagcugu    |agcugguaaaauggaaccaaau     |agcugguaaaauggaaccaaaucaacuguugaauggauuugguccccuucaaccagcugu                                |17:69274673..69274733:-     |"
```

```r
write.table(md, file = "../3_extracted_mirdeep2_core_mature_mirna_detected.csv", sep = "\t", col.names=TRUE)
```

Now I can open the `novel_mirna_predicted.csv` file and filter by various thresholds

### 1. miRDeep2 score of 10 or more


```r
novelmirna<-read.table("../2_extracted_mirdeep2_core_predicted_novel_mirna.csv")

colnames(novelmirna)
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
dim(novelmirna)	
```

```
## [1] 1199   17
```

```r
head(kable(novelmirna))
```

```
## [1] "|provisional.id   | miRDeep2.score|estimated.probability.that.the.miRNA.candidate.is.a.true.positive |rfam.alert | total.read.count| mature.read.count| loop.read.count| star.read.count|significant.randfold.p.value |miRBase.miRNA |example.miRBase.miRNA.with.the.same.seed |UCSC.browser |NCBI.blastn |consensus.mature.sequence |consensus.star.sequence        |consensus.precursor.sequence                                                                    |precursor.coordinate        |"
## [2] "|:----------------|--------------:|:-----------------------------------------------------------------|:----------|----------------:|-----------------:|---------------:|---------------:|:----------------------------|:-------------|:----------------------------------------|:------------|:-----------|:-------------------------|:------------------------------|:-----------------------------------------------------------------------------------------------|:---------------------------|"
## [3] "|GL892871.2_43622 |       572746.3|91 +/- 2%                                                         |-          |          1123410|           1119829|               0|            3581|yes                          |-             |hsa-miR-26a-5p                           |-            |-           |uucaaguaauucaggauagguu    |ccuguucuccauuacuuggcuc         |uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc                                       |GL892871.2:56019..56076:+   |"
## [4] "|X_39130          |       136780.7|91 +/- 2%                                                         |-          |           268282|            268080|               0|             202|yes                          |-             |hsa-miR-660-5p                           |-            |-           |uacccauugcauaucggaguug    |accuccuaugugcaugguuuac         |uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac                                      |X:48640826..48640884:+      |"
## [5] "|GL896425.1_44884 |       128796.6|91 +/- 2%                                                         |-          |           252627|            249026|               0|            3601|yes                          |-             |-                                        |-            |-           |uggugccugacgucuuggcagu    |agccagggcugcaggcacugaca        |agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu                                   |GL896425.1:1750..1811:-     |"
## [6] "|16_36270         |       113215.6|91 +/- 2%                                                         |-          |           222068|            222028|               0|              40|no                           |-             |hsa-let-7a-5p                            |-            |-           |ugagguaguaggcugugugg      |aauagcaacaacaacaaca            |aauagcaacaacaacaacaacaaaagaaugagguaguaggcugugugg                                                |16:6629237..6629285:-       |"
```

```r
score10novelmirna<-novelmirna[novelmirna$miRDeep2.score >= 10, ]

dim(score10novelmirna)
```

```
## [1] 352  17
```

```r
nrow(score10novelmirna)
```

```
## [1] 352
```

So, there are 352 predicted miRNA precursors with a miRDeep2 score > or = 10

Estimated false positives is 32 +/- 6 (obtained from `1_extracted_mirdeep2_core_output_stats.csv`)


```r
32 - 6
```

```
## [1] 26
```

```r
32 + 6
```

```
## [1] 38
```

```r
26/352
```

```
## [1] 0.07386364
```

```r
38/352
```

```
## [1] 0.1079545
```

### 2. Significant Randfold p-value

Now, subset this again into those that had a significant Randfold p value, indicating ability of secondary structure formation


```r
head(score10novelmirna$significant.randfold.p.value)
```

```
## [1] yes yes yes no  yes yes
## Levels: no yes
```

```r
sum(score10novelmirna$significant.randfold.p.value=="yes")
```

```
## [1] 271
```

```r
randfoldsigpval<-score10novelmirna[score10novelmirna$significant.randfold.p.value == "yes", ]
dim(randfoldsigpval)
```

```
## [1] 271  17
```

Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences


```r
rfam<-randfoldsigpval$rfam.alert
sum(rfam=="-")
```

```
## [1] 271
```

```r
sum(rfam !="-")
```

```
## [1] 0
```

### 3. Do the putative novel sequences have a homologous human miRNA seed sequence?


```r
sum(randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-")
```

```
## [1] 100
```

This indicates that 100 of the sequences have a homologous human seed sequence


```r
homologseed<-randfoldsigpval[randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-",]
```

Subset of the 100 sequences


```r
dim(homologseed)
```

```
## [1] 100  17
```

```r
write.table(homologseed, "../4_extracted_predicted_novel_sigranfold_homologseed.txt", sep = "\t", col.names=TRUE)
```

