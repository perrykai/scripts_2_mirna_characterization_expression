**Script:** ``

**Directory of Code:**  ``

**Date:**  

**Input File Directory:**  ``

**Input File(s):** ``

**Output File Directory:** ``

**Output File(s):** ``

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives
## Install libraries
## Load data
## Analysis
## Visualize
## Save data


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/")

library(plyr)

system.time(
rfamtblout<-readLines("../infernal_rfam_query_output/174_split_collapsed_reads_51_tblout.txt")
)
```

```
##    user  system elapsed 
##   0.307   0.009   0.837
```

```r
length(rfamtblout)
```

```
## [1] 48276
```

Eliminate row the last 9 lines (summary info, not useful)


```r
rfamtblout<-rfamtblout[-c(1,2,((length(rfamtblout)-9):length(rfamtblout)))]

length(rfamtblout)
```

```
## [1] 48264
```

```r
str(rfamtblout)
```

```
##  chr [1:48264] "b55                  RF01783   seq_232446615_x1     -          cm        1       12        7       18      +    3'    3 0.75   "| __truncated__ ...
```

```r
head(rfamtblout)
```

```
## [1] "b55                  RF01783   seq_232446615_x1     -          cm        1       12        7       18      +    3'    3 0.75   0.0    9.3       1.3 ?   -"
## [2] "TB11Cs2H1            RF01537   seq_232446622_x1     -          cm       51       70        1       20      +    5'    2 0.55   0.0    8.8       2.1 ?   -"
## [3] "sR40                 RF01297   seq_232446624_x1     -          cm       40       61        1       21      +    5'    2 0.52   0.0    9.2      0.55 ?   -"
## [4] "SAM-SAH              RF01727   seq_232446628_x1     -          cm       14       21        5       12      +    no    4 0.62   0.0    5.7       2.7 ?   -"
## [5] "mir-1                RF00103   seq_232446634_x1     -          cm       51       70        1       20      +    5'    4 0.25   0.0   17.3     0.013 ?   -"
## [6] "MIR1846              RF02000   seq_232446637_x1     -          cm       24       47        1       24      +    3'    4 0.75   0.0    9.8       1.8 ?   -"
```

```r
tail(rfamtblout)
```

```
## [1] "psbNH                RF01753   seq_232573387_x1     -          cm        1       18       12       29      +    3'    3 0.39   0.0    6.3       2.3 ?   -"
## [2] "UPSK                 RF00390   seq_232573390_x1     -          cm        1       11       11       21      +    3'    3 0.45   0.0    6.5       6.9 ?   -"
## [3] "CRISPR-DR38          RF01348   seq_232573391_x1     -          cm       22       35        7       20      +    5'    2 0.36   0.0    7.4       5.2 ?   -"
## [4] "CRISPR-DR36          RF01346   seq_232573391_x1     -          cm       14       36        3       20      +    5'    2 0.39   0.0    4.7       7.1 ?   -"
## [5] "Gl_U4                RF02493   seq_232573393_x1     -          cm        1       18        3       20      +    3'    3 0.61   0.0    7.9       2.4 ?   -"
## [6] "neisseria_FSE        RF01843   seq_232573394_x1     -          cm       29       40        1       12      +    5'    2 0.50   0.0    8.1       4.7 ?   -"
```

```r
system.time(
rfamsummary<-ldply(rfamtblout, function(x) (read.table(textConnection(x), colClasses = c(rep("character", 5), rep("numeric", 4), rep("character", 2), rep("numeric", 5), "character"))[1:17]))
)
```

```
##    user  system elapsed 
##  54.167   0.100  54.226
```

```r
colnames(rfamsummary) <- c("target_name", "accession_target", "query_name", "accession_query", "mdl", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc", "pass", "gc", "bias", "score", "E_value", "inc")

head(rfamsummary)
```

```
##   target_name accession_target       query_name accession_query mdl
## 1         b55          RF01783 seq_232446615_x1               -  cm
## 2   TB11Cs2H1          RF01537 seq_232446622_x1               -  cm
## 3        sR40          RF01297 seq_232446624_x1               -  cm
## 4     SAM-SAH          RF01727 seq_232446628_x1               -  cm
## 5       mir-1          RF00103 seq_232446634_x1               -  cm
## 6     MIR1846          RF02000 seq_232446637_x1               -  cm
##   mdl_from mdl_to seq_from seq_to strand trunc pass   gc bias score
## 1        1     12        7     18      +    3'    3 0.75    0   9.3
## 2       51     70        1     20      +    5'    2 0.55    0   8.8
## 3       40     61        1     21      +    5'    2 0.52    0   9.2
## 4       14     21        5     12      +    no    4 0.62    0   5.7
## 5       51     70        1     20      +    5'    4 0.25    0  17.3
## 6       24     47        1     24      +    3'    4 0.75    0   9.8
##   E_value inc
## 1   1.300   ?
## 2   2.100   ?
## 3   0.550   ?
## 4   2.700   ?
## 5   0.013   ?
## 6   1.800   ?
```

```r
tail(rfamsummary)
```

```
##         target_name accession_target       query_name accession_query mdl
## 48259         psbNH          RF01753 seq_232573387_x1               -  cm
## 48260          UPSK          RF00390 seq_232573390_x1               -  cm
## 48261   CRISPR-DR38          RF01348 seq_232573391_x1               -  cm
## 48262   CRISPR-DR36          RF01346 seq_232573391_x1               -  cm
## 48263         Gl_U4          RF02493 seq_232573393_x1               -  cm
## 48264 neisseria_FSE          RF01843 seq_232573394_x1               -  cm
##       mdl_from mdl_to seq_from seq_to strand trunc pass   gc bias score
## 48259        1     18       12     29      +    3'    3 0.39    0   6.3
## 48260        1     11       11     21      +    3'    3 0.45    0   6.5
## 48261       22     35        7     20      +    5'    2 0.36    0   7.4
## 48262       14     36        3     20      +    5'    2 0.39    0   4.7
## 48263        1     18        3     20      +    3'    3 0.61    0   7.9
## 48264       29     40        1     12      +    5'    2 0.50    0   8.1
##       E_value inc
## 48259     2.3   ?
## 48260     6.9   ?
## 48261     5.2   ?
## 48262     7.1   ?
## 48263     2.4   ?
## 48264     4.7   ?
```

```r
str(rfamsummary)
```

```
## 'data.frame':	48264 obs. of  17 variables:
##  $ target_name     : chr  "b55" "TB11Cs2H1" "sR40" "SAM-SAH" ...
##  $ accession_target: chr  "RF01783" "RF01537" "RF01297" "RF01727" ...
##  $ query_name      : chr  "seq_232446615_x1" "seq_232446622_x1" "seq_232446624_x1" "seq_232446628_x1" ...
##  $ accession_query : chr  "-" "-" "-" "-" ...
##  $ mdl             : chr  "cm" "cm" "cm" "cm" ...
##  $ mdl_from        : num  1 51 40 14 51 24 83 1 1 1 ...
##  $ mdl_to          : num  12 70 61 21 70 47 101 34 12 15 ...
##  $ seq_from        : num  7 1 1 5 1 1 12 7 16 3 ...
##  $ seq_to          : num  18 20 21 12 20 24 30 29 27 17 ...
##  $ strand          : chr  "+" "+" "+" "+" ...
##  $ trunc           : chr  "3'" "5'" "5'" "no" ...
##  $ pass            : num  3 2 2 4 4 4 4 3 3 3 ...
##  $ gc              : num  0.75 0.55 0.52 0.62 0.25 0.75 0.37 0.61 0.75 0.6 ...
##  $ bias            : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ score           : num  9.3 8.8 9.2 5.7 17.3 9.8 5 9.6 6.5 5.6 ...
##  $ E_value         : num  1.3 2.1 0.55 2.7 0.013 1.8 2.5 1.5 7.7 4.9 ...
##  $ inc             : chr  "?" "?" "?" "?" ...
```

```r
save(rfamsummary, file = "../summary_infernal_rfam/test_rfam_summary_table.Rdata")
```

