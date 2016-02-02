**Script:** `1_create_mirna_expression_df.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/scripts`

**Date:**  1/26/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output`

**Input File(s):** `miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna`

**Output File(s):** `1_mean_mature_mirna_expression_unfiltered.txt`

                    `1_mean_mature_mirna_expression_unfiltered.Rdata`

                    `2_mean_mature_mirna_expression_filtered.txt`

                    `2_mean_mature_mirna_expression_filtered.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives
The objective of this script is to create a matrix of miRNA expression data (read counts) for incorporation into the miRNA eQTL analysis, beginning with the known miRNA expression data output from the miRDeep2 core module.
Later, expression data of the putative novel miRNA candidates can also be included, also output from the miRDeep2 core module. 
To acheive one read count per miRNA per animal, I will need to take the average of the mature read counts in instances of multiple precursor sequences.

The result of this script will be two data frames containing one average read count per miRNA per animal: one unfiltered, and one filtered for miRNAs expressed greater than the number of animals in the population (174) and transposed.

So, what I need to do:

1. Extract the columns of the mature miRNA read counts for each animal
2. Use the 'by' function to apply the function colmeans to the entire data frame of read counts (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that group of miRNA. The result of this will be a list containing the average read counts for each miRNA for each animal)
3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function
4. Filter the data for expression threshold: The read count for the miRNA needs to be greater than the number of animals in the population
5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR
6. Transpose the data.frame to make the animals the rows and the miRNAs the columns

## Install libraries


```r
library(plyr)
```

## Load data
###1. Read in the csv of the miRNA expression profiles, using check.names default to correct column names with non-ASCII characters


```r
rc<-read.csv("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv", sep = "\t", header = TRUE, row.names=NULL)
```

Set the name of the first column to "miRNA":


```r
colnames(rc)[[1]]<-"miRNA"
```

Remove the "X" character from the beginning of each column:


```r
colnames(rc)<-gsub("X", "", colnames(rc))
colnames(rc)
```

```
##   [1] "miRNA"      "read_count" "precursor"  "total"      "001"       
##   [6] "002"        "003"        "004"        "005"        "006"       
##  [11] "007"        "008"        "009"        "010"        "011"       
##  [16] "012"        "013"        "014"        "015"        "016"       
##  [21] "017"        "018"        "019"        "020"        "021"       
##  [26] "022"        "023"        "024"        "025"        "026"       
##  [31] "027"        "028"        "029"        "030"        "031"       
##  [36] "032"        "033"        "034"        "035"        "036"       
##  [41] "037"        "038"        "039"        "040"        "041"       
##  [46] "042"        "043"        "044"        "045"        "046"       
##  [51] "047"        "048"        "049"        "050"        "051"       
##  [56] "052"        "053"        "054"        "055"        "056"       
##  [61] "057"        "058"        "059"        "060"        "061"       
##  [66] "062"        "063"        "064"        "065"        "066"       
##  [71] "067"        "068"        "069"        "070"        "071"       
##  [76] "072"        "073"        "074"        "075"        "076"       
##  [81] "077"        "078"        "079"        "080"        "081"       
##  [86] "082"        "083"        "084"        "085"        "086"       
##  [91] "087"        "088"        "089"        "090"        "091"       
##  [96] "092"        "093"        "094"        "095"        "096"       
## [101] "097"        "098"        "099"        "100"        "101"       
## [106] "102"        "103"        "104"        "105"        "106"       
## [111] "107"        "108"        "109"        "110"        "111"       
## [116] "112"        "113"        "114"        "115"        "116"       
## [121] "117"        "118"        "119"        "120"        "121"       
## [126] "122"        "123"        "124"        "125"        "126"       
## [131] "127"        "128"        "129"        "130"        "131"       
## [136] "132"        "133"        "134"        "135"        "136"       
## [141] "137"        "138"        "139"        "140"        "141"       
## [146] "142"        "143"        "144"        "145"        "146"       
## [151] "147"        "148"        "149"        "150"        "151"       
## [156] "152"        "153"        "154"        "155"        "156"       
## [161] "157"        "158"        "159"        "160"        "161"       
## [166] "162"        "163"        "164"        "165"        "166"       
## [171] "167"        "168"        "169"        "170"        "171"       
## [176] "172"        "173"        "174"        "001.norm."  "002.norm." 
## [181] "003.norm."  "004.norm."  "005.norm."  "006.norm."  "007.norm." 
## [186] "008.norm."  "009.norm."  "010.norm."  "011.norm."  "012.norm." 
## [191] "013.norm."  "014.norm."  "015.norm."  "016.norm."  "017.norm." 
## [196] "018.norm."  "019.norm."  "020.norm."  "021.norm."  "022.norm." 
## [201] "023.norm."  "024.norm."  "025.norm."  "026.norm."  "027.norm." 
## [206] "028.norm."  "029.norm."  "030.norm."  "031.norm."  "032.norm." 
## [211] "033.norm."  "034.norm."  "035.norm."  "036.norm."  "037.norm." 
## [216] "038.norm."  "039.norm."  "040.norm."  "041.norm."  "042.norm." 
## [221] "043.norm."  "044.norm."  "045.norm."  "046.norm."  "047.norm." 
## [226] "048.norm."  "049.norm."  "050.norm."  "051.norm."  "052.norm." 
## [231] "053.norm."  "054.norm."  "055.norm."  "056.norm."  "057.norm." 
## [236] "058.norm."  "059.norm."  "060.norm."  "061.norm."  "062.norm." 
## [241] "063.norm."  "064.norm."  "065.norm."  "066.norm."  "067.norm." 
## [246] "068.norm."  "069.norm."  "070.norm."  "071.norm."  "072.norm." 
## [251] "073.norm."  "074.norm."  "075.norm."  "076.norm."  "077.norm." 
## [256] "078.norm."  "079.norm."  "080.norm."  "081.norm."  "082.norm." 
## [261] "083.norm."  "084.norm."  "085.norm."  "086.norm."  "087.norm." 
## [266] "088.norm."  "089.norm."  "090.norm."  "091.norm."  "092.norm." 
## [271] "093.norm."  "094.norm."  "095.norm."  "096.norm."  "097.norm." 
## [276] "098.norm."  "099.norm."  "100.norm."  "101.norm."  "102.norm." 
## [281] "103.norm."  "104.norm."  "105.norm."  "106.norm."  "107.norm." 
## [286] "108.norm."  "109.norm."  "110.norm."  "111.norm."  "112.norm." 
## [291] "113.norm."  "114.norm."  "115.norm."  "116.norm."  "117.norm." 
## [296] "118.norm."  "119.norm."  "120.norm."  "121.norm."  "122.norm." 
## [301] "123.norm."  "124.norm."  "125.norm."  "126.norm."  "127.norm." 
## [306] "128.norm."  "129.norm."  "130.norm."  "131.norm."  "132.norm." 
## [311] "133.norm."  "134.norm."  "135.norm."  "136.norm."  "137.norm." 
## [316] "138.norm."  "139.norm."  "140.norm."  "141.norm."  "142.norm." 
## [321] "143.norm."  "144.norm."  "145.norm."  "146.norm."  "147.norm." 
## [326] "148.norm."  "149.norm."  "150.norm."  "151.norm."  "152.norm." 
## [331] "153.norm."  "154.norm."  "155.norm."  "156.norm."  "157.norm." 
## [336] "158.norm."  "159.norm."  "160.norm."  "161.norm."  "162.norm." 
## [341] "163.norm."  "164.norm."  "165.norm."  "166.norm."  "167.norm." 
## [346] "168.norm."  "169.norm."  "170.norm."  "171.norm."  "172.norm." 
## [351] "173.norm."  "174.norm."
```

View the data set:


```r
head(rc[1:8])
```

```
##           miRNA read_count    precursor   total   001   002   003   004
## 1    ssc-let-7a    5066155 ssc-let-7a-1 5066155 48170 23457 28447 29899
## 2    ssc-let-7a    5055934 ssc-let-7a-2 5055934 48095 23397 28448 29821
## 3    ssc-let-7c    3238382   ssc-let-7c 3238382 32745 14987 18144 18681
## 4 ssc-let-7d-5p     490052   ssc-let-7d  490052  4925  1938  2511  3076
## 5 ssc-let-7d-3p      51589   ssc-let-7d   51589   381   192   198   269
## 6    ssc-let-7e     274145   ssc-let-7e  274145  2811  1302  1463  1690
```

###2. Read in the config file, maintaining the characters in the 3-digit code names


```r
configfile<-read.table("../../../1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/config_for_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)
```

```
##                V1  V2
## 1 1034_express.fa 001
## 2 1036_express.fa 002
## 3 1041_express.fa 003
## 4 1049_express.fa 004
## 5 1058_express.fa 005
## 6 1060_express.fa 006
```

Remove the "_express.fa" from each file name to leave the pig id:


```r
configfile$V1<-gsub("_express.fa", "", configfile$V1)
```

Make filenames more informative:


```r
colnames(configfile)<-c("pigid","code")
colnames(configfile)
```

```
## [1] "pigid" "code"
```

## Analysis
###1. Extract the columns of mature read counts for each miRNA for each animal


```r
mirquant<-rc[,c(1,5:178)]
colnames(mirquant)
```

```
##   [1] "miRNA" "001"   "002"   "003"   "004"   "005"   "006"   "007"  
##   [9] "008"   "009"   "010"   "011"   "012"   "013"   "014"   "015"  
##  [17] "016"   "017"   "018"   "019"   "020"   "021"   "022"   "023"  
##  [25] "024"   "025"   "026"   "027"   "028"   "029"   "030"   "031"  
##  [33] "032"   "033"   "034"   "035"   "036"   "037"   "038"   "039"  
##  [41] "040"   "041"   "042"   "043"   "044"   "045"   "046"   "047"  
##  [49] "048"   "049"   "050"   "051"   "052"   "053"   "054"   "055"  
##  [57] "056"   "057"   "058"   "059"   "060"   "061"   "062"   "063"  
##  [65] "064"   "065"   "066"   "067"   "068"   "069"   "070"   "071"  
##  [73] "072"   "073"   "074"   "075"   "076"   "077"   "078"   "079"  
##  [81] "080"   "081"   "082"   "083"   "084"   "085"   "086"   "087"  
##  [89] "088"   "089"   "090"   "091"   "092"   "093"   "094"   "095"  
##  [97] "096"   "097"   "098"   "099"   "100"   "101"   "102"   "103"  
## [105] "104"   "105"   "106"   "107"   "108"   "109"   "110"   "111"  
## [113] "112"   "113"   "114"   "115"   "116"   "117"   "118"   "119"  
## [121] "120"   "121"   "122"   "123"   "124"   "125"   "126"   "127"  
## [129] "128"   "129"   "130"   "131"   "132"   "133"   "134"   "135"  
## [137] "136"   "137"   "138"   "139"   "140"   "141"   "142"   "143"  
## [145] "144"   "145"   "146"   "147"   "148"   "149"   "150"   "151"  
## [153] "152"   "153"   "154"   "155"   "156"   "157"   "158"   "159"  
## [161] "160"   "161"   "162"   "163"   "164"   "165"   "166"   "167"  
## [169] "168"   "169"   "170"   "171"   "172"   "173"   "174"
```

```r
dim(mirquant)
```

```
## [1] 492 175
```

```r
head(mirquant[1:8])
```

```
##           miRNA   001   002   003   004   005   006   007
## 1    ssc-let-7a 48170 23457 28447 29899 40817 29774 38804
## 2    ssc-let-7a 48095 23397 28448 29821 40698 29725 38794
## 3    ssc-let-7c 32745 14987 18144 18681 34313 15294 28022
## 4 ssc-let-7d-5p  4925  1938  2511  3076  3472  3276  3705
## 5 ssc-let-7d-3p   381   192   198   269   778   239   774
## 6    ssc-let-7e  2811  1302  1463  1690  2512  1410  2229
```

Take a subset of this data.frame for testing:


```r
test<-mirquant[1:20,1:8]
head(test)
```

```
##           miRNA   001   002   003   004   005   006   007
## 1    ssc-let-7a 48170 23457 28447 29899 40817 29774 38804
## 2    ssc-let-7a 48095 23397 28448 29821 40698 29725 38794
## 3    ssc-let-7c 32745 14987 18144 18681 34313 15294 28022
## 4 ssc-let-7d-5p  4925  1938  2511  3076  3472  3276  3705
## 5 ssc-let-7d-3p   381   192   198   269   778   239   774
## 6    ssc-let-7e  2811  1302  1463  1690  2512  1410  2229
```

###2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts
(What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
The result of this will be a list containing the average read counts for each miRNA for each animal

Example: by(data, index, function)


```r
head(by(test[,2:ncol(test)], test[,1], colMeans))
```

```
## $`ssc-let-7a`
##     001     002     003     004     005     006     007 
## 48132.5 23427.0 28447.5 29860.0 40757.5 29749.5 38799.0 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007 
## 32745 14987 18144 18681 34313 15294 28022 
## 
## $`ssc-let-7d-3p`
## 001 002 003 004 005 006 007 
## 381 192 198 269 778 239 774 
## 
## $`ssc-let-7d-5p`
##  001  002  003  004  005  006  007 
## 4925 1938 2511 3076 3472 3276 3705 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007 
## 2811 1302 1463 1690 2512 1410 2229 
## 
## $`ssc-let-7f`
##     001     002     003     004     005     006     007 
## 35432.5 16422.0 20161.5 22985.5 16512.5 24475.5 19175.0
```

Apply the by function to the full dataframe:


```r
meanrc<-by(mirquant[,2:ncol(mirquant)], mirquant[,1], colMeans)
```

This should be 411, the number of mature pig miRNAs in miRBase:


```r
length(meanrc)
```

```
## [1] 411
```

```r
head(meanrc)
```

```
## $`ssc-let-7a`
##     001     002     003     004     005     006     007     008     009 
## 48132.5 23427.0 28447.5 29860.0 40757.5 29749.5 38799.0 31459.5 29857.0 
##     010     011     012     013     014     015     016     017     018 
## 21333.0 35977.0 21304.0 43336.5 46963.0 41826.5 41174.0 36678.5 45540.0 
##     019     020     021     022     023     024     025     026     027 
## 59997.0 34675.0 36866.0 34397.5 42057.0 49982.0 35410.5  2366.0 25973.5 
##     028     029     030     031     032     033     034     035     036 
## 41518.5 42680.0 30864.0 44394.0 43819.0 38598.5 50949.5 46225.0 47411.0 
##     037     038     039     040     041     042     043     044     045 
## 35912.5 37383.0 24522.0 35894.5 15876.5 42837.0 36219.0 30672.0 28829.0 
##     046     047     048     049     050     051     052     053     054 
## 28805.0 24042.0 30291.5 25858.0 31424.0 10909.0 33168.0 23342.0 28218.0 
##     055     056     057     058     059     060     061     062     063 
## 17887.5 15043.5 16239.5 22696.0 13045.5 28524.0 17197.0 32303.5 28680.0 
##     064     065     066     067     068     069     070     071     072 
## 37439.0 47521.5 42747.0 41810.5 20253.0 28069.0 47363.0 29313.0 26772.0 
##     073     074     075     076     077     078     079     080     081 
## 35965.5 27357.5 42022.0 44591.0 48267.5 35515.5 30704.5 37564.5 38083.5 
##     082     083     084     085     086     087     088     089     090 
## 36868.5 40084.0 35862.5 23067.5 37179.0 38866.5 29532.5 28213.0 39304.5 
##     091     092     093     094     095     096     097     098     099 
## 37340.5 33901.5 25400.0 28632.0 26621.0 23326.0 27332.5 29679.5 24334.5 
##     100     101     102     103     104     105     106     107     108 
## 25708.0 29901.0 30335.0 23033.5 23954.0 17311.0 19595.5 28878.5 22962.5 
##     109     110     111     112     113     114     115     116     117 
## 18231.0 25629.5 17310.5 14579.0 13451.5 14827.5 35094.0 22924.0 22506.5 
##     118     119     120     121     122     123     124     125     126 
## 21439.0 19489.5 22782.0 24263.5 21087.0 23661.5 25093.0 25866.0 20312.5 
##     127     128     129     130     131     132     133     134     135 
## 23853.0 28593.0 28171.0 27868.5 23409.5 21416.0 18956.5 29712.0 20679.0 
##     136     137     138     139     140     141     142     143     144 
## 22174.0 23523.5 26256.5 23873.5 23470.0 16744.5 23586.5 32558.0 23146.5 
##     145     146     147     148     149     150     151     152     153 
## 13914.0 15831.0 30359.5 21163.5 17806.5 32457.5 28573.0 26978.5 15400.5 
##     154     155     156     157     158     159     160     161     162 
## 20649.0 20217.0 34496.0 27724.0 32625.0 27000.0 31787.0 25896.0 24638.0 
##     163     164     165     166     167     168     169     170     171 
## 22857.0 23807.0 19823.5 28145.5 34188.5 29473.5 16803.0 24594.5 21326.0 
##     172     173     174 
## 19700.5 20084.0 18977.0 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007   008   009   010   011   012 
## 32745 14987 18144 18681 34313 15294 28022 19512 16741 12522 20772 12449 
##   013   014   015   016   017   018   019   020   021   022   023   024 
## 26354 28968 28788 29017 29022 31957 45232 23212 23731 21229 29962 37824 
##   025   026   027   028   029   030   031   032   033   034   035   036 
## 23684  3077 15389 28463 30503 17962 31221 31310 25291 37902 32444 32430 
##   037   038   039   040   041   042   043   044   045   046   047   048 
## 24015 25399 14477 25063 14181 31924 24089 20474 21753 18148 14466 18737 
##   049   050   051   052   053   054   055   056   057   058   059   060 
## 16355 19779  5973 20260 14482 18147 10798  8526  9964 15666  7201 17704 
##   061   062   063   064   065   066   067   068   069   070   071   072 
##  9342 24491 20807 26015 28812 26738 26126 17415 16928 34465 17004 15809 
##   073   074   075   076   077   078   079   080   081   082   083   084 
## 22883 16277 27810 31228 34068 25721 18956 24904 24512 22762 25993 24232 
##   085   086   087   088   089   090   091   092   093   094   095   096 
## 14980 19671 24475 18429 18530 25909 23004 20516 14806 15282 13519 13682 
##   097   098   099   100   101   102   103   104   105   106   107   108 
## 15386 17262 13800 16006 18661 20238 13522 14096 10464 10629 16443 14089 
##   109   110   111   112   113   114   115   116   117   118   119   120 
## 10157 15627  8287  7647  6847  8201 22353 13792 12973 12208 11024 13482 
##   121   122   123   124   125   126   127   128   129   130   131   132 
## 12707 12936 14476 15825 14335 11533 15121 18384 18005 15034 13522 12939 
##   133   134   135   136   137   138   139   140   141   142   143   144 
## 10933 21979 13147 13539 14349 18145 15222 14097 10471 15412 21981 13851 
##   145   146   147   148   149   150   151   152   153   154   155   156 
##  8184  8413 19455 13359  9886 19750 17841 16159  8110 12113 12460 24578 
##   157   158   159   160   161   162   163   164   165   166   167   168 
## 18455 22931 16695 20916 14315 14228 13955 13703 10785 16482 23031 17416 
##   169   170   171   172   173   174 
##  9464 14097 11644 11737 13957 10290 
## 
## $`ssc-let-7d-3p`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
##  381  192  198  269  778  239  774  355  243  166  440  182  389  503  514 
##  016  017  018  019  020  021  022  023  024  025  026  027  028  029  030 
##  449  710  460 1071  411  356  330  319  750  449  525  339  378  393  256 
##  031  032  033  034  035  036  037  038  039  040  041  042  043  044  045 
##  438  407  352  498  569  527  269  326  199  265 1033  365  261  268  652 
##  046  047  048  049  050  051  052  053  054  055  056  057  058  059  060 
##  195  209  203  183  236   84  170  167  212  149   84  101  186   69  147 
##  061  062  063  064  065  066  067  068  069  070  071  072  073  074  075 
##  116  566  308  379  383  536  323  979  264  506  254  161  310  180  487 
##  076  077  078  079  080  081  082  083  084  085  086  087  088  089  090 
##  507  460  492  252  345  369  334  502  292  214  271  306  317  249  417 
##  091  092  093  094  095  096  097  098  099  100  101  102  103  104  105 
##  405  265  295  170  293  245  327  212  246  302  278  293  183  225  130 
##  106  107  108  109  110  111  112  113  114  115  116  117  118  119  120 
##  192  245  211  162  292  150   96   89  113  381  281  224  227  174  152 
##  121  122  123  124  125  126  127  128  129  130  131  132  133  134  135 
##  169  245  276  233  360  255  172  249  224  241  258  211  206  239  184 
##  136  137  138  139  140  141  142  143  144  145  146  147  148  149  150 
##  236  202  233  178  212  180  195  255  219  242   95  254  221  220  439 
##  151  152  153  154  155  156  157  158  159  160  161  162  163  164  165 
##  263  205   99  204  157  333  196  375  265  250  239  203  184  168  178 
##  166  167  168  169  170  171  172  173  174 
##  268  183  196  195  245  204  141  202  153 
## 
## $`ssc-let-7d-5p`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
## 4925 1938 2511 3076 3472 3276 3705 3301 3094 2068 3259 1873 5003 5616 4916 
##  016  017  018  019  020  021  022  023  024  025  026  027  028  029  030 
## 4470 3083 5035 6661 3659 3756 3558 4141 5225 3845  175 2610 4251 4195 2588 
##  031  032  033  034  035  036  037  038  039  040  041  042  043  044  045 
## 4543 4655 3962 5754 5456 5398 3037 3840 2029 3507 1356 4167 3663 3019 2545 
##  046  047  048  049  050  051  052  053  054  055  056  057  058  059  060 
## 2906 2304 2973 2427 3006  826 3037 2314 2058 1582 1233 1421 1378 1004 2518 
##  061  062  063  064  065  066  067  068  069  070  071  072  073  074  075 
## 1635 2799 2618 3885 4754 4920 4554 1831 2821 4780 2946 2586 3729 2315 5000 
##  076  077  078  079  080  081  082  083  084  085  086  087  088  089  090 
## 5502 5521 3304 2918 3832 4164 3616 4835 3364 2023 3814 3651 2961 2957 3820 
##  091  092  093  094  095  096  097  098  099  100  101  102  103  104  105 
## 4027 3215 2537 2731 2630 2136 2927 2860 2266 2406 2753 3117 2171 2170 1505 
##  106  107  108  109  110  111  112  113  114  115  116  117  118  119  120 
## 1951 2850 1942 1478 2574 1361 1299 1134 1208 3604 2010 2039 1824 1689 2046 
##  121  122  123  124  125  126  127  128  129  130  131  132  133  134  135 
## 2228 1847 2206 2194 2288 1822 2150 2506 2392 2583 2178 1856 1639 2861 1609 
##  136  137  138  139  140  141  142  143  144  145  146  147  148  149  150 
## 2020 2380 2388 2079 1946 1376 1955 2807 2051 1063 1138 2979 1606 1515 3108 
##  151  152  153  154  155  156  157  158  159  160  161  162  163  164  165 
## 2557 2661 1426 1710 1558 3197 2455 2962 2619 3118 2474 2587 1960 2102 1726 
##  166  167  168  169  170  171  172  173  174 
## 2873 3267 2413 1278 2121 2183 1609 1567 1741 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
## 2811 1302 1463 1690 2512 1410 2229 1577 1523  993 1898  954 3130 3081 2900 
##  016  017  018  019  020  021  022  023  024  025  026  027  028  029  030 
## 2163 2069 2654 3942 2326 2373 2205 2829 2481 2471   75 1345 2686 2324 1716 
##  031  032  033  034  035  036  037  038  039  040  041  042  043  044  045 
## 2746 2749 2181 3820 3060 3045 1752 2440 1034 2036  557 2730 2396 1674 1383 
##  046  047  048  049  050  051  052  053  054  055  056  057  058  059  060 
## 1900 1114 1508 1524 1782  349 1968 1254 1413  977  665  728  935  494 1225 
##  061  062  063  064  065  066  067  068  069  070  071  072  073  074  075 
##  768 1786 1261 2065 3047 2534 2852  996 1365 2816 1496 1127 2030 1532 2639 
##  076  077  078  079  080  081  082  083  084  085  086  087  088  089  090 
## 2478 3076 1880 1578 2163 1887 2096 2469 2071 1160 2017 2711 1471 1670 2262 
##  091  092  093  094  095  096  097  098  099  100  101  102  103  104  105 
## 2314 1888 1126 1320 1291 1209 1435 1684 1272 1079 1554 2120 1090 1148  711 
##  106  107  108  109  110  111  112  113  114  115  116  117  118  119  120 
## 1063 1532 1039  760 1402  787  552  474  596 2063 1165  956  977  900 1101 
##  121  122  123  124  125  126  127  128  129  130  131  132  133  134  135 
## 1049  915 1254 1097 1112  775 1069 1734 1319 1430 1111  845  803 1799  910 
##  136  137  138  139  140  141  142  143  144  145  146  147  148  149  150 
## 1000 1193 1540 1079 1021  743 1234 1800 1057  487  555 1786 1020  695 1886 
##  151  152  153  154  155  156  157  158  159  160  161  162  163  164  165 
## 1443 1418  643  850 1019 1993 1403 1891 1319 1532 1267 1212  984 1077  852 
##  166  167  168  169  170  171  172  173  174 
## 1678 1683 1422  641 1090  925  910 1254  839 
## 
## $`ssc-let-7f`
##     001     002     003     004     005     006     007     008     009 
## 35432.5 16422.0 20161.5 22985.5 16512.5 24475.5 19175.0 22094.5 22425.5 
##     010     011     012     013     014     015     016     017     018 
## 16483.0 29076.5 15185.0 32137.0 36462.5 31451.5 29666.0 13446.5 34278.0 
##     019     020     021     022     023     024     025     026     027 
## 34402.5 25967.0 28076.5 27274.5 29413.5 29411.0 26171.5   576.0 20826.5 
##     028     029     030     031     032     033     034     035     036 
## 29465.5 29266.0 20800.0 31999.0 32016.0 28727.0 36303.5 34956.5 34056.0 
##     037     038     039     040     041     042     043     044     045 
## 27197.0 25904.0 20112.0 24842.0 10064.0 29654.5 25786.5 20981.5 18572.5 
##     046     047     048     049     050     051     052     053     054 
## 20946.5 18135.5 22221.5 19081.0 23174.5  7665.0 24685.5 16910.5 20088.5 
##     055     056     057     058     059     060     061     062     063 
## 11720.5  9891.0 10624.0 16741.5  8510.0 20538.5 12141.5 12570.0 21874.0 
##     064     065     066     067     068     069     070     071     072 
## 27558.5 35355.5 32477.0 31694.0 10026.0 21929.0 35707.0 22301.5 19183.5 
##     073     074     075     076     077     078     079     080     081 
## 26026.5 20536.5 31470.0 34575.5 33807.0 19231.0 23504.5 27732.0 29918.0 
##     082     083     084     085     086     087     088     089     090 
## 27981.5 32045.5 24423.5 17376.0 30385.0 27878.5 22126.0 20542.0 29883.0 
##     091     092     093     094     095     096     097     098     099 
## 28315.5 25277.5 19161.0 21219.0 21407.5 17459.5 22679.0 23619.5 18688.0 
##     100     101     102     103     104     105     106     107     108 
## 18800.0 19519.0 21239.0 16502.5 18830.0 13480.5 15197.0 22217.5 16950.5 
##     109     110     111     112     113     114     115     116     117 
## 14487.5 21120.5 14236.0 12651.0 11789.0 11817.0 27063.0 16380.0 18726.0 
##     118     119     120     121     122     123     124     125     126 
## 18098.0 14844.5 17390.5 19551.0 16101.0 16648.0 18694.5 19417.0 16761.0 
##     127     128     129     130     131     132     133     134     135 
## 17132.5 22043.0 20539.5 21527.0 17039.5 16169.5 15178.5 19578.5 14941.5 
##     136     137     138     139     140     141     142     143     144 
## 17284.5 18355.0 19715.0 17365.0 17702.0 12597.5 16530.5 21863.0 17525.5 
##     145     146     147     148     149     150     151     152     153 
## 10830.0 12593.5 21235.5 15161.0 14276.5 24075.0 21701.0 20042.5 12906.5 
##     154     155     156     157     158     159     160     161     162 
## 15595.0 15125.5 23704.0 20668.5 23188.0 21315.0 23538.0 19776.5 20946.0 
##     163     164     165     166     167     168     169     170     171 
## 17024.5 18061.0 15794.0 21232.5 25268.5 22609.5 12634.5 18261.5 17327.0 
##     172     173     174 
## 14902.0 14506.5 15007.5
```


###3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function
Example: ldply(.data, .fun, .id)

id = name of the index column (used if data is a named list). Pass NULL to avoid creation
     of the index column. For compatibility, omit this argument or pass NA to avoid converting the index column
     to a factor; in this case, ".id" is used as column name.


```r
dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(dfmeanrc[1:8])
```

```
##             .id     001   002     003     004     005     006   007
## 1    ssc-let-7a 48132.5 23427 28447.5 29860.0 40757.5 29749.5 38799
## 2    ssc-let-7c 32745.0 14987 18144.0 18681.0 34313.0 15294.0 28022
## 3 ssc-let-7d-3p   381.0   192   198.0   269.0   778.0   239.0   774
## 4 ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0  3276.0  3705
## 5    ssc-let-7e  2811.0  1302  1463.0  1690.0  2512.0  1410.0  2229
## 6    ssc-let-7f 35432.5 16422 20161.5 22985.5 16512.5 24475.5 19175
```

```r
dim(dfmeanrc)
```

```
## [1] 411 175
```

These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase
And there are 174 animals in the analysis, plus the miRNA column.

Check that the correct miRNA name went with the correct data:


```r
sum(names(meanrc)==dfmeanrc[,1])
```

```
## [1] 411
```

```r
identical(names(meanrc), dfmeanrc[,1])
```

```
## [1] TRUE
```

```r
colnames(dfmeanrc)[[1]]<-"miRNA"

sum(colnames(dfmeanrc)==colnames(mirquant))
```

```
## [1] 175
```

```r
sum(colnames(dfmeanrc)!=colnames(mirquant))
```

```
## [1] 0
```

```r
head(dfmeanrc[,1:10])
```

```
##           miRNA     001   002     003     004     005     006   007
## 1    ssc-let-7a 48132.5 23427 28447.5 29860.0 40757.5 29749.5 38799
## 2    ssc-let-7c 32745.0 14987 18144.0 18681.0 34313.0 15294.0 28022
## 3 ssc-let-7d-3p   381.0   192   198.0   269.0   778.0   239.0   774
## 4 ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0  3276.0  3705
## 5    ssc-let-7e  2811.0  1302  1463.0  1690.0  2512.0  1410.0  2229
## 6    ssc-let-7f 35432.5 16422 20161.5 22985.5 16512.5 24475.5 19175
##       008     009
## 1 31459.5 29857.0
## 2 19512.0 16741.0
## 3   355.0   243.0
## 4  3301.0  3094.0
## 5  1577.0  1523.0
## 6 22094.5 22425.5
```

Check that each position of the let-7a element of the list matches the let-7a row of the dataframe:


```r
sum((meanrc$'ssc-let-7a')-(dfmeanrc[1,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

Check a second miRNA in the same way:


```r
sum((meanrc$'ssc-miR-369')-(dfmeanrc[206,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

And a third miRNA in the same way:


```r
sum((meanrc$'ssc-miR-9')-(dfmeanrc[321,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

###4. Filter the data for expression threshold: The mean read count for the miRNA needs to be greater than the number of animals in the population


```r
rowSums(dfmeanrc[1:3,2:ncol(dfmeanrc)])
```

```
##       1       2       3 
## 5061044 3238382   51589
```

```r
sum(dfmeanrc[1,2:ncol(dfmeanrc)])
```

```
## [1] 5061044
```

```r
sum(dfmeanrc[2,2:ncol(dfmeanrc)])
```

```
## [1] 3238382
```

```r
sum(dfmeanrc[3,2:ncol(dfmeanrc)])
```

```
## [1] 51589
```

```r
sum(rowSums(dfmeanrc[,2:ncol(dfmeanrc)]) > 174)
```

```
## [1] 285
```

So, 285 miRNAs have expression greater than 174 (number of animals in population).

Create a subset of the large dataframe containing only the miRNAs expressed > 174 times:


```r
filtermeanrc<-dfmeanrc[which(rowSums(dfmeanrc[,2:ncol(dfmeanrc)]) > 174),]
dim(filtermeanrc)
```

```
## [1] 285 175
```

```r
filtermeanrc[1:10,1]
```

```
##  [1] "ssc-let-7a"    "ssc-let-7c"    "ssc-let-7d-3p" "ssc-let-7d-5p"
##  [5] "ssc-let-7e"    "ssc-let-7f"    "ssc-let-7g"    "ssc-let-7i"   
##  [9] "ssc-miR-1"     "ssc-miR-100"
```

Check that the column order did not switch in the merge:


```r
sum(colnames(filtermeanrc) == colnames(dfmeanrc))
```

```
## [1] 175
```

All the column names match!

Identify which rows contain 0s and which do not:


```r
rowsums.zero<-rowSums(filtermeanrc[,2:ncol(filtermeanrc)]==0)
rowsums.zero
```

```
##   1   2   3   4   5   6   7   8   9  10  11  12  15  16  17  18  19  20 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0 
##  22  24  25  26  27  28  29  30  31  32  34  35  36  37  38  39  40  41 
##   0   0   0   0   0   0   2   7   0   0   0  32   3   3   0   0   0   0 
##  42  43  44  45  46  47  48  49  51  52  53  54  55  56  57  58  59  60 
##   0   0   0   0   0   0  27   1  23   4   0   0   0   0   0   0   0   0 
##  61  62  63  65  66  67  68  69  70  71  72  73  74  75  76  78  79  80 
##   0   0   0   0   0   0   0   0   0   0   0   0   0   0  82   0   0   0 
##  81  82  83  84  85  86  88  89  90  92  93  94  95  97  98  99 100 101 
##   0   0   0   0   0   0   0   0  27  18   3   0   0   0   3   0   1   0 
## 102 103 104 105 108 109 110 111 112 113 114 115 116 117 119 120 121 122 
##   0   0   0   0   0   0   5   0   0   0   0   0   1   0   5   0   5   0 
## 123 124 125 126 127 128 129 132 133 135 136 139 140 141 142 144 145 146 
##   0   0   0   0   0  75   0  52   0   0   0   0  23   0   0   0   0   0 
## 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 
##  57   0   0  60   0   0   0   3   0   0   0   0   0   0   0   0   0   0 
## 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 
##   0   0   0   0   0   2   0   0   1   0   0   0   0  51   0   0   0   0 
## 184 185 186 187 188 189 190 191 192 193 194 195 196 197 199 200 201 202 
##   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0 
## 203 204 206 207 209 210 211 212 213 214 216 217 218 219 220 222 223 225 
##   0   0   0   0   0   0   1   0   0   1   2   0   0  10   0   0   0   0 
## 226 227 228 229 230 231 233 234 235 236 242 243 244 245 246 247 248 249 
##   0   0   0   0   0  75   0  30   0   0   0   0   0   0   0   0  13   0 
## 250 252 253 254 255 256 258 259 260 262 263 264 265 266 267 268 269 270 
##   0   0  46   0  20  20   0   0   0   0   0   0   0   0   0   0   0   0 
## 271 272 273 274 276 277 278 279 280 281 282 283 284 285 287 288 289 290 
##   0   0   0   0   0   0  10   0  23   0   0   8   0   0   0   0   1   0 
## 291 292 293 294 297 298 300 302 307 313 314 315 316 318 319 320 321 322 
##   0   0   7  68   6   0  53  71   0   0   0   0   0   0   0   0   1   1 
## 323 324 325 327 328 329 331 346 357 388 390 398 407 410 411 
##   1   0   0   0   0  57  44   0  37   1   0   6   4   0   0
```

Notice that some miRNAs have multiple 0 read counts across animals, while some only have one or two 0 read counts.

How many miRNAs have 0s?


```r
sum(rowsums.zero!=0)
```

```
## [1] 58
```

So, 58 miRNA profiles expressed greater than 174 times still contain zeros and need to be adjusted prior log-transformation via voom function.

Create a logical vector containing true if a row contains a 0:


```r
rowcontains.zeros.logical<-rowsums.zero!=0

head(rowcontains.zeros.logical)
```

```
##     1     2     3     4     5     6 
## FALSE FALSE FALSE FALSE FALSE FALSE
```




```r
sum(rowcontains.zeros.logical==TRUE)
```

```
## [1] 58
```

58 miRNAs contain 0s.


```r
sum(rowcontains.zeros.logical!=TRUE)
```

```
## [1] 227
```

227 miRNAs contain no 0s.
###5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR


```r
head(configfile)
```

```
##   pigid code
## 1  1034  001
## 2  1036  002
## 3  1041  003
## 4  1049  004
## 5  1058  005
## 6  1060  006
```

Check that the 3 miRNAs maintained the same positions in filtermeanrc versus dfmeanrc:

ssc-let-7a:


```r
sum((filtermeanrc[1,2:ncol(filtermeanrc)])-(dfmeanrc[1,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

ssc-miR-4332:


```r
sum((filtermeanrc[181,2:ncol(filtermeanrc)])-(dfmeanrc[203,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

ssc-miR-99a:


```r
sum((filtermeanrc[284,2:ncol(filtermeanrc)])-(dfmeanrc[410,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

Now I need to substitute the pid ids for the 3-digit code, ensuring the names stay in the correct order:

filtermeanrc: matrix of average read counts filtered for expression > 174

Set first column of filtermeanrc (miRNA ids) as the row.names:


```r
rownames(filtermeanrc)<-filtermeanrc$miRNA
```

Eliminate column of row names:


```r
filtermeanrc$miRNA<-NULL
```

Use match function to find positional vector and match column names:


```r
configfile[match(configfile$code,colnames(filtermeanrc)),1]
```

```
##   [1] "1034" "1036" "1041" "1049" "1058" "1060" "1080" "1082" "1085" "1091"
##  [11] "1096" "1100" "1107" "1110" "1111" "1113" "1116" "1123" "1134" "1136"
##  [21] "1145" "1147" "1152" "1154" "1158" "1170" "1177" "1179" "1192" "1194"
##  [31] "1197" "1199" "1205" "1207" "1237" "1239" "1240" "1242" "1265" "1267"
##  [41] "1278" "1282" "1291" "1295" "1300" "1304" "1321" "1323" "1423" "1424"
##  [51] "1425" "1426" "1431" "1434" "1435" "1444" "1445" "1449" "1456" "1458"
##  [61] "1482" "1484" "1491" "1493" "1502" "1504" "1510" "1512" "1517" "1523"
##  [71] "1529" "1532" "1533" "1534" "1537" "1543" "1578" "1580" "1589" "1592"
##  [81] "1593" "1594" "1625" "1627" "1638" "1640" "1644" "1646" "1652" "1662"
##  [91] "1669" "1677" "1685" "1687" "1695" "1697" "1746" "1758" "1760" "1776"
## [101] "1778" "1782" "1784" "1785" "1789" "1793" "1798" "1800" "1818" "1819"
## [111] "1820" "1833" "1836" "1839" "1843" "1844" "1879" "1881" "1884" "1886"
## [121] "1889" "1891" "1903" "1904" "1907" "1910" "1914" "1916" "1928" "1930"
## [131] "1965" "1971" "1976" "1980" "1989" "1991" "1999" "2003" "2018" "2020"
## [141] "2022" "2024" "2026" "2027" "2029" "2030" "2064" "2071" "2073" "2076"
## [151] "2094" "2100" "2118" "2119" "2120" "2123" "2131" "2135" "2141" "2143"
## [161] "2152" "2154" "2164" "2168" "2195" "2197" "2229" "2231" "2261" "2263"
## [171] "2297" "2303" "2311" "2317"
```

Assign the column names using match:


```r
colnames(filtermeanrc)<- configfile[match(configfile$code,colnames(filtermeanrc)),1]
head(filtermeanrc[1:5])
```

```
##                  1034  1036    1041    1049    1058
## ssc-let-7a    48132.5 23427 28447.5 29860.0 40757.5
## ssc-let-7c    32745.0 14987 18144.0 18681.0 34313.0
## ssc-let-7d-3p   381.0   192   198.0   269.0   778.0
## ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0
## ssc-let-7e     2811.0  1302  1463.0  1690.0  2512.0
## ssc-let-7f    35432.5 16422 20161.5 22985.5 16512.5
```

Do the same positional check for the three miRNAs:

ssc-let-7a


```r
sum((filtermeanrc[1,])-(dfmeanrc[1,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

ssc-miR-363


```r
sum((filtermeanrc[181,])-(dfmeanrc[203,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

ssc-miR-99a


```r
sum((filtermeanrc[284,])-(dfmeanrc[410,2:ncol(dfmeanrc)]))
```

```
## [1] 0
```

###6. Transpose the data.frame to make the animals take the rows and the miRNAs the columns


```r
transposefiltermeanrc<-t(filtermeanrc)

dim(transposefiltermeanrc)
```

```
## [1] 174 285
```

```r
head(transposefiltermeanrc[,1:5])
```

```
##      ssc-let-7a ssc-let-7c ssc-let-7d-3p ssc-let-7d-5p ssc-let-7e
## 1034    48132.5      32745           381          4925       2811
## 1036    23427.0      14987           192          1938       1302
## 1041    28447.5      18144           198          2511       1463
## 1049    29860.0      18681           269          3076       1690
## 1058    40757.5      34313           778          3472       2512
## 1060    29749.5      15294           239          3276       1410
```

```r
is.numeric(transposefiltermeanrc)
```

```
## [1] TRUE
```

## Save data
What I am saving here is the unfiltered average read counts in one data.frame, and the filtered, transposed average read counts as another data.frame.

###1. Save the unfiltered average read counts as a data.frame and an Rdata object


```r
write.table(dfmeanrc, file = "../1_mean_mature_mirna_expression_unfiltered.txt", quote = FALSE, sep = "\t", col.names = TRUE)
save(dfmeanrc, file = "../1_mean_mature_mirna_expression_unfiltered.Rdata")
```

###2. Save the filtered, transposed average read counts as a data.frame and an Rdata object


```r
write.table(transposefiltermeanrc, file = "../2_mean_mature_mirna_expression_filtered.txt", quote = FALSE, sep = "\t ", col.names = TRUE)
save(transposefiltermeanrc, file = "../2_mean_mature_mirna_expression_filtered.Rdata")
```

