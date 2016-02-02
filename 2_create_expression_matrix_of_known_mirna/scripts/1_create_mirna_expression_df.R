#' **Script:** `1_create_mirna_expression_df.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/scripts`
#' 
#' **Date:**  1/26/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output`
#' 
#' **Input File(s):** `miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna`
#' 
#' **Output File(s):** `1_mean_mature_mirna_expression_unfiltered.txt`
#' 
#'                     `1_mean_mature_mirna_expression_unfiltered.Rdata`
#' 
#'                     `2_mean_mature_mirna_expression_filtered.txt`
#' 
#'                     `2_mean_mature_mirna_expression_filtered.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)
#' 
#' ## Objectives
#' The objective of this script is to create a matrix of miRNA expression data (read counts) for incorporation into the miRNA eQTL analysis, beginning with the known miRNA expression data output from the miRDeep2 core module.
#' Later, expression data of the putative novel miRNA candidates can also be included, also output from the miRDeep2 core module. 
#' To acheive one read count per miRNA per animal, I will need to take the average of the mature read counts in instances of multiple precursor sequences.
#' 
#' The result of this script will be two data frames containing one average read count per miRNA per animal: one unfiltered, and one filtered for miRNAs expressed greater than the number of animals in the population (174) and transposed.
#' 
#' So, what I need to do:
#' 
#' 1. Extract the columns of the mature miRNA read counts for each animal
#' 2. Use the 'by' function to apply the function colmeans to the entire data frame of read counts (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that group of miRNA. The result of this will be a list containing the average read counts for each miRNA for each animal)
#' 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function
#' 4. Filter the data for expression threshold: The read count for the miRNA needs to be greater than the number of animals in the population
#' 5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR
#' 6. Transpose the data.frame to make the animals the rows and the miRNAs the columns
#' 

#' ## Install libraries
library(plyr)

#' ## Load data

#' ###1. Read in the csv of the miRNA expression profiles, using check.names default to correct column names with non-ASCII characters
rc<-read.csv("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv", sep = "\t", header = TRUE, row.names=NULL)
#' Set the name of the first column to "miRNA":
colnames(rc)[[1]]<-"miRNA"
#' Remove the "X" character from the beginning of each column:
colnames(rc)<-gsub("X", "", colnames(rc))
colnames(rc)
#' View the data set:
head(rc[1:8])

#' ###2. Read in the config file, maintaining the characters in the 3-digit code names
configfile<-read.table("../../../1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/config_for_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)

#' Remove the "_express.fa" from each file name to leave the pig id:
configfile$V1<-gsub("_express.fa", "", configfile$V1)

#' Make filenames more informative:
colnames(configfile)<-c("pigid","code")
colnames(configfile)

#' ## Analysis

#' ###1. Extract the columns of mature read counts for each miRNA for each animal
mirquant<-rc[,c(1,5:178)]
colnames(mirquant)

dim(mirquant)

head(mirquant[1:8])

#' Take a subset of this data.frame for testing:
test<-mirquant[1:20,1:8]
head(test)

#' ###2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts
#' (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
#' The result of this will be a list containing the average read counts for each miRNA for each animal
#' 
#' Example: by(data, index, function)
head(by(test[,2:ncol(test)], test[,1], colMeans))

#' Apply the by function to the full dataframe:
meanrc<-by(mirquant[,2:ncol(mirquant)], mirquant[,1], colMeans)

#' This should be 411, the number of mature pig miRNAs in miRBase:
length(meanrc)

head(meanrc)

#' 
#' ###3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function

#' Example: ldply(.data, .fun, .id)
#' 
#' id = name of the index column (used if data is a named list). Pass NULL to avoid creation
#'      of the index column. For compatibility, omit this argument or pass NA to avoid converting the index column
#'      to a factor; in this case, ".id" is used as column name.

dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(dfmeanrc[1:8])

dim(dfmeanrc)
#' These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase
#' And there are 174 animals in the analysis, plus the miRNA column.

#' 
#' Check that the correct miRNA name went with the correct data:
sum(names(meanrc)==dfmeanrc[,1])

identical(names(meanrc), dfmeanrc[,1])

colnames(dfmeanrc)[[1]]<-"miRNA"

sum(colnames(dfmeanrc)==colnames(mirquant))

sum(colnames(dfmeanrc)!=colnames(mirquant))

head(dfmeanrc[,1:10])

#' Check that each position of the let-7a element of the list matches the let-7a row of the dataframe:
sum((meanrc$'ssc-let-7a')-(dfmeanrc[1,2:ncol(dfmeanrc)]))

#' Check a second miRNA in the same way:
sum((meanrc$'ssc-miR-369')-(dfmeanrc[206,2:ncol(dfmeanrc)]))

#' And a third miRNA in the same way:
sum((meanrc$'ssc-miR-9')-(dfmeanrc[321,2:ncol(dfmeanrc)]))


#' ###4. Filter the data for expression threshold: The mean read count for the miRNA needs to be greater than the number of animals in the population

rowSums(dfmeanrc[1:3,2:ncol(dfmeanrc)])

sum(dfmeanrc[1,2:ncol(dfmeanrc)])

sum(dfmeanrc[2,2:ncol(dfmeanrc)])

sum(dfmeanrc[3,2:ncol(dfmeanrc)])

sum(rowSums(dfmeanrc[,2:ncol(dfmeanrc)]) > 174)
#' So, 285 miRNAs have expression greater than 174 (number of animals in population).

#' 
#' Create a subset of the large dataframe containing only the miRNAs expressed > 174 times:
filtermeanrc<-dfmeanrc[which(rowSums(dfmeanrc[,2:ncol(dfmeanrc)]) > 174),]
dim(filtermeanrc)

filtermeanrc[1:10,1]

#' Check that the column order did not switch in the merge:
sum(colnames(filtermeanrc) == colnames(dfmeanrc))
#' All the column names match!

#' 
#' Identify which rows contain 0s and which do not:
rowsums.zero<-rowSums(filtermeanrc[,2:ncol(filtermeanrc)]==0)
rowsums.zero
#' Notice that some miRNAs have multiple 0 read counts across animals, while some only have one or two 0 read counts.
#' 

#' How many miRNAs have 0s?
sum(rowsums.zero!=0)
#' So, 58 miRNA profiles expressed greater than 174 times still contain zeros and need to be adjusted prior log-transformation via voom function.

#' 
#' Create a logical vector containing true if a row contains a 0:
rowcontains.zeros.logical<-rowsums.zero!=0

head(rowcontains.zeros.logical)

#' 
sum(rowcontains.zeros.logical==TRUE)
#' 58 miRNAs contain 0s.

sum(rowcontains.zeros.logical!=TRUE)
#' 227 miRNAs contain no 0s.


#' ###5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR

head(configfile)

#' Check that the 3 miRNAs maintained the same positions in filtermeanrc versus dfmeanrc:

#' 
#' ssc-let-7a:
sum((filtermeanrc[1,2:ncol(filtermeanrc)])-(dfmeanrc[1,2:ncol(dfmeanrc)]))

#' ssc-miR-4332:
sum((filtermeanrc[181,2:ncol(filtermeanrc)])-(dfmeanrc[203,2:ncol(dfmeanrc)]))

#' ssc-miR-99a:
sum((filtermeanrc[284,2:ncol(filtermeanrc)])-(dfmeanrc[410,2:ncol(dfmeanrc)]))


#' Now I need to substitute the pid ids for the 3-digit code, ensuring the names stay in the correct order:
#' 
#' filtermeanrc: matrix of average read counts filtered for expression > 174
#' 
#' Set first column of filtermeanrc (miRNA ids) as the row.names:
rownames(filtermeanrc)<-filtermeanrc$miRNA

#' Eliminate column of row names:
filtermeanrc$miRNA<-NULL

#' Use match function to find positional vector and match column names:
configfile[match(configfile$code,colnames(filtermeanrc)),1]

#' Assign the column names using match:
colnames(filtermeanrc)<- configfile[match(configfile$code,colnames(filtermeanrc)),1]
head(filtermeanrc[1:5])


#' Do the same positional check for the three miRNAs:
#' 
#' ssc-let-7a
sum((filtermeanrc[1,])-(dfmeanrc[1,2:ncol(dfmeanrc)]))

#' ssc-miR-363
sum((filtermeanrc[181,])-(dfmeanrc[203,2:ncol(dfmeanrc)]))

#' ssc-miR-99a
sum((filtermeanrc[284,])-(dfmeanrc[410,2:ncol(dfmeanrc)]))


#' ###6. Transpose the data.frame to make the animals take the rows and the miRNAs the columns
transposefiltermeanrc<-t(filtermeanrc)

dim(transposefiltermeanrc)

head(transposefiltermeanrc[,1:5])

is.numeric(transposefiltermeanrc)

#' ## Save data
#' What I am saving here is the unfiltered average read counts in one data.frame, and the filtered, transposed average read counts as another data.frame.
#'

#' ###1. Save the unfiltered average read counts as a data.frame and an Rdata object
write.table(dfmeanrc, file = "../1_mean_mature_mirna_expression_unfiltered.txt", quote = FALSE, sep = "\t", col.names = TRUE)
save(dfmeanrc, file = "../1_mean_mature_mirna_expression_unfiltered.Rdata")

#' ###2. Save the filtered, transposed average read counts as a data.frame and an Rdata object
write.table(transposefiltermeanrc, file = "../2_mean_mature_mirna_expression_filtered.txt", quote = FALSE, sep = "\t ", col.names = TRUE)
save(transposefiltermeanrc, file = "../2_mean_mature_mirna_expression_filtered.Rdata")

