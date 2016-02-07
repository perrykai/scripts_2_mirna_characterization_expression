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

#' Do a for loop to check the class of each column in the data.frame

colclass<-NULL

for (i in colnames(rc)) {
   colclass<-c(colclass,class(rc[,i]))
}

table(colclass)
colclass


#' Notice here that the mature miRNA names and the miRNA precursor names are considered factors (columns 1 and 3).

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

test[1:10,]

#' ###2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts
#' (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
#' The result of this will be a list containing the average read counts for each miRNA for each animal.
#' 
#' Example: by(data, index, function)
bytst<-by(test[,2:ncol(test)], test[,1], colMeans)
bytst[1:25]
#' Notice here that the rest of the miRNA names remain since the miRNA name is a factor, but since there is no data for them they are filled with NULL.
#' 

#' Apply the by function to the full dataframe:
meanrc<-by(mirquant[,2:ncol(mirquant)], mirquant[,1], colMeans)

#' This should be 411 (the number of mature pig miRNAs in miRBase), meaning we have one expression profile for each mature miRNA:
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
#' and there are 174 animals in the analysis, plus the miRNA column.

#' 
#' Check that the correct miRNA name went with the correct data:
if (sum(names(meanrc)!=dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(dfmeanrc)!=colnames(mirquant)) != 0) stop ("animal order not the same")

head(dfmeanrc[,1:10])

#' ###4. Filter the data for expression threshold: The mean read count for the miRNA needs to be greater than the number of animals in the population
#' First, a quick check that the rows didn't get switched around somehow:
if (rowSums(dfmeanrc[1,2:ncol(dfmeanrc)]) != sum(dfmeanrc[1,2:ncol(dfmeanrc)])) stop ("rowsums not the same")

if (rowSums(dfmeanrc[2,2:ncol(dfmeanrc)]) != sum(dfmeanrc[2,2:ncol(dfmeanrc)])) stop ("rowsums not the same")

if (rowSums(dfmeanrc[3,2:ncol(dfmeanrc)]) != sum(dfmeanrc[3,2:ncol(dfmeanrc)])) stop ("rowsums not the same")


#' 
#' Identify which rows contain 0s and which do not:
rowsums.zero<-rowSums(filtermeanrc[,2:ncol(filtermeanrc)]==0)
rowsums.zero
#' Notice that some miRNAs have multiple 0 read counts across animals, while some only have one or two 0 read counts.
#' 

#' How many miRNAs have 0s?
#' 
#' Create a logical vector containing true if a row contains a 0:
rowcontains.zeros.logical<-rowsums.zero!=0

head(rowcontains.zeros.logical)

table(rowcontains.zeros.logical)
#' So, 58 miRNA profiles contain 0 read counts
#' 
#' ###5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR

head(configfile)

#' Now I need to substitute the pid ids for the 3-digit code, ensuring the names stay in the correct order:
#' 
#' filtermeanrc: matrix of average read counts filtered for expression > 174
#' 
#' Set first column of filtermeanrc (miRNA ids) as the row.names:
rownames(filtermeanrc)<-filtermeanrc$miRNA

#' Eliminate column of row names:
filtermeanrc<-filtermeanrc[,-c(1)]
head(filtermeanrc[,1:10])
dim(filtermeanrc)


#' Use match function to find positional index and match column names:
#' 
#' The object filtermeanrc has column names that need to be re-named. I have the config file which contains
#' the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
#' based on where the config file "code" column matches the position of the filtermeanrc object's column names, then having it return the corresponding value in column "pigid". 
#' 
#' So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 
#' 
#' "Where does [vector] match in [matrix]?" or "Match the column names of filtermeanrc to the configfile "code" column."
configfile[match(colnames(filtermeanrc),configfile$code),"pigid"]

#' Assign the column names using match:
colnames(filtermeanrc)<- configfile[match(colnames(filtermeanrc),configfile$code),"pigid"]
head(filtermeanrc[1:5])


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

