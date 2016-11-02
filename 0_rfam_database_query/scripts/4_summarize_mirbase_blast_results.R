#' **Script:** `4_summarize_mirbase_blast_results.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/`
#' 
#' **Date:**  `7/11/16`
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Input File(s):** `0_susscrofa_mirbase_blastn_output_e5.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Output File(s):** 
#' 
#' 1. `0a_susscrofa_filtered_mirbase_blastn_results.txt`
#' 2. `0b_susscrofa_filtered_uniqseq_miRBase_blastn_results.csv`
#' 3. `0c_susscrofa_filtered_totseq_miRBase_blastn_results.csv`
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
#' 
#' The objective of this script is to filter and summarize the output from the first of the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
#' This script filters and summarizes the blast results from the Sus scrofa miRBase database query.
#' 
#' First, if there is only one hit for a sequence, return that hit.
#' If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
#' Then, retain matches with > 96% sequence identity.
#' Finally, the first filter returns sequences with the maximum bitscore.
#' 
#' The second filter returns one BLAST hit if the sequence with greatest percent match to blast hit.
#' Finally, if there are still remaining sequences with multiple hits, retain only the first hit. 
#' 
#' ## Install libraries
library(parallel)
library(stringr)
rm(list=ls())

#' ## Load data
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts")

system.time({
ssc<-read.table("../0_susscrofa_mirbase_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})

dim(ssc)
head(ssc)
str(ssc)

#' ## Analysis
ssc$dbseq_id<-as.character(ssc$dbseq_id)
ssc$query_id<-as.character(ssc$query_id)

#' Create a subset data frame containing the pertinent information for filtering
mirbasesubset<-ssc[,c("query_id", "dbseq_id", "perc_identical", "evalue", "bitscore")]
dim(mirbasesubset)
head(mirbasesubset)

#' The number of unique sequences with susscrofa miRBase hits in this dataset
length(unique(mirbasesubset$query_id))

#' The number of unique miRBase miRNAs identified by this BLAST query
length(unique(mirbasesubset$dbseq_id))

######################################################################
#' Obtained the following code from Deborah Velez, how she filtered her BLAST results for gene annotation (/mnt/research/ernstc_lab/fat_eqtl/BLAST/code/annotation_pigoligoarray.R)
######################################################################
#' Format the miRBase blast output data.frame to a list to filter the sequences with multiple miRBase hits
# Number of unique sequences
n <- length(unique(mirbasesubset$query_id))
n

# Create the sequence list to filter
system.time({
idx1 <- mclapply(as.character(unique(mirbasesubset$query_id))[1:round(n/5)], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx2 <- mclapply(as.character(unique(mirbasesubset$query_id))[(round(n/5) + 1):(2*round(n/5))], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx3 <- mclapply(as.character(unique(mirbasesubset$query_id))[((2*round(n/5)) + 1):(3*round(n/5))], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx4 <- mclapply(as.character(unique(mirbasesubset$query_id))[((3*round(n/5)) + 1):(4*round(n/5))], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx5 <- mclapply(as.character(unique(mirbasesubset$query_id))[((4*round(n/5)) + 1):n], function(x)
        mirbasesubset[as.character(mirbasesubset$query_id) == x,], mc.cores=10)
idx <- c(idx1,idx2,idx3,idx4,idx5)
})

length(idx)

#' Function to filter multiple miRBase blast hits per sequence
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
        bit <- seqblast[seqblast$bitscore == max(seqblast$bitscore),]
                        if (nrow(bit) == 1){
                        return(bit)
        }

        return(bit)
}

#' Apply filter to the sequence list
system.time({
p1 <- mclapply(idx[1:round(n/5)], filter, mc.cores=10)
p2 <- mclapply(idx[(round(n/5) + 1):(2*round(n/5))], filter, mc.cores=10)
p3 <- mclapply(idx[((2*round(n/5)) + 1):(3*round(n/5))], filter, mc.cores=10)
p4 <- mclapply(idx[((3*round(n/5)) + 1):(4*round(n/5))], filter, mc.cores=10)
p5 <- mclapply(idx[((4*round(n/5)) + 1):n], filter, mc.cores=10)
stp1 <- c(p1,p2,p3,p4,p5)
})

length(stp1)
head(stp1)
tail(stp1)

#' How many sequences have more than one hit after the filtering
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
table(unlist(lapply(stp2,nrow)))

filter2 <- filter <- function(seqblast){
# Sequence has only one miRBase blast hit
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

#' Apply a second filter to the filtered sequence list
system.time({
q1 <- mclapply(stp1[1:round(n/5)], filter2, mc.cores=10)
q2 <- mclapply(stp1[(round(n/5) + 1):(2*round(n/5))], filter2, mc.cores=10)
q3 <- mclapply(stp1[((2*round(n/5)) + 1):(3*round(n/5))], filter2, mc.cores=10)
q4 <- mclapply(stp1[((3*round(n/5)) + 1):(4*round(n/5))], filter2, mc.cores=10)
q5 <- mclapply(stp1[((4*round(n/5)) + 1):n], filter2, mc.cores=10)
stp3 <- c(q1,q2,q3,q4,q5)
})

#' This command sums the number of sequence IDs appearing more than once in the dataset.
length(stp3[unlist(lapply(stp3,nrow)) > 1])

stp4 <- stp3[unlist(lapply(stp3,nrow)) > 1]
length(stp4)
table(unlist(lapply(stp4,nrow)))

#' Summary file of small RNA sequence miRBase blast results
sumblast <- do.call(rbind,stp3)
head(sumblast)
dim(sumblast)
length(unique(sumblast$query_id))

table(sumblast$dbseq_id)

uniqsumblast<-as.data.frame(table(sumblast$dbseq_id))
colnames(uniqsumblast)<-c("miRNA", "Freq")
uniqsumblast

#' Add the column of sequence count to the sumblast data frame
sumblast$seq_count<-as.numeric(str_split_fixed(sumblast$query_id, "_x", 2)[,2])
head(sumblast)

#' Use the "by" function to sum the sequence counts by their gene biotypes
totalsumbiotype<-as.matrix(by(sumblast$seq_count, sumblast$dbseq_id, sum))
totalsumbiotype

if (sum(rownames(totalsumbiotype) != uniqsumblast$dbseq_id)) stop ("miRNA names not equal")

#' As a check, manually sum the ssc-miR-1 and the ssc-miR-339 counts:
sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-1"])

sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-339"])

if (sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-1"]) != totalsumbiotype["ssc-miR-1",]) stop ("ssc-miR-1 counts not equal")
if (sum(sumblast$seq_count[sumblast$dbseq_id == "ssc-miR-339"]) != totalsumbiotype["ssc-miR-339",]) stop ("ssc-miR-339 counts not equal")

#' ## Save data
save(sumblast, uniqsumblast, totalsumbiotype, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/0a_susscrofa_filtered_miRBase_blastn_results.Rdata"))
write.csv(uniqsumblast, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/0b_susscrofa_filtered_uniqseq_miRBase_blastn_results.csv")
write.csv(totalsumbiotype, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/0c_susscrofa_filtered_totseq_miRBase_blastn_results.csv", row.names=TRUE)