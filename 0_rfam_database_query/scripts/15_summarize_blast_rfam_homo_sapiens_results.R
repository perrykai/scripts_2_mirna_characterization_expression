#' **Script:** `15_summarize_blast_rfam_homo_sapiens_results.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  7/13/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Input File(s):** `3_homosapiens_rfam_blastn_output_e5.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Output File(s):** 
#' 
#' 1. `3a_homosapiens_filtered_rfam_blastn_results.Rdata`
#' 2. `3b_homosapiens_filtered_uniqseq_rfam_blastn_results.csv`
#' 3. `3c_homosapiens_filtered_totseq_rfam_blastn_results.csv`
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
#' The objective of this script is to filter and summarize the output from the sequential BLAST queries in order to characterize the small RNA classes present in the small RNA seq data.
#' This script filters and summarizes the blast results from the Homo sapiens Rfam database query.
#' 
#' First, if there is only one hit for a sequence, return that hit.
#' If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
#' Then, retain matches with > 96% sequence identity.
#' Finally, the first filter returns sequences with the maximum bitscore.
#' 
#' The second filter returns one BLAST hit with the maximum percent identical between redundant hits.
#' Finally, if there are still remaining sequences with multiple hits, retain only the first hit. 
#' 
#' ## Install libraries
library(parallel)
library(stringr)
rm(list=ls())

#' ## Load data
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts")

system.time({
hsa<-read.table("../3_homosapiens_rfam_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})

dim(hsa)
head(hsa)
str(hsa)

#' ## Analysis
hsa$dbseq_id<-as.character(hsa$dbseq_id)
hsa$query_id<-as.character(hsa$query_id)

#' Summarize the Rfam database hits by separating the components of the "dbseq_id" column at the semicolons.
hsa$rfam_accession <- as.character(lapply(strsplit(as.character(hsa$dbseq_id), ';', fixed = TRUE), "[", 1))
hsa$rfam_biotype <- as.character(lapply(strsplit(as.character(hsa$dbseq_id), ';', fixed = TRUE), "[", 2))
hsa$rfam_seqnamestartend <- as.character(lapply(strsplit(as.character(hsa$dbseq_id), ';', fixed = TRUE), "[", 3))
head(hsa)

#' The number of unique sequences with homo sapiens Rfam hits in this dataset
length(unique(hsa$query_id))

#' Create a subset data frame containing the pertinent information for filtering
rfamsubset<-hsa[,c("query_id", "perc_identical", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
head(rfamsubset)

#' Create an index of the sequences to filter through
system.time({
idx <- mclapply(unique(rfamsubset$query_id), function(x)
        rfamsubset[as.character(rfamsubset$query_id) == x,], mc.cores=10)
})

length(idx)

#' Function to filter multiple Rfam blast hits per sequence
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


#' Apply filter to the sequence list
system.time({
stp1 <- mclapply(idx, filter, mc.cores=10)
})

length(stp1)
head(stp1)
tail(stp1)

#' How many sequences have more than one annotation after the filtering
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
table(unlist(lapply(stp2,nrow)))

filter2 <- filter <- function(seqblast){
# Sequence has only one Rfam blast hit
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
stp3 <- mclapply(stp1, filter2, mc.cores=10)
})

#' How many sequences have more than one annotation after the second filtering
stp4 <- stp3[unlist(lapply(stp3,nrow)) > 1]
length(stp4)

#' Summary file of small RNA sequence ensembl blast results
sumblast4 <- do.call(rbind,stp3)
head(sumblast4)
dim(sumblast4)
length(unique(sumblast4$query_id))

table(sumblast4$rfam_biotype)

uniqsumblast4<-as.data.frame(table(sumblast4$rfam_biotype))
colnames(uniqsumblast4)<-c("Gene_Biotype", "Freq")
uniqsumblast4

#' Add the column of sequence count to the sumblast data frame
sumblast4$seq_count<-as.numeric(str_split_fixed(sumblast4$query_id, "_x", 2)[,2])
head(sumblast4)

#' Use the "by" function to sum the sequence counts by their gene biotypes
totalsumbiotype4<-as.matrix(by(sumblast4$seq_count, sumblast4$rfam_biotype, sum))
totalsumbiotype4

if (sum(rownames(totalsumbiotype4) != uniqsumblast4$Gene_Biotype)) stop ("Gene_Biotypes not equal")

#' As a check, manually sum the 5s_rRNAs and the tRNA fields:
sum(sumblast4$seq_count[sumblast4$rfam_biotype == "5S_rRNA"])

sum(sumblast4$seq_count[sumblast4$rfam_biotype == "tRNA"])

if (sum(sumblast4$seq_count[sumblast4$rfam_biotype == "5S_rRNA"]) != totalsumbiotype4["5S_rRNA",]) stop ("5S_rRNA counts not equal")
if (sum(sumblast4$seq_count[sumblast4$rfam_biotype == "tRNA"]) != totalsumbiotype4["tRNA",]) stop ("tRNA counts not equal")

#' ## Save data
save(sumblast4, uniqsumblast4, totalsumbiotype4, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/3a_homosapiens_filtered_rfam_blastn_results.Rdata"))
write.csv(uniqsumblast4, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/3b_homosapiens_filtered_uniqseq_rfam_blastn_results.csv")
write.csv(totalsumbiotype4, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/3c_homosapiens_filtered_totseq_rfam_blastn_results.csv", row.names=TRUE)