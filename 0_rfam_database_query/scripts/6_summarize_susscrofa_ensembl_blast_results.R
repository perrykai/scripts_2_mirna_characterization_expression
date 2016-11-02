#' **Script:** `6_summarize_susscrofa_ensembl_blast_results.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Date:**  `6/14/16 UPDATED 7/5/16`
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Input File(s):** `1_susscrofa_ensembl_blastn_output_e5.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Output File(s):** 
#' 
#' 1. `1a_susscrofa_filtered_ensembl_blastn_results.Rdata`
#' 2. `1b_susscrofa_filtered_uniqseq_ensembl_blastn_results.csv`
#' 3. `1c_susscrofa_filtered_totseq_ensembl_blastn_results.csv`
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
#' This script filters and summarizes the blast results from the Sus scrofa Ensembl database query.
#' 
#' First, if there is only one hit for a sequence, return that hit.
#' If there is more than one hit for a sequence, retain those BLAST hits that are 100% identical to the query sequence.
#' Then, retain matches with > 96% sequence identity.
#' Finally, the first filter returns sequences with the maximum bitscore.
#' 
#' The second filter returns the sequence with greatest percent match to blast hit.
#' Finally, if there are still remaining sequences with multiple hits, retain only the first hit. 
#' 
#' ## Install libraries
#' 
library(biomaRt)
library(parallel)
library(stringr)
rm(list=ls())

#' ## Load data
#' 

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts")

system.time({
blastresults<-read.table("../1_susscrofa_ensembl_blastn_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
})

dim(blastresults)
head(blastresults)
str(blastresults)

#' ## Analysis
#' 
blastresults$dbseq_id<-as.character(blastresults$dbseq_id)
blastresults$query_id<-as.character(blastresults$query_id)

#' Extract the ensembl id to another column:
#' 
#' Move the version numbers (".1", ".2" or ".3") after the ensembl ids to another column for each sequence to make it compatible with biomaRt
#' 
#' See website http://useast.ensembl.org/Help/View?id=181 for details on the stability of ensembl IDs
blastresults$ensid<-as.character(blastresults$dbseq_id)
blastresults$ensid<-gsub("\\..*","",blastresults$ensid)
blastresults$ensid_version <- as.character(lapply(strsplit(as.character(blastresults$dbseq_id), '.', fixed = TRUE), "[", 2))

head(blastresults)
head(table(blastresults$ensid))

#' Check how many unique small RNA sequences had ensembl hits:
length(unique(blastresults$query_id))
#' Check how many unique ensembl database sequences had hits:
length(unique(blastresults$dbseq_id))

#' Use biomaRt package (v/2.20.0) to obtain annotation information for Ensembl dataset BLAST hits:
mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="dec2015.archive.ensembl.org", dataset="sscrofa_gene_ensembl")

system.time({
rst <- getBM(
attributes = c("ensembl_transcript_id","chromosome_name", "ensembl_gene_id", "rfam", "rfam_transcript_name", "version", "transcript_version", "transcript_source", "status", "transcript_status", "gene_biotype"),
filters = "ensembl_transcript_id", 
values=blastresults$ensid, 
mart = mart)
})

dim(rst)
head(rst)
tail(rst)
table(rst$rfam)
table(rst$gene_biotype)

sum(rst$rfam == "")
sum(rst$gene_biotype == "")

#' Use match function to obtain the gene_biotype and gene status information for each ensembl match and add it as a column to blastresults:
blastresults$gene_biotype<-as.character(rst[match(blastresults$ensid, rst$ensembl_transcript_id), "gene_biotype"])
blastresults$status<-rst[match(blastresults$ensid, rst$ensembl_transcript_id), "status"]

dim(blastresults)
head(blastresults)
tail(blastresults)

table(blastresults$gene_biotype)

######################################################################
#' Obtained the following code from Deborah Velez, how she filtered her BLAST results for gene annotation (/mnt/research/ernstc_lab/fat_eqtl/BLAST/code/annotation_pigoligoarray.R)
######################################################################
#' Format the ensembl blast output data.frame to a list to filter the sequences with multiple ensembl hits
# Number of unique sequences
n <- length(unique(blastresults$query_id))
n

# Create the sequence list to filter
system.time({
idx1 <- mclapply(as.character(unique(blastresults$query_id))[1:round(n/5)], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx2 <- mclapply(as.character(unique(blastresults$query_id))[(round(n/5) + 1):(2*round(n/5))], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx3 <- mclapply(as.character(unique(blastresults$query_id))[((2*round(n/5)) + 1):(3*round(n/5))], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx4 <- mclapply(as.character(unique(blastresults$query_id))[((3*round(n/5)) + 1):(4*round(n/5))], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx5 <- mclapply(as.character(unique(blastresults$query_id))[((4*round(n/5)) + 1):n], function(x)
        blastresults[as.character(blastresults$query_id) == x,], mc.cores=10)
idx <- c(idx1,idx2,idx3,idx4,idx5)
})

length(idx)

#' Function to filter multiple ensembl blast hits per sequence
filter <- function(seqblast){
# Sequence has only one ensembl blast hit
        if (nrow(seqblast) == 1){
                return(seqblast)
        }
# Sequence with 100% match to ensembl blast hit
        ident <- seqblast[sapply(seqblast$perc_identical, function(x) x == 100),]
                if (nrow(ident) == 1){
                        return(ident)
        }
# Sequence with greater than 96% match to ensembl blast hit
        ident2 <- seqblast[sapply(seqblast$perc_identical, function(x) x > 96),]
                if (nrow(ident2) == 1){
                        return(ident2)
        }
# Select ensembl blast hit with the greatest bitscore for a sequence
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
# Sequence has only one ensembl blast hit
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

#' Summary file of small RNA sequence ensembl blast results
sumblast2 <- do.call(rbind,stp3)
head(sumblast2)
dim(sumblast2)
length(unique(sumblast2$query_id))

table(sumblast2$gene_biotype)

uniqsumblast2<-as.data.frame(table(sumblast2$gene_biotype))
colnames(uniqsumblast2)<-c("Gene_Biotype", "Freq")
uniqsumblast2

#' Add the column of sequence count to the sumblast2 data frame
sumblast2$seq_count<-as.numeric(str_split_fixed(sumblast2$query_id, "_x", 2)[,2])
head(sumblast2)

#' Use the "by" function to sum the sequence counts by their gene biotypes
totalsumbiotype2<-as.matrix(by(sumblast2$seq_count, sumblast2$gene_biotype, sum))
totalsumbiotype2

if (sum(rownames(totalsumbiotype2) != uniqsumblast2$Gene_Biotype)) stop ("Gene_Biotypes not equal")

#' As a check, manually sum the miRNAs and the processed_transcript fields:
sum(sumblast2$seq_count[sumblast2$gene_biotype == "miRNA"])

sum(sumblast2$seq_count[sumblast2$gene_biotype == "processed_transcript"])

if (sum(sumblast2$seq_count[sumblast2$gene_biotype == "miRNA"]) != totalsumbiotype2["miRNA",]) stop ("miRNA counts not equal")
if (sum(sumblast2$seq_count[sumblast2$gene_biotype == "processed_transcript"]) != totalsumbiotype2["processed_transcript",]) stop ("processed_transcript counts not equal")

#' ## Save data
save(sumblast2, uniqsumblast2, totalsumbiotype2, file=("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1a_susscrofa_filtered_ensembl_blastn_results.Rdata"))
write.csv(uniqsumblast2, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1b_susscrofa_filtered_uniqseq_ensembl_blastn_results.csv")
write.csv(totalsumbiotype2, file="/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1c_susscrofa_filtered_totseq_ensembl_blastn_results.csv", row.names=TRUE)