#' **Script:** `4_filter_candidate_novel_miRNA_blast_results.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts/`
#' 
#' **Date:**  `10/24/16`
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/`
#' 
#' **Input File(s):** 
#' 
#' 1. `1_candidate_novel_miRNA_precursor_blastn_human_output_e5.txt`
#' 2. `1_candidate_novel_miRNA_mature_blastn_human_output_e5.txt`
#' 3. `2_candidate_novel_miRNA_precursor_blastn_human_output_eval1.txt`
#' 4. `2_candidate_novel_miRNA_mature_blastn_human_output_eval1.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/`
#' 
#' **Output File(s):** 
#' 
#' 1. `3_filtered_candidate_novel_miRNA_precursor_abundance_e5.txt`
#' 2. `4_filtered_candidate_novel_precursor_eval1.txt`
#' 3. `5_filtered_candidate_novel_mature_eval1.txt`
#' 4. `6_filtered_candidate_novel_miRNA_precursor_ids_e5.txt`
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
#' The objective of this script is to summarize the output from the BLAST query in order to characterize the candidate novel miRNAs present in the small RNA seq data.
#' Additionally, I want to filter the miRNA BLAST results with the e-value = 1.0, to remove some of the redundant hits and make the results more managable. 
#' 
#' So, need to load the precursor BLAST dataset and the full dataset of candidate novel miRNA and compare the names of sequences to see if the most abundant miRNA had BLAST results.
#' 
#' ## Install libraries
rm(list=ls())

#' ## Load data
setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts")

hsa.blast.precursor<-read.table("../8_blastn_candidate_novel_miRNA_output/1_candidate_novel_miRNA_precursor_blastn_human_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
mature.hsa.blast<-read.table("../8_blastn_candidate_novel_miRNA_output/1_candidate_novel_miRNA_mature_blastn_human_output_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
eval1.hsa.blast.precursor<-read.table("../8_blastn_candidate_novel_miRNA_output/2_candidate_novel_miRNA_precursor_blastn_human_output_eval1.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
eval1.hsa.blast.mature<-read.table("../8_blastn_candidate_novel_miRNA_output/2_candidate_novel_miRNA_mature_blastn_human_output_eval1.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
#' --------------------
dim(hsa.blast.precursor)
head(hsa.blast.precursor)
str(hsa.blast.precursor)

hsa.blast.precursor$dbseq_id<-as.character(hsa.blast.precursor$dbseq_id)
hsa.blast.precursor$query_id<-as.character(hsa.blast.precursor$query_id)
#' --------------------
dim(mature.hsa.blast)
head(mature.hsa.blast)
str(mature.hsa.blast)

mature.hsa.blast$dbseq_id<-as.character(mature.hsa.blast$dbseq_id)
mature.hsa.blast$query_id<-as.character(mature.hsa.blast$query_id)
#' --------------------
dim(eval1.hsa.blast.precursor)
head(eval1.hsa.blast.precursor)
str(eval1.hsa.blast.precursor)

eval1.hsa.blast.precursor$dbseq_id<-as.character(eval1.hsa.blast.precursor$dbseq_id)
eval1.hsa.blast.precursor$query_id<-as.character(eval1.hsa.blast.precursor$query_id)
#' --------------------
dim(eval1.hsa.blast.mature)
head(eval1.hsa.blast.mature)
str(eval1.hsa.blast.mature)

eval1.hsa.blast.mature$dbseq_id<-as.character(eval1.hsa.blast.mature$dbseq_id)
eval1.hsa.blast.mature$query_id<-as.character(eval1.hsa.blast.mature$query_id)
#' --------------------
load("../5_putative_novel_miRNA_filtered_candidates.Rdata")

dim(novelmir10sigrandfoldmincounts)
names(novelmir10sigrandfoldmincounts)
#' Make the consensus.mature.sequence column into a character vector to count the length of the strings
novelmir10sigrandfoldmincounts$consensus.mature.sequence<-as.character(novelmir10sigrandfoldmincounts$consensus.mature.sequence)
rownames(novelmir10sigrandfoldmincounts)<- novelmir10sigrandfoldmincounts$provisional.id
head(novelmir10sigrandfoldmincounts)
novelmir10sigrandfoldmincounts$consensus.mature.sequence.length<-nchar(novelmir10sigrandfoldmincounts$consensus.mature.sequence)
head(novelmir10sigrandfoldmincounts)
#' ## Analysis
#' 
#' ### 1. The goal is to compare the candidate novel miRNAs with BLAST results to the most abundant candidate novel miRNAs.
#' 
#' First, return the name of the sequences to the way they were before BLASTing at e-value = 1x10^-5.

hsa.blast.precursor$seqname<-sapply(strsplit(hsa.blast.precursor$query_id, "|", fixed=TRUE),'[',2)
head(hsa.blast.precursor)

#' Identify the unique number of sequences with BLAST hits at eval = 1x10^-5
seqids<-unique(hsa.blast.precursor$seqname)
length(seqids)
seqids

#' Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset
match(seqids, rownames(novelmir10sigrandfoldmincounts))
rownames(novelmir10sigrandfoldmincounts)%in%seqids

#' Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.
novelmir.abundance<-novelmir10sigrandfoldmincounts[seqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count", "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]

#' Is this object ordered by miRDeep2 score or by total.read.count?
sum(rownames(novelmir.abundance[order(novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(novelmir.abundance))
sum(rownames(novelmir.abundance[order(novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(novelmir.abundance))

#' Combine the information; are the predicted miRNA precursors with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?
#' 
#' This object (candidate novel precursors with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.
novelmir.abundance<-cbind(novelmir.abundance[match(hsa.blast.precursor$seqname, rownames(novelmir.abundance)),], hsa.blast.precursor)
sum(novelmir.abundance$seqname != novelmir.abundance$provisional.id)
head(novelmir.abundance)

novelmir.abundance<-novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(novelmir.abundance)
head(novelmir.abundance)
novelmir.provisional.ids<-unique(novelmir.abundance$provisional.id)

#' ---------------------------------------
#' 
#' Now to repeat the same analysis on the candidate novel mature miRNAs with BLAST results at an e-value of 1x10-5
#' 
#' First, return the name of the mature sequences to the way they were before BLASTing at e-value = 1x10^-5.

mature.hsa.blast$seqname<-sapply(strsplit(mature.hsa.blast$query_id, "|", fixed=TRUE),'[',2)
head(mature.hsa.blast)

#' Identify the unique number of sequences with BLAST hits at eval = 1x10^-5
matureseqids<-unique(mature.hsa.blast$seqname)
length(matureseqids)
matureseqids

#' Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset
match(matureseqids, rownames(novelmir10sigrandfoldmincounts))
rownames(novelmir10sigrandfoldmincounts)%in%matureseqids

#' Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.
mature.novelmir.abundance<-novelmir10sigrandfoldmincounts[matureseqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]

#' Is this object ordered by miRDeep2 score or by total.read.count?
sum(rownames(mature.novelmir.abundance[order(mature.novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(mature.novelmir.abundance))
sum(rownames(mature.novelmir.abundance[order(mature.novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(mature.novelmir.abundance))

#' Combine the information; are the predicted miRNAs with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?
#' 
#' This object (candidate novels with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.
mature.novelmir.abundance<-cbind(mature.novelmir.abundance[match(mature.hsa.blast$seqname, rownames(mature.novelmir.abundance)),], mature.hsa.blast)
sum(mature.novelmir.abundance$seqname != mature.novelmir.abundance$provisional.id)
head(mature.novelmir.abundance)

mature.novelmir.abundance<-mature.novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(mature.novelmir.abundance)
head(mature.novelmir.abundance)
mature.novelmir.provisional.ids<-unique(mature.novelmir.abundance$provisional.id)

#' ---------------------------------------

#' ### 2. Create a subset data frame containing the pertinent information for filtering the miRNA precursor sequences BLASTed at eval = 1.0
eval1.precursor.subset<-eval1.hsa.blast.precursor[,c("query_id", "dbseq_id", "perc_identical", "evalue", "bitscore")]
dim(eval1.precursor.subset)
head(eval1.precursor.subset)

#' The number of unique sequences with homo sapien miRBase hits in this dataset at eval = 1.0
length(unique(eval1.precursor.subset$query_id))

#' The number of unique homo sapiens miRBase miRNAs identified by this BLAST query at eval = 1.0
length(unique(eval1.precursor.subset$dbseq_id))

######################################################################
#' Obtained the following code from Deborah Velez, how she filtered her BLAST results for gene annotation (/mnt/research/ernstc_lab/fat_eqtl/BLAST/code/annotation_pigoligoarray.R)
######################################################################
#' Format the miRBase blast output data.frame to a list to filter the sequences with multiple miRBase hits
# Number of unique sequences
length(unique(eval1.precursor.subset$query_id))

# Create the sequence list to filter
idx <- lapply(as.character(unique(eval1.precursor.subset$query_id)), function(x)
        eval1.precursor.subset[as.character(eval1.precursor.subset$query_id) == x,])

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
#        bit <- seqblast[seqblast$bitscore == max(seqblast$bitscore),]
#                        if (nrow(bit) == 1){
#                        return(bit)
#        }

        return(ident2)
}

#' Apply filter to the candidate novel precursor sequence list at eval = 1.0
stp1 <- lapply(idx, filter)

length(stp1)
head(stp1)
tail(stp1)

#' How many candidate novel precursor sequences have more than one hit after the filtering?
stp2 <- stp1[unlist(lapply(stp1,nrow)) > 1]
length(stp2)
table(unlist(lapply(stp2,nrow)))
head(stp2)

#' Summary file of candidate novel miRNA human miRBase blast results at eval = 1.0
sumblast <- do.call(rbind,stp1)
head(sumblast)
dim(sumblast)
#' How many unqiue candidate novel miRNA precursors had miRBase hits at eval = 1.0
length(unique(sumblast$query_id))
#' How many unique human miRBase miRNAs are identified in this dataset at eval = 1.0?
length(unique(sumblast$dbseq_id))
head(table(sumblast$query_id))
head(table(sumblast$dbseq_id))


#' ---------------------------------------
#' 
#' ### 3. Repeat the same filtering step on the candidate novel mature miRNA at eval = 1.0

#' Create a subset data frame containing the pertinent information for filtering the novel mature miRNA sequences BLASTed at eval = 1.0
length(unique(eval1.hsa.blast.mature$query_id))
eval1.mature.subset<-eval1.hsa.blast.mature[,c("query_id", "dbseq_id", "perc_identical", "evalue", "bitscore")]
dim(eval1.mature.subset)
head(eval1.mature.subset)

#' The number of unique candidate novel mature miRNA sequences with human miRBase hits in this dataset at eval = 1.0
length(unique(eval1.mature.subset$query_id))

#' The number of unique human miRBase miRNAs identified by this BLAST query at eval = 1.0
length(unique(eval1.mature.subset$dbseq_id))

#' Format the miRBase blast output data.frame to a list to filter the sequences with multiple miRBase hits
# Number of unique candidate mature miRNA sequences at eval = 1.0
length(unique(eval1.mature.subset$query_id))

# Create the unique candidate mature miRNA sequence list to filter at eval = 1.0
idx2 <- lapply(as.character(unique(eval1.mature.subset$query_id)), function(x)
        eval1.mature.subset[as.character(eval1.mature.subset$query_id) == x,])

length(idx2)

#' Apply filter to the unique candidate mature miRNA sequence list to filter at eval = 1.0
stp1m <- lapply(idx2, filter)

length(stp1m)
head(stp1m)
tail(stp1m)

#' How many sequences have more than one hit after the filtering
stp2m <- stp1m[unlist(lapply(stp1m,nrow)) > 1]
length(stp2m)
table(unlist(lapply(stp2m,nrow)))

#' Summary file of small RNA sequence miRBase blast results
sumblast.mature <- do.call(rbind,stp1m)
head(sumblast.mature)
dim(sumblast.mature)
#' How many unique homo sapiens miRNAs were identified in this BLAST query at eval = 1.0
length(unique(sumblast.mature$dbseq_id))
#' How many unique candidate novel mature miRNAs were identified in this BLAST query at eval = 1.0
length(unique(sumblast.mature$query_id))
#' Notice that there are two fewer candidate novel miRNAs than prior to filtering;
#' 
#' This is because there are two BLAST hits here that are less than 96% identical and were removed in the filtering process.
#' 
#' These include the sequence 13_28851 and 14_31131
#' 
#' Just to check:
eval1.mature.subset[9:11,]
eval1.mature.subset[111:113,]

#' ## Save data
write.table(novelmir.abundance, file="../8_blastn_candidate_novel_miRNA_output/3_candidate_novel_miRNA_precursor_abundance_e5.txt")
write.table(mature.novelmir.abundance, file="../8_blastn_candidate_novel_miRNA_output/3_candidate_novel_miRNA_mature_abundance_e5.txt")
write.table(sumblast, file="../8_blastn_candidate_novel_miRNA_output/4_filtered_candidate_novel_precursor_eval1.txt")
write.table(sumblast.mature, file="../8_blastn_candidate_novel_miRNA_output/5_filtered_candidate_novel_mature_eval1.txt")
#' Create a list of the pertinent precursor provisional.ids to extract the correct pdfs for a supplemental figure
write.table(novelmir.provisional.ids, file="../8_blastn_candidate_novel_miRNA_output/6_filtered_candidate_novel_miRNA_precursor_ids_e5.txt", row.names=FALSE, col.names=FALSE)
write.table(mature.novelmir.provisional.ids, file="../8_blastn_candidate_novel_miRNA_output/6_filtered_candidate_novel_miRNA_mature_ids_e5.txt", row.names=FALSE, col.names=FALSE)