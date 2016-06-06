q#' **Script:** `6_summarize_blast_results.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Date:**  `5/18/16`
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Input File(s):** ``
#' 
#' **Output File Directory:** ``
#' 
#' **Output File(s):** ``
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' 
#' ## Install libraries
#' 
library("biomaRt")

rm(list=ls())

#' ## Load data
#' 
system.time(
blastresults<-read.table("../1_susscrofa_ensembl_blastn_output_e6.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
)
dim(blastresults)
head(blastresults)
str(blastresults)

#' ## Analysis
#' 
blastresults$dbseq_id<-as.character(blastresults$dbseq_id)

#' Extract the ensembl id to another column:
#' 
#' Move the version numbers (".1", ".2" or ".3") after the ensembl ids to another column for each sequence to make it compatible with biomaRt
#' 
#' See website http://useast.ensembl.org/Help/View?id=181 for details on the stability of ensembl IDs
blastresults$ensid<-as.character(blastresults$dbseq_id)
blastresults$ensid<-gsub("\\..*","",blastresults$ensid)
blastresults$ensid_version <- lapply(strsplit(as.character(blastresults$dbseq_id), '.', fixed = TRUE), "[", 2)

head(blastresults)
dim(table(blastresults$ensid))
head(table(blastresults$ensid))

#' Check how many unique small RNA sequences had ensembl hits:
length(unique(blastresults$query_id))
#' Check how many unique ensembl database sequences had hits:
length(unique(blastresults$dbseq_id))

#' Use biomaRt package (v/2.20.0) to obtain annotation information for Ensembl dataset BLAST hits:
mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="dec2015.archive.ensembl.org", dataset="sscrofa_gene_ensembl")

system.time(
rst <- getBM(
attributes = c("ensembl_transcript_id","chromosome_name", "ensembl_gene_id", "rfam", "rfam_transcript_name", "version", "transcript_version", "transcript_source", "status", "transcript_status", "gene_biotype", "transcript_biotype"),
filters = "ensembl_transcript_id", 
values=blastresults$ensid, 
mart = mart)
)

dim(rst)
head(rst)
tail(rst)
table(rst$rfam)
table(rst$gene_biotype)

sum(rst$rfam == "")
sum(rst$gene_biotype == "")
sum(rst$status != rst$transcript_status)
sum(rst$gene_biotype != rst$transcript_biotype)
table(rst$transcript_biotype)
table(rst$gene_biotype)

#' Use match function to obtain the gene_biotype, transcript_biotype, gene status, and transcript_status information for each ensembl match and add it as a column to blastresults:
blastresults$gene_biotype<-rst[match(blastresults$ensid, rst$ensembl_transcript_id), "gene_biotype"]
blastresults$transcript_biotype<-rst[match(blastresults$ensid, rst$ensembl_transcript_id), "transcript_biotype"]
blastresults$status<-rst[match(blastresults$ensid, rst$ensembl_transcript_id), "status"]
blastresults$transcript_status<-rst[match(blastresults$ensid, rst$ensembl_transcript_id), "transcript_status"]

dim(blastresults)
head(blastresults)
tail(blastresults)

table(blastresults$gene_biotype)
table(blastresults$transcript_biotype)

sum(blastresults$status != blastresults$transcript_status)

amborf<-blastresults[(blastresults$transcript_biotype == "ambiguous_orf"),c(1,16)]
dim(amborf)
length(unique(amborf$query_id))

retint<-blastresults[(blastresults$transcript_biotype == "retained_intron"),c(1,16)]
dim(retint)
length(unique(retint$query_id))

protran<-blastresults[(blastresults$transcript_biotype == "processed_transcript"),c(1,16)]
dim(protran)
length(unique(protran$query_id))

genebio<-blastresults[(blastresults$gene_biotype == "processed_transcript"),c(1,16)]
dim(genebio)

genebio$protranmatch<-match(genebio$query_id, protran$query_id)
genebio$retintmatch<-match(genebio$query_id, retint$query_id)
genebio$amborfmatch<-match(genebio$query_id, amborf$query_id)
head(genebio)
genebio
genebio[1:50,]


length(unique(genebio$query_id))


#' Summarize the Rfam database hits by separating the components of the "dbseq_id" column at the semicolons.
rfamresults$rfam_accession <- lapply(strsplit(as.character(rfamresults$dbseq_id), ';', fixed = TRUE), "[", 1)
rfamresults$rfam_biotype <- lapply(strsplit(as.character(rfamresults$dbseq_id), ';', fixed = TRUE), "[", 2)
#rfamresults$rfam_seqnamestartend <- lapply(strsplit(as.character(rfamresults$dbseq_id), ';', fixed = TRUE), "[", 3)

head(rfamresults)

rfamsubset<-rfamresults[,c("query_id", "perc_identical", "mismatch", "gapopen", "evalue", "bitscore", "rfam_accession", "rfam_biotype")]
dim(rfamsubset)
head(rfamsubset)

#' Subset the biotype column only to see if any duplicates remain when unique is run on the dataframe
rfambiotype<-rfamsubset[,c("query_id", "rfam_biotype")]
dim(rfambiotype)
head(rfambiotype)
rfambiouniq<-unique(rfambiotype)
dim(rfambiouniq)
head(rfambiouniq)

#' 
#' Use unique to create the list of sequences with ensembl results, then create a positional index for the rfam results:
ensseq<-unique(blastresults$query_id)
length(ensseq)
head(ensseq)

#' Match the positions of the 
test<-ensseq %in% rfamsubset$query_id
head(test)
length(test)
nrow(rfamsubset)

rfamnoens<-subset(rfamsubset, !(test))
dim(rfamnoens)
head(rfamnoens)
tail(rfamnoens)

#' 
#' Order the rfamresults by the seq identifier then the evalue
#system.time(
#rfamordered<- rfamresults[with(rfamresults, order(query_id, evalue)), ]
#)
#dim(rfamordered)
#head(rfamordered)

#' 
#' Next, subset the rfam results dataframe, essentially removing any sequence ids that match the sequence id in the ensembl database

#' ## Visualize
#' 

#' ## Save data
