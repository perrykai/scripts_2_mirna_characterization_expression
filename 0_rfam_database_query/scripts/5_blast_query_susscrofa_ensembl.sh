#!/bin/sh  -login
#PBS -l nodes=1:ppn=9,walltime=02:00:00,mem=5Gb
#PBS -N blastn_sscrofa_ensembl_query
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/
#PBS -m abe
#PBS -M perrykai@msu.edu

#' **Script:** `5_blast_query_susscrofa_ensembl.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  `5/31/16 (full dataset run 5/20/16)`
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Input File(s):** `174_library_uniq_exp.fa`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/`
#' 
#' **Output File(s):** `1_susscrofa_ensembl_blastn_output_e6.txt`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 
#' ## Objectives
#' 
#' The objective of this script is to utilize BLAST+ to characterize the ncRNA species (tRNA, sn/snoRNA, miRNA, etc) present in the small RNA seq dataset.
#' The databases used for this query were created from the Ensembl Sus Scrofa (10.2) ncRNA database, and the Rfam database version 11.0.
#' See script /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/4_download_ensembl_rfam_build_blastdb.sh for code building those databases.
#' First, the sequences will be blasted against the Sus scrofa Ensembl noncoding sequences, then move sequentially to the Sus scrofa, Homo sapiens, and Mus musculus Rfam sequences. 
#' 
#' ## Install libraries
module load BLAST+/2.2.30

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/

#' 
#' BLASTing the small RNAseq unique sequences (~6.8 million sequences) against the Sus scrofa Rfam and Ensembl databases:
#' 
#' Parameters used: (See websites http://www.ncbi.nlm.nih.gov/BLAST/Why.shtml, http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ, and documentation http://www.ncbi.nlm.nih.gov/books/NBK279690/ for details)
#' 
#' "blastn" option used for nucleotide BLAST:
#' 
#' word-size: 7 (length of initial exact match); 
#' reward: 1 (reward for a nt match);
#' penalty: -3 (penalty for a nt mismatch);
#' 
#' "query" defines the input dataset (in this case, the small RNA seq data)
#' 
#' "task" option defines which type of blastn to run; using blastn-short option which is optimized for sequences shorter than 50 bases
#' 
#' "db" defines the databases used for the BLAST query: multiple databases can be utilized by including space-separated database names in quotes
#' 
#' "out" defines the name of the output file
#' 
#' "evalue" defines the 'expect value' used in the analysis, and can be used as a measure for the quality of the match between query and database.
#' The "expect value" is a parameter describing the number of hits one could expect "by chance" when searching a database of a particular size. Essentially, it describes the random background noise.
#' 
#' "For example, an E value of 1 assigned to a hit can be interpreted as meaning that in a database of the current size one might expect to see 1 match with a similar score simply by chance."
#' So, for this analysis, the e-value threshold is set to 0.000001, meaning we can expect that a hit with that e-value would have 0.000001 match with a similar score by chance.
#' The lower the e-value, the more "significant" the match is. 

#' 
#' Tabular results:
blastn -query /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/174_library_uniq_exp.fa -task blastn-short -db /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/"Sus_scrofa.Sscrofa10.2.ncrna.fa" -out /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/1_susscrofa_ensembl_blastn_output_e6.txt -evalue 0.000001 -outfmt 6 -num_threads 8

qstat -f ${PBS_JOBID}

