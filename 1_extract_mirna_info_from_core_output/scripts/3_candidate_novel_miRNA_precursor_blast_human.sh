#!/bin/sh  -login
#PBS -l nodes=1:ppn=1,walltime=00:05:00,mem=1Gb
#PBS -N 1_candidate_novel_miRNA_precursor_blastn_mirbase_human_precursor
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/5_novel_miRNA_characterization
#PBS -m abe
#PBS -M perrykai@msu.edu

#' **Script:** `3_candidate_novel_miRNA_precursor_blast_human.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts`
#' 
#' **Date:**  10/12/16
#' 
#' **Input File Directory:**  
#'
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/`
#' 
#' **Input File(s):** 
#' 
#' 1. `6_candidate_novel_mature_mir.fa`
#' 2. `7_candidate_novel_precursor_mir.fa`
#' 3. `hsa_mature_mir.fa`
#' 4. `hsa_hairpin_mir.fa`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_candidate_novel_miRNA_precursor_blastn_human_output_e5.txt`
#' 2. `1_candidate_novel_miRNA_mature_blastn_human_output_e5.txt`
#' 3. ``
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Analysis](#analysis)
#' 
#' ## Objectives
#' 
#' The objective of this script is to BLAST the candidate novel precursor miRNA sequences against the known human miRNA precursor sequences.
#' The table output format will be utilized for filtering the output in R.
#' 
#' ## Install libraries
#' 
module load BLAST+/2.2.30

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/

#'
#' BLASTing the candidate novel miRNA precursor sequences against the human miRBase precursor databases:
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
#' So, for this analysis, the e-value threshold is set to 0.00001, meaning we can expect that a hit with that e-value would have 0.00001 match with a similar score by chance.
#' The lower the e-value, the more "significant" the match is.
#'
#' Tabular results:

blastn -query /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/7_candidate_novel_precursor_mir.fa -task blastn-short -db /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/"hsa_hairpin_mir.fa" -out /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/1_candidate_novel_miRNA_precursor_blastn_human_output_e5.txt -evalue 0.00001 -outfmt 6

blastn -query /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/6_candidate_novel_mature_mir.fa -task blastn-short -db /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/"hsa_mature_mir.fa" -out /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/1_candidate_novel_miRNA_mature_blastn_human_output_e5.txt -evalue 0.00001 -outfmt 6

#' Repeat the same analysis but at e-value of 1

blastn -query /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/7_candidate_novel_precursor_mir.fa -task blastn-short -db /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/"hsa_hairpin_mir.fa" -out /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/2_candidate_novel_miRNA_precursor_blastn_human_output_eval1.txt -evalue 1.0 -outfmt 6

blastn -query /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/6_candidate_novel_mature_mir.fa -task blastn-short -db /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/"hsa_mature_mir.fa" -out /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/2_candidate_novel_miRNA_mature_blastn_human_output_eval1.txt -evalue 1.0 -outfmt 6


qstat -f ${PBS_JOBID}
