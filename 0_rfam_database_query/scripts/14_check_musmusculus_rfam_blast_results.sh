#' **Script:** `14_check_musmusculus_rfam_blast_results.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  7/13/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Input File(s):** 
#' 
#' 1. `4_musmusculus_rfam_blastn_output_e5.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Output File(s):** 
#' 
#' 1. `musmusculus_rfam_blastn_output_e5_ids.txt`
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
#' The objective of this analysis is to utilize BLAST+ to characterize the ncRNA species (tRNA, sn/snoRNA, miRNA, etc) present in the small RNA seq dataset.
#' The databases used for this query were created from the Ensembl homo sapiens (10.2) ncRNA database, and the Rfam database version 11.0.
#' See script /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/4_download_ensembl_rfam_build_blastdb.sh for code building those databases.
#' First, the sequences will be blasted against the homo sapiens Ensembl noncoding sequences, then move sequentially to the homo sapiens, Homo sapiens, and Mus musculus Rfam sequences. 
#' This script is used to check the output of the mus musculus Rfam blast results  for redundancies between the previous analyses (check that no sequences in this output appear in previous output)
#' 
#' ## Install libraries

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#' The index of the sequence ids from the unique sequences input into the Rfam analysis already exists, as the file nonensembl_rfam_susscrofa_input_ids.txt:
#' 
#' Subset the ensembl blast result sequence IDs from the blast output
cat 4_musmusculus_rfam_blastn_output_e5.txt |cut -f1 | sort | uniq > musmusculus_rfam_blastn_output_e5_ids.txt

#' Check that the mus musculus sequence ids do not match previous output sequence ids.
#'
#' grep -F option: Fixed strings: interpret the pattern as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
#' grep -v option: invert match
#' grep -f option: file
#' 
#' As a check, make sure that no sequence ids match between the the previous mirbase, ensembl and rfam-hit ID files
grep -F -f musmusculus_rfam_blastn_output_e5_ids.txt homosapiens_rfam_blastn_output_e5_ids.txt

grep -F -f musmusculus_rfam_blastn_output_e5_ids.txt susscrofa_rfam_blastn_output_e5_ids.txt

grep -F -f musmusculus_rfam_blastn_output_e5_ids.txt susscrofa_ensembl_blastn_output_e5_ids.txt

grep -F -f musmusculus_rfam_blastn_output_e5_ids.txt susscrofa_mirbase_blastn_output_e5_ids.txt

#' Additionally, count how many total reads are present in the dataset (sum all read counts in 174_library_uniq_exp.fa)
grep '_x' 174_library_uniq_exp.fa | cut -d x -f2 | paste -sd+ - | bc
# 227425784


