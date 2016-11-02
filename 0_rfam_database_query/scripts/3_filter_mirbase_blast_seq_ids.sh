#' **Script:** `3_filter_mirbase_blast_seq_ids.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  7/11/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Input File(s):** 
#' 
#' 1. `174_library_uniq_exp.fa`
#' 2. `0_susscrofa_mirbase_blastn_output_e5.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Output File(s):** 
#' 
#' 1. `174_library_uniq_exp_ids.txt`
#' 2. `susscrofa_mirbase_blastn_output_e5_ids.txt`
#' 3. `nonmirbase_susscrofa_ensembl_input_ids.txt`
#' 4. `174_library_nonmirbase_uniq_seq.fa`
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
#' The databases used for this query were created from the Sus scrofa miRBase, the Ensembl Sus Scrofa (10.2) ncRNA database, and the Rfam database version 11.0.
#' See script /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/1_download_ensembl_rfam_build_blastdb_uniq_exp.sh for code building those databases.
#' First, the sequences will be blasted against the Sus scrofa miRBase database, then move sequentially to the Sus scrofa Ensembl, and Sus scrofa, Homo sapiens, and Mus musculus Rfam sequences. 
#' This script is used to filter the output of the mirbase blast results to filter out the sequences that did have mirbase hits. 
#' 
#' 
#' 
#' ## Install libraries

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#' Create an index of the sequence ids from the unique sequences
grep '^>seq_*' 174_library_uniq_exp.fa |cut -c2- | sort | uniq > 174_library_uniq_exp_ids.txt

#' Subset the mirbase blast result sequence IDs from the blast output
cat 0_susscrofa_mirbase_blastn_output_e5.txt |cut -f1 | sort | uniq > susscrofa_mirbase_blastn_output_e5_ids.txt

#' Subset the non-mirbase blast results from the 174_library_uniq_exp.fa (will be input for Sus scrofa Ensembl blast)
#'
#' grep -F option: Fixed strings: interpret the pattern as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
#' grep -v option: invert match
#' grep -f option: file
grep -F -v -f susscrofa_mirbase_blastn_output_e5_ids.txt 174_library_uniq_exp_ids.txt > nonmirbase_susscrofa_ensembl_input_ids.txt

#' As a check, make sure that no sequence ids match between the nonmiRBase-hit ID file and the miRBase hit ID file
grep -F -f nonmirbase_susscrofa_ensembl_input_ids.txt susscrofa_mirbase_blastn_output_e5_ids.txt
# No results

#' Additionally, the line count from the miRBase-hit file and the nonmiRBase-hit files should add up to the total number of sequence IDs.
wc -l nonmirbase_susscrofa_ensembl_input_ids.txt
# 1297756 nonmirbase_susscrofa_ensembl_input_ids.txt
wc -l susscrofa_mirbase_blastn_output_e5_ids.txt
# 147657 susscrofa_mirbase_blastn_output_e5_ids.txt
wc -l 174_library_uniq_exp_ids.txt
# 1445413 174_library_uniq_exp_ids.txt

#' Then, use the awk script obained from Biostars to extract the Sus scrofa Ensembl input sequences to reduce the previously used sequences, leaving the non-miRBase sequences:
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' 174_library_uniq_exp.fa | awk -F"\t" 'BEGIN{while((getline k < "nonmirbase_susscrofa_ensembl_input_ids.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 174_library_ensembl_blast_input_uniq_seq.fa

#' The length of the 174_library_nonensembl_uniq_seq.fa should be twice the length of the nonensembl ID file
wc -l 174_library_ensembl_blast_input_uniq_seq.fa
# 2595512 174_library_ensembl_blast_input_uniq_seq.fa

wc -l 174_library_uniq_exp.fa
# 2890826 174_library_uniq_exp.fa

#' The number of lines in the complete fasta file (174_library_uniq_exp.fa) minus the number of lines in the ensembl input fasta file (174_library_ensembl_blast_input_uniq_seq.fa)
#' should be equivalent to the number of sequences that had miRBase blast hits (number of IDs then would be that divided by two)
# 2890826 - 2595512 = 295314
# 295314/2 = 147657 (Number of IDs in the miRBase blast hit file)

#1445413 - 1297756 = 147657