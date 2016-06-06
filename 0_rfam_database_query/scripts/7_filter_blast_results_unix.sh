#' **Script:** `7_filter_blast_results_unix.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  5/31/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Input File(s):** 
#' 
#' 1. `174_library_uniq_exp.fa`
#' 2. `1_susscrofa_ensembl_blastn_output_e6.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Output File(s):** 
#' 
#' 1. `174_library_uniq_exp_ids.txt`
#' 2. `susscrofa_ensembl_blastn_output_e6_ids.txt`
#' 3. `nonensembl_rfam_susscrofa_input_ids.txt`
#' 4. `174_library_nonensembl_uniq_seq.fa`
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
#' The databases used for this query were created from the Ensembl Sus Scrofa (10.2) ncRNA database, and the Rfam database version 11.0.
#' See script /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/4_download_ensembl_rfam_build_blastdb.sh for code building those databases.
#' First, the sequences will be blasted against the Sus scrofa Ensembl noncoding sequences, then move sequentially to the Sus scrofa, Homo sapiens, and Mus musculus Rfam sequences. 
#' This script is used to filter the output of the ensembl blast to filter out the sequences that did have ensembl hits. 
#' 
#' 
#' 
#' ## Install libraries

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#' Create an index of the sequence ids from the unique sequences
grep '^>seq_*' 174_library_uniq_exp.fa |cut -c2- | sort | uniq > 174_library_uniq_exp_ids.txt

#' Subset the ensembl blast result sequence IDs from the blast output
cat 1_susscrofa_ensembl_blastn_output_e6.txt |cut -f1 | sort | uniq > susscrofa_ensembl_blastn_output_e6_ids.txt

#' Subset the non-ensembl blast results from the 174_library_uniq_exp.fa (will be input for Rfam sus scrofa blast)
#'
#' grep -F option: Fixed strings: interpret the pattern as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
#' grep -v option: invert match
#' grep -f option: file
grep -F -v -f susscrofa_ensembl_blastn_output_e6_ids.txt 174_library_uniq_exp_ids.txt > nonensembl_rfam_susscrofa_input_ids.txt

#' As a check, make sure that no sequence ids match between the nonensembl-hit ID file and the ensembl hit ID file
grep -F -f nonensembl_rfam_susscrofa_input_ids.txt susscrofa_ensembl_blastn_output_e6_ids.txt

#' Additionally, the line count from the ensembl-hit file and the nonensembl-hit files should add up to the total number of sequence IDs.
wc -l nonensembl_rfam_susscrofa_input_ids.txt
# 1317528 nonensembl_rfam_susscrofa_input_ids.txt
wc -l susscrofa_ensembl_blastn_output_e6_ids.txt
# 127885 susscrofa_ensembl_blastn_output_e6_ids.txt
wc -l 174_library_uniq_exp_ids.txt
# 1445413 174_library_uniq_exp_ids.txt

#' Then, use the same awk script as previously used to extract the Rfam Sus scrofa sequences to reduce the previously used sequences, leaving the non-ensembl sequences:
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' 174_library_uniq_exp.fa | awk -F"\t" 'BEGIN{while((getline k < "nonensembl_rfam_susscrofa_input_ids.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 174_library_nonensembl_uniq_seq.fa

#' The length of the 174_library_nonensembl_uniq_seq.fa should be twice the length of the nonensembl ID file
wc -l 174_library_nonensembl_uniq_seq.fa
# 2635056 174_library_nonensembl_uniq_seq.fa

wc -l nonensembl_rfam_susscrofa_input_ids.txt
# 1317528 nonensembl_rfam_susscrofa_input_ids.txt

wc -l 174_library_uniq_exp.fa
# 2890826 174_library_uniq_exp.fa
# 2890826 - 2635056 = 255770
# 255770/2 = 127885

wc -l susscrofa_ensembl_blastn_output_e6_ids.txt
# 127885 susscrofa_ensembl_blastn_output_e6_ids.txt


#' 
# We are first interested in the ensembl blast matches; so, delete entries from Rfam results if ensembl matches exist for that seq:
# id=(`cat blast_seq_ids.txt`)

# for ((i=0; i<${#id[@]} ; i++)) do
# if grep --quiet ${id[$i]} ensembl_ids.txt; then
# grep ${id[$i]} ensembl_blast_results.txt >> reduced_blast.txt
# else
# grep ${id[$i]} rfam_blast_results.txt >> reduced_blast.txt
# fi
# done
