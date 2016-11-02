#' **Script:** `9_filter_sscrofa_rfam_blast_results.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  7/12/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Input File(s):** 
#' 
#' 1. `2_susscrofa_rfam_blastn_output_e5.txt`
#' 2. `174_library_susscrofa_rfam_blast_input_uniq_seq.fa`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Output File(s):** 
#' 
#' 1. `susscrofa_rfam_blastn_output_e5_ids.txt`
#' 2. `nonsusscrofa_homosapiens_rfam_input_ids.txt`
#' 3. `174_library_homosapiens_rfam_blast_input_uniq_seq.fa`
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
#' The databases used for this query were created from the Sus Scrofa miRBase (release 21), Ensembl Sus Scrofa (10.2) ncRNA database, and the Rfam database version 11.0.
#' See script /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/4_download_ensembl_rfam_build_blastdb.sh for code building those databases.
#' First, the sequences will be blasted against the Sus scrofa miRBase and Ensembl noncoding sequences, then move sequentially to the Sus scrofa, Homo sapiens, and Mus musculus Rfam sequences. 
#' This script is used to filter the output of the sus scrofa Rfam blast results to filter out the sequences that did have sus scrofa rfam hits and create the input for the homo sapiens BLAST query
#' 
#' 
#' 
#' ## Install libraries

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#' The index of the sequence ids from the unique sequences input into the Sus scrofa Rfam analysis already exists, as the file nonensembl_susscrofa_rfam_input_ids.txt:
#' 
#' Subset the sus scrofa rfam blast result sequence IDs from the blast output
cat 2_susscrofa_rfam_blastn_output_e5.txt |cut -f1 | sort | uniq > susscrofa_rfam_blastn_output_e5_ids.txt

#' Subset the non-susscrofa rfam blast results from the 174_library_susscrofa_rfam_blast_input_uniq_seq.fa (will be input for Rfam homo sapiens blast)
#'
#' grep -F option: Fixed strings: interpret the pattern as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
#' grep -v option: invert match
#' grep -f option: file
#' 
#' First, extract the IDs that were NOT sus scrofa Rfam blast hits by doing the inverse match of the sus scrofa Rfam blast hits and the input sequence ids:
grep -F -v -f susscrofa_rfam_blastn_output_e5_ids.txt nonensembl_susscrofa_rfam_input_ids.txt > nonsusscrofa_homosapiens_rfam_input_ids.txt

#' As a check, make sure that no sequence ids match between the nonsusscrofa-rfam-hit ID file and the ensembl and rfam-hit ID files
grep -F -f nonsusscrofa_homosapiens_rfam_input_ids.txt susscrofa_rfam_blastn_output_e5_ids.txt

grep -F -f nonsusscrofa_homosapiens_rfam_input_ids.txt susscrofa_ensembl_blastn_output_e5_ids.txt

#' Additionally, the line count from the rfam-hit file and the nonsusscrofa-rfam-hit files should add up to the total number of sequence IDs in nonensembl_susscrofa_rfam_input_ids.txt.
wc -l nonmirbase_susscrofa_ensembl_input_ids.txt
# 1297756 nonmirbase_susscrofa_ensembl_input_ids.txt
wc -l nonensembl_susscrofa_rfam_input_ids.txt
# 1165617 nonensembl_susscrofa_rfam_input_ids.txt
wc -l susscrofa_rfam_blastn_output_e5_ids.txt
# 6552 susscrofa_rfam_blastn_output_e5_ids.txt
wc -l nonsusscrofa_homosapiens_rfam_input_ids.txt
# 1159065 nonsusscrofa_homosapiens_rfam_input_ids.txt

# 1159065 + 6552 = 1165617

#' Then, use the same awk script as previously used to extract the non-sus scrofa rfam hit sequences to reduce the previously used sequences, leaving the homo sapiens rfam input sequences:
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' 174_library_susscrofa_rfam_blast_input_uniq_seq.fa | awk -F"\t" 'BEGIN{while((getline k < "nonsusscrofa_homosapiens_rfam_input_ids.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 174_library_homosapiens_rfam_blast_input_uniq_seq.fa

#' The length of the 174_library_homosapiens_rfam_blast_input_uniq_seq.fa should be twice the length of the nonsusscrofa ID file
wc -l 174_library_homosapiens_rfam_blast_input_uniq_seq.fa
# 2318130 174_library_homosapiens_rfam_blast_input_uniq_seq.fa

wc -l nonsusscrofa_homosapiens_rfam_input_ids.txt
# 1159065 nonsusscrofa_homosapiens_rfam_input_ids.txt

wc -l 174_library_uniq_exp.fa
# 2890826 174_library_uniq_exp.fa

wc -l 174_library_susscrofa_rfam_blast_input_uniq_seq.fa
# 2331234 174_library_susscrofa_rfam_blast_input_uniq_seq.fa

#' Some quick math to summarize the results and check that file lengths are correct:
#' 
#' The length of the original input .fa file minus the ensembl input .fa file should be equivalent to 2*(mirbase hits)
# 2890826 - 2595512 = 295314 (this should be equivalent to 2*number of miRBase hits)
# 295314/2 = 147657

#' The length of the ensembl input .fa file minus the sus scrofa rfam input .fa file should be equivalent to 2*(ensembl hits)
# 2595512 - 2331234 = 264278 (this should be equivalent to 2*number of ensembl hits)
# 264278/2 = 132139

wc -l susscrofa_ensembl_blastn_output_e5_ids.txt
# 132139 susscrofa_ensembl_blastn_output_e5_ids.txt

#' The length of the sus scrofa rfam input .fa file minus the homo sapiens rfam input .fa file should be equivalent to 2*(sus scrofa rfam hits)
# 2331234 - 2318130 = 13104 (this should be equivalent to 2*number of sus scrofa rfam hits)
# 13104/2 = 6552

# Check that both of these processes worked correctly: Subtract the number of mirbase hits, ensembl hits, and sus scrofa rfam hits from the total number of unique sequences input originally
# 2890826 - 295314 - 264278 - 13104 = 2318130 (this should be equivalent to the input .fa file for homo sapiens rfam blast)

#' Also, check that no sequences match between the Sus scrofa Rfam blast hits and the fasta file of input sequences
grep -F -f susscrofa_rfam_blastn_output_e5_ids.txt 174_library_homosapiens_rfam_blast_input_uniq_seq.fa | wc -l
# 0

#' So, the sequences in the 174_library_homosapiens_rfam_blast_input_uniq_seq.fa will be input into a BLAST query against the homo sapiens Rfam database sequences.
#' 
