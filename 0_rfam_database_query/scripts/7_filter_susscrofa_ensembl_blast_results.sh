#' **Script:** `7_filter_susscrofa_ensembl_blast_results.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  `5/31/16 UPDATED 7/11/16`
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Input File(s):** 
#' 
#' 1. `174_library_ensembl_blast_input_uniq_seq.fa`
#' 2. `1_susscrofa_ensembl_blastn_output_e5.txt`
#' 3. `nonmirbase_susscrofa_ensembl_input_ids.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Output File(s):** 
#' 
#' 1. `susscrofa_ensembl_blastn_output_e5_ids.txt`
#' 2. `nonensembl_susscrofa_rfam_input_ids.txt`
#' 3. `174_library_susscrofa_rfam_blast_input_uniq_seq.fa`
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
#' The databases used for this query were created from the Sus scrofa miRBase (release 21), the Ensembl Sus Scrofa (10.2) ncRNA database, and the Sus scrofa, Homo sapiens, and Mus musculus Rfam database version 11.0.
#' See script /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/4_download_ensembl_rfam_build_blastdb.sh for code building those databases.
#' First, the sequences will be blasted against the Sus scrofa miRBase, the Sus scrofa Ensembl noncoding sequences, then move sequentially to the Sus scrofa, Homo sapiens, and Mus musculus Rfam sequences. 
#' This script is used to filter the output of the ensembl blast to filter out the sequences that did have ensembl hits. 
#' 
#' 
#' 
#' ## Install libraries

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#' The index of the sequence ids from the unique sequences input into the ensembl analysis already exists, as the file nonmirbase_susscrofa_ensembl_input_ids.txt:
#' 
#' Subset the ensembl blast result sequence IDs from the blast output
cat 1_susscrofa_ensembl_blastn_output_e5.txt |cut -f1 | sort | uniq > susscrofa_ensembl_blastn_output_e5_ids.txt

#' Subset the non-ensembl blast results from the 174_library_ensembl_blast_input_uniq_seq.fa (will be input for Sus scrofa Rfam blast)
#'
#' grep -F option: Fixed strings: interpret the pattern as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
#' grep -v option: invert match
#' grep -f option: file
#' 
#' First, extract the IDs that were NOT sus scrofa ensembl blast hits by doing the inverse match of the sus scrofa ensembl blast hits and the ensembl blast input sequence ids:
grep -F -v -f susscrofa_ensembl_blastn_output_e5_ids.txt nonmirbase_susscrofa_ensembl_input_ids.txt > nonensembl_susscrofa_rfam_input_ids.txt

#' As a check, make sure that no sus scrofa rfam input sequence ids match the ensembl-hit ID files
grep -F -f nonensembl_susscrofa_rfam_input_ids.txt susscrofa_ensembl_blastn_output_e5_ids.txt

#' Additionally, the line count from the ensembl-hit file and the susscrofa rfam input files should add up to the total number of sequence IDs in nonmirbase_susscrofa_ensembl_input_ids.txt.
wc -l nonensembl_susscrofa_rfam_input_ids.txt
# 1165617 nonensembl_susscrofa_rfam_input_ids.txt
wc -l susscrofa_ensembl_blastn_output_e5_ids.txt
# 132139 susscrofa_ensembl_blastn_output_e5_ids.txt
wc -l nonmirbase_susscrofa_ensembl_input_ids.txt
# 1297756 nonmirbase_susscrofa_ensembl_input_ids.txt

# 1165617 + 132139 = 1297756

#' Then, use the same awk script as previously used to extract the non-ensembl hit sequences to reduce the previously used sequences, leaving the sus scrofa rfam input sequences:
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' 174_library_ensembl_blast_input_uniq_seq.fa | awk -F"\t" 'BEGIN{while((getline k < "nonensembl_susscrofa_rfam_input_ids.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 174_library_susscrofa_rfam_blast_input_uniq_seq.fa

#' The length of the 174_library_susscrofa_rfam_blast_input_uniq_seq.fa should be twice the length of the nonensembl ID file
wc -l 174_library_susscrofa_rfam_blast_input_uniq_seq.fa
# 2331234 174_library_susscrofa_rfam_blast_input_uniq_seq.fa

wc -l nonensembl_susscrofa_rfam_input_ids.txt
# 1165617 nonensembl_susscrofa_rfam_input_ids.txt

# 1165617 * 2 = 2331234

wc -l 174_library_uniq_exp.fa
# 2890826 174_library_uniq_exp.fa

wc -l 174_library_ensembl_blast_input_uniq_seq.fa 
# 2595512 174_library_ensembl_blast_input_uniq_seq.fa

#' Some quick math to summarize the results and check that file lengths are correct:
#' 
#' The length of the original input fa file minus the ensembl input fa file should be equivalent to 2*(mirbase hits)
# 2890826 - 2595512 = 295314 (this should be equivalent to 2*number of miRBase hits)
# 295314/2 = 147657

#' The length of the ensembl input fa file minus the sus scrofa rfam input fa file should be equivalent to 2*(ensembl hits)
# 2595512 - 2331234 = 264278 (this should be equivalent to 2*number of ensembl hits)
# 264278/2 = 132139

# Check that both of these processes worked correctly: Subtract the number of mirbase hits and ensembl hits from the total number of unique sequences input originally
# 2890826 - 295314 - 264278 = 2331234 (this should be equivalent to the input .fa file for sus scrofa rfam blast)

wc -l susscrofa_ensembl_blastn_output_e5_ids.txt
# 132139 susscrofa_ensembl_blastn_output_e5_ids.txt

#' Also, check that no sequences match between the ensembl blast hits and the fasta file of sus scrofa Rfam input sequences
grep -F -f susscrofa_ensembl_blastn_output_e5_ids.txt 174_library_susscrofa_rfam_blast_input_uniq_seq.fa | wc -l
# 0


#' So, the sequences in the 174_library_susscrofa_rfam_blast_input_uniq_seq.fa will be input into a BLAST query against the sus scrofa Rfam database sequences.