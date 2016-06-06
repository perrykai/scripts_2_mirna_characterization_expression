#' **Script:** `11_filter_hsapiens_rfam_blast_results.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts`
#' 
#' **Date:**  6/2/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_homosapiens_rfam_blastn_output_e6.txt`
#' 2. `174_library_homosapiens_rfam_input_uniq_seq.fa`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Output File(s):** 
#' 
#' 1. `homosapiens_rfam_blastn_output_e6_ids.txt`
#' 2. `nonhomosapiens_rfam_musmusculus_input_ids.txt`
#' 3. `174_library_musmusculus_rfam_input_uniq_seq.fa`
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
#' This script is used to filter the output of the homo sapiens Rfam blast results to filter out the sequences that did have homo sapiens rfam hits and create the input for the mus musculus BLAST query
#' 
#' 
#' 
#' ## Install libraries

#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#' The index of the sequence ids from the unique sequences input into the Rfam analysis already exists, as the file nonensembl_rfam_susscrofa_input_ids.txt:
#' 
#' Subset the ensembl blast result sequence IDs from the blast output
cat 3_homosapiens_rfam_blastn_output_e6.txt |cut -f1 | sort | uniq > homosapiens_rfam_blastn_output_e6_ids.txt

#' Subset the non-susscrofa rfam blast results from the 174_library_nonensembl_uniq_seq.fa (will be input for Rfam homo sapiens blast)
#'
#' grep -F option: Fixed strings: interpret the pattern as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
#' grep -v option: invert match
#' grep -f option: file
#' 
#' First, extract the IDs that were NOT homo sapiens Rfam blast hits by doing the inverse match of the homo sapiens Rfam blast hits and the input sequence ids:
grep -F -v -f homosapiens_rfam_blastn_output_e6_ids.txt nonsusscrofa_rfam_homosapiens_input_ids.txt > nonhomosapiens_rfam_musmusculus_input_ids.txt

#' As a check, make sure that no sequence ids match between the nonhomosapiens-rfam-hit ID file and the previous ensembl and rfam-hit ID files
grep -F -f nonhomosapiens_rfam_musmusculus_input_ids.txt homosapiens_rfam_blastn_output_e6_ids.txt

grep -F -f nonhomosapiens_rfam_musmusculus_input_ids.txt susscrofa_rfam_blastn_output_e6_ids.txt

grep -F -f nonhomosapiens_rfam_musmusculus_input_ids.txt susscrofa_ensembl_blastn_output_e6_ids.txt


#' Additionally, the line count from the susscrofa rfam-hit file and the nonsusscrofa-rfam-hit files should add up to the total number of sequence IDs in nonsusscrofa_rfam_homosapiens_input_ids.txt.
wc -l nonsusscrofa_rfam_homosapiens_input_ids.txt
# 1311674 nonsusscrofa_rfam_homosapiens_input_ids.txt
wc -l homosapiens_rfam_blastn_output_e6_ids.txt
# 13328 homosapiens_rfam_blastn_output_e6_ids.txt
wc -l nonhomosapiens_rfam_musmusculus_input_ids.txt
# 1298346 nonhomosapiens_rfam_musmusculus_input_ids.txt

# 1298346 + 13328 = 1311674

#' Then, use the same awk script as previously used to extract the Rfam homo sapiens sequences to reduce the previously used sequences, leaving the non-homosapiens rfam sequences:
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' 174_library_homosapiens_rfam_input_uniq_seq.fa | awk -F"\t" 'BEGIN{while((getline k < "nonhomosapiens_rfam_musmusculus_input_ids.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 174_library_musmusculus_rfam_input_uniq_seq.fa

#' Also, check that no sequences match between the homosapiens Rfam blast hits and the fasta file of mus musculus rfam input sequences
grep -F -f homosapiens_rfam_blastn_output_e6_ids.txt 174_library_musmusculus_rfam_input_uniq_seq.fa | wc -l

#' The length of the 174_library_musmusculus_rfam_input_uniq_seq.fa should be twice the length of the nonhomosapiens ID file
wc -l 174_library_musmusculus_rfam_input_uniq_seq.fa
# 2596692 174_library_musmusculus_rfam_input_uniq_seq.fa

wc -l nonhomosapiens_rfam_musmusculus_input_ids.txt
# 1298346 nonhomosapiens_rfam_musmusculus_input_ids.txt

# 2596692/2 = 1298346


#' As another check, subtract the number of lines in each of the fasta files and make sure that they line up with the number of IDs extracted from each BLAST query.
wc -l 174_library_uniq_exp.fa
# 2890826 174_library_uniq_exp.fa

wc -l 174_library_nonensembl_uniq_seq.fa
# 2635056 174_library_nonensembl_uniq_seq.fa

wc -l 174_library_homosapiens_rfam_input_uniq_seq.fa
# 2623348 174_library_homosapiens_rfam_input_uniq_seq.fa

wc -l 174_library_musmusculus_rfam_input_uniq_seq.fa
# 2596692 174_library_musmusculus_rfam_input_uniq_seq.fa

# 2890826 - 2635056 = 255770 (this should be equivalent to 2*number of sus scrofa ensembl hits)
# 255770/2 = 127885

# 2635056 - 2623348 = 11708 (this should be equivalent to 2*number of sscrofa rfam hits)
# 11708/2 = 5854

# 2623348 - 2596692 = 26656 (this should be equivalent to 2*number of homosapiens rfam hits)
# 26656/2 = 13328



# Check that both of these processes worked correctly: Subtract the number of ensembl hits and sus scrofa rfam hits from the total number of unique sequences input originally
# 2890826 - 255770 - 11708 = 2623348 (this should be equivalent to the input .fa file for homo sapiens rfam blast)

# Check that both of these processes worked correctly: Subtract the number of ensembl hits, sus scrofa, and homo sapiens rfam hits from the total number of unique sequences input originally
# 2890826 - 255770 - 11708 - 26656 = 2596692 (this should be equivalent to the input .fa file for mus musculus rfam blast)

#' So, the sequences in the 174_library_musmusculus_rfam_input_uniq_seq.fa will be input into a BLAST query against the mus musculus Rfam database sequences.
#' 


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
