#!/bin/sh  -login
#PBS -l nodes=1:ppn=1,walltime=00:05:00,mem=1Gb
#PBS -N 5_copy_filtered_candidate_novel_precursor_ids_e5
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/5_novel_miRNA_characterization
#PBS -m a
#PBS -M perrykai@msu.edu

#' **Script:** `5_copy_filtered_candidate_novel_precursor_ids_e5.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts`
#' 
#' **Date:**  11/2/16
#' 
#' **Input File Directory:**  
#'
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/pdfs_21_01_2016_t_20_01_10
#' 
#' **Input File(s):** 
#' 
#' 1. `6_filtered_candidate_novel_miRNA_precursor_ids_e5.txt`
#' 2. One predicted secondary structure pdf file for each of the 27 unique candidate novel precursors. 
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output`
#' 
#' **Output File(s):** 
#' 
#' 1. The predicted secondary structure pdf files for each of the 27 unique candidate novel precursors will be copied to the output directory from 
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Analysis](#analysis)
#' 3. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to copy the predicted secondary structure pdf files for each of the 27 unique candidate novel precursors in order to build a figure containing graphics of the secondary structures of potential novel miRNAs.
#' 
#' ## Analysis
#' 
#' 
#' ## Save data

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output

cut -d '"' -f2 6_filtered_candidate_novel_miRNA_precursor_ids_e5.txt > 7_cut_filtered_ids.txt
sed -e 's/$/.pdf/' -i 7_cut_filtered_ids.txt

f1=`cat /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/7_cut_filtered_ids.txt`
cd /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/pdfs_21_01_2016_t_20_01_10

for file in `cat ../../../2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/8_blastn_candidate_novel_miRNA_output/7_cut_filtered_ids.txt`;
do cp "$file" ../../../2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/9_candidate_novel_miRNA_precursor_pdf_e5/;
done

qstat -f ${PBS_JOBID}