# Script: 0_collapse_processed_reads_mirdeep2.sh

# Directory of Code: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts

# Date:  2/01/2016

# Input File Directory: /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/9_mirdeep2_genome_mapper_output

# Input File(s): 174_library_reads_processed.fa
 
# Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Output File(s): 174_library_collapsed_reads.fa

#  1. [Objectives]
#      The objective of this analysis is to identify the different species of small RNA present in the small RNA seq dataset.
#      This is accomplished by querying the small RNA sequencing reads against the Rfam database
#      using the cmscan command of the Infernal software. To decrease computation time, the reads processed by the miRDeep2 module 
#      will be collapsed into unique sequences with total read counts using the miRDeep2 mapper module -m option.
#      Additionally, the "toponly" option of the cmscan command will be used. 
#      This script accomplishes the collapsing of the small RNA seq reads, then splits the file into 54 equal parts for parallel processing.

#  2. [Modules Used]
#     miRDeep2/0.0.5


# ===============================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel14,walltime=00:20:00,mem=1Gb
#PBS -N collapse_processed_reads_for_cmscan
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/
#PBS -m a
#PBS -M perrykai@msu.edu

module load miRDeep2/0.0.5

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query

# Collapse reads processed (output file from miRDeep2 genome mapper step) using -m option of miRDeep2/0.0.5
# -c indicates input is in fasta format
# -m collapses the reads
# -s writes the output file to the designated file
mapper.pl ../../1_preprocess_fastq_files/9_mirdeep2_genome_mapper_output/174_library_reads_processed.fa -c -m -s 174_library_collapsed_reads.fa

# Split the files so each one has 253578 lines, meaning each file has 126789 sequences and can run efficiently through Infernal
split -d -l 253578 174_library_collapsed_reads.fa ./split_174_library_collapsed_files/174_split_collapsed_reads_ 

qstat -f ${PBS_JOBID}