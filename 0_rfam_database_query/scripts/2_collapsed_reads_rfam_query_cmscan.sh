# Script: 1_collapse_processed_reads_run_cmscan.sh

# Directory of Code: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts

# Date:  2/01/2016

# Input File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Input File(s): 174_library_collapsed_reads.fa
 
# Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Output File(s): 174_collapsed_libraries_cmscan_tblout.txt
#                 174_collapsed_libraries_cmscan_output.txt

#  1. [Objectives]
#      The objective of this analysis is to identify the different species of small RNA present in the small RNA seq dataset.
#      This is accomplished by querying the small RNA sequencing reads against the Rfam database
#      using the cmscan command of the Infernal software. To decrease computation time, the reads processed by the miRDeep2 module 
#      will be collapsed into unique sequences with total read counts using the miRDeep2 mapper module -m option.
#      Additionally, the "toponly" option of the cmscan command will be used. 

#  2. [Modules Used]
#     infernal/1.1rc1

# ===============================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=5:intel14,walltime=3:00:00,mem=5Gb
#PBS -N 5cpu_174_collapsed_libraries_rfam_cmscan
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/
#PBS -m abe
#PBS -M perrykai@msu.edu

module load Infernal/1.1rc1

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Query the small RNA sequencing data using cmscan command of infernal, using 5 cpus and the --toponly option. 
cmscan --tblout 174_collapsed_libraries_cmscan_tblout.txt --cpu 5 --toponly -o 174_collapsed_libraries_cmscan_output.txt Rfam.cm 174_library_collapsed_reads.fa

qstat -f ${PBS_JOBID}
