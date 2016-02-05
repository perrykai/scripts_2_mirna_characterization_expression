# Script: 2_collapsed_reads_rfam_query_cmscan.sh

# Directory of Code: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts

# Date:  2/05/2016

# Input File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Input File(s): 174_split_collapsed_reads_* (54 files in total)
 
# Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/infernal_rfam_query_output

# Output File(s): 174_split_collapsed_reads_*_tblout.txt
#                 174_split_collapsed_reads_*_output.txt

#  1. [Objectives]
#      The objective of this analysis is to identify the different species of small RNA present in the small RNA seq dataset.
#      This is accomplished by querying the small RNA sequencing reads against the Rfam database
#      using the cmscan command of the Infernal software. To decrease computation time, the reads processed by the miRDeep2 module 
#      will be collapsed into unique sequences with total read counts using the miRDeep2 mapper module -m option.
#      Additionally, the "toponly" option of the cmscan command will be used. 

#  2. [Modules Used]
#     infernal/1.1rc1

# ===============================

f1=(`ls /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/split_174_library_collapsed_files/`)

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/

for (( i = 0 ; i < ${#f1[@]} ; i++ )) do
	echo '#!/bin/sh  -login' > base.sh
	echo '#PBS -l nodes=2:ppn=8:intel14,walltime=6:00:00,mem=5Gb' >> base.sh
	echo '#PBS -N '${f1[$i]}'_rfamscan' >> base.sh
	echo '#PBS -j oe' >> base.sh
	echo '#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/' >> base.sh
	echo '#PBS -m a' >> base.sh
	echo '#PBS -M perrykai@msu.edu' >> base.sh

	echo 'module load infernal/1.1rc1' >> base.sh

	# Query the small RNA sequencing data using cmscan command of infernal, using 5 cpus and the --toponly option. 
	echo 'cmscan --tblout /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/infernal_rfam_query_output/'${f1[$i]}'_tblout.txt --cpu 16 --toponly -o /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/infernal_rfam_query_output/'${f1[$i]}'_output.txt /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/Rfam.cm /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/split_174_library_collapsed_files/'${f1[$i]} >> base.sh

	echo 'qstat -f ${PBS_JOBID}' >> base.sh

	qsub base.sh

done