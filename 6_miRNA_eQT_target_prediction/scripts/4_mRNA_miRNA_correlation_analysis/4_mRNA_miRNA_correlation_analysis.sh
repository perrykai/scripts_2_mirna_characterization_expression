
#=====================================================================
# This script runs: 4_mRNA_miRNA_correlation_analysis.R
# Submited on: Mon Jul 3 13:51:20 EDT 2017
#=====================================================================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=11,walltime=00:15:00,mem=12G
#PBS -N 4_mRNA_miRNA_correlation_analysis
#PBS -j oe
#PBS -m abe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/scripts/4_mRNA_miRNA_correlation_analysis

cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction/scripts/4_mRNA_miRNA_correlation_analysis

R -e 'library("knitr");knitr::spin ("4_mRNA_miRNA_correlation_analysis.R")'

qstat -f ${PBS_JOBID}
