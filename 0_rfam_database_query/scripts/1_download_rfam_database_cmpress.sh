# Script: 1_download_rfam_database_cmpress.sh

# Directory of Code: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts

# Date:  2/01/2016

# Input File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Input File(s): 174_library_collapsed_reads.fa
 
# Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

# Output File(s): 174_collapsed_libraries_cmscan_tblout.txt
#                 174_collapsed_libraries_cmscan_output.txt

#  1. [Objectives]
#      The objective of this project  is to identify the different species of small RNA present in the small RNA seq dataset.
#      This is accomplished by querying the small RNA sequencing reads against the Rfam database
#      using the cmscan command of the Infernal software. To decrease computation time, the reads processed by the miRDeep2 module 
#      will be collapsed into unique sequences with total read counts using the miRDeep2 mapper module -m option.
#      Additionally, the "toponly" option of the cmscan command will be used. 
#	   This script downloads the Rfam.cm file from the Rfam ftp site, unzips it, and uses cmpress from Infernal 
#      to prepare it for use in the cmscan query. 


#  2. [Modules Used]
#     infernal/1.1rc1
#     wget

# ===============================
module load infernal/1.1rc1


cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/

#Download Rfam 12.0 from the ftp site
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz
	# --2016-02-01 16:36:31--  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz
	#            => “Rfam.cm.gz”
	# Resolving ftp.ebi.ac.uk... 193.62.192.4
	# Connecting to ftp.ebi.ac.uk|193.62.192.4|:21... connected.
	# Logging in as anonymous ... Logged in!
	# ==> SYST ... done.    ==> PWD ... done.
	# ==> TYPE I ... done.  ==> CWD (1) /pub/databases/Rfam/12.0 ... done.
	# ==> SIZE Rfam.cm.gz ... 28161375
	# ==> PASV ... done.    ==> RETR Rfam.cm.gz ... done.
	# Length: 28161375 (27M) (unauthoritative)

	# 100%[===================================================================================================================================>] 28,161,375  3.72M/s   in 6.7s

	# 2016-02-01 16:36:40 (4.02 MB/s) - “Rfam.cm.gz” saved [28161375]


#Unzip the file
gunzip Rfam.cm.gz

cmpress Rfam.cm
	# Working...    done.
	# Pressed and indexed 2450 CMs and p7 HMM filters (2450 names and 2450 accessions).
	# Covariance models and p7 filters pressed into binary file:  Rfam.cm.i1m
	# SSI index for binary covariance model file:                 Rfam.cm.i1i
	# Optimized p7 filter profiles (MSV part)  pressed into:      Rfam.cm.i1f
	# Optimized p7 filter profiles (remainder) pressed into:      Rfam.cm.i1p