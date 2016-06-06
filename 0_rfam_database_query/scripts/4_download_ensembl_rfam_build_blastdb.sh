#' **Script:** `4_download_ensembl_rfam_build_blastdb.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query`
#' 
#' **Date:**  `5/12/16 MODIFIED 5/29/16`
#' 
#' **Input File Directory:**  
#' 
#' 1. `ftp://ftp.ensembl.org/pub/release-84/fasta/sus_scrofa/ncrna/`
#' 2. `ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/`
#' 
#' **Input File(s):** 
#' 
#' 1. `Sus_scrofa.Sscrofa10.2.ncrna.fa.gz`
#' 2. `Rfam.fasta.gz`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences`
#' 
#' **Output File(s):** 
#'
#' 1. `Sus_scrofa.Sscrofa10.2.ncrna.fa.gz`
#' 2. `Rfam.fasta.gz`
#' 3. `Susscrofa_Rfam_seq.fa`
#' 4. `Homosapiens_Rfam_seq.fa`
#' 5. `Musmusculus_Rfam_seq.fa`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 2. [Analysis](#analysis)
#' 
#' ## Objectives

#' The objectives of this script are to download the Sus scrofa ncRNA sequences from Ensembl 
#' and the Rfam database (version 11.0) for use in characterizing the smallRNA species present in the 
#' small RNA seq data obtained from the 174 F2 MSUPRP samples. Additionally, these two databases will be built into
#' BLAST databases for use in the BLAST query. 

#' 
#' ## Install libraries
module load BLAST+/2.2.30


#' ## Analysis
#' 
#' Download the Sus_scrofa.Sscrofa10.2.ncrna.fa.gz file from the ensembl ftp site:
cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences

wget ftp://ftp.ensembl.org/pub/release-84/fasta/sus_scrofa/ncrna/Sus_scrofa.Sscrofa10.2.ncrna.fa.gz
# --2016-05-12 11:05:10--  ftp://ftp.ensembl.org/pub/release-84/fasta/sus_scrofa/ncrna/Sus_scrofa.Sscrofa10.2.ncrna.fa.gz
#            => “Sus_scrofa.Sscrofa10.2.ncrna.fa.gz”
# Resolving ftp.ensembl.org... 193.62.203.85
# Connecting to ftp.ensembl.org|193.62.203.85|:21... connected.
# Logging in as anonymous ... Logged in!
# ==> SYST ... done.    ==> PWD ... done.
# ==> TYPE I ... done.  ==> CWD (1) /pub/release-84/fasta/sus_scrofa/ncrna ... done.
# ==> SIZE Sus_scrofa.Sscrofa10.2.ncrna.fa.gz ... 255606
# ==> PASV ... done.    ==> RETR Sus_scrofa.Sscrofa10.2.ncrna.fa.gz ... done.
# Length: 255606 (250K) (unauthoritative)

# 100%[===========================================================================================================================================>] 255,606      578K/s   in 0.4s

# 2016-05-12 11:05:12 (578 KB/s) - “Sus_scrofa.Sscrofa10.2.ncrna.fa.gz” saved [255606]

#' 
#' Download the Rfam ncRNA sequence database version 11.0
cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences

wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.fasta.gz
# Rfam/11.0/Rfam.fasta.gz
# --2016-05-13 14:08:48--  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.fasta.gz
#            => “Rfam.fasta.gz”
# Resolving ftp.ebi.ac.uk... 193.62.194.182
# Connecting to ftp.ebi.ac.uk|193.62.194.182|:21... connected.
# Logging in as anonymous ... Logged in!
# ==> SYST ... done.    ==> PWD ... done.
# ==> TYPE I ... done.  ==> CWD (1) /pub/databases/Rfam/11.0 ... done.
# ==> SIZE Rfam.fasta.gz ... 19463617
# ==> PASV ... done.    ==> RETR Rfam.fasta.gz ... done.
# Length: 19463617 (19M) (unauthoritative)

# 100%[=============================================>] 19,463,617  4.77M/s   in 5.7s

# 2016-05-13 14:08:56 (3.28 MB/s) - “Rfam.fasta.gz” saved [19463617]

#' 
#' Unzip the two database files and build the BLAST database
gunzip Sus_scrofa.Sscrofa10.2.ncrna.fa.gz
gunzip Rfam.fasta.gz

makeblastdb -in Sus_scrofa.Sscrofa10.2.ncrna.fa -input_type fasta -dbtype nucl -title "Sus Scrofa 10.2 ncRNA"

# Building a new DB, current time: 05/13/2016 14:52:02
# New DB name:   Sus_scrofa.Sscrofa10.2.ncrna.fa
# New DB title:  Sus Scrofa 10.2 ncRNA
# Sequence type: Nucleotide
# Keep Linkouts: T
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 3215 sequences in 0.13694 seconds.

#' Updated 5/29/16: 
#'

#' Prior to BLAST analysis, first filter both the unique sequences and the Rfam databases to reduce computation time.
cd /mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query

#' First, from the unique sequences, remove all sequences only expressed one time (no biologically influential inference can be drawn from a BLAST hit from one molecule of sequence output)
sed '/_x1$/,+1 d' 174_library_collapsed_reads.fa > 174_library_uniq_exp.fa

#' As a check, grep the same expression to be sure it doesn't appear in the output: 
grep '_x1$' 174_library_uniq_exp.fa


#' Now, extract the pig, human, and mouse Rfam database entries to create smaller Rfam databases to BLAST against:
cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/

#' First, extract the Sus scrofa, Homo sapiens, and Mus musculus sequences from Rfam:
grep '9823:Sus scrofa (pig)$' Rfam.fasta | cut -c2- > Susscrofa_RfamIDs.txt

grep '9606:Homo sapiens (human)$' Rfam.fasta | cut -c2- > Homosapiens_RfamIDs.txt

grep '10090:Mus musculus (house mouse)$' Rfam.fasta | cut -c2- > Musmusculus_RfamIDs.txt

#' Then, using the awk script found on BioStars: https://www.biostars.org/p/127141/
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' Rfam.fasta | awk -F"\t" 'BEGIN{while((getline k < "Susscrofa_RfamIDs.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > Susscrofa_Rfam_seq.fa

head Susscrofa_Rfam_seq.fa

awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' Rfam.fasta | awk -F"\t" 'BEGIN{while((getline k < "Homosapiens_RfamIDs.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > Homosapiens_Rfam_seq.fa

head Homosapiens_Rfam_seq.fa

awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' Rfam.fasta | awk -F"\t" 'BEGIN{while((getline k < "Musmusculus_RfamIDs.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > Musmusculus_Rfam_seq.fa

head Musmusculus_Rfam_seq.fa

#' Then, make BLAST+ databases out of the Rfam species-specific IDs:
#' 
#' First, Sus scrofa:
makeblastdb -in Susscrofa_Rfam_seq.fa -input_type fasta -dbtype nucl -title "Sus Scrofa Rfam Database"

# Building a new DB, current time: 05/29/2016 15:54:16
# New DB name:   Susscrofa_Rfam_seq.fa
# New DB title:  Sus Scrofa Rfam Database
# Sequence type: Nucleotide
# Keep Linkouts: T
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 592 sequences in 0.021997 seconds.

#' 
#' Second, Homo sapiens:
makeblastdb -in Homosapiens_Rfam_seq.fa -input_type fasta -dbtype nucl -title "Homo Sapiens Rfam Database"

# Building a new DB, current time: 05/29/2016 15:57:46
# New DB name:   Homosapiens_Rfam_seq.fa
# New DB title:  Homo Sapiens Rfam Database
# Sequence type: Nucleotide
# Keep Linkouts: T
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 2838 sequences in 0.094084 seconds.

#' 
#' Finally, Mus musculus:
makeblastdb -in Musmusculus_Rfam_seq.fa -input_type fasta -dbtype nucl -title "Mus Musculus Rfam Database"

# Building a new DB, current time: 05/29/2016 15:57:51
# New DB name:   Musmusculus_Rfam_seq.fa
# New DB title:  Mus Musculus Rfam Database
# Sequence type: Nucleotide
# Keep Linkouts: T
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 5935 sequences in 0.181584 seconds.

#' Test that the format of the Rfam databases built with the awk script will work with the makeblastdb and isn't missing anything in its blast against the small RNA seq data:
#' First, take a few sequences out of the Rfam Sus scrofa database that had blast hits:

# seq_161536921_x29165    RF00001;5S_rRNA;CT009542.22/163378-163269       100.00  20      0       0       1       20      72      91      6e-07   40.1
# seq_164820205_x24792    RF00708;mir-450;AEMK01130021.1/1834-1924        100.00  21      0       0       1       21      18      38      1e-07   42.1
# seq_167301252_x22072    RF00737;mir-322;CU468475.7/117556-117473        100.00  22      0       0       1       22      11      32      4e-08   44.1
# seq_171050786_x17914    RF00708;mir-450;AEMK01130021.1/1834-1924        100.00  20      0       0       1       20      18      37      6e-07   40.1
# seq_176615697_x13351    RF00737;mir-322;CU468475.7/117556-117473        100.00  20      0       0       1       20      11      30      6e-07   40.1
# seq_178776422_x11814    RF00737;mir-322;CU468475.7/117556-117473        100.00  22      0       0       1       22      11      32      4e-08   44.1

makeblastdb -in Susscrofa_Rfam_seq.fa -input_type fasta -dbtype nucl -title "Sus Scrofa Rfam Database"