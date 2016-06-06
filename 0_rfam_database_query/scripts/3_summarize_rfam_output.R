#' **Script:** ``
#' 
#' **Directory of Code:**  ``
#' 
#' **Date:**  
#' 
#' **Input File Directory:**  ``
#' 
#' **Input File(s):** ``
#' 
#' **Output File Directory:** ``
#' 
#' **Output File(s):** ``
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
#' ## Install libraries
#' ## Load data
#' ## Analysis
#' ## Visualize
#' ## Save data

setwd("/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/0_rfam_database_query/scripts/")

library(plyr)

system.time(
rfamtblout<-readLines("../infernal_rfam_query_output/174_split_collapsed_reads_51_tblout.txt")
)

length(rfamtblout)

#' Eliminate row the last 9 lines (summary info, not useful)
rfamtblout<-rfamtblout[-c(1,2,((length(rfamtblout)-9):length(rfamtblout)))]

length(rfamtblout)
str(rfamtblout)
head(rfamtblout)
tail(rfamtblout)

system.time(
rfamsummary<-ldply(rfamtblout, function(x) (read.table(textConnection(x), colClasses = c(rep("character", 5), rep("numeric", 4), rep("character", 2), rep("numeric", 5), "character"))[1:17]))
)

colnames(rfamsummary) <- c("target_name", "accession_target", "query_name", "accession_query", "mdl", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc", "pass", "gc", "bias", "score", "E_value", "inc")

head(rfamsummary)
tail(rfamsummary)
str(rfamsummary)

save(rfamsummary, file = "../summary_infernal_rfam/test_rfam_summary_table.Rdata")

#126789