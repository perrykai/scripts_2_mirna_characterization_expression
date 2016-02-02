#' **File:** `1_extract_mirna_info_from_core_output.R`
#' 
#' **Directory Code:**  `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_extract_mirna_info_from_core_output/scripts/`
#' 
#' **Date:**  1/26/16
#' 
#' **Description:**  This code extracts the miRDeep2 core module output from the `result_21_01_2016_t_20_01_10.csv` file
#'                 This .csv file contains the miRDeep2 score distribution for the novel miRNA prediction step, 
#'                 the novel miRNAs predicted by miRDeep2, the mature miRBase miRNAs detected by miRDeep2, and the 
#'                 miRBase miRNAs not detected by miRDeep2. The first three of these items are extracted from this .csv file using this script.
#'                 The objective here is to characterize the known and novel miRNAs present in this dataset, isolate the sequences meeting miRDeep2 score of 10 or greater,
#'                 to isolate the sequences at that score cutoff having homologous seed sequence with a human miRBase miRNA,
#'                 and to estimate the false discovery rate of the miRDeep2 prediction step at miRDeep2 score 10. 
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output`
#' 
#' **Input File(s):**  `result_21_01_2016_t_20_01_10.csv`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/1_mirna_info_extracted_core_output/`
#' 
#' **Output File(s):** 
#' 
#' -  `1_extracted_mirdeep2_core_output_stats.csv`
#' -  `2_extracted_mirdeep2_core_predicted_novel_mirna.csv`
#' -  `3_extracted_mirdeep2_core_mature_mirna_detected.csv`
#' 

#' To isolate the first section of the csv containing the miRDeep2 distribution scores:
sts<-read.table("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/result_21_01_2016_t_20_01_10.csv", sep = "\t", nrows=21, header = TRUE)

head(kable(sts))

write.table(sts, file="../1_extracted_mirdeep2_core_output_stats.csv", sep = "\t", col.names = TRUE)

#' To isolate the second section of the csv containing the novel predicted miRNAs:
cv<-read.table("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/result_21_01_2016_t_20_01_10.csv", sep = "\t", skip=26, nrow=1199,  header = TRUE, fill=TRUE)

head(kable(cv))

write.table(cv, file="../2_extracted_mirdeep2_core_predicted_novel_mirna.csv", sep = "\t", col.names=TRUE)

#' To isolate the third section of the csv containing the miRBase miRNAs detected by miRDeep2:
md<-read.table("../../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/result_21_01_2016_t_20_01_10.csv", sep = "\t", skip = 1230, nrow = 306, header = TRUE, fill = TRUE)

head(kable(md))

write.table(md, file = "../3_extracted_mirdeep2_core_mature_mirna_detected.csv", sep = "\t", col.names=TRUE)

#' Now I can open the `novel_mirna_predicted.csv` file and filter by various thresholds
#' 
#' ### 1. miRDeep2 score of 10 or more

novelmirna<-read.table("../2_extracted_mirdeep2_core_predicted_novel_mirna.csv")

colnames(novelmirna)

dim(novelmirna)	

head(kable(novelmirna))

score10novelmirna<-novelmirna[novelmirna$miRDeep2.score >= 10, ]

dim(score10novelmirna)

nrow(score10novelmirna)

#' So, there are 352 predicted miRNA precursors with a miRDeep2 score > or = 10
#' 
#' Estimated false positives is 32 +/- 6 (obtained from `1_extracted_mirdeep2_core_output_stats.csv`)

32 - 6

32 + 6

26/352

38/352


#' ### 2. Significant Randfold p-value
#' 
#' Now, subset this again into those that had a significant Randfold p value, indicating ability of secondary structure formation
head(score10novelmirna$significant.randfold.p.value)

sum(score10novelmirna$significant.randfold.p.value=="yes")
randfoldsigpval<-score10novelmirna[score10novelmirna$significant.randfold.p.value == "yes", ]
dim(randfoldsigpval)

#' Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences
rfam<-randfoldsigpval$rfam.alert
sum(rfam=="-")
sum(rfam !="-")


#' ### 3. Do the putative novel sequences have a homologous human miRNA seed sequence?
sum(randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-")

#' This indicates that 100 of the sequences have a homologous human seed sequence
homologseed<-randfoldsigpval[randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-",]

#' Subset of the 100 sequences
dim(homologseed)

write.table(homologseed, "../4_extracted_predicted_novel_sigranfold_homologseed.txt", sep = "\t", col.names=TRUE)
