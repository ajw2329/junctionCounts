library(polyester)

library(tidyverse)



filename <- "abundance.tsv"

import_fpkms <- function(samples, base_path, filename) {
  
  count_df_list = list()
  fpk_df_list = list()
  tpm_df_list = list()
  
  for (sample in samples) {
    
    path = paste0(base_path, "/", sample, "/", filename)
    temp_df <- read.table(path, 
                          header = TRUE, 
                          sep = "\t",
                          stringsAsFactors = FALSE) %>%
      rename(transcript_id = target_id)
    
    sample_tpm <- paste0(sample, "_tpm")
    
    tpm_df <- temp_df %>%
      select(transcript_id, tpm) %>%
      rename(!!sample_tpm := tpm)
    
    tpm_df_list <- c(tpm_df_list, list(tpm_df))
    
    count_df <- temp_df %>%
      select(transcript_id, est_counts) %>%
      rename(!!sample := est_counts)
    
    count_df_list <- c(count_df_list, list(count_df))
    
    fpk_df <- temp_df %>%
      transmute(transcript_id = transcript_id, !!sample := est_counts/(length/1000))
    
    fpk_df_list <- c(fpk_df_list, list(fpk_df))
    
  }
  
  tpm_df <- tpm_df_list %>%
    purrr::reduce(inner_join, by = "transcript_id")
  
  count_df <- count_df_list %>%
    purrr::reduce(inner_join, by = "transcript_id")
  
  fpk_df <- fpk_df_list %>%
    purrr::reduce(inner_join, by = "transcript_id")
  
  fragment_totals <- colSums(count_df[,2:ncol(count_df)])
  
  print(fragment_totals)
  
  per_million_factors <- fragment_totals/1000000
  
  fpkm_df <- fpk_df
  
  fpkm_df[,2:ncol(fpkm_df)] <- fpkm_df[,2:ncol(fpkm_df)]/per_million_factors
  
  
  return(list(fpkm_df, fragment_totals, tpm_df))
  
}

print("Function defined . . . ")

###experiment 1 

###condition a file paths

exp1_transcript_fasta <- "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/gencode.v29lift37.basic.annotation.one_col.fa"

exp1_base <- "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/kallisto/"

exp1a_samples <- c("SRR536342","SRR536344","SRR536346")

exp1b_samples <- c("SRR536348","SRR536350","SRR536352")

###exp1a mean tpms

exp1a_res <- import_fpkms(exp1a_samples, exp1_base, filename)
exp1b_res <- import_fpkms(exp1b_samples, exp1_base, filename)

exp1a_tpm <- exp1a_res[[3]]
exp1b_tpm <- exp1b_res[[3]]

exp1a_fpkm <- exp1a_res[[1]]
exp1b_fpkm <- exp1b_res[[1]]

exp1a_depth <- exp1a_res[[2]]
exp1b_depth <- exp1b_res[[2]]

print(exp1a_depth)
print(exp1b_depth)
print(mean(c(exp1a_depth, exp1b_depth))/2)

exp1_fpkm_tpm <- exp1a_fpkm %>% 
  inner_join(exp1b_fpkm, by = "transcript_id") %>%
  inner_join(exp1a_tpm, by = "transcript_id") %>%
  inner_join(exp1b_tpm, by = "transcript_id")

transcript_id <- exp1_fpkm_tpm$transcript_id

exp1_fpkm_tpm <- cbind(transcript_id,
                       sample_frac(exp1_fpkm_tpm[,c(2,3,4,8,9,10,5,6,7,11,12,13)]))


# randomize order ot transcript names to proliferate changing splicing events

exp1_fpkm_tpm %>%
  write.table(file = "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/simulations/shuffled_tx_ids/fpkm_tpm_table.tsv",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)



exp1_fpkm_mat <- as.matrix(exp1_fpkm_tpm[,2:7])
rownames(exp1_fpkm_mat) <- exp1_fpkm_tpm$transcript_id

simulate_experiment_empirical(fpkmMat = exp1_fpkm_mat, 
							  grouplabels=c("a","a","a","b","b","b"), 
							  fasta = exp1_transcript_fasta,
   							  mean_rps=mean(c(exp1a_depth, exp1b_depth))/2, 
   							  seed=1247,
   							  outdir="/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/simulations/shuffled_tx_ids/fasta/",
   							  bias = "rnaf",
   							  strand_specific = TRUE,
   							  gzip = TRUE,
   							  paired = TRUE)

