library(polyester)

library(tidyverse)



filename <- "abundance.tsv"

import_fpkms <- function(samples, base_path, filename) {

	count_df_list = list()
	fpk_df_list = list()

	for (sample in samples) {

		path = paste0(base_path, "/", sample, "/", filename)
		temp_df <- read.table(path, 
			header = TRUE, 
			sep = "\t",
			stringsAsFactors = FALSE) %>%
			rename(transcript_id = target_id)

		count_df <- temp_df %>%
			select(transcript_id, est_counts) %>%
			rename(!!sample := est_counts)

		count_df_list <- c(count_df_list, list(count_df))

		fpk_df <- temp_df %>%
			transmute(transcript_id = transcript_id, !!sample := est_counts/(length/1000))

		fpk_df_list <- c(fpk_df_list, list(fpk_df))

	}

	count_df <- count_df_list %>%
		purrr::reduce(inner_join, by = "transcript_id")

	fpk_df <- fpk_df_list %>%
		purrr::reduce(inner_join, by = "transcript_id")

	fragment_totals <- colSums(count_df[,2:ncol(count_df)])

	print(fragment_totals)

	per_million_factors <- fragment_totals/1000000

	fpkm_df <- fpk_df

	fpkm_df[,2:ncol(fpkm_df)] <- fpkm_df[,2:ncol(fpkm_df)]/per_million_factors


	return(list(fpkm_df, fragment_totals))

}

print("Function defined . . . ")

###experiment 1 

###condition a file paths

exp1_transcript_fasta <- "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014"

exp1_base <- "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/kallisto/"

exp1a_samples <- c("SRR1158578","SRR1158580","SRR1158582")

exp1b_samples <- c("SRR1158546","SRR1158548","SRR1158550")

###exp1a mean tpms

exp1a_res <- import_fpkms(exp1a_samples, exp1_base, filename)
exp1b_res <- import_fpkms(exp1b_samples, exp1_base, filename)

exp1a_fpkm <- exp1a_res[[1]]
exp1b_fpkm <- exp1b_res[[1]]

exp1a_depth <- exp1a_res[[2]]
exp1b_depth <- exp1b_res[[2]]

print(exp1a_depth)
print(exp1b_depth)
print(mean(c(exp1a_depth, exp1b_depth))/2)

exp1_fpkm <- exp1a_fpkm %>% 
	inner_join(exp1b_fpkm, by = "transcript_id")

exp1_fpkm_mat <- as.matrix(exp1_fpkm[,2:ncol(exp1_fpkm)])
rownames(exp1_fpkm_mat) <- exp1_fpkm$transcript_id

simulate_experiment_empirical(fpkmMat = exp1_fpkm_mat, 
							  grouplabels=c("a","a","a","b","b","b"), 
							  fasta = exp1_transcript_fasta,
   							  mean_rps=mean(c(exp1a_depth, exp1b_depth))*2, 
   							  seed=1247,
   							  outdir="/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/simulations/",
   							  bias = "rnaf",
   							  strand_specific = TRUE,
   							  gzip = TRUE,
   							  paired = TRUE)

