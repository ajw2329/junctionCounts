library(polyester)
library(truncnorm)
library(aster)
library(dplyr)
library(Biostrings)


fasta_file <- '/hive/users/anjowall/projects/polyester_exploration/gencode.v26.annotation.50000_tx.fa'
new_fasta_file <- '/hive/users/anjowall/projects/polyester_exploration/gencode.v26.annotation.50000_tx_subset.fa'
kallisto_abundance <- "/hive/users/anjowall/projects/polyester_exploration/h9esc_cyto_jd004_kallisto_out/abundance.tsv"


filename <- "abundance.tsv"

import_fpkms <- function(samples) {

	df_list = list()

	for (sample in samples) {

		path = paste0(exp1_base, sample, "/", filename)
		temp_df <- read.table(path, 
			header = TRUE, 
			sep = "\t",
			stringsAsFactors = FALSE) %>% 
		select(target_id, est_counts) %>%
		mutate(fpk = est_counts/eff_length) %>%
		rename(transcript_id = target_id, !!sample := fpk)
		df_list <- c(df_list, list(temp_df))

	}

	fpk_df <- df_list %>%
		reduce(inner_join, by = "transcript_id")

	fragment_totals <- colSums(fpk_df[,2:ncol(fpk_df)])

	per_million_factors <- fragment_totals/1000000

	fpkm_df <- fpk_df/per_million_factors

	return(list(fpkm_df, fragment_totals))

}


###experiment 1 

###condition a file paths

exp1_base <- "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/kallisto/"

exp1a_samples <- c("SRR536342","SRR536344","SRR536346")

exp1b_samples <- c("SRR536348","SRR536350","SRR536352")

###exp1a mean tpms

exp1a_res <- import_fpkms(exp1a_samples)
exp1b_res <- import_fpkms(exp1b_samples)

exp1a_fpkm <- exp1a_res[[1]]
exp1b_fpkm <- exp1b_res[[1]]

exp1a_depth <- exp1a_res[[2]]
exp1b_depth <- exp1b_res[[2]]

exp1_fpkm <- exp1a_fpkm %>% 
	inner_join(exp1b_fpkm, by = "transcript_id")

exp1_fpkm_mat <- as.matrix(exp1_fpkm[,2:ncol(exp1_fpkm)])
rownames(exp1_fpkm_mat) <- exp1_fpkm$transcript_id

simulate_experiment_empirical(fpkmMat = exp1_fpkm_mat, 
							  grouplabels=pData(bg)$group, 
							  gtf=gtf,
    seqpath=chr22seq, mean_rps=mean(c(exp1a_depth, exp1b_depth)), outdir='empirical_reads', seed=1247)





simulate_experiment(new_fasta_file, reads_per_transcript=readspertx, 
    num_reps=c(3,3), fold_changes=fc_mat, outdir='/hive/users/anjowall/projects/polyester_exploration/polyester_out/', paired = TRUE, bias = "rnaf", strand_specific = TRUE, gzip = TRUE)