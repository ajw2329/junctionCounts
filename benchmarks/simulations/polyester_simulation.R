library(polyester)
library(truncnorm)
library(aster)
library(dplyr)
library(Biostrings)


fasta_file <- '/hive/users/anjowall/projects/polyester_exploration/gencode.v26.annotation.50000_tx.fa'
new_fasta_file <- '/hive/users/anjowall/projects/polyester_exploration/gencode.v26.annotation.50000_tx_subset.fa'
kallisto_abundance <- "/hive/users/anjowall/projects/polyester_exploration/h9esc_cyto_jd004_kallisto_out/abundance.tsv"


filename <- "abundance.tsv"

import_tpms <- function(samples) {

	df_list = list()

	for (sample in samples) {

		path = paste0(exp1_base, sample, "/", filename)
		temp_df <- read.table(path, 
			header = TRUE, 
			sep = "\t",
			stringsAsFactors = FALSE) %>% 
		mutate(sample_id = sample) %>%
		rename(transcript_id = target_id)
		exp1a_df_list <- c(df_list, list(temp_df))

	}

	mean_tpm_df <- do.call(rbind, df_list) %>%
		group_by(transcript_id) %>% 
		summarize(mean_tpm = mean(tpm), 
			mean_est_counts = mean(est_counts))

	return(mean_tpm_df)

}


###experiment 1 

###condition a file paths

exp1_base <- "/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/kallisto/"

exp1a_samples <- c("SRR536342","SRR536344","SRR536346")

exp1b_samples <- c("SRR536348","SRR536350","SRR536352")

###exp1a mean tpms

exp1a_mean_tpm <- import_tpms(exp1a_samples)
exp1b_mean_tpm <- import_tpms(exp1b_samples)

exp1a_depth <- exp1a_mean_tpm %>%
	summarize(depth = sum(mean_est_counts)) %>%
	pull(depth)

exp1b_depth <- exp1b_mean_tpm %>%
	summarize(depth = sum(mean_est_counts)) %>%
	pull(depth)


## get fold changes






frac_diff <- 0.2

fasta <- readDNAStringSet(fasta_file)

names(fasta) <- lapply(strsplit(names(fasta), "\\s"), function(x) x[1]) %>% unlist()


a <- exp(1)^c(rtruncnorm(mean = 1.5, sd = 1, a = 1, n = 10000), rep(0, 40027))
b <- exp(1)^c(rep(0, 40027), rtruncnorm(mean = 1.5, sd = 1, a = 1, n = 10000))

fc_df <- data.frame(a = a, b = b) %>% sample_n(size = 50027)
fc_mat <- as.matrix(fc_df)
head(fc_mat)

#readspertx <- rktnb(n = 50000, mu = 1000, size = 0.1, k = 0)

tpm <- read.table(kallisto_abundance, header = TRUE, sep = "\t") %>% dplyr::filter(est_counts >= 1)##read in kallisto abundance.tsv
head(tpm)
desired_depth <- 40000000
print(desired_depth)
actual_depth <- sum(tpm$est_counts)
print(actual_depth)
readspertx <- tpm %>% mutate(modified_counts = desired_depth*est_counts/actual_depth) %>% pull(modified_counts) %>% round()
head(readspertx)
length(readspertx)

tx_count <- nrow(tpm)
diff_count <- round(frac_diff*tx_count)
nondiff_count <- tx_count - diff_count

tx_count
diff_count
nondiff_count

fasta_subset <- fasta[tpm$target_id]
writeXStringSet(fasta_subset, new_fasta_file)

a <- exp(1)^c(rtruncnorm(mean = 1.5, sd = 1, a = 1, n = diff_count), rep(0, nondiff_count))
b <- exp(1)^c(rep(0, nondiff_count), rtruncnorm(mean = 1.5, sd = 1, a = 1, n = diff_count))

fc_df <- data.frame(a = a, b = b) %>% sample_n(size = tx_count)
fc_mat <- as.matrix(fc_df)
head(fc_mat)


dim(fc_mat)

head(fc_mat*readspertx)

simulate_experiment(new_fasta_file, reads_per_transcript=readspertx, 
    num_reps=c(3,3), fold_changes=fc_mat, outdir='/hive/users/anjowall/projects/polyester_exploration/polyester_out/', paired = TRUE, bias = "rnaf", strand_specific = TRUE, gzip = TRUE)