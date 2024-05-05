#!/usr/bin/env Rscript

if (!require("pacman")) install.packages("pacman")
pacman::p_load(DEXSeq, dplyr, tidyr, optparse)

option_list = list(
  make_option(c("-e", "--events"), action = "store", type="character", default=NULL, 
              help="path to infer_pairwise_events directory", metavar="infer_pairwise_events_dir"),
  make_option(c("-p", "--psi_files"), action = "store", type="character", default=NULL, 
              help="path to junctionCounts *psi_outfile.tsv files", metavar="junctionCounts_psi_tsv_dir"),
  make_option(c("-s", "--sample_table"), action = "store", type="character", default=NULL, 
              help="path to sample table: sample,condition,file (.csv)", metavar="sample_table"),
  make_option("--min_jc", action = "store", type="double", default=15.0, 
              help="minimum number of total mean junction read support across conditions to report an event (default: 15.0)", metavar="minimum_jc"),
  make_option("--min_psi", action = "store", type="double", default=0.1, 
              help="minimum max(mean(PSI)) across conditions to report an event  (default: 0.1)", metavar="minimum_psi"),
  make_option("--ri_span", action = "store", type="double", default=0.03, 
              help="maximum max(mean(span_PSI)) across conditions to report a RI/MR event.\n
              span_PSI is the difference in PSI between the flanking junctions;\n
              the smaller the value, the more stringent  (default: 0.03).", metavar="ri_span_psi")) 

opt <- parse_args(OptionParser(option_list=option_list))

splice_lib_events <- opt$events
dataDir <- opt$psi_files
sampleTable <- read.csv(opt$sample_table) %>%
  mutate(count_file = paste0(dataDir, "/", gsub('.{4}$', '', file), ".txt"))

# Collect event data.
coords <- read.csv(file.path(splice_lib_events, 'splice_lib_events.bed'), sep='\t', header=F)[c(1, 2, 3, 4, 6)]
genes <- read.csv(file.path(splice_lib_events, 'splice_lib_events.ioe'), sep='\t', header=T)[2:3] %>%
  separate(event_id, c("gene", "event"), ";") %>% 
  left_join(coords, by=c("event"="V4")) %>%
  dplyr::select(-2) %>%
  setNames(., c("gene", "event_id", "chr", "start", "end", "strand"))

# Generate the flattened event exonic part file for DEXSeq.
gtf <- read.csv(file.path(splice_lib_events, 'splice_lib_events.bed'), sep='\t', header=F)
aggregate <- gtf %>% mutate(source = 'splice_lib_event', feature = 'aggregate_gene', name = V4, score = '.', 
                            attribute = paste0('gene_id "', V4, '"'))
included <- gtf %>% mutate(source = 'splice_lib_event', feature = 'exonic_part', name = V4, score = '.',
                           attribute = paste0('transcripts "', V4, 
                                              '_included"; exonic_part_number "001"; gene_id "', V4, '"'))
excluded <- gtf %>% mutate(source = 'splice_lib_event', feature = 'exonic_part', name = V4, score = '.',
                           attribute = paste0('transcripts "', V4, 
                                              '_excluded"; exonic_part_number "002"; gene_id "', V4, '"'))
gff <- rbind(aggregate, included, excluded) %>% dplyr::select(1, 7, 8, 2, 3, 9, 6, 10, 11) %>% arrange(name, feature)
write.table(gff, file.path(splice_lib_events, 'splice_lib_events_dexseq.gff'), row.names=F, col.names=F, quote=F, sep="\t")

# Generate separate junction read count files for the included and excluded form of events.
for (sample in sampleTable$file) {
  name <- paste0(gsub('.{4}$', '', sample), ".txt")
  x <- read.csv(file.path(dataDir, sample), sep="\t")
  ijc <- gff %>% filter(grepl('exonic_part_number "001"', attribute)) %>%
    left_join(x[c("event_id", "avg_ijc")], by=c("name"="event_id")) %>% 
    rename(counts = avg_ijc) %>% 
    mutate(name = paste0(name, ":001")) %>% dplyr::select(6, 10)
  ejc <- gff %>% filter(grepl('exonic_part_number "002"', attribute)) %>%
    left_join(x[c("event_id", "avg_ejc")], by=c("name"="event_id")) %>% 
    rename(counts = avg_ejc) %>% 
    mutate(name = paste0(name, ":002")) %>% dplyr::select(6, 10)
  counts <- rbind(ijc, ejc) %>% arrange(name) %>%
    mutate_at(2, ~replace_na(.,0)) %>%
    mutate(counts = floor(0.5 + counts))
  write.table(counts, file.path(dataDir, name), row.names=F, col.names=F, quote=F, sep="\t")
}

# Compare the dispersion of included and excluded junction read counts using DEXSeq.
print("Initializing DEXSeq. Petit à petit, l'oiseau fait son nid...")
counts <- sampleTable$count_file
flattenedFile <- file.path(splice_lib_events, 'splice_lib_events_dexseq.gff')
dxd <- DEXSeqDataSetFromHTSeq(
     counts,
     sampleData=sampleTable,
     design= ~ sample + exon + condition:exon,
     flattenedfile=flattenedFile)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
dxr <- DEXSeqResults(dxd)
event_qval <- as.data.frame(perGeneQValue(dxr)) %>%
  tibble::rownames_to_column("event_id") %>% 
  setNames(., c("event_id", "event_qval")) 
  
# Calculate mean PSI, ijc and ejc per condition, then collate data into final output.
collectData <- function(cond) {
  events <- read.csv(file.path(dataDir, sampleTable$file[1]), sep="\t")[2] %>% arrange(event_id)
  meanPsi <- events
  spanPsi <- events
  ijc <- events
  ejc <- events
  
  for (file in subset(sampleTable, sampleTable$condition == cond)$file) {
    x <- read.csv(file.path(dataDir, file), sep="\t")[c(2, 6, 12, 10, 11)]
    meanPsi <- meanPsi %>% left_join(x[1:2], by="event_id")
    spanPsi <- spanPsi %>% left_join(x[c(1, 3)], by="event_id")
    ijc <- ijc %>% left_join(x[c(1, 4)], by="event_id")
    ejc <- ejc %>% left_join(x[c(1, 5)], by="event_id")
  }
  
  meanPsi <- meanPsi %>% replace(is.na(.), 0) %>% mutate(mean_psi = select(., -c(1)) %>% rowMeans(na.rm = TRUE))
  spanPsi <- spanPsi %>% replace(is.na(.), 0) %>% mutate(mean_spanPsi = select(., -c(1)) %>% rowMeans(na.rm = TRUE))
  ijc <- ijc %>% replace(is.na(.), 0) %>% mutate(mean_ijc = select(., -c(1)) %>% rowMeans(na.rm = TRUE))
  ejc <- ejc %>% replace(is.na(.), 0) %>% mutate(mean_ejc = select(., -c(1)) %>% rowMeans(na.rm = TRUE))
  quantData <- cbind(ijc[c("event_id", "mean_ijc")], ejc["mean_ejc"],
                   meanPsi["mean_psi"], spanPsi["mean_spanPsi"])
  names(quantData)[2] <- paste0(cond, "_mean_ijc")
  names(quantData)[3] <- paste0(cond, "_mean_ejc")
  names(quantData)[4] <- paste0(cond, "_mean_psi")
  names(quantData)[5] <- paste0(cond, "_mean_spanPsi")
  
  return(quantData)
}

quantA <- collectData(unique(sampleTable$condition)[1])
quantB <- collectData(unique(sampleTable$condition)[2])
psi <- cbind(quantA[c(1, 4)], quantB[4]) %>% setNames(., c("event_id", "a", "b")) %>% 
  filter(pmax(a, b, na.rm = TRUE) >= opt$min_psi)
spanPsi <- cbind(quantA[c(1, 5)], quantB[5]) %>% setNames(., c("event_id", "a", "b")) %>% 
  filter(grepl('R', event_id), pmax(a, b, na.rm = TRUE) > opt$ri_span)

merged <- cbind(quantA[1:4], quantB[2:4]) %>%
  filter(rowSums(dplyr::select(., c(colnames(.)[c(2, 3, 5, 6)])), na.rm = TRUE) >= opt$min_jc,
         event_id %in% psi$event_id, !(event_id %in% spanPsi$event_id)) %>%
  mutate(dpsi = get(names(.)[7]) - get(names(.)[4])) %>%
  left_join(genes, by="event_id") %>%
  left_join(event_qval, by="event_id") %>%
  mutate_at(14, ~replace_na(., 1)) %>%
  mutate(sig = ifelse((dpsi >= 0.1) & (event_qval <= 0.05), 1, 0),
         event_type = substr(event_id, start = 1, stop = 2)) %>%
  mutate(sig = ifelse((dpsi <= -0.1) & (event_qval <= 0.05), -1, sig)) %>%
  dplyr::select(9, 1, 16, 10, 11, 12, 13, 2, 3, 4, 5, 6, 7, 8, 14, 15)

outName <- paste0(unique(sampleTable$condition)[1], "_", unique(sampleTable$condition)[2], "_dpsi.csv")
write.csv(merged, file.path(dataDir, outName), row.names=F)

print("Statistical testing complete. C’est la fin des haricots.")
