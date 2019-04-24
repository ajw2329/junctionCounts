library("argparser")
library("tidyverse")

parser <- arg_parser(description = "Calculate dPSI from junctionCounts output")
parser <- parser %>% 
	add_argument(
		"--all_junctioncounts", 
		type = "character", 
		help = "Path to concatenated junctionCounts count file")
parser <- parser %>% 
	add_argument(
		"--sample_info",
		type = "character",
		help = "Path to file containing sample info.  Should be a tsv with a header where the first column 'sample_name' contains the sample ids that were passed to junctionCounts, and the remaining columns contain grouping variables of interest")
parser <- parser %>% 
	add_argument(
		"--outdir",
		type = "character",
		help = "Path to output directory")

args <- parse_args(parser)


print("Reading data")


sample_info <- read.table(args$sample_info,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          header = TRUE)

## add check that first column name is sample_name

all_counts <- read.table(args$all_junctioncounts,
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         header = TRUE) %>% 
  left_join(sample_info, by = "sample_name") %>% 
  filter(!grepl("CF|CL", event_id))



grouping_variables <- colnames(sample_info)[2:length(sample_info)]

print("Calculating per-condition event span")


events_full_span <- all_counts %>% 
                      mutate(min_jc_row = pmap_dbl(list(x = min_IJC, y = min_SJC), 
                                                   function(x,y) max(x,y)),
                             sum_jc = min_IJC + min_SJC) %>%
                      group_by_at(c("event_id",
                               grouping_variables)) %>% 
                      filter(!all(is.na(PSI_lo)) & !all(is.na(PSI_hi))) %>%
                      summarize(min_PSI_lo = min(PSI_lo), 
                                max_PSI_hi = max(PSI_hi), 
                                min_jc_condition = min(min_jc_row),
                                min_sum_jc = min(sum_jc)) %>% 
                      ungroup()

print("Writing event span to output file")

events_full_span %>%
  write.table(file = paste0(args$outdir,"/events_full_span.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)

rm(all_counts)

print("Defining dPSI calculation functions")


conditional_filter <- function(df, 
                               column_name,
                               must_equal, 
                               filter_bool) {
  if (filter_bool) {
    #print(must_equal)
    return(filter(df, UQ(as.name(column_name)) == must_equal))
  } else {
    return(df)
  }
}


calc_dpsi <- function(events_full_span, 
                      comparisons, 
                      comparisons_factor, 
                      constant_factor = NA, 
                      constant_factor_value = NA, 
                      use_constant_factor = FALSE) {
all_condition_dpsi_df_list <- list()

 for (comparison in comparisons) {
   
   condition_1 <- events_full_span %>% 
   filter(UQ(as.name(comparisons_factor)) == comparison[1]) %>%
   conditional_filter(constant_factor, constant_factor_value, use_constant_factor) %>%
   transmute(event_id = event_id, cond1_min_PSI_lo = min_PSI_lo, cond1_max_PSI_hi = max_PSI_hi, min_jc_cond1 = min_jc_condition)
   
   condition_2 <- events_full_span %>% 
   filter(UQ(as.name(comparisons_factor)) == comparison[2]) %>%
   conditional_filter(constant_factor, constant_factor_value, use_constant_factor) %>%
   transmute(event_id = event_id, cond2_min_PSI_lo = min_PSI_lo, cond2_max_PSI_hi = max_PSI_hi, min_jc_cond2 = min_jc_condition)
   
   if (use_constant_factor) {
     #print("Using constant factor")
   dpsi <- condition_1 %>%
           left_join(condition_2, by = c("event_id")) %>%
           mutate(min_dpsi = cond2_min_PSI_lo - cond1_max_PSI_hi, 
                  max_dpsi = cond2_max_PSI_hi - cond1_min_PSI_lo, 
                  mid_dpsi = (min_dpsi + max_dpsi)/2, 
                  span_dpsi = max_dpsi - min_dpsi, 
                  comparison = paste0(comparison[1], 
                                      "_", 
                                      comparison[2]), 
                  UQ(as.name(constant_factor)) := constant_factor_value) %>% 
           select(event_id, 
                  span_dpsi, 
                  max_dpsi, 
                  min_dpsi, 
                  comparison, 
                  min_jc_cond1, 
                  min_jc_cond2, 
                  UQ(as.name(constant_factor))) %>%
           mutate(min_jc_comparison = pmap_dbl(list(x = min_jc_cond1, y = min_jc_cond2), 
                                               function(x,y) min(x,y)), 
                  inner_dpsi = pmap_dbl(list(x = min_dpsi, y = max_dpsi), 
                                        function(x,y) ifelse(sign(x) == sign(y), 
                                                             sign(x)*min(abs(x), abs(y)), 0)))
   
   } else {

   dpsi <- condition_1 %>%
           left_join(condition_2, by = c("event_id")) %>%
           mutate(min_dpsi = cond2_min_PSI_lo - cond1_max_PSI_hi, 
                  max_dpsi = cond2_max_PSI_hi - cond1_min_PSI_lo, 
                  mid_dpsi = (min_dpsi + max_dpsi)/2, 
                  span_dpsi = max_dpsi - min_dpsi, 
                  comparison = paste0(comparison[1], 
                                      "_", 
                                      comparison[2])) %>% 
           select(event_id, 
                  span_dpsi, 
                  max_dpsi, 
                  min_dpsi, 
                  comparison, 
                  min_jc_cond1, 
                  min_jc_cond2) %>%
           mutate(min_jc_comparison = pmap_dbl(list(x = min_jc_cond1, y = min_jc_cond2), 
                                               function(x,y) min(x,y)), 
                  inner_dpsi = pmap_dbl(list(x = min_dpsi, y = max_dpsi), 
                                        function(x,y) ifelse(sign(x) == sign(y), 
                                                             sign(x)*min(abs(x), abs(y)), 0)))
     
   }
   
   all_condition_dpsi_df_list <- c(all_condition_dpsi_df_list, list(dpsi))
   
 }

   all_condition_dpsi_df_actual <- do.call(rbind, all_condition_dpsi_df_list)
   all_condition_dpsi_df_actual <- all_condition_dpsi_df_actual %>%
      mutate(mid_dpsi = (max_dpsi + min_dpsi)/2)
   
   return(all_condition_dpsi_df_actual)
   

 }


 print("Calculating dPSI")


 ### Add fn for ddPSI calculation


ribo_fraction_levels = c("untr","mono","poll","polm","polh")
all_fraction_levels = c("cell","cyto","untr","mono","poll","polm","polh")


all_fraction_comparisons = list()

for (i in 1:(length(all_fraction_levels)-1)) {
  for (j in (i+1):length(all_fraction_levels)) {
    all_fraction_comparisons <- c(all_fraction_comparisons, list(c(all_fraction_levels[i], all_fraction_levels[j])))
  }
}

upf3b_dn_fraction_comparisons <- calc_dpsi(events_full_span, 
                                           all_fraction_comparisons, 
                                           "fraction", 
                                           "condition", 
                                           "upf3b_dn", 
                                           use_constant_factor = TRUE)

upf3b_di_fraction_comparisons <- calc_dpsi(events_full_span, 
                                           all_fraction_comparisons, 
                                           "fraction", 
                                           "condition", 
                                           "upf3b_di", 
                                           use_constant_factor = TRUE)


upf1x_dn_fraction_comparisons <- calc_dpsi(events_full_span, 
                                           all_fraction_comparisons, 
                                           "fraction", 
                                           "condition", 
                                           "upf1x_dn", 
                                           use_constant_factor = TRUE)

upf1x_di_fraction_comparisons <- calc_dpsi(events_full_span, 
                                           all_fraction_comparisons, 
                                           "fraction", 
                                           "condition", 
                                           "upf1x_di", 
                                           use_constant_factor = TRUE)



esc_fraction_comparisons <- calc_dpsi(events_full_span, 
                                           all_fraction_comparisons, 
                                           "fraction", 
                                           "condition", 
                                           "h9esc", 
                                           use_constant_factor = TRUE)

npc_fraction_comparisons <- calc_dpsi(events_full_span, 
                                           all_fraction_comparisons, 
                                           "fraction", 
                                           "condition", 
                                           "h9npc", 
                                           use_constant_factor = TRUE)

condition_comparisons_list <- list()

for (fraction in all_fraction_levels) {
  
  temp_df <- calc_dpsi(events_full_span,
            list(c("upf3b_dn","upf3b_di"), 
            	c("h9esc","h9npc"),
            	c("upf1x_dn","upf1x_di")),
            "condition",
            "fraction",
            fraction,
            use_constant_factor = TRUE)
  
  condition_comparisons_list <- c(condition_comparisons_list, list(temp_df))
}

condition_comparisons <- do.call(rbind, condition_comparisons_list)
rm(condition_comparisons_list)
rm(temp_df)

fraction_comparisons <- rbind(upf3b_di_fraction_comparisons,
                              upf3b_dn_fraction_comparisons,
                              upf1x_di_fraction_comparisons,
                              upf1x_dn_fraction_comparisons,
                              esc_fraction_comparisons,
                              npc_fraction_comparisons)

rm(upf3b_di_fraction_comparisons,
  upf3b_dn_fraction_comparisons,
  esc_fraction_comparisons,
  npc_fraction_comparisons,
  upf1x_di_fraction_comparisons,
  upf1x_dn_fraction_comparisons)


#rm(fraction_comparisons)
gc()

print("Writing dPSI values to output files")

fraction_comparisons %>%
  write.table(file = paste0(args$outdir,"/fraction_comparisons.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)


condition_comparisons %>%
  write.table(file = paste0(args$outdir,"/condition_comparisons.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)

print("Ceci n'est pas un algorithme bioinformatique.")