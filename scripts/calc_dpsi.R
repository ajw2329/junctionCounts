library("argparser")
library("tidyverse")
library("furrr")


version = "0.1.0"

pipeline_message <- function(x, message) {
  print(message)
  return(x)
}



calc_distance <- function(all_dpsi) {

  distances <- 
    all_dpsi %>% 
    group_by(
      event_id,
      event_type,
      comparison, 
      constant_vars) %>%
    filter(n() > 1) %>%
    mutate(min_mid_dpsi = min(mid_dpsi),
           max_mid_dpsi = max(mid_dpsi),
           min_inner_dpsi = min(inner_dpsi),
           max_inner_dpsi = max(inner_dpsi),
           minmax_mid_dpsi = 
            (mid_dpsi - min_mid_dpsi)/
            (max_mid_dpsi - min_mid_dpsi),
           minmax_inner_dpsi = 
            (inner_dpsi - min_inner_dpsi)/
            (max_inner_dpsi - min_inner_dpsi)) %>%
    summarize(
      sum_span_dpsi = sum(span_dpsi),
      dist = sum(abs(mid_dpsi)),
      inner_dist = sum(abs(inner_dpsi)),
      signed_dist = abs(sum(mid_dpsi)),
      signed_inner_dist = abs(sum(inner_dpsi)),
      minmax_dist = sum(abs(minmax_mid_dpsi)),
      minmax_inner_dist = sum(abs(minmax_inner_dpsi)),
      signed_minmax_dist = abs(sum(minmax_mid_dpsi)),
      signed_minmax_inner_dist = abs(sum(minmax_inner_dpsi)))

  return(distances)

}

get_all_ddpsi <- function(
  all_dpsi_df,
  level_list) {

  ddpsi_list <- list()

  print("Identifying input for secondary comparisons")

  print("Making nested source df and unlisting")

  scndry_comparison_source_list <-
    all_dpsi_df %>% 
    filter(constants != "none") %>%
    select(comparison, constants) %>%
    unique() %>%  
    nest(constants) %>%
    pull(data) %>%
    map(unlist) %>%
    unique()

  print(scndry_comparison_source_list)


  for (i in scndry_comparison_source_list) {

    broken = FALSE

    for (j in level_list) {
      if (all(i %in% j)) {

        print(
          paste0("Secondary comparison levels found ",
                 "in user-specified levels. Ordering appropriately"))

        scndry_comparisons <- identify_comparisons(j)

        print("Calculating ddPSI")

        ddpsi <- calc_ddpsi(all_dpsi_df,
                            scndry_comparisons)

        ddpsi_list <- c(ddpsi_list, list(ddpsi))

        broken = TRUE
        break
      }

    }

    if (!broken) {
      scndry_comparisons <- identify_comparisons(i)

      print(
        paste0("Secondary comparison levels not found ",
               "in user-specified levels. Ordering as found."))      

      print("Calculating ddPSI")

      ddpsi <- calc_ddpsi(all_dpsi_df,
                            scndry_comparisons)

      ddpsi_list <- c(ddpsi_list, list(ddpsi))

    } 

  }

  print("Binding ddPSI results into one data frame")

  all_ddpsi_df <- do.call(rbind, ddpsi_list)

  return(all_ddpsi_df)

}


calc_ddpsi <- function(
  all_dpsi_df,
  scndry_comparisons) {

  ddpsi_df_list <- list()

  for (scndry_comparison in scndry_comparisons) {

    condition_1 <- all_dpsi_df %>% 

       filter(
        constants == 
        scndry_comparison[1]) %>%

       mutate( 
        comp1_max_dpsi = max_dpsi, 
        comp1_min_dpsi = min_dpsi, 
        comp1_mid_dpsi = mid_dpsi,
        comp1_span_dpsi = span_dpsi,
        comp1_min_jc_comp = min_jc_comparison) %>%

       select(
        c(
          "event_id",
          "event_type",
          "comp1_max_dpsi",
          "comp1_min_dpsi",
          "comp1_mid_dpsi",
          "comp1_span_dpsi",
          "comp1_min_jc_comp",
          "comparison",
          "constant_vars"))


       condition_2 <- all_dpsi_df %>% 

       filter(
        constants == 
        scndry_comparison[2]) %>%

       mutate( 
        comp2_max_dpsi = max_dpsi, 
        comp2_min_dpsi = min_dpsi, 
        comp2_mid_dpsi = mid_dpsi,
        comp2_span_dpsi = span_dpsi,
        comp2_min_jc_comp = min_jc_comparison) %>%

       select(
        c(
          "event_id",
          "event_type",
          "comp2_max_dpsi",
          "comp2_min_dpsi",
          "comp2_mid_dpsi",
          "comp2_span_dpsi",
          "comp2_min_jc_comp",
          "comparison",
          "constant_vars"))   

    ddpsi <- condition_1 %>%

       left_join(condition_2, 
                 by = c("event_id",
                        "event_type",
                        "comparison",
                        "constant_vars")) %>%

       mutate(
              min_ddpsi = 
        comp2_min_dpsi - comp1_max_dpsi, 
              max_ddpsi = 
        comp2_max_dpsi - comp1_min_dpsi, 
              mid_ddpsi = 
        (min_ddpsi + max_ddpsi)/2, 
              span_ddpsi = 
        max_ddpsi - min_ddpsi, 
              scndry_comparison = 
                paste0(scndry_comparison[1], 
                       "_", 
                       scndry_comparison[2])) %>% 

       select(c("event_id",
                "event_type", 
                "span_ddpsi", 
                "max_ddpsi", 
                "min_ddpsi",
                "mid_ddpsi", 
                "scndry_comparison", 
                "comp1_min_jc_comp", 
                "comp2_min_jc_comp",
                "comparison",
                "constant_vars"
                )) %>%

       mutate(
        min_jc_comparison = future_pmap_int(
          list(x = comp1_min_jc_comp, y = comp2_min_jc_comp), 
          function(x,y) min(x,y), .progress = TRUE),
        comp1_min_jc_comp = as.integer(comp1_min_jc_comp),
        comp2_min_jc_comp = as.integer(comp2_min_jc_comp), 
        abs_mid_ddpsi = abs(mid_ddpsi),
        inner_ddpsi = future_pmap_dbl(
          list(x = min_ddpsi, y = max_ddpsi), 
          function(x,y) ifelse(sign(x) == sign(y), 
          sign(x)*min(abs(x), abs(y)), 0), .progress = TRUE),
        abs_inner_ddpsi = abs(inner_ddpsi))

      ddpsi_df_list <- c(
        ddpsi_df_list, 
        list(ddpsi))
  }

   ddpsi_df <- 
    do.call(
      rbind, 
      ddpsi_df_list)
     
  return(ddpsi_df)

}


get_all_dpsi <- function(
  all_counts,
  comparisons_constants,
  sample_info,
  levels_arg,
  level_list = NULL,
  outdir = NULL
  ) {

    dpsi_list <- list()

    print("Establishing requested comparisons")

    for (a in comparisons_constants) {

     for (i in names(a)) {
      print(i)
      if (!is.na(levels_arg)) {

        if (i %in% names(level_list)) {

          all_levels <- level_list[i][[1]]
          } else {

            all_levels <- sample_info %>%
                          pull(UQ(as.name(i))) %>%
                          unique()
          }
      } else {

        all_levels <- sample_info %>%
                      pull(UQ(as.name(i))) %>%
                      unique()
      }

      ## get comparisons

      all_comparisons <- 
        identify_comparisons(all_levels)

      print("Comparisons established.")

      ## get constant levels

      constants <- a[i][[1]]

      print("Calculating span")

      if (!is.na(constants)) {
        span <- calculate_span(
                  all_counts,
                  c(i, constants),
                  outdir)
      } else {
          span <- calculate_span(
                    all_counts,
                    i,
                    outdir)
      }

      print("Calculating dPSI")

      if (!is.na(constants)) {
        dpsi <- calc_dpsi(span, 
                          all_comparisons, 
                          i, 
                          constants)
        } else {

        dpsi <- calc_dpsi(span, 
                          all_comparisons, 
                          i)

        }

      dpsi_list <- c(dpsi_list, list(dpsi))

     }
    }

    print("Concatenating all dPSI dfs")

    dpsi_df <- do.call(rbind, dpsi_list)


}



calc_dpsi <- function(events_full_span, 
                      comparisons, 
                      comparisons_factor, 
                      constant_grouping_vars = c()) {

  dpsi_df_list <- list()

   for (comparison in comparisons) {

    print("Calculating dPSI for ")
    print(comparison)


    print("Filtering span for condition 1")
    print(comparisons_factor)
     
     condition_1 <- events_full_span %>% 

     filter(
      UQ(as.name(comparisons_factor)) == 
         comparison[1]) %>% 
     mutate( 
      cond1_min_min_psi = min_min_psi, 
      cond1_max_max_psi = max_max_psi, 
      cond1_min_jc = min_jc_condition) %>%
     select(
      c(
        "event_id",
        "event_type",
        "cond1_min_min_psi",
        "cond1_max_max_psi",
        "cond1_min_jc",
        constant_grouping_vars))

    print("Filtering span for condition 2")

     condition_2 <- events_full_span %>% 
     filter(
      UQ(as.name(comparisons_factor)) == 
         comparison[2]) %>%

     mutate(
      cond2_min_min_psi = min_min_psi, 
      cond2_max_max_psi = max_max_psi, 
      cond2_min_jc = min_jc_condition) %>%

     select(
      c(
        "event_id",
        "event_type",
        "cond2_min_min_psi",
        "cond2_max_max_psi",
        "cond2_min_jc",
        constant_grouping_vars))


     dpsi <- condition_1 %>%

             left_join(condition_2, 
                       by = c("event_id",
                              "event_type",
                              constant_grouping_vars)) %>%

             mutate(
                    min_dpsi = 
              cond2_min_min_psi - cond1_max_max_psi, 
                    max_dpsi = 
              cond2_max_max_psi - cond1_min_min_psi, 
                    mid_dpsi = 
              (min_dpsi + max_dpsi)/2, 
                    span_dpsi = 
              max_dpsi - min_dpsi, 
                    comparison = 
                      paste0(comparison[1], 
                             "_", 
                             comparison[2]),
              constant_vars = if_else(length(constant_grouping_vars) > 0, 
                                      paste(constant_grouping_vars,
                                            collapse = ","),
                                      "none")) %>% 

             select(c(
                    "event_id",
                    "event_type", 
                    "span_dpsi", 
                    "max_dpsi", 
                    "min_dpsi",
                    "mid_dpsi", 
                    "comparison", 
                    "cond1_min_jc", 
                    "cond2_min_jc",
                    constant_grouping_vars,
                    "constant_vars"
                    )) %>%

             mutate(
              min_jc_comparison = future_pmap_int(
                list(x = cond1_min_jc, y = cond2_min_jc), 
                function(x,y) min(x,y), .progress = TRUE),
              abs_mid_dpsi = abs(mid_dpsi), 
              inner_dpsi = future_pmap_dbl(
                list(x = min_dpsi, y = max_dpsi), 
                function(x,y) ifelse(sign(x) == sign(y), 
                sign(x)*min(abs(x), abs(y)), 0), .progress = TRUE),
              abs_inner_dpsi = abs(inner_dpsi),
              cond1_min_jc = as.integer(cond1_min_jc),
              cond2_min_jc = as.integer(cond2_min_jc)
              )

      dpsi_df_list <- c(
        dpsi_df_list, 
        list(dpsi))
     
     }
     

     dpsi_df <- 
      do.call(
        rbind, 
        dpsi_df_list)

     if(length(constant_grouping_vars) > 0) {

       dpsi_df <- 
        dpsi_df %>%
        mutate(
          comparison_factor = comparisons_factor) %>%
        unite(
          "constants", 
          constant_grouping_vars, 
          sep = ",") } else {

       dpsi_df <- 
        dpsi_df %>%
        mutate(
          comparison_factor = comparisons_factor,
          constants = "none")

        }
     
  return(dpsi_df)
   
 }


calculate_span <- function(
                    input_counts,
                    grouping_variables,
                    outdir = NULL) {

  print("Starting span calculation")

  span <- 

    input_counts %>%

    mutate(
      min_jc_row = future_pmap_dbl(
        list(x = min_ijc, y = min_ejc), 
        function(x,y) max(x,y), .progress = TRUE),
      sum_jc = min_ijc + min_ejc) %>%

    pipeline_message("Starting group_by_at in calculate_span") %>%

    group_by_at(
      c("event_id",
        "event_type",
        grouping_variables)) %>% 

    pipeline_message("Finished group_by_at. Starting summarize") %>%

    summarize(
           min_na = all(is.na(min_psi)),
           max_na = all(is.na(max_psi)),
           min_min_psi = min(min_psi),
           max_max_psi = max(max_psi),
           min_jc_condition = as.integer(min(min_jc_row)),
           min_sum_jc = as.integer(min(sum_jc))) %>%

    pipeline_message("Finished summarize. Starting filter") %>%

    select(-min_na, -max_na) %>%

    pipeline_message("Finished select. Starting ungroup") %>%

    ungroup()


  if(!is.null(outdir)) {

    filename = paste0(
        outdir,
        "/",
        paste(grouping_variables, collapse = "-"),
        "-",
        "event_span.tsv")

    span %>% 
    unite(
      "constants", 
      grouping_variables, 
      sep = ",") %>%
    mutate_if(is.double, ~sprintf("%.4f", .x)) %>%
    write.table(
      file = filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE)
  }

  return(span)

}


identify_comparisons <- function(input_levels) {

  all_comparisons <- list()


  for (i in 1:(length(input_levels)-1)) {
    for (j in (i+1):length(input_levels)) {
      all_comparisons <- c(
        all_comparisons, 
        list(c(input_levels[i], input_levels[j])))
    }
  }

  return(all_comparisons)
}


check_levels <- function(
  level_list,
  sample_info_names,
  sample_info) {

  ## check that level_list names match sample_info column names
  ## check also that level_list values match sample_info column vals

  for (i in names(level_list)) {

    stopifnot(i %in% sample_info_names)

    true_vals <- 
      sample_info %>%

      pull(UQ(as.name(i))) %>%
      
      unique()

    stopifnot(all(level_list[i][[1]] %in% true_vals))
  }
}


check_comparisons_constants <- function(
  comparisons_constants,
  sample_info_names) {
  ## check that provided factors for dPSI calculation
  ## are actually present in sample_info

  for (i in comparisons_constants) {
    for (j in names(i)) {
      stopifnot(j %in% sample_info_names)
      for (k in i[j]) {
        if ( !is.na(k)) {
          stopifnot(k %in% sample_info_names)
        }
      }
    }
  }
}


check_sample_info <- function(sample_info_names) {

  ## check that sample_name is first column in sample name

  stopifnot(sample_info_names[1] == "sample_name")

}

parse_levels <- function(levels_arg) {

  ### parse levels

  level_list <- list()

      if (!is.na(levels_arg)) {

        levels <- 
          stringr::str_split(
            levels_arg, 
            pattern = ",", 
            simplify = TRUE) %>%
          as.character()

        for (var_lev in levels) {
          var_lev_split <-
            stringr::str_split(
              var_lev,
              pattern = ":",
              simplify = TRUE) %>%
            as.character()

          var = var_lev_split[1]

          lev <- 
            stringr::str_split(
              var_lev_split[2],
              pattern = "-",
              simplify = TRUE) %>%
            as.character()

          level_list[var] = list(lev)
        }
    }

  return(level_list)

}


parse_comparisons <- function(
  comparisons_arg) {

  ### parse comparisons

    comparisons_constants <- 
      list()

    comparisons <- 
      stringr::str_split(
        comparisons_arg, 
        pattern = ",", 
        simplify = TRUE) %>%
      as.character()

    ## now store comparison factor, constant 
    ## factor as key-value pairs in list


    for (i in 
         1:length(comparisons)) {

      comparisons_constants[[i]] <- list()

      t <- 
        stringr::str_split(
          comparisons[i], 
            pattern = ":",
            simplify = TRUE) %>% 
        as.character()

      f <- t[1] #primary factor name

      if (length(t) > 1) {

      c <- 
        t[2] %>%
        stringr::str_split(
          pattern = "-",
          simplify = TRUE) %>%
        as.character()

      ## Note that the reason for creating lists 
      ## lists is to allow the same comparison factor(s)
      ## to be used with multiple different constant factors

      comparisons_constants[[i]][f] <- list(c)

      } else {
        comparisons_constants[[i]][f] <- NA_character_
      }

    }

  return(comparisons_constants)      
}




call_argparser <- function(
  arg_input = commandArgs(trailingOnly = TRUE)) {


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
        help = paste0("Path to file containing sample info. ",
                      "Should be a tsv with a header where the ",
                      "first column 'sample_name' contains the ",
                      "sample ids that were passed to junctionCounts, ",
                      "and the remaining columns contain grouping ",
                      "variables of interest"))

    parser <- parser %>% 
      add_argument(
        "--outdir",
        type = "character",
        help = "Path to output directory")

    parser <- parser %>%
      add_argument(
        "--comparisons",
        type = "character",
        help = paste0("Comma-separated list of factors between",
                      " whose levels dPSI should be calculated. ",
                      "Can optionally specify other factors across ",
                       "which the comparisons should be made ",
                       "individually as a hyphen-separated list ",
                       "separated from the primary factor with a ",
                       "colon.  E.g. providing ",
                       "condition:time,time:condition would result ",
                       "in 1) the calculation of dPSI between all ",
                       "pairs of condition values at each time value ",
                       "(timepoint) and 2) the calculation of dPSI ",
                       "between all pairs of timepoints at each ",
                       "condition value.  For a simple A vs B ",
                       "condition comparison, just provide the ",
                       "name of the condition variable in the ",
                       "provided sample info tsv."))

    parser <- parser %>%
      add_argument(
        "--levels",
        type = "character",
        help = paste0("Comma-separated list of lists of ",
                      "hyphen-separated factor levels to ",
                      "allow user-specified directionality of ",
                      "dPSI values.  E.g. ",
                      "time:week1-week2-week3,condition:no_treatment-treatment"),
        default = NA_character_)

    parser <- parser %>%
      add_argument(
        "--num_threads",
        type = "integer",
        help = paste0("Number of threads for multiprocessing"),
        default = 1)

    parser <- parser %>%
      add_argument(
        "--calc_ddpsi",
        help = "Calculate ddPSI (delta delta) values",
        flag = TRUE)

    parser <- parser %>%
      add_argument(
        "--calc_distance",
        help = paste0("Calculate distance (assorted ",
                      "aggregation of dPSI values across groups ",
                      "derived from --comparisons. ",
                      "E.g. in case of condition:fraction passed ",
                      "to --comparisons, the 'distance' across ",
                      "values of 'condition' would be calculated ",
                      "by summing the absolute values of dPSI ",
                      "for each pairwise condition comparison ",
                      "across different fractions."),
        flag = TRUE)

    parser <- parser %>%
      add_argument(
        "--version",
        help = "Print version and exit.",
        flag = TRUE)


    args <- parse_args(parser)

    if (args$version) {
      print(paste0("calc_dpsi.R is part of junctionCounts version ", version))

      quit(save = "no")
    }

    return(args)
}



main <- function() {

  print("Parsing CL arguments")
  args <- call_argparser()

  options(future.globals.maxSize = 2097152000)

  plan(multicore, workers = args$num_threads)

  print("Parsing provided comparisons")
  comparisons_constants <- 
    parse_comparisons(args$comparisons)


  any_na <- c()

  for (i in comparisons_constants) {
    for(j in i) {
      any_na <- c(all(is.na(j)), any_na)
    }
  }

  all_na_constants <- all(any_na)

  print("Parsing provided levels (if any)")
  level_list <- 
    parse_levels(args$levels)

  print("Reading in sample info")
  sample_info <- read.table(
    args$sample_info,
    sep = "\t",
    stringsAsFactors = FALSE,
    header = TRUE)

  sample_info_names <- 
    colnames(sample_info)

  print("Checking comparisons and levels against sample info")

  check_sample_info(sample_info_names)

  check_comparisons_constants(
    comparisons_constants,
    sample_info_names)

  check_levels(
    level_list,
    sample_info_names,
    sample_info)

  print("Importing counts table.  This may take a while . . . ")

  all_counts <- read.table(
      args$all_junctioncounts,
      sep = "\t",
      stringsAsFactors = FALSE,
      header = TRUE) %>% 
    left_join(sample_info, by = "sample_name")

  for ( i in sample_info_names[2:length(sample_info_names)] ) {
    all_counts[[i]] <- as.factor(as.character(all_counts[[i]]))
  }

  print("Calculating span and dPSI for all comparisons.  This may take a while . . . ")

  all_dpsi <- get_all_dpsi(
    all_counts,
    comparisons_constants,
    sample_info,
    args$levels,
    level_list = level_list,
    outdir = args$outdir)

  print("Writing dPSI output")

  all_dpsi %>% 
    mutate_if(is.double, ~sprintf("%.4f", .x)) %>%
    write.table(
      file = paste0(
        args$outdir,
        "/all_dpsi.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE)

  if (args$calc_ddpsi) {

      if (!all_na_constants) {

        print("ddPSI calculation indicated. Calculating ddPSI values")

        all_ddpsi <- get_all_ddpsi(
          all_dpsi,
          level_list)

        print("Writing ddPSI output")

        all_ddpsi %>%
          mutate_if(is.double, ~sprintf("%.4f", .x)) %>%
          write.table(
            file = paste0(
              args$outdir,
              "/all_ddpsi.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE) } else {

            print("Warning: ddPSI calculation indicated but no levels across which it could be calculated. Skipping . . . ")
          }
    }

  if (args$calc_distance) {

    if(!all_na_constants) {

      print("Distance calculation indicated. Calculating distance values.")

      distances <- calc_distance(all_dpsi)

      print("Writing distance output")

      distances %>%
        mutate_if(is.double, ~sprintf("%.4f", .x)) %>%
        write.table(
          file = paste0(
            args$outdir,
            "/all_distances.tsv"),
          sep = "\t",
          quote = FALSE,
          row.names = FALSE,
          col.names = TRUE)  

    } else {

        print("Warning: distance calculation indicated but no levels across which it could be calculated. Skipping . . . ")
      }
  }

  print("calc_dpsi.R finished.")

  print("Ceci n'est pas un algorithme bioinformatique.")
}



if (sys.nframe() == 0){
  main()
}
