library("argparser")
library("tidyverse")
library("furrr")
library("splines")


version = "0.1.0"


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
			"--full_model",
			type = "character",
			help = "Full model to be fit by quasibinomial GLM")

	parser <- parser %>%
		add_argument(
			"--reduced_model",
			type = "character",
			help = "reduced model to be fit by GLM")

    parser <- parser %>%
		add_argument(
		"--num_rounds",
		type = "integer",
		help = paste0("Number of times to sample from min, max values and repeat GLM"),
		default = 10)

    parser <- parser %>%
		add_argument(
		"--frac_agreement",
		type = "numeric",
		help = paste0("Minimum fraction of rounds that must meet q < alpha"),
		default = 0.95)

    parser <- parser %>%
		add_argument(
		"--alpha",
		type = "numeric",
		help = paste0("Desired FDR"),
		default = 0.05)

    parser <- parser %>%
      add_argument(
        "--num_threads",
        type = "integer",
        help = paste0("Number of threads for multiprocessing"),
        default = 1)


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

	print("This script doesn't do anything yet. ")
}

if (sys.nframe() == 0){
  main()
}
