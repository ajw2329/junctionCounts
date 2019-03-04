import splice_lib
import argparse
import subprocess


def parse_majiq_tsv(majiq_tsv, event_dict, junction_indexed_event_dict, outdir):

	majiq_file = open(majiq_tsv, 'r')

	outfile = open(outdir + "/majiq_file_with_events.tsv", 'w')

	header = majiq_file.readline().split()

	chrom_index = header.index("Chromosome")
	inc_junc_index = header.index("IncJunc")
	exc_junc_index = header.index("ExcJunc")
	strand_index = header.index("Strand")

	new_header = "\t".join(header + ["events"])
	outfile.write(new_header + "\n")

	for line in majiq_file:

		entry = line.split()
		chrom = entry[chrom_index]
		inc_junc = entry[inc_junc_index]
		exc_junc = entry[exc_junc_index]
		strand = entry[strand_index]

		inc_search_key = chrom + "_" + "_".join(inc_junc.split("-")) + "_" + strand
		exc_search_key = chrom + "_" + "_".join(exc_junc.split("-")) + "_" + strand

		putative_included_forms = set()
		putative_excluded_forms = set()

		if inc_search_key in junction_indexed_event_dict:

			putative_included_forms = set(["_".join(i.split("_")[0:-1]) for i in junction_indexed_event_dict[inc_search_key] if "included" in i])

		if exc_search_key in junction_indexed_event_dict:

			putative_excluded_forms = set(["_".join(i.split("_")[0:-1]) for i in junction_indexed_event_dict[exc_search_key] if "excluded" in i])

		putatively_matching_events = putative_included_forms & putative_excluded_forms
		if len(putatively_matching_events) > 0:
			events = ",".join(putatively_matching_events)
		else:
			events = "NA"

		outfile.write("\t".join(entry + [events]) + "\n")

	majiq_file.close()
	outfile.close()




def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--event_gtf", type = str, help = "Path to splice_lib event_gtf", required = True)
	parser.add_argument("--majiq_tsv", type = str, help = "Path to majiq tsv file", required = True)
	parser.add_argument("--outdir", type = str, help = "Path to output directory", required = True)
	args = parser.parse_args()

	subprocess.call("mkdir -p " + args.outdir, shell = True)

	event_dict = splice_lib.generate_standard_event_dict(args.event_gtf)
	splice_lib.complete_event_dict(event_dict)

	junction_indexed_event_dict = splice_lib.generate_junction_indexed_event_dict(event_dict)

	parse_majiq_tsv(args.majiq_tsv, event_dict, junction_indexed_event_dict, args.outdir)



if __name__ == "__main__":

	main()