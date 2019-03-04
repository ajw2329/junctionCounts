import argparse
import splice_lib
import subprocess

def associate_cassettes(event_dict, rmats_tsv, outdir):

	cassette_indexed_event_dict = {}

	rmats_file_events = open(outdir + "/rmats_file_events.tsv", 'w')

	for event in event_dict:

		if event_dict[event]["event_type"] == "SE":

			key = event_dict[event]["chrom"] + "_" + str(event_dict[event]["included_exons"][1][0]) + "_" + str(event_dict[event]["included_exons"][1][1]) + "_" + event_dict[event]["strand"]

			cassette_indexed_event_dict.setdefault(key, set()).add(event)

	rmats_file = open(rmats_tsv, 'r')

	header = rmats_file.readline().strip().split("\t")
	new_header = header + ["events"]

	print new_header
	rmats_file_events.write("\t".join(new_header) + "\n")

	chrom_index = header.index("Chromosome")
	strand_index = header.index("Strand")
	exon_start_index = header.index("Exon Start")
	exon_end_index = header.index("Exon End")

	for line in rmats_file:

		entry = line.strip().split("\t")

		chrom = entry[chrom_index]
		strand = entry[strand_index]
		exon_start = str(int(entry[exon_start_index]) + 1)
		exon_end = entry[exon_end_index]

		search_key = chrom + "_" + exon_start + "_" + exon_end + "_" + strand

		if search_key in cassette_indexed_event_dict:

			event_set = ",".join(cassette_indexed_event_dict[search_key])

		else:

			event_set = "NA"

		rmats_file_events.write("\t".join(entry + [event_set]) + "\n")


	rmats_file.close()
	rmats_file_events.close()


def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--event_gtf", type = str, help = "Path to event gtf", required = True)
	parser.add_argument("--rmats_tsv", type = str, help = "Path to rMATs rt-PCR tsv file", required = True)
	parser.add_argument("--outdir", type = str, help = "Path to output directory", required = True)
	args = parser.parse_args()

	subprocess.call("mkdir -p " + args.outdir, shell = True)

	event_dict = splice_lib.generate_standard_event_dict(args.event_gtf)
	splice_lib.complete_event_dict(event_dict)
	associate_cassettes(event_dict, args.rmats_tsv, args.outdir)


if __name__ == "__main__":

	main()