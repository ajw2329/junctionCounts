import sys
import re
import argparse
import splice_lib
import copy
import subprocess
import pickle
import gzip
from datetime import datetime
import time
import numpy as np
from ncls import NCLS
from random import randint
import pysam
from collections import defaultdict



def create_eij_ncls_dict(standard_event_dict):

	eij_by_chrom_strand = {}
	eij_indexed_event_dict = {}
	eij_only_count_dict = {}
	ncls_by_chrom_strand = {}

	for chrom in standard_event_dict:

		for strand in standard_event_dict[chrom]:

			for event in standard_event_dict[chrom][strand]:

				for eij in standard_event_dict[chrom][strand][event]["included_ei_junctions"]:

					eij_by_chrom_strand.setdefault(chrom, {}).setdefault(strand, set()).add(int(eij))

					eij_index = chrom + "_" + str(eij) + "_" + strand

					eij_indexed_event_dict.setdefault(chrom, {}).setdefault(strand, {}).setdefault(eij_index, set()).add(event + "_included")

					eij_only_count_dict.setdefault(eij_index, 0)

				for eij in standard_event_dict[chrom][strand][event]["excluded_ei_junctions"]:

					eij_by_chrom_strand.setdefault(chrom, {}).setdefault(strand, set()).add(int(eij))

					eij_index = chrom + "_" + str(eij) + "_" + strand

					eij_indexed_event_dict.setdefault(chrom, {}).setdefault(strand, {}).setdefault(eij_index, set()).add(event + "_excluded")

					eij_only_count_dict.setdefault(eij_index, 0)


	for chrom in eij_by_chrom_strand:

		ncls_by_chrom_strand[chrom] = {}

		for strand in eij_by_chrom_strand[chrom]:

			starts = np.array(list(eij_by_chrom_strand[chrom][strand])) - 1 ## to make 0-based
			ends = starts
			ids = starts

			ncls_by_chrom_strand[chrom][strand] = NCLS(starts,ends,ids)

	return ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict



def read_pair_generator(bam, region_string = None):
	"""
    Slight modification of function written by biostars user gizone1: https://www.biostars.org/p/306041/
	Generate read pairs in a BAM file or within a region string.
	Reads are added to read_dict until a pair is found.
	"""
	read_dict = defaultdict(lambda: [None, None])
	for read in bam.fetch(region = region_string):
		if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
			continue
		qname = read.query_name
		if qname not in read_dict:
			if read.is_read1:
				read_dict[qname][0] = read
			else:
				read_dict[qname][1] = read
		else:
			if read.is_read1:
				yield read, read_dict[qname][1]
			else:
				yield read_dict[qname][0], read
			del read_dict[qname]


def parse_reads(read_list, forward_read):
	'''
		Takes as input a read or read pair in a list of pysam AlignedSegments, returns a dictionary containing a list of junctions (where a junction is described by the chromosome, the start/end position, and the strand, all joined as a string), as well as (separately) the chrom and strand information.  Some code derived from find_introns_slow() in pysam.
	'''

	chrom_list = list(set([i.reference_name for i in read_list]))

	if len(chrom_list) > 1: ##implies some sort of chimeric read - not currently handled

		return None

	else:

		chrom = re.sub("_","&",chrom_list[0].strip())

	exons = []
	junctions = []


	if forward_read == "R1":

		if read_list[0].is_reverse:

			strand = "-"

		else:

			strand = "+"

	elif forward_read == "R2":

		if read_list[1].is_reverse:

			strand = "-"

		else:

			strand = "+"

	elif forward_read == "unstranded":

		strand = ""

	else:

		sys.exit("Improperly designated forward read in parse_read().  Exiting . . . ")


	#############derived from pysam find_introns()########################

	match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
	BAM_CREF_SKIP = 3

	for read in read_list:

		exons += [list(i) for i in read.get_blocks()]

		#############derived from pysam find_introns_slow()########################

		#if 'N' in read.cigarstring:
		#	last_read_pos = False
		#	for read_loc, genome_loc in read.get_aligned_pairs():
		#		if read_loc is None and last_read_pos:
		#			start = genome_loc
		#		elif read_loc and last_read_pos is None:
		#			stop = genome_loc  # we are right exclusive ,so this is correct
		#			junctions.append(chrom + "_" + str(start) + "_" + str(stop + 1) + "_" + strand)
		#			del start
		#			del stop
		#		last_read_pos = read_loc

		base_position = read.pos

		for op, nt in read.cigartuples:
			if op in match_or_deletion: 
				base_position += nt
			elif op == BAM_CREF_SKIP: 
				junc_start = base_position
				base_position += nt
				junctions.append(chrom + "_" + str(junc_start) + "_" + str(base_position + 1) + "_" + strand)

		############################################################################


	if len(exons) > 1:
		
		outdict = {"junctions": list(set(junctions)), "exons": exons, "strand": strand, "chrom": chrom}

	else:

		outdict = {"junctions": [], "exons": exons, "strand": strand, "chrom": chrom}

	return outdict

				
def intersection_collapse(input_list):

	'''
		Takes as input a list of lists of event IDs matched by junctions in a read.  Attempts to eliminate ambiguous assignments by recursively identifying pairwise intersections of the sublists.
		Returns a final minimal list of lists of event IDs

	'''

	input_set = map(set, input_list)

	intersections = []

	#print input_set

	condition = True
	while condition:

		for i in range(0,len(input_set)-1):

			for j in range(i+1,len(input_set)):

				intersect = input_set[i] & input_set[j]

				if len(intersect) > 0:

					del input_set[i], input_set[j-1]

					intersections.append(list(intersect))

					break
					condition = False

		else:
			condition = False

	if len(intersections) > 0:

		input_set.extend(intersections)

		return intersection_collapse(input_set)

	elif len(intersections) == 0:

		#print "No more recursion!"
		#print input_set
		minimal_list = map(list, input_set)

		return minimal_list


def process_reads(bam_filename, junction_indexed_event_dict, junction_only_count_dict, standard_event_dict, ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict, forward_read = "R2", single_end = False, bootstraps = False):

	'''
		Primary read assignment functino - takes as input a junction read filename (i.e. a BED12 file - can be gzipped), parses the read with "parse_junction_read", collapses matching events with "intersection_collapse", and adds 1 count to all surviving event matches in the event ID-indexed dictionary.
	'''


	all_read_info = []


	bam = pysam.AlignmentFile(bam_filename, 'rb')

	if single_end:

		for read in bam:

			read_properties = parse_reads([read], forward_read)
			read_info = assign_reads(read_properties, junction_indexed_event_dict, junction_only_count_dict, standard_event_dict, ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict, forward_read, bootstraps)
			if bootstraps:
				all_read_info.append(read_info)

	else:

		for read1, read2 in read_pair_generator(bam):

			read_properties = parse_reads([read1, read2], forward_read)
			read_info = assign_reads(read_properties, junction_indexed_event_dict, junction_only_count_dict, standard_event_dict, ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict, forward_read, bootstraps)
			if bootstraps:
				all_read_info.append(read_info)

	return all_read_info


def assign_reads(read_properties, junction_indexed_event_dict, junction_only_count_dict, standard_event_dict, ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict, forward_read, bootstraps):

		exons = read_properties["exons"]
		junctions = read_properties["junctions"]
		strand = read_properties["strand"]
		chrom = read_properties["chrom"]

		overlapping_eij = set()

		candidate_isoforms = []

		event_junction_dict = {}
		event_eij_dict = {}
		possible_strands = set() ### for use when forward_read == "unstranded"


		if forward_read != "unstranded":

			for junction in junctions:

				if chrom in junction_indexed_event_dict:

					if strand in junction_indexed_event_dict[chrom]: 

						if junction in junction_indexed_event_dict[chrom][strand]:

							#junction_only_count_dict[chrom][strand][junction] += 1

							matching_events = junction_indexed_event_dict[chrom][strand][junction][:]

							if matching_events not in candidate_isoforms:

								candidate_isoforms.append(matching_events)

							for matching_event in matching_events:

								event_junction_dict.setdefault(matching_event, set()).add(junction)

			for exon in exons:

				if chrom in ncls_by_chrom_strand:

					if strand in ncls_by_chrom_strand[chrom]:

						[overlapping_eij.add(chrom + "_" + str(i[0] + 1) + "_" + strand) for i in ncls_by_chrom_strand[chrom][strand].find_overlap(exon[0], exon[1])]


			for eij in overlapping_eij:

				if chrom in eij_indexed_event_dict:

					if strand in eij_indexed_event_dict[chrom]:

						if eij in eij_indexed_event_dict[chrom][strand]:

							#eij_only_count_dict[eij] += 1

							matching_events = list(eij_indexed_event_dict[chrom][strand][eij])

							if matching_events not in candidate_isoforms:

								candidate_isoforms.append(matching_events)

							for matching_event in matching_events:

								event_eij_dict.setdefault(matching_event, set()).add(eij)

		else:

			junctions = [i + "+" for i in junctions] + [i + "-" for i in junctions]

			for junction in junctions:

				if chrom in junction_indexed_event_dict:

					for local_strand in junction_indexed_event_dict[chrom]: 

						if junction in junction_indexed_event_dict[chrom][local_strand]:

							possible_strands.add(local_strand)

							#junction_only_count_dict[chrom][strand][junction] += 1

							matching_events = junction_indexed_event_dict[chrom][local_strand][junction]

							if matching_events not in candidate_isoforms:

								candidate_isoforms.append(matching_events)	

							for matching_event in matching_events:

								event_junction_dict.setdefault(matching_event, set()).add(junction)

			for exon in exons:

				if chrom in ncls_by_chrom_strand:

					for local_strand in ncls_by_chrom_strand[chrom]:

						[overlapping_eij.add(chrom + "_" + str(i[0] + 1) + "_" + local_strand) for i in ncls_by_chrom_strand[chrom][local_strand].find_overlap(exon[0], exon[1])]


			for eij in overlapping_eij:

				if chrom in eij_indexed_event_dict:

					test_strand = eij.split("_")[-1]

					if test_strand in eij_indexed_event_dict[chrom]:

						if eij in eij_indexed_event_dict[chrom][test_strand]:

							possible_strands.add(test_strand)

							#eij_only_count_dict[eij] += 1
							matching_events = list(eij_indexed_event_dict[chrom][test_strand][eij])

							if matching_events not in candidate_isoforms:

								candidate_isoforms.append(matching_events)

							for matching_event in matching_events:

								event_eij_dict.setdefault(matching_event, set()).add(eij)


		if len(candidate_isoforms) > 0:

			minimal_candidate_isoform_list = intersection_collapse(candidate_isoforms)

			#print minimal_candidate_isoform_list
			#print event_junction_dict
			#print event_eij_dict

			for i in minimal_candidate_isoform_list:

				for j in i:

					identifier = j.split("_")
					
					event_id = "_".join(identifier[0:-1])
					event_form = identifier[-1]

					#standard_event_dict[chrom][strand][event_id][event_form + "_" + "count"] += 1 ## total count

					if forward_read != "unstranded":

						if len(junctions) > 0:

							if j in event_junction_dict:

								for junction in event_junction_dict[j]:

									standard_event_dict[chrom][strand][event_id][event_form + "_junction_counts"][junction] += 1

						if len(overlapping_eij) > 0:

							if j in event_eij_dict:

								for eij in event_eij_dict[j]:

									standard_event_dict[chrom][strand][event_id][event_form + "_junction_counts"][eij] += 1

					else:

						if len(junctions) > 0:

							if j in event_junction_dict:

								for junction in event_junction_dict[j]:

									for test_strand in possible_strands:

										if event_id in standard_event_dict[chrom][test_strand]:

											standard_event_dict[chrom][test_strand][event_id][event_form + "_junction_counts"][junction] += 1

						if len(overlapping_eij) > 0:

							if j in event_eij_dict:

								for eij in event_eij_dict[j]:

									for test_strand in possible_strands:

										if event_id in standard_event_dict[chrom][test_strand]:

											standard_event_dict[chrom][test_strand][event_id][event_form + "_junction_counts"][eij] += 1


		else:

			minimal_candidate_isoform_list = []

		if bootstraps:

			read_info = {"minimal_candidate_isoform_list": minimal_candidate_isoform_list,
						 "junctions": event_junction_dict,
						 "eij": event_eij_dict,
						 "chrom": chrom,
						 "strand": strand,
						 "possible_strands": possible_strands}

			return read_info


def bootstrap_junction_counts(all_read_info, junction_only_count_dict, standard_event_dict, eij_only_count_dict, junction_indexed_event_dict, eij_indexed_event_dict, n_reads, forward_read):

	for n in range(0,n_reads):

		read = all_read_info[randint(0,n_reads - 1)]

		minimal_candidate_isoform_list = read["minimal_candidate_isoform_list"]
		event_junction_dict = read["junctions"]
		event_eij_dict = read["eij"]
		chrom = read["chrom"]
		strand = read["strand"]
		possible_strands = read["possible_strands"]

		#for junction in junctions:

		#	if junction in junction_indexed_event_dict:

		#		junction_only_count_dict[junction] += 1

		#for eij in overlapping_eij:

		#	if eij in eij_indexed_event_dict:

		#		eij_only_count_dict[eij] += 1


		for i in minimal_candidate_isoform_list:

			for j in i:

				identifier = j.split("_")
				
				event_id = "_".join(identifier[0:-1])
				event_form = identifier[-1]

				if forward_read != "unstranded":

					if j in event_junction_dict:

						for junction in event_junction_dict[j]:

							standard_event_dict[chrom][strand][event_id][event_form + "_junction_counts"][junction] += 1


					if j in event_eij_dict:

						for eij in event_eij_dict[j]:

							standard_event_dict[chrom][strand][event_id][event_form + "_junction_counts"][eij] += 1

				else:

					if j in event_junction_dict:

						for junction in event_junction_dict[j]:

							for test_strand in possible_strands:

								if event_id in standard_event_dict[chrom][test_strand]:

									standard_event_dict[chrom][test_strand][event_id][event_form + "_junction_counts"][junction] += 1

					if j in event_eij_dict:

						for eij in event_eij_dict[j]:

							for test_strand in possible_strands:

								if event_id in standard_event_dict[chrom][test_strand]:

									standard_event_dict[chrom][test_strand][event_id][event_form + "_junction_counts"][eij] += 1




def get_exon_edge_counts(junction_only_count_dict, junction_indexed_event_dict, standard_event_dict):

	'''
		Creates dictionary of edge counts (edge is one-half of junction - an edge's counts is the sum of all of the junction counts that include that edge) 
		as well as a dictionary indexing even forms (included/excluded) by edges that are unique to that particular form.

		This information will be used to better inform the span_PSI in certain event types e.g. RI.  In the case of RI, false positives can be created by 
		(among other issues) changing SE events with a relatively high baseline intron coverage.  The emergence of a cassette exon above the background of intron
		looks (if considering only the junction involving the flanking exon and the cassette exon) exactly like an intron whose retention is decreasing.  Only by 
		considering the other junction (flanking exon to flanking exon) does it become clear that this is not the case.  In practical terms, this information will
		contribute to a greater span_PSI so these instances will be filtered. 

	'''

	edge_count_dict = {}
	edge_indexed_event_dict = {}

	for junction in junction_only_count_dict:

		junction_list = junction.split("_")

		edgel = "_".join(junction_list[0:2] + [junction_list[3]])
		edger = "_".join([junction_list[0]] + junction_list[2:4])

		edge_count_dict.setdefault(edgel, 0)
		edge_count_dict.setdefault(edger, 0)

		edge_indexed_event_dict.setdefault(edgel, []).extend(junction_indexed_event_dict[junction])
		edge_indexed_event_dict.setdefault(edger, []).extend(junction_indexed_event_dict[junction])

		edge_count_dict[edgel] += junction_only_count_dict[junction]
		edge_count_dict[edger] += junction_only_count_dict[junction]

	for edge in edge_indexed_event_dict:

		event_forms = set(edge_indexed_event_dict[edge])

		final_event_forms = []

		for event_form in event_forms:

			edge_coord = int(edge.split("_")[1])

			identifier = event_form.split("_")

			event_id = "_".join(identifier[0:-1])
			form = identifier[-1]

			if edge in standard_event_dict[event_id][form + "_junction_counts"]:

				standard_event_dict[event_id][form + "_junction_counts"][edge] = edge_count_dict[edge]



def max_jnc_gene_dict(event_ioe, standard_event_dict):
	'''
		Recovers gene-event associations from ioe file, then generates lists of all junction counts in a particular gene in a dictionary index by gene symbol.  Returns said dictionary.
	'''


	###Get gene-event associations from ioe file

	if event_ioe.endswith(".gz"):

		ioe = gzip.open(event_ioe, 'rb')

	else:

		ioe = open(event_ioe, 'r')

	ioe.readline()

	gene_event_dict = {}

	for line in ioe:

		entry = line.strip().split()
		chrom = entry[0].strip()
		gene = entry[1].strip()
		event = entry[2].strip().split(";")[1]

		gene_event_dict.setdefault(gene, {"chrom": [], "events": []})

		if event not in gene_event_dict[gene]["events"]:
			gene_event_dict[gene]["events"].append(event)

		if chrom not in gene_event_dict[gene]["chrom"]:
			gene_event_dict[gene]["chrom"].append(chrom)

	###Collect junction counts from event_dict

	gene_jc_dict = {}

	for gene in gene_event_dict:

		for event in gene_event_dict[gene]["events"]:

			assigned = False

			for c in gene_event_dict[gene]["chrom"]:

				for s in standard_event_dict[c]:

					if event in standard_event_dict[c][s]:

						chrom = c
						strand = s
						assigned = True

			if not assigned:

				sys.exit("Can't find event in event dict during max_jnc_gene_dict")

			if "gene" not in standard_event_dict[chrom][strand][event]:
				standard_event_dict[chrom][strand][event]["gene"] = [gene]
			else:
				if gene not in standard_event_dict[chrom][strand][event]["gene"]:
					standard_event_dict[chrom][strand][event]["gene"].append(gene)

			if gene not in gene_jc_dict:

				gene_jc_dict[gene] = []

			for form in ["included", "excluded"]:

				for junction in standard_event_dict[chrom][strand][event][form + "_junction_counts"]:

					gene_jc_dict[gene].append(standard_event_dict[chrom][strand][event][form + "_junction_counts"][junction])

	return gene_jc_dict





def calc_psi(standard_event_dict, outdir, sample_name, gzipped, gene_jc_dict, suppress_output, bootstrap_num = "NA", filename_addendum = "", file_write_mode = "w", header = True):

	'''
		Calculates PSI values for all events, writes tab-separated outfile containing:

		event_id	event_type	included_count	excluded_count	included_junction_count	excluded_junction_count	psi

		Note that included and excluded junction count are normalization factors (e.g. the included form of SE event (skipped exon) has 2 junctions while the excluded form has 1).  Normalization for retained introns is done internally so the factors for these events are just set to "1") - these factors are for use with rMATS-STAT.

	'''

	if not suppress_output:

		if gzipped:
			count_psi_outfile = gzip.open(outdir + "/" + sample_name + "_" + "count_psi_outfile" + filename_addendum + ".tsv.gz", file_write_mode + 'b')	
		else:
			count_psi_outfile = open(outdir + "/" + sample_name  + "_" + "count_psi_outfile" + filename_addendum + ".tsv", file_write_mode)

		if header:

			count_psi_outfile.write("sample_name" + "\t" + "event_id" + "\t" + "event_type" + "\t" + "min_IJC" + "\t" + "min_SJC" + "\t" + "IncFormLen" + "\t" + "SkipFormLen" + "\t" + "PSI" + "\t" + "max_gene_frac" + "\t" + "all_IJC" + "\t" + "all_SJC" + "\t" + "span_PSI" + "\t" + "PSI_lo" + "\t" + "PSI_hi" + "\t" + "bootstrap_num" + "\n")


	for chrom in standard_event_dict:
		for strand in standard_event_dict[chrom]:
			for event in standard_event_dict[chrom][strand]:

				included_counts = [standard_event_dict[chrom][strand][event]["included_junction_counts"][i] for i in standard_event_dict[chrom][strand][event]["included_junction_counts"]]
				min_included = min(included_counts)
				avg_included = float(sum(included_counts))/float(len(included_counts))

				excluded_counts = [standard_event_dict[chrom][strand][event]["excluded_junction_counts"][i] for i in standard_event_dict[chrom][strand][event]["excluded_junction_counts"]]
				min_excluded = min(excluded_counts)
				avg_excluded = float(sum(excluded_counts))/float(len(excluded_counts))

				if min_included + min_excluded > 0:

					min_psi = round(float(min_included)/(float(min_included) + float(min_excluded)),4)

					all_psi_values = []

					for i in included_counts:

						for j in excluded_counts:

							temp_PSI = float(i)/(float(i) + float(j))
							all_psi_values.append(temp_PSI)

					span_psi = round(max(all_psi_values) - min(all_psi_values), 4)
					psi_lo = round(min(all_psi_values),4)
					psi_hi = round(max(all_psi_values),4)

				else:

					min_psi = "NA"
					span_psi = "NA"
					psi_lo = "NA"
					psi_hi = "NA"


				if avg_included + avg_excluded > 0:

					avg_psi = float(avg_included)/(float(avg_included) + float(avg_excluded))

				else:
					avg_psi = "NA"

				if gene_jc_dict != None:
					if "gene" in standard_event_dict[chrom][strand][event]:

						genes_max_jn = []

						for gene in standard_event_dict[chrom][strand][event]["gene"]: ### looped in case event "belongs" to multiple genes

							genes_max_jn.append(max(gene_jc_dict[gene]))

						gene_max_jn = max(genes_max_jn)

						if gene_max_jn > 0:

							event_max_frac = round(max((float(min_included)/gene_max_jn), float(min_excluded)/gene_max_jn),4)

						else:

							event_max_frac = "NA"
					else:
						event_max_frac = "NA"
				else:
					event_max_frac = "NA"

				standard_event_dict[chrom][strand][event]["avg_psi"] = avg_psi
				standard_event_dict[chrom][strand][event]["min_psi"] = min_psi
				standard_event_dict[chrom][strand][event]["span_psi"] = span_psi
				standard_event_dict[chrom][strand][event]["psi_lo"] = psi_lo
				standard_event_dict[chrom][strand][event]["psi_hi"] = psi_hi
				standard_event_dict[chrom][strand][event]["event_max_frac"] = event_max_frac
				standard_event_dict[chrom][strand][event]["min_included"] = min_included
				standard_event_dict[chrom][strand][event]["min_excluded"] = min_excluded
				standard_event_dict[chrom][strand][event]["avg_included"] = avg_included
				standard_event_dict[chrom][strand][event]["avg_excluded"] = avg_excluded
				standard_event_dict[chrom][strand][event]["all_included_counts"] = ",".join(map(str, included_counts))
				standard_event_dict[chrom][strand][event]["all_excluded_counts"] = ",".join(map(str, excluded_counts))


				if not suppress_output:
					###below "1" is now used for included and excluded numbers of junctions. The normalization for junction number is not needed as only the minimum count value is used now.
					count_psi_outfile.write(sample_name + "\t" + event + "\t" + standard_event_dict[chrom][strand][event]["event_type"] + "\t" + str(min_included) + "\t" + str(min_excluded) + "\t" + "1" + "\t" + "1" + "\t" + str(avg_psi) + "\t" + str(event_max_frac) + "\t" + ",".join(map(str, included_counts)) + "\t" + ",".join(map(str, excluded_counts)) + "\t" + str(span_psi) + "\t" + str(psi_lo) + "\t" + str(psi_hi) + "\t" + str(bootstrap_num) + "\n")




def main(args, event_dict = None):

	start_time = time.time()

	print "{0}: {1} seconds elapsed. Starting junctionCounts.  Now parsing arguments.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	parser = argparse.ArgumentParser()
	parser.add_argument("--event_gtf", type = str, help = "GTF describing pairwise events")
	parser.add_argument("--bam", type = str, help = "BAM read file for counting junction reads", required = True)
	parser.add_argument("--se", action = "store_true", default = False, help = "BAM file is single-ended.  Default assumes paired-end reads.")
	parser.add_argument("--forward_read", type = str, help = "Specify read that matches the 'forward' strand.  Options are 'R1','R2' or 'unstranded' if the library is unstranded.  Unstranded use is not currently recommended. default = 'R2'", default = "R2")
	parser.add_argument("--outdir", type = str, help = "Path for output files", required = True)
	parser.add_argument("--sample_name", type = str, help = "Sample name", required = True)
	parser.add_argument("--dump_pkl_dict", action = "store_true", help = "Dumps a pkl file of the splice event dict with all quantifications")
	parser.add_argument("--gzipped", action = "store_true", help = "Output files will be gzipped if set")
	parser.add_argument("--event_ioe", type = str, help = "Event ioe file - used to recover event-gene association")
	parser.add_argument("--calc_gene_frac", action = "store_true", help = "Requires IOE file to be passed to --event_ioe.  If set, a maximal event fraction of gene expression will be estimated by max(min_excluded, min_included)/max(gene_junctions)")
	parser.add_argument("--suppress_output", action = "store_true", help = "Suppresses output files if set")
	parser.add_argument("--enable_edge_use", action = "store_false", default = True, help = "Use individual exon edges that are unique to form in quantification")
	parser.add_argument("--turn_off_no_ends", action = "store_false", default = True, help = "Disable the exclusion of transcript-termi from isoform-specific exon edge quantification (default is to exclude ends)")
	parser.add_argument("--suppress_eij_use", action = "store_true", default = False, help = "Don't use exon-overlapping exon-intron junctions for quantification except for RI events.")
	parser.add_argument("--disable_ri_extrapolation", action = "store_false", default = True, help = "Disable the use of RI/MR included forms to inform quantification of other events that contain these exons")
	parser.add_argument("--n_bootstraps", type = int, help = "Number of bootstraps. default = 0", default = 0)

	args = parser.parse_args(args)


	input_gtf = args.event_gtf
	bam_filename = args.bam
	se = args.se
	forward_read = args.forward_read
	output_directory = args.outdir
	sample_name = args.sample_name
	dump_pkl_dict = args.dump_pkl_dict
	gzipped = args.gzipped
	event_ioe = args.event_ioe
	calc_gene_frac = args.calc_gene_frac
	suppress_output = args.suppress_output
	n_bootstraps = args.n_bootstraps
	bootstraps = n_bootstraps > 0

	if forward_read not in ["R1","R2","unstranded"]:

		sys.exit("Only values R1, R2, or unstranded supported for --forward_read, not " + forward_read + ".")

	if se:

		forward_read = "R1"


	if event_dict is None and input_gtf is None:

		sys.exit("No input events! Must supply --event_gtf option OR pass event_dict argument to main function.  junctionCounts exiting  . . . ")

	elif event_dict is not None and input_gtf is not None:

		sys.exit("Both event_dict and event_gtf supplied - need to supply one or the other.  junctionCounts exiting  . . . ")



	if calc_gene_frac:

		if event_ioe == None:

			sys.exit("--calc_gene_frac set without passing IOE file to --event_ioe.  Exiting . . . ")

	subprocess.call("mkdir -p " + output_directory, shell = True)

	print "{0}: {1} seconds elapsed. Arguments parsed. Now importing event GTF file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	##Generate event dict

	if input_gtf is not None:
		
		standard_event_dict = splice_lib.generate_standard_event_dict_chrom_strand(input_gtf)

		splice_lib.complete_event_dict_chrom_strand(standard_event_dict, args.enable_edge_use, args.suppress_eij_use, args.turn_off_no_ends, args.disable_ri_extrapolation)

		if bootstraps:

			standard_event_dict_pristine = copy.deepcopy(standard_event_dict) ## maintains a copy with 0 counts for use in bootstrapping

	else:

		standard_event_dict = event_dict

	event_type_counts = splice_lib.assess_event_types_chrom_strand(standard_event_dict)

	#"Events indexed by junction. SE:nnn RI:nnn ..

	#SHOULD BE: 2017-12-11 13:12:24  12.5 seconds elapsed. Did this. Going to do that

	print "{0}: {1} seconds elapsed. Imported event dict.  Found {2} total events with {3} MS, {4} SE, {5} A3, {6} A5, {7} AF, {8} AL, {9} MF, {10} ML, {11} CF, {12} CL, {13} UF, {14} UL, {15} AT, {16} AP, {17} RI, {18} MR, {19} MX, {20} CO and {21} AB events.  Now indexing events by junctions . . . ".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)), str(event_type_counts["total"]), str(event_type_counts["MS"]), str(event_type_counts["SE"]), str(event_type_counts["A3"]), str(event_type_counts["A5"]), str(event_type_counts["AF"]), str(event_type_counts["AL"]), str(event_type_counts["MF"]), str(event_type_counts["ML"]), str(event_type_counts["CF"]), str(event_type_counts["CL"]), str(event_type_counts["UF"]), str(event_type_counts["UL"]), str(event_type_counts["AT"]), str(event_type_counts["AP"]),  str(event_type_counts["RI"]), str(event_type_counts["MR"]), str(event_type_counts["MX"]),  str(event_type_counts["CO"]), str(event_type_counts["AB"]))


	junction_indexed_event_dict = splice_lib.generate_junction_indexed_event_dict_chrom_strand(standard_event_dict)

	junction_only_count_dict = {}
	#junction_only_count_dict = {key:0 for key in junction_indexed_event_dict}

	if bootstraps:

		junction_only_count_dict_pristine = copy.deepcopy(junction_only_count_dict)


	ri_present = event_type_counts["RI"] > 0 or event_type_counts["MR"] > 0

	if not args.suppress_eij_use or ri_present:

		print "{0}: {1} seconds elapsed. Events indexed by junction. Now generating exon-intron junction-indexed event dict and nested containment lists for exon-intron junctions.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

		ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict = create_eij_ncls_dict(standard_event_dict)

		print "{0}: {1} seconds elapsed. Exon-intron junction-indexed event dict and nested containment lists built.  Counting splice and exon-intron junction reads for quantification. This may take a while".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	else:

		print "{0}: {1} seconds elapsed. Events indexed by junction.  Skipping exon-intron junction assessment as no intron retention events (RI or MR) events present and exon-intron junction use for other event types not requested. Now counting junction reads to quantify junction-containing isoforms.  This may take a while.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

		ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict = {},{},{}

	if bootstraps:

		eij_only_count_dict_pristine = copy.deepcopy(eij_only_count_dict)


	all_read_info = process_reads(bam_filename, junction_indexed_event_dict, junction_only_count_dict, standard_event_dict, ncls_by_chrom_strand, eij_indexed_event_dict, eij_only_count_dict, forward_read = forward_read, single_end = se, bootstraps = bootstraps)

	###get counts of each exon edge (i.e. the sum of all junctions involving that exon)

	if not args.enable_edge_use:
		get_exon_edge_counts(junction_only_count_dict, junction_indexed_event_dict, standard_event_dict)
	

	if calc_gene_frac:
		print "{0}: {1} seconds elapsed. Junction counting complete.  Now calculating gene fraction.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
		gene_jc_dict = max_jnc_gene_dict(event_ioe, standard_event_dict)
		print "{0}: {1} seconds elapsed. Gene fraction calculation complete.  Now calculating PSI values and writing output file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
	else:
		gene_jc_dict = None
		print "{0}: {1} seconds elapsed Junction counting complete.  Skipping gene fraction calculation. Calculating PSI values and writing output file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	calc_psi(standard_event_dict, output_directory, sample_name, gzipped, gene_jc_dict, suppress_output)

	print "{0}: {1} seconds elapsed. File output complete.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	if dump_pkl_dict:

		print "{0}: {1} seconds elapsed. Dumping event dict pkl file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
		pickle.dump(standard_event_dict, open(output_directory + "/" + sample_name + "_junctionCounts_standard_event_dict.pkl", 'wb'), pickle.HIGHEST_PROTOCOL)


	if bootstraps:

		print "{0}: {1} seconds elapsed. Beginning bootstrapping.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

		n_reads = len(all_read_info)

		for i in range(0, n_bootstraps):

			junction_only_count_dict_bootstrap = copy.deepcopy(junction_only_count_dict_pristine)
			standard_event_dict_bootstrap = copy.deepcopy(standard_event_dict_pristine)
			eij_only_count_dict_bootstrap = copy.deepcopy(eij_only_count_dict_pristine)

			bootstrap_junction_counts(all_read_info, junction_only_count_dict_bootstrap, standard_event_dict_bootstrap, eij_only_count_dict_bootstrap, junction_indexed_event_dict, eij_indexed_event_dict,  n_reads, forward_read)

			if not args.enable_edge_use:
				get_exon_edge_counts(junction_only_count_dict_bootstrap, junction_indexed_event_dict, standard_event_dict_bootstrap)

			calc_psi(standard_event_dict_bootstrap, output_directory, sample_name, gzipped, gene_jc_dict, suppress_output, bootstrap_num = i, filename_addendum = "_bootstraps", file_write_mode = "a", header = True if i == 0 else False)



	print "{0}: {1} seconds elapsed. junctionCounts complete.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
	print "Ceci n'est pas un algorithme bioinformatique."

	return standard_event_dict

if __name__ == 	'__main__':

	main(sys.argv[1:])

	#minimal_list = intersection_collapse([[1],[1,2,3],[2,3],[5,6],[1,7]])
	#print minimal_list

	# Ultimately will have something in the list like [[1,2,3],2,5], which should reduce to [2,5]

	# [[1,2,3],[2,3],5] reduces to [[2,3],5]