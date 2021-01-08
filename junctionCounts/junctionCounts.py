import sys
import re
import argparse
from splice_lib import splice_lib
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

__version__ = "0.1.0"


def create_eij_ncls_dict(standard_event_dict):

	eij_by_chrom_strand = {}
	eij_indexed_event_dict = {}
	eij_only_count_dict = {}
	ncls_by_chrom_strand = {}

	for event, event_val in standard_event_dict.iteritems():

		strand = event_val["strand"]
		chrom = event_val["chrom"]

		for eij in event_val["included_ei_junctions"]:

			((eij_by_chrom_strand
				).setdefault(chrom, {})
					).setdefault(strand, set()).add(int(eij))

			eij_index = chrom + "_" + str(eij) + "_" + strand

			eij_indexed_event_dict.setdefault(
				eij_index, set()).add(event + "_included")

			eij_only_count_dict.setdefault(eij_index, 0)

		for eij in event_val["excluded_ei_junctions"]:

			((eij_by_chrom_strand
				).setdefault(chrom, {})
					).setdefault(strand, set()).add(int(eij))

			eij_index = chrom + "_" + str(eij) + "_" + strand

			eij_indexed_event_dict.setdefault(
				eij_index, set()).add(event + "_excluded")

			eij_only_count_dict.setdefault(eij_index, 0)


	for chrom, chrom_dict in eij_by_chrom_strand.iteritems():

		ncls_by_chrom_strand[chrom] = {}

		for strand, strand_dict in chrom_dict.iteritems():

			starts = np.array(list(strand_dict)) - 1 ## to make 0-based
			ends = starts
			ids = starts

			ncls_by_chrom_strand[chrom][strand] = NCLS(starts,ends,ids)


	return (ncls_by_chrom_strand,
		    eij_indexed_event_dict, 
		    eij_only_count_dict)



def read_pair_generator(bam):
	"""
    	Slight modification of function written by biostars user gizone1: https://www.biostars.org/p/306041/
	Generate read pairs in a BAM file or within a region string.
	Reads are added to read_dict until a pair is found.
	"""

	read_dict = defaultdict(lambda: [None, None])

	for read in bam.fetch():

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


def parse_reads(read_list, strand_list):
	'''
		Takes as input a read or read pair in a list of pysam AlignedSegments, 
		returns a dictionary containing a list of junctions 
		(where a junction is described by the chromosome, 
		the start/end position, and the strand, all joined as a string), 
		as well as (separately) the chrom and strand information.  
		Some code derived from find_introns() in pysam.
	'''

	chrom_list = list(set([i.reference_name for i in read_list]))


	if len(chrom_list) > 1: ##implies some sort of chimeric read - not currently handled

		return None

	else:

		chrom = re.sub("_","&",chrom_list[0].strip())


	exons = []
	junctions = []


	if read_list[0].is_reverse:

		strand = strand_list[0]

	else:

		strand = strand_list[1]


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
		
		outdict = {"junctions": list(set(junctions)), 
				   "exons": exons, 
				   "strand": strand, 
				   "chrom": chrom}

	else:

		outdict = {"junctions": [], 
				   "exons": exons, 
				   "strand": strand, 
				   "chrom": chrom}

	return outdict

				
def intersection_collapse(input_list):

	'''
		Takes as input a list of lists of event IDs matched by junctions in a read.  Attempts to eliminate ambiguous assignments by recursively identifying pairwise intersections of the sublists.
		Returns a final minimal list of lists of event IDs

	'''

	input_set = [set(i) for i in input_list]

	intersections = []

	#print input_set

	condition = True
	while condition:

		for i in range(0,len(input_set)-1):

			for j in range(i+1,len(input_set)):

				intersect = input_set[i] & input_set[j]

				if len(intersect) > 0:

					del input_set[i], input_set[j-1]

					intersections.append(intersect)
					break


		if intersections :

			input_set.extend(intersections)

			intersections = []

		else:
			minimal_list = [list(k) for k in input_set]
			condition = False

	return minimal_list



def process_reads(bam_filename, 
				  junction_indexed_event_dict, 
				  junction_only_count_dict, 
				  standard_event_dict, 
				  ncls_by_chrom_strand,
				  eij_indexed_event_dict, 
				  eij_only_count_dict, 
				  output_directory, 
				  sample_name, 
				  forward_read = "R2", 
				  single_end = False, 
				  bootstraps = False):

	'''
		Primary read assignment function - takes as input a bam file, 
		parses the read with "parse_junction_read", 
		collapses matching events with "intersection_collapse", 
		and adds 1 count to all surviving event matches in the event ID-indexed dictionary.
	'''

	if forward_read == "R1":

		strand_list = ["-", "+"]
		stranded_unstranded_fn = assign_reads_stranded

	elif forward_read == "R2":

		strand_list = ["+", "-"]
		stranded_unstranded_fn = assign_reads_stranded

	elif forward_read == "unstranded":

		strand_list = ["", ""]
		stranded_unstranded_fn = assign_reads_unstranded

	else:

		sys.exit("Improperly designated forward read in parse_read().  Exiting . . . ")	
	


	bam = pysam.AlignmentFile(bam_filename, 'rb')


	idxstats = pysam.idxstats(bam_filename).split("\n")

	read_count = 0

	for i in idxstats:

		entry = i.split()
		if len(entry) == 4:
			read_count += int(entry[2])

	size = read_count if single_end else read_count/2

	all_read_info = [None]*size



	if bootstraps:

		bootstrap_function = store_for_bootstrapping

	else:

		bootstrap_function = do_nothing

	list_index_counter = 0


	if single_end:

		for read in bam.fetch():

			read_properties = parse_reads([read], 
										  strand_list)

			if not read_properties:
				continue

			read_info = assign_reads(
				 read_properties, 
				 junction_indexed_event_dict, 
				 junction_only_count_dict, 
				 standard_event_dict,
				 ncls_by_chrom_strand,
				 eij_indexed_event_dict, 
				 eij_only_count_dict,
				 stranded_unstranded_fn = stranded_unstranded_fn, 
				 bootstrap_function = bootstrap_function)

			all_read_info[list_index_counter] = read_info
			list_index_counter += 1


	else:

		for read1, read2 in read_pair_generator(bam):
			
			if read1 is not None and read2 is not None:

				read_properties = parse_reads([read1, read2], strand_list)

				if not read_properties:
					continue

				read_info = assign_reads(read_properties, 
					 junction_indexed_event_dict, 
					 junction_only_count_dict, 
					 standard_event_dict, 
					 ncls_by_chrom_strand,
					 eij_indexed_event_dict, 
					 eij_only_count_dict,
					 stranded_unstranded_fn = stranded_unstranded_fn, 
					 bootstrap_function = bootstrap_function)

				all_read_info[list_index_counter] = read_info
				list_index_counter += 1


	bam.close()

	return all_read_info, size




def append_matching(
		candidate_isoforms,
		matching_events, 
		junction_dict, 
		junction):		

	if matching_events not in candidate_isoforms:

		candidate_isoforms.append(matching_events)	

	for matching_event in matching_events:

		junction_dict.setdefault(matching_event, set()).add(junction)



def reduce_candidates(
		candidate_isoforms,
		event_junction_dict,
		event_eij_dict):

	minimal_candidate_isoform_list = intersection_collapse(candidate_isoforms)

	#print minimal_candidate_isoform_list
	#print event_junction_dict
	#print event_eij_dict

	minimal_candidate_isoform_list_flat = [i for j in minimal_candidate_isoform_list for i in j]

	## discard events from event_junction_dict and event_eij_dict that are eliminated by intersection_collapse

	event_junction_dict = {k:v for k,v in event_junction_dict.items() 
						   if k in minimal_candidate_isoform_list_flat}


	event_eij_dict = {k:v for k,v in event_eij_dict.items() 
					  if k in minimal_candidate_isoform_list_flat}

	return event_junction_dict, event_eij_dict


def store_for_bootstrapping(
		event_junction_dict,
		event_eij_dict,
		chrom,
		strand):

	event_junction_dict_list = "=".join([event + ":" + ",".join(junctions) 
												 for event, junctions in 
												 event_junction_dict.items()])

	event_eij_dict_list = "=".join([event + ":" + ",".join(eij) 
		 							for event, eij in 
		 							event_eij_dict.items()])

	#read_info = {"minimal_candidate_isoform_list": minimal_candidate_isoform_list,
	#			 "junctions": event_junction_dict,
	#			 "eij": event_eij_dict,
	#			 "chrom": chrom,
	#			 "strand": strand,
	#			 "possible_strands": possible_strands}

	return "&".join([event_junction_dict_list, 
					 event_eij_dict_list, 
					 chrom, 
					 strand])


def do_nothing(*args):

	pass



def assign_reads(
			 read_properties, 
			 junction_indexed_event_dict, 
			 junction_only_count_dict, 
			 standard_event_dict, 
			 ncls_by_chrom_strand,
			 eij_indexed_event_dict, 
			 eij_only_count_dict,
			 stranded_unstranded_fn, 
			 bootstrap_function):

	exons = read_properties["exons"]
	junctions = read_properties["junctions"]
	strand = read_properties["strand"]
	chrom = read_properties["chrom"]

	overlapping_eij = set()

	candidate_isoforms = []

	event_junction_dict = {}
	event_eij_dict = {}
	possible_strands = set() 
	### for use when forward_read == "unstranded"

	stranded_unstranded_fn(
		exons,
		junctions,
		strand,
		chrom,
		overlapping_eij,
		candidate_isoforms,
		ncls_by_chrom_strand,
		event_junction_dict,
		event_eij_dict,
		possible_strands,
		standard_event_dict,
		junction_indexed_event_dict,
		eij_indexed_event_dict,
		junction_only_count_dict,
		eij_only_count_dict)

	return bootstrap_function(
				event_junction_dict,
				event_eij_dict,
				chrom,
				strand)



def assign_reads_stranded(
		exons,
		junctions,
		strand,
		chrom,
		overlapping_eij,
		candidate_isoforms,
		ncls_by_chrom_strand,
		event_junction_dict,
		event_eij_dict,
		possible_strands,
		standard_event_dict,
		junction_indexed_event_dict,
		eij_indexed_event_dict,
		junction_only_count_dict,
		eij_only_count_dict):

	for junction in junctions:

		if junction in junction_indexed_event_dict:

			#junction_only_count_dict[chrom][strand][junction] += 1

			matching_events = junction_indexed_event_dict[junction][:]

			append_matching(
				candidate_isoforms,
				matching_events, 
				event_junction_dict, 
				junction)

	if chrom in ncls_by_chrom_strand:

		if strand in ncls_by_chrom_strand[chrom]:

			for exon in exons:

				[overlapping_eij.add(chrom + "_" + str(i[0] + 1) + "_" + strand) 
				 for i in 
				 ncls_by_chrom_strand[chrom][strand].find_overlap(exon[0], exon[1])]


	for eij in overlapping_eij:

		if eij in eij_indexed_event_dict:

			#eij_only_count_dict[eij] += 1

			matching_events = list(eij_indexed_event_dict[eij])

			append_matching(
				candidate_isoforms,
				matching_events, 
				event_eij_dict, 
				eij)


	if len(candidate_isoforms) > 0:

		event_junction_dict, event_eij_dict = reduce_candidates(
												candidate_isoforms,
												event_junction_dict,
												event_eij_dict)


	for j in event_junction_dict:

		identifier = j.split("_")
		event_id = "_".join(identifier[0:-1])
		event_form = identifier[-1]

		for junction in event_junction_dict[j]:

			standard_event_dict[event_id][event_form + "_junction_counts"][junction] += 1

	for j in event_eij_dict:

		identifier = j.split("_")
		event_id = "_".join(identifier[0:-1])
		event_form = identifier[-1]

		for eij in event_eij_dict[j]:

			standard_event_dict[event_id][event_form + "_junction_counts"][eij] += 1



def assign_reads_unstranded(
		exons,
		junctions,
		strand,
		chrom,
		overlapping_eij,
		candidate_isoforms,
		ncls_by_chrom_strand,
		event_junction_dict,
		event_eij_dict,
		possible_strands,
		standard_event_dict,
		junction_indexed_event_dict,
		eij_indexed_event_dict,
		junction_only_count_dict,
		eij_only_count_dict):

	junctions = [i + "+" for i in junctions] + [i + "-" for i in junctions]

	for junction in junctions:

		matching_events = junction_indexed_event_dict.get(junction)

		if matching_events:

			local_strand = i.split("_")[-1]

			possible_strands.add(local_strand)

			#junction_only_count_dict[chrom][strand][junction] += 1

			append_matching(
				matching_events, 
				event_junction_dict, 
				junction)

	for exon in exons:

		ncls_by_chrom_strand_chrom_entry = ncls_by_chrom_strand.get(chrom)

		if ncls_by_chrom_strand_chrom_entry:

			for local_strand in ncls_by_chrom_strand_chrom_entry:

				[overlapping_eij.add(chrom + "_" + str(i[0] + 1) + "_" + local_strand) 
				 for i in 
				 ncls_by_chrom_strand[chrom][local_strand].find_overlap(exon[0], exon[1])]


	for eij in overlapping_eij:

		matching_events_set = eij_indexed_event_dict.get(eij)

		if matching_events_set:

			test_strand = eij.split("_")[-1]

			possible_strands.add(test_strand)

			#eij_only_count_dict[eij] += 1
			matching_events = list(matching_events_set)

			append_matching(
				matching_events, 
				event_eij_dict, 
				eij)


	if len(candidate_isoforms) > 0:

		event_junction_dict, event_eij_dict = reduce_candidates(
												candidate_isoforms,
												event_junction_dict,
												event_eij_dict)

		for j in event_junction_dict:

			for junction in event_junction_dict[j]:

				standard_event_dict[event_id][event_form + "_junction_counts"][junction] += 1


		for j in event_eij_dict:

			for eij in event_eij_dict[j]:

				standard_event_dict[event_id][event_form + "_junction_counts"][eij] += 1



def bootstrap_junction_counts(junction_only_count_dict, 
							  standard_event_dict, 
						  	  eij_only_count_dict, 
						  	  junction_indexed_event_dict, 
						  	  eij_indexed_event_dict, 
						  	  n_reads, 
						  	  all_read_info):


	for n in range(0,int(n_reads)):

		index = randint(0,n_reads - 1)

		read = all_read_info[index].split("&")

		event_junction_dict_list = [[i.split(":")[0], i.split(":")[1].split(",")] 
		 							for i in 
		 							read[0].split("=")] if read[0] != "" else []

		event_eij_dict_list = [[i.split(":")[0], i.split(":")[1].split(",")] 
							   for i in 
							   read[1].split("=")] if read[1] != "" else []

		chrom = read[2]
		strand = read[3]

		for j in event_junction_dict_list:

			identifier = j[0].split("_")				
			event_id = "_".join(identifier[0:-1])
			event_form = identifier[-1]

			for junction in j[1]:

				standard_event_dict[event_id][event_form + "_junction_counts"][junction] += 1


		for j in event_eij_dict_list:

			identifier = j[0].split("_")				
			event_id = "_".join(identifier[0:-1])
			event_form = identifier[-1]

			for eij in j[1]:

				standard_event_dict[event_id][event_form + "_junction_counts"][eij] += 1



def get_exon_edge_counts(junction_only_count_dict, 
						 junction_indexed_event_dict, 
					 	 standard_event_dict):

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

	for junction, junction_val in junction_only_count_dict.iteritems():

		junction_list = junction.split("_")

		edgel = "_".join(junction_list[0:2] + [junction_list[3]])
		edger = "_".join([junction_list[0]] + junction_list[2:4])

		edge_count_dict.setdefault(edgel, 0)
		edge_count_dict.setdefault(edger, 0)

		edge_indexed_event_dict.setdefault(edgel, []).extend(junction_indexed_event_dict[junction])
		edge_indexed_event_dict.setdefault(edger, []).extend(junction_indexed_event_dict[junction])

		edge_count_dict[edgel] += junction_val
		edge_count_dict[edger] += junction_val

	for edge, edge_val in edge_indexed_event_dict.iteritems():

		event_forms = set(edge_val)

		final_event_forms = []

		for event_form in event_forms:

			edge_coord = int(edge.split("_")[1])

			identifier = event_form.split("_")

			event_id = "_".join(identifier[0:-1])
			form = identifier[-1]

			if edge in standard_event_dict[event_id][form + "_junction_counts"]:

				standard_event_dict[event_id][form + "_junction_counts"][edge] = edge_count_dict[edge]



def max_jnc_gene_dict(event_ioe, 
					  standard_event_dict):
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

		gene_event_dict.setdefault(gene, []).append(event)

	###Collect junction counts from event_dict

	gene_jc_dict = {}

	for gene, event_list in gene_event_dict.iteritems():

		for event in event_list:

			event_entry = standard_event_dict.get(event)

			if event_entry:

				event_entry.setdefault("gene", set()).add(gene)

				gene_jc_dict.setdefault(gene, [])

				gene_jc_entry = gene_jc_dict[gene]

				for form in ["included", "excluded"]:

					for junc_val in event_entry[form + "_junction_counts"].itervalues():

						gene_jc_entry.append(junc_val)

	ioe.close()

	return gene_jc_dict





def calc_psi(standard_event_dict, 
			 outdir, 
			 sample_name, 
			 gzipped, 
			 gene_jc_dict, 
			 suppress_output, 
			 bootstrap_num = "NA", 
			 filename_addendum = "", 
			 file_write_mode = "w", 
			 header = True):

	'''
		Calculates PSI values for all events, writes tab-separated outfile containing:

		event_id	event_type	included_count	excluded_count	included_junction_count	excluded_junction_count	psi

		Note that included and excluded junction count are normalization factors (e.g. the included form of SE event (skipped exon) has 2 junctions while the excluded form has 1).  Normalization for retained introns is done internally so the factors for these events are just set to "1") - these factors are for use with rMATS-STAT.

	'''

	if not suppress_output:

		if gzipped:
			count_psi_outfile = gzip.open(outdir + 
											"/" + 
											sample_name + 
											"_" + 
											"count_psi_outfile" + 
											filename_addendum + 
											".tsv.gz", 
											file_write_mode + 
											'b')	

		else:
			count_psi_outfile = open(outdir + 
									"/" + 
									sample_name  + 
									"_" + 
									"count_psi_outfile" + 
									filename_addendum + 
									".tsv", 
									file_write_mode)

		if header:

			header_string = "\t".join(["sample_name", 
									   "event_id", 
									   "event_type", 
									   "min_ijc", 
									   "min_ejc", 
									   "avg_psi", 
									   "max_gene_frac", 
									   "all_ijc", 
									   "all_ejc",
									   "avg_ijc",
									   "avg_ejc", 
									   "span_psi", 
									   "min_psi",
									   "ijc_min_psi",
									   "ejc_min_psi",
									   "max_psi", 
									   "ijc_max_psi",
									   "ejc_max_psi",
									   "mid_psi",
									   "bootstrap_num"])

			count_psi_outfile.write(header_string + "\n")


	for event, event_entry in standard_event_dict.iteritems():

		included_counts = [ event_entry["included_junction_counts"][i] 
							for i in 
							event_entry["included_junction_counts"] ]


		min_included = min(included_counts)

		avg_included = float(sum(included_counts))/float(len(included_counts))

		excluded_counts = [ event_entry["excluded_junction_counts"][i] 
							for i in 
							event_entry["excluded_junction_counts"] ]


		min_excluded = min(excluded_counts)
		avg_excluded = float(sum(excluded_counts))/float(len(excluded_counts))

		if min_included + min_excluded > 0:

			psi_min_counts = float(min_included)/(float(min_included) + float(min_excluded))

			all_psi_values = []
			all_psi_values_count_pairs = []

			for i in included_counts:

				for j in excluded_counts:

					temp_PSI = float(i)/(float(i) + float(j))
					all_psi_values.append(temp_PSI)
					all_psi_values_count_pairs.append((i,j))

			psi_span = max(all_psi_values) - min(all_psi_values)

			psi_avg = (sum(all_psi_values)/len(all_psi_values))

			psi_lo_full = min(all_psi_values)
			psi_lo = psi_lo_full
			psi_lo_count_pair = all_psi_values_count_pairs[ all_psi_values.index(psi_lo_full) ]
			psi_lo_inc_counts = psi_lo_count_pair[0]
			psi_lo_exc_counts = psi_lo_count_pair[1]

			psi_hi_full = max(all_psi_values)
			psi_hi = psi_hi_full
			psi_hi_count_pair = all_psi_values_count_pairs[ all_psi_values.index(psi_hi_full) ]
			psi_hi_inc_counts = psi_hi_count_pair[0]
			psi_hi_exc_counts = psi_hi_count_pair[1]

			psi_mid = (psi_hi_full + psi_lo_full)/2

		else:

			psi_avg = "NA"
			psi_min_counts = "NA"
			psi_span = "NA"
			psi_lo = "NA"
			psi_lo_inc_counts = "NA"
			psi_lo_exc_counts = "NA"
			psi_hi = "NA"
			psi_hi_inc_counts = "NA"
			psi_hi_exc_counts = "NA"
			psi_mid = "NA"


		if avg_included + avg_excluded > 0:

			psi_avg_counts = float(avg_included)/(float(avg_included) + float(avg_excluded))

		else:
			psi_avg_counts = "NA"

		if gene_jc_dict is not None:
			if "gene" in event_entry:

				genes_max_jn = []

				for gene in event_entry["gene"]: ### looped in case event "belongs" to multiple genes

					genes_max_jn.append(max(gene_jc_dict[gene]))

				gene_max_jn = max(genes_max_jn)

				if gene_max_jn > 0:

					event_max_frac = max(
										(float(min_included)/gene_max_jn), 
										(float(min_excluded)/gene_max_jn)
										)

				else:

					event_max_frac = "NA"
			else:
				event_max_frac = "NA"
		else:
			event_max_frac = "NA"

		event_entry["avg_counts_psi"] = psi_avg_counts
		event_entry["psi_avg"] = psi_avg
		event_entry["psi_min_counts"] = psi_min_counts
		event_entry["psi_span"] = psi_span
		event_entry["psi_lo"] = psi_lo
		event_entry["psi_lo_inc_counts"] = psi_lo_inc_counts
		event_entry["psi_lo_exc_counts"] = psi_lo_exc_counts
		event_entry["psi_hi_inc_counts"] = psi_hi_inc_counts
		event_entry["psi_hi_exc_counts"] = psi_hi_exc_counts
		event_entry["psi_hi"] = psi_hi
		event_entry["psi_mid"] = psi_mid
		event_entry["event_max_frac"] = event_max_frac
		event_entry["min_included"] = min_included
		event_entry["min_excluded"] = min_excluded
		event_entry["avg_included"] = avg_included
		event_entry["avg_excluded"] = avg_excluded
		event_entry["all_included_counts"] = ",".join(map(str, included_counts))
		event_entry["all_excluded_counts"] = ",".join(map(str, excluded_counts))


		if not suppress_output:

			###below "1" is now used for included and excluded numbers of junctions. 
			###The normalization for junction number is not needed as only a 
			###single count value is used for numerator and denominator now.

			output_entry = "\t".join([sample_name,
									  event,
									  event_entry["event_type"],
									  str(min_included),
									  str(min_excluded),
									  "%.4f" % psi_avg if psi_avg != "NA" else "NA",
									  "%.4f" % event_max_frac if event_max_frac != "NA" else "NA",
									  ",".join(map(str, included_counts)),
									  ",".join(map(str, excluded_counts)),
									  "%.1f" % avg_included,
									  "%.1f" % avg_excluded,
									  "%.4f" % psi_span if psi_span != "NA" else "NA",
									  "%.4f" % psi_lo if psi_lo != "NA" else "NA",
									  str(psi_lo_inc_counts),
									  str(psi_lo_exc_counts),
									  "%.4f" % psi_hi if psi_hi != "NA" else "NA",
									  str(psi_hi_inc_counts),
									  str(psi_hi_exc_counts),
									  "%.4f" % psi_mid if psi_mid != "NA" else "NA",
									  str(bootstrap_num)])

			count_psi_outfile.write(output_entry + "\n")




def main(args, event_dict = None):

	start_time = time.time()

	print ("{0}: {1} seconds elapsed. Starting junctionCounts." +
		   "Now parsing arguments."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		   			str(round(time.time() - start_time, 1)))

	parser = argparse.ArgumentParser()
	parser.add_argument("--event_gtf", 
						type = str, 
						help = "GTF describing pairwise events")

	parser.add_argument("--bam", 
						type = str, 
						help = "BAM read file for counting junction reads", 
						required = True)

	parser.add_argument("--se", 
						action = "store_true", 
						default = False, 
						help = "BAM file is single-ended.  Default assumes paired-end reads.")

	parser.add_argument("--forward_read", 
						type = str, 
						help = "Specify read that matches the 'forward' strand. " +
						       "Options are 'R1','R2' or 'unstranded' if the " +
						       "library is unstranded.  Unstranded use is not " +
						       "currently recommended. default = 'R2'", 
						default = "R2")

	parser.add_argument("--outdir", 
						type = str, 
						help = "Path for output files", 
						required = True)

	parser.add_argument("--sample_name", 
						type = str, 
						help = "Sample name", 
						required = True)

	parser.add_argument("--dump_pkl_dict", 
						action = "store_true", 
						help = "Dumps a pkl file of the splice event dict with " +
						 	   "all quantifications")

	parser.add_argument("--gzipped", 
						action = "store_true",
						help = "Output files will be gzipped if set")

	parser.add_argument("--event_ioe", 
						type = str, 
						help = "Event ioe file - used to recover event-gene association")

	parser.add_argument("--calc_gene_frac", 
						action = "store_true", 
						help = "Requires IOE file to be passed to --event_ioe. " + 
						       "If set, a maximal event fraction of gene expression " + 
						       "will be estimated by max(min_excluded, " + 
						       "min_included)/max(gene_junctions)")

	parser.add_argument("--suppress_output", 
						action = "store_true", 
						help = "Suppresses output files if set")

	parser.add_argument("--enable_edge_use", 
						action = "store_false", 
						default = True, 
						help = "Use individual exon edges that are unique to form " + 
							   "in quantification")

	parser.add_argument("--turn_off_no_ends", 
						action = "store_false", 
						default = True, 
						help = "Disable the exclusion of transcript-termini from " +
							   "isoform-specific exon edge quantification " +
							   "(default is to exclude ends)")

	parser.add_argument("--suppress_eij_use", 
						action = "store_true", 
						default = False, 
						help = "Don't use exon-overlapping exon-intron junctions " + 
							   "for quantification except for RI events.")

	parser.add_argument("--disable_ri_extrapolation", 
						action = "store_false", 
						default = True, 
						help = "Disable the use of RI/MR included forms to inform " + 
						       "quantification of other events that contain these exons")

	parser.add_argument("--n_bootstraps", 
						type = int, 
						help = "Number of bootstraps. default = 0", 
						default = 0)

	parser.add_argument("--version",
						action = "version",
						version = "junctionCounts.py is part of junctionCounts version " + __version__)

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

		sys.exit("Only values R1, R2, or unstranded supported for " + 
				 "--forward_read, not " + 
				 forward_read + 
				 ".")

	if se:

		forward_read = "R1"


	if event_dict is None and input_gtf is None:

		sys.exit("No input events! Must supply --event_gtf option " + 
				 "OR pass event_dict argument to main function. " + 
				 "junctionCounts exiting  . . . ")

	elif event_dict is not None and input_gtf is not None:

		sys.exit("Both event_dict and event_gtf supplied - " + 
				 "need to supply one or the other.  junctionCounts exiting  . . . ")



	if calc_gene_frac:

		if event_ioe is None:

			sys.exit("--calc_gene_frac set without passing IOE file to " + 
				     "--event_ioe.  Exiting . . . ")

	subprocess.call("mkdir -p " + output_directory, shell = True)

	print ("{0}: {1} seconds elapsed. Arguments parsed." + 
		   "Now importing event GTF file."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		            str(round(time.time() - start_time, 1)))

	##Generate event dict

	if input_gtf is not None:
		
		standard_event_dict = splice_lib.generate_standard_event_dict(input_gtf)

		splice_lib.complete_event_dict(standard_event_dict, 
										args.enable_edge_use, 
										args.suppress_eij_use, 
										args.turn_off_no_ends, 
										args.disable_ri_extrapolation)

		if bootstraps:

			standard_event_dict_pristine = copy.deepcopy(standard_event_dict) ## maintains a copy with 0 counts for use in bootstrapping


	else: 

		standard_event_dict = event_dict


	event_type_counts = splice_lib.assess_event_types(standard_event_dict)

	#"Events indexed by junction. SE:nnn RI:nnn ..

	#SHOULD BE: 2017-12-11 13:12:24  12.5 seconds elapsed. Did this. Going to do that

	print ("{0}: {1} seconds elapsed. Imported event dict. " + 
		   "Found {2} total events with {3} MS, {4} SE, " + 
		   "{5} A3, {6} A5, {7} AF, {8} AL, {9} MF, {10} ML, " + 
		   "{11} CF, {12} CL, {13} UF, {14} UL, {15} AT, " + 
		   "{16} AP, {17} RI, {18} MR, {19} MX, {20} CO and " + 
		   "{21} AB events.  Now indexing events by junctions . . . "
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		            str(round(time.time() - start_time, 1)), 
		            str(event_type_counts["total"]), 
		            str(event_type_counts["MS"]), 
		            str(event_type_counts["SE"]), 
		            str(event_type_counts["A3"]), 
		            str(event_type_counts["A5"]), 
		            str(event_type_counts["AF"]), 
		            str(event_type_counts["AL"]), 
		            str(event_type_counts["MF"]), 
		            str(event_type_counts["ML"]), 
		            str(event_type_counts["CF"]), 
		            str(event_type_counts["CL"]), 
		            str(event_type_counts["UF"]), 
		            str(event_type_counts["UL"]), 
		            str(event_type_counts["AT"]), 
		            str(event_type_counts["AP"]),  
		            str(event_type_counts["RI"]), 
		            str(event_type_counts["MR"]), 
		            str(event_type_counts["MX"]),  
		            str(event_type_counts["CO"]), 
		            str(event_type_counts["AB"]))


	junction_indexed_event_dict = splice_lib.generate_junction_indexed_event_dict(standard_event_dict)

	junction_only_count_dict = {}
	#junction_only_count_dict = {key:0 for key in junction_indexed_event_dict}

	if bootstraps:

		junction_only_count_dict_pristine = copy.deepcopy(junction_only_count_dict)


	ri_present = event_type_counts["RI"] > 0 or event_type_counts["MR"] > 0

	if not args.suppress_eij_use or ri_present:

		print ("{0}: {1} seconds elapsed. Events indexed by junction. " +
			   "Now generating exon-intron junction-indexed event dict " + 
			   "and nested containment lists for exon-intron junctions."
			   ).format(str(datetime.now().replace(microsecond = 0)), 
			            str(round(time.time() - start_time, 1)))

		(ncls_by_chrom_strand, 
		 eij_indexed_event_dict, 
		 eij_only_count_dict) = create_eij_ncls_dict(standard_event_dict)

		print ("{0}: {1} seconds elapsed. Exon-intron " +
		       "junction-indexed event dict and nested containment " + 
		       "lists built.  Counting splice and exon-intron junction " +
		       "reads for quantification. This may take a while"
		       ).format(str(datetime.now().replace(microsecond = 0)), 
		                str(round(time.time() - start_time, 1)))

	else:

		print ("{0}: {1} seconds elapsed. Events indexed by junction. " +
		       "Skipping exon-intron junction assessment as no intron " + 
		       "retention events (RI or MR) events present and " + 
		       "exon-intron junction use for other event types not " + 
		       "requested. Now counting junction reads to quantify " + 
		       "junction-containing isoforms.  This may take a while."
		       ).format(str(datetime.now().replace(microsecond = 0)), 
		                str(round(time.time() - start_time, 1)))

		(ncls_by_chrom_strand, 
		 eij_indexed_event_dict, 
		 eij_only_count_dict) = {},{},{}

	if bootstraps:

		eij_only_count_dict_pristine = copy.deepcopy(eij_only_count_dict)


	all_read_info, n_reads = process_reads(bam_filename, 
										   junction_indexed_event_dict, 
										   junction_only_count_dict, 
										   standard_event_dict, 
										   ncls_by_chrom_strand, 
										   eij_indexed_event_dict, 
										   eij_only_count_dict, 
										   output_directory, 
										   sample_name, 
										   forward_read = forward_read, 
										   single_end = se, 
										   bootstraps = bootstraps)

	###get counts of each exon edge (i.e. the sum of all junctions involving that exon)

	if not args.enable_edge_use:
		get_exon_edge_counts(junction_only_count_dict, 
							 junction_indexed_event_dict, 
							 standard_event_dict)
	

	if calc_gene_frac:

		print ("{0}: {1} seconds elapsed. Junction counting complete. " + 
			   "Now calculating gene fraction."
			   ).format(str(datetime.now().replace(microsecond = 0)), 
			            str(round(time.time() - start_time, 1)))

		gene_jc_dict = max_jnc_gene_dict(event_ioe, standard_event_dict)


		print ("{0}: {1} seconds elapsed. Gene fraction calculation " +
			   "complete.  Now calculating PSI values and writing " + 
			   "output file."
			   ).format(str(datetime.now().replace(microsecond = 0)), 
			            str(round(time.time() - start_time, 1)))


	else:

		gene_jc_dict = None

		print ("{0}: {1} seconds elapsed Junction counting complete. " + 
			   "Skipping gene fraction calculation. Calculating PSI " + 
			   "values and writing output file."
			   ).format(str(datetime.now().replace(microsecond = 0)), 
			            str(round(time.time() - start_time, 1)))

	calc_psi(standard_event_dict, 
			 output_directory, 
			 sample_name, 
			 gzipped, 
			 gene_jc_dict, 
			 suppress_output)

	print ("{0}: {1} seconds elapsed. File output complete."
		  ).format(str(datetime.now().replace(microsecond = 0)), 
		           str(round(time.time() - start_time, 1)))

	if dump_pkl_dict:

		print ("{0}: {1} seconds elapsed. " + 
			   "Dumping event dict pkl file."
			   ).format(str(datetime.now().replace(microsecond = 0)), 
			            str(round(time.time() - start_time, 1)))


		pickle.dump(standard_event_dict, 
					open(output_directory + 
						 "/" + 
						 sample_name + 
						 "_junctionCounts_standard_event_dict.pkl", 
						 'wb'), 
					pickle.HIGHEST_PROTOCOL)


	if bootstraps:


		print ("{0}: {1} seconds elapsed. Beginning bootstrapping."
			).format(str(datetime.now().replace(microsecond = 0)), 
			         str(round(time.time() - start_time, 1)))


		for i in range(0, n_bootstraps):

			print ("{0}: {1} seconds elapsed. " + 
				   "Beginning bootstrapping round {2}."
				   ).format(str(datetime.now().replace(microsecond = 0)), 
				     str(round(time.time() - start_time, 1)), 
				     str(i))

			junction_only_count_dict_bootstrap = copy.deepcopy(junction_only_count_dict_pristine)
			standard_event_dict_bootstrap = copy.deepcopy(standard_event_dict_pristine)
			eij_only_count_dict_bootstrap = copy.deepcopy(eij_only_count_dict_pristine)

			bootstrap_junction_counts(junction_only_count_dict_bootstrap, 
									  standard_event_dict_bootstrap, 
									  eij_only_count_dict_bootstrap, 
									  junction_indexed_event_dict, 
									  eij_indexed_event_dict,  
									  n_reads, 
									  all_read_info)


			if not args.enable_edge_use:
				get_exon_edge_counts(junction_only_count_dict_bootstrap, 
									 junction_indexed_event_dict, 
									 standard_event_dict_bootstrap)

			calc_psi(standard_event_dict_bootstrap, 
					 output_directory, 
					 sample_name, gzipped, 
					 gene_jc_dict, 
					 suppress_output, 
					 bootstrap_num = i, 
					 filename_addendum = "_bootstraps", 
					 file_write_mode = "a", 
					 header = True if i == 0 else False)

			print ("{0}: {1} seconds elapsed. Bootstrapping round {2} complete."
				   ).format(str(datetime.now().replace(microsecond = 0)), 
				            str(round(time.time() - start_time, 1)), str(i))

		print ("{0}: {1} seconds elapsed. Bootstrapping complete."
			   ).format(str(datetime.now().replace(microsecond = 0)), 
			            str(round(time.time() - start_time, 1)))


	print ("{0}: {1} seconds elapsed. junctionCounts complete."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		            str(round(time.time() - start_time, 1)))

	print "Ceci n'est pas un algorithme bioinformatique."

	return standard_event_dict

if __name__ == 	'__main__':

	main(sys.argv[1:])

	#minimal_list = intersection_collapse([[1],[1,2,3],[2,3],[5,6],[1,7]])
	#print minimal_list

	# Ultimately will have something in the list like [[1,2,3],2,5], which should reduce to [2,5]

	# [[1,2,3],[2,3],5] reduces to [[2,3],5]
