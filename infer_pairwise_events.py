import splice_lib
import sys
import argparse
import copy
import subprocess
import cPickle as pickle
import time
import re
from datetime import datetime


####TODO: improve by writing the following functions:
########## 1) Select flanking exon
########## 2) Add event to dict
########## 3)  


####TODO: Fix alt events


def filter_standard_transcript_dict(full_transcript_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length):
	'''
		Filters transcripts for nonsensical exons/introns. E.g. exons with start > end position, introns with length < 20, exons with length < 3.
	'''

	transcripts_to_delete = []
	
	for transcript in full_transcript_dict:

		i_should_continue = False

		if i_should_continue:

			continue

		for exon in full_transcript_dict[transcript]["exons"]:

			if (exon[1] - exon[0] + 1) < min_exon_length or (exon[1] - exon[0] + 1) > max_exon_length:

				transcripts_to_delete.append(transcript)
				i_should_continue = True
				break

		if i_should_continue:

			continue

		for junction in full_transcript_dict[transcript]["junctions"]:

			if (junction[1] - junction[0] - 1) < min_intron_length or (junction[1] - junction[0] - 1) > max_intron_length:

				transcripts_to_delete.append(transcript)
				break

	filtered_transcript_log = open(outdir + "/filtered_transcript_log.txt", 'w')
	for transcript in transcripts_to_delete:

		del full_transcript_dict[transcript]
		filtered_transcript_log.write(transcript + "\n")

	filtered_transcript_log.close()

def filter_event_dict(standard_event_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length):

	events_to_delete = []

	for event in standard_event_dict:

		i_should_continue = False

		if i_should_continue:

			continue

		for exon in standard_event_dict[event]["included_exons"] + standard_event_dict[event]["excluded_exons"]:

			if (int(exon[1]) - int(exon[0]) + 1) < min_exon_length or (int(exon[1]) - int(exon[0]) + 1) > max_exon_length:

				events_to_delete.append(event)
				i_should_continue = True
				break

		if i_should_continue:

			continue

		if standard_event_dict[event]["event_type"] == "RI":

			for junction in standard_event_dict[event]["excluded_junctions"]:

				if (int(junction[1]) - int(junction[0]) - 1) < min_intron_length or (int(junction[1]) - int(junction[0]) - 1) > max_intron_length: 

					events_to_delete.append(event)
					break

		else:
			for junction in standard_event_dict[event]["included_junctions"] + standard_event_dict[event]["excluded_junctions"]:

				if (int(junction[1]) - int(junction[0]) - 1) < min_intron_length or (int(junction[1]) - int(junction[0]) - 1) > max_intron_length:

					events_to_delete.append(event)
					break

	filtered_event_log = open(outdir + "/filtered_event_log.txt", 'w')
	for event in events_to_delete:

		#print standard_event_dict[event]
		del standard_event_dict[event]
		filtered_event_log.write(event + "\n")

	del events_to_delete


def cassette_events(standard_transcript_dict, standard_junction_indexed_transcript_dict, standard_donor_acceptor_indexed_transcript_dict):

	SE_counter = 0
	MSE_counter = 0

	standard_event_dict = {}

	for junction in standard_junction_indexed_transcript_dict:

		putative_excluded_form_transcripts = standard_junction_indexed_transcript_dict[junction]

		chrom = junction.split("_")[0]
		strand = junction.split("_")[3]
		left_junction = junction.split("_")[1]
		right_junction = junction.split("_")[2]

		left_junction_key = chrom + "_" + left_junction + "_" + strand
		right_junction_key = chrom + "_" + right_junction + "_" + strand

		if left_junction_key in standard_donor_acceptor_indexed_transcript_dict:

			left_junction_transcript_set = set(standard_donor_acceptor_indexed_transcript_dict[left_junction_key])

		if right_junction_key in standard_donor_acceptor_indexed_transcript_dict:

			right_junction_transcript_set = set(standard_donor_acceptor_indexed_transcript_dict[right_junction_key])

		putative_included_form_transcripts = list((left_junction_transcript_set & right_junction_transcript_set) - set(putative_excluded_form_transcripts))

		left_flanking_exons = [] ##Could be multiple alternative distal coordinates - must pick the shortest
		right_flanking_exons = []

		####Add all possible excluded form flanking exons to list
		for transcript in putative_excluded_form_transcripts:

			for exon in standard_transcript_dict[transcript]["exons"]:

				if int(left_junction) == int(exon[1]):

					if exon not in left_flanking_exons:

						left_flanking_exons.append(exon)

				elif int(right_junction) == int(exon[0]):

					if exon not in right_flanking_exons:
					
						right_flanking_exons.append(exon)

		if len(putative_included_form_transcripts) > 0:

			#print "Putative SE or MSE event found!"
			#print "Extracting the cassette exons . . . "

			cassette_exon_sets = [] #Each set corresponds to a different event - allows for MSE to overlap with SE etc

			cassette_set_transcript_dict = {}


			####Pull out all possible flanking exons, all possible (sets of) cassette exons
			for transcript in putative_included_form_transcripts:

				cassette_exons = []
				i_should_collect_cassette = False

				for exon in standard_transcript_dict[transcript]["exons"]:

					if int(exon[1]) == int(left_junction):

						if exon not in left_flanking_exons:
							left_flanking_exons.append(exon)

						i_should_collect_cassette = True 

						continue

					elif i_should_collect_cassette and int(exon[0]) != int(right_junction):

						cassette_exons.append(exon)

					elif i_should_collect_cassette and int(exon[0]) == int(right_junction):

						if exon not in right_flanking_exons:
							right_flanking_exons.append(exon)

						if cassette_exons not in cassette_exon_sets:
							cassette_exon_sets.append(cassette_exons)

						cassette_exons_key = "_".join(["_".join(map(str, i)) for i in cassette_exons])

						if cassette_exons_key not in cassette_set_transcript_dict:

							cassette_set_transcript_dict[cassette_exons_key] = []

						if transcript not in cassette_set_transcript_dict[cassette_exons_key]:
							cassette_set_transcript_dict[cassette_exons_key].append(transcript)

						break

				if len(left_flanking_exons) == 0 or len(right_flanking_exons) == 0 or len(cassette_exon_sets) == 0:

					#print "False positive - possible instance of donor + acceptor with same coordinates"
					continue


			#####Identify minimal flanking exons

			left_flanking_exon_lengths = []

			for exon in left_flanking_exons:

				left_flanking_exon_lengths.append(splice_lib.calc_length_exon_list([exon]))

			final_left_flanking_exon = left_flanking_exons[left_flanking_exon_lengths.index(min(left_flanking_exon_lengths))]

			right_flanking_exon_lengths = []

			for exon in right_flanking_exons:

				right_flanking_exon_lengths.append(splice_lib.calc_length_exon_list([exon]))

			final_right_flanking_exon = right_flanking_exons[right_flanking_exon_lengths.index(min(right_flanking_exon_lengths))]


			#print "There are (is?)", str(len(cassette_exon_sets)), "included forms associated with this excluded form"
			for cassette_set in cassette_exon_sets:

				if len(cassette_set) == 1:

					event_type = "SE"
					SE_counter += 1
					event_id = "SE." + str(SE_counter)

				elif len(cassette_set) > 1:

					event_type = "MSE"
					MSE_counter += 1
					event_id = "MSE." + str(MSE_counter)

				else:

					sys.exit("Cassette exon set of length 0 found in cassette_events")

				included_exons = [final_left_flanking_exon] + cassette_set + [final_right_flanking_exon]
				excluded_exons = [final_left_flanking_exon] + [final_right_flanking_exon]

				cassette_set_key = "_".join(["_".join(map(str, i)) for i in cassette_set])
				included_form_transcripts = cassette_set_transcript_dict[cassette_set_key]

				if event_id not in standard_event_dict:

					standard_event_dict[event_id] = {
							"event_type": event_type,
							"included_exons": copy.deepcopy(included_exons),
							"excluded_exons": copy.deepcopy(excluded_exons),
							"chrom": chrom,
							"strand": strand,
							"included_count": 0,
							"excluded_count": 0,
							"included_form_transcripts": copy.deepcopy(included_form_transcripts),
							"excluded_form_transcripts": copy.deepcopy(putative_excluded_form_transcripts)
						}
				else:
					sys.exit("Event id already found in dict - problem with naming strategy in cassette_events")


	return standard_event_dict




def alt_events(exon_pair_dict, standard_event_dict):

	A3_counter = 0
	A5_counter = 0

	for junction_that_is_missing in exon_pair_dict: ##Inner left missing or inner right missing

		which_junction = junction_that_is_missing.split("_")[1] ##right or left - combined with strand this identifies event as A3 or A5

		for key in exon_pair_dict[junction_that_is_missing]: 

			chrom = key.split("_")[0]
			strand = key.split("_")[4]

			forms_dict = {}  ###collect all unique forms

			for index,pair in enumerate(exon_pair_dict[junction_that_is_missing][key]["exon_pairs"]):

				if len(pair) > 2:

					sys.exit("Failure in generating exon pair dict - pair length > 2 found.")

				full_pair_key = chrom + "_" + "_".join(map(str, pair[0])) + "_" + "_".join(map(str, pair[1])) + "_" + strand

				transcript = exon_pair_dict[junction_that_is_missing][key]["transcripts"][index][:]

				if full_pair_key not in forms_dict:

					forms_dict[full_pair_key] = []

				if transcript not in forms_dict[full_pair_key]:

					forms_dict[full_pair_key].append(transcript)

			else:

				###Inner loop complete - all forms collected.  Infer events.

				if strand == "+" and which_junction == "left":

					event_type = "A5"

				elif strand == "+":

					event_type = "A3"

				elif strand == "-" and which_junction == "left":

					event_type = "A3"

				elif strand == "-":

					event_type = "A5"

				else:

					sys.exit("Error in alt_events - event type could not be classified. Exiting . . . ")

				forms_dict_keys = copy.deepcopy(forms_dict.keys())

				forms_pair_dict = {}


				for index, forms_key in enumerate(forms_dict_keys):

					if index < (len(forms_dict_keys) - 1):

						for sub_index in range(index + 1, len(forms_dict_keys)):

							forms_pair_key = forms_key + "_" + forms_dict_keys[sub_index]

							if forms_pair_key not in forms_pair_dict:

								forms_pair_dict[forms_pair_key] = {forms_key[:]: copy.deepcopy(forms_dict[forms_key]), forms_dict_keys[sub_index][:]: copy.deepcopy(forms_dict[forms_dict_keys[sub_index]])}

				else:

					###Identify included form, excluded form

					for event in forms_pair_dict:

						if len(forms_pair_dict[event]) != 2:

							sys.exit("Incorrect event size in form_pair_dict (alt_events). Exiting . . . ")


						event_form_keys = copy.deepcopy(forms_pair_dict[event].keys())

						if (int(event_form_keys[0].split("_")[3]) - int(event_form_keys[0].split("_")[2])) < (int(event_form_keys[1].split("_")[3]) - int(event_form_keys[1].split("_")[2])):

							included_form_index = 0
							excluded_form_index = 1

						elif (int(event_form_keys[0].split("_")[3]) - int(event_form_keys[0].split("_")[2])) > (int(event_form_keys[1].split("_")[3]) - int(event_form_keys[1].split("_")[2])):

							included_form_index = 1
							excluded_form_index = 0

						else:

							print event_form_keys
							print forms_pair_dict
							print event

							sys.exit("Included form seems to be same length as excluded form in alt_events. Exiting . . . ")

						if event_type == "A3":

							A3_counter += 1
							event_id = "A3." + str(A3_counter)

						elif event_type == "A5":

							A5_counter += 1
							event_id = "A5." + str(A5_counter)

						included_exons = [[int(event_form_keys[included_form_index].split("_")[1]), int(event_form_keys[included_form_index].split("_")[2])], [int(event_form_keys[included_form_index].split("_")[3]), int(event_form_keys[included_form_index].split("_")[4])]]

						included_form_transcripts = copy.deepcopy(forms_pair_dict[event][event_form_keys[included_form_index]])


						excluded_exons = [[int(event_form_keys[excluded_form_index].split("_")[1]), int(event_form_keys[excluded_form_index].split("_")[2])], [int(event_form_keys[excluded_form_index].split("_")[3]), int(event_form_keys[excluded_form_index].split("_")[4])]]

						excluded_form_transcripts = copy.deepcopy(forms_pair_dict[event][event_form_keys[excluded_form_index]])


						if event_id not in standard_event_dict:

							standard_event_dict[event_id] = {
									"event_type": event_type,
									"included_exons": copy.deepcopy(included_exons),
									"excluded_exons": copy.deepcopy(excluded_exons),
									"chrom": chrom[:],
									"strand": strand[:],
									"included_count": 0,
									"excluded_count": 0,
									"included_form_transcripts": copy.deepcopy(included_form_transcripts),
									"excluded_form_transcripts": copy.deepcopy(excluded_form_transcripts)
								}
						else:
							sys.exit("Event id already found in dict - problem with naming strategy in alt_events")



def write_intron_bedfile(standard_junction_indexed_transcript_dict, outdir):

	intron_bedfile = open(outdir + "/RI_introns.bed", 'w')
	intron_bedfile.truncate()

	for junction in standard_junction_indexed_transcript_dict:

		chrom = re.sub("&","_",junction.split("_")[0])

		strand = junction.split("_")[-1]



		putative_excluded_form_transcripts = ",".join(standard_junction_indexed_transcript_dict[junction])

		intron_start = str(int(junction.split("_")[1]) - 1) ####added - 1
		intron_end = str(int(junction.split("_")[2])) ###removed - 1

		if (int(intron_end) - int(intron_start)) > 19:

			intron_bedfile.write("\t".join([chrom, intron_start, intron_end, junction + "|" + putative_excluded_form_transcripts, "1000", strand]) + "\n")
		else:

			print "Warning: Intron between", junction, "is either too short, or does not have a positive, nonzero length. Skipping . . . "


def write_exon_bedfile(standard_transcript_dict, outdir):

	exon_bedfile = open(outdir + "/transcript_exon_bedfile.bed", 'w')
	exon_bedfile.truncate()

	for transcript in standard_transcript_dict:

		chrom = re.sub("&","_",standard_transcript_dict[transcript]["chrom"])
		strand = standard_transcript_dict[transcript]["strand"]

		for exon in standard_transcript_dict[transcript]["exons"]:

			start = str(exon[0] - 1)
			end = str(exon[1])

			exon_bedfile.write("\t".join([chrom, start, end, transcript, "1000", strand]) + "\n")

	exon_bedfile.close()



def call_bedtools_intersect(outdir):

	bedtools_command = "bedtools intersect -a " + outdir + "/RI_introns.bed -b " + outdir + "/transcript_exon_bedfile.bed -wa -wb -s -f 1.00 > " + outdir + "/RI_introns_exon_overlaps.bed"


	try:
		subprocess.check_output(bedtools_command, shell = True, stderr = subprocess.STDOUT)
	except subprocess.CalledProcessError as e:
		print e.output
		sys.exit("Bedtools intersect failed in 'call_bedtools_intersect'.  RI inference failed. Exiting . . . ")

	#subprocess.call(bedtools_command, shell = True)


def retained_intron_events(outdir, standard_transcript_dict, standard_event_dict):

	event_type = "RI"

	with open(outdir + "/RI_introns_exon_overlaps.bed", 'r') as file:

		putative_ri_events = {}

		for line in file:

			entry = line.strip().split()

			junction = entry[3].split("|")[0]
			chrom = junction.split("_")[0]
			strand = junction.split("_")[-1]
			full_junction = [int(junction.split("_")[1]), int(junction.split("_")[2])]

			putative_excluded_form_transcripts = entry[3].split("|")[1].split(",")

			putative_included_form_transcript = entry[9]

			putative_included_form_exon = [int(entry[7]) + 1, int(entry[8])] ####check + 1

			event_id = junction

			if event_id not in putative_ri_events: 

				putative_ri_events[event_id] = {"left_exons": [], "right_exons": [], "chrom": chrom, "strand": strand, "included_exons": [], "included_form_transcripts": [], "excluded_form_transcripts": copy.deepcopy(putative_excluded_form_transcripts), "junction": copy.deepcopy(full_junction)}

			if putative_included_form_transcript not in putative_ri_events[event_id]["included_form_transcripts"]:

				putative_ri_events[event_id]["included_form_transcripts"].append(putative_included_form_transcript)

			if putative_included_form_exon not in putative_ri_events[event_id]["included_exons"]:

				putative_ri_events[event_id]["included_exons"].append(copy.deepcopy(putative_included_form_exon))

			for transcript in putative_excluded_form_transcripts:

				for exon in standard_transcript_dict[transcript]["exons"]:

					if exon[1] == full_junction[0]:

						if exon not in putative_ri_events[event_id]["left_exons"]:

							putative_ri_events[event_id]["left_exons"].append(copy.deepcopy(exon))

					elif exon[0] == full_junction[1]:

						if exon not in putative_ri_events[event_id]["right_exons"]:

							putative_ri_events[event_id]["right_exons"].append(copy.deepcopy(exon))

	RI_counter = 0

	for event_id in putative_ri_events:

		RI_counter += 1
		new_id = "RI." + str(RI_counter)

		chrom = putative_ri_events[event_id]["chrom"]
		strand = putative_ri_events[event_id]["strand"]
		included_form_transcripts = putative_ri_events[event_id]["included_form_transcripts"]
		excluded_form_transcripts = putative_ri_events[event_id]["excluded_form_transcripts"]

		left_coord_list = []

		for exon in putative_ri_events[event_id]["left_exons"] + putative_ri_events[event_id]["included_exons"]:

			left_coord_list.append(exon[0])

		final_left_coord = max(left_coord_list)


		right_coord_list = []

		for exon in putative_ri_events[event_id]["right_exons"] + putative_ri_events[event_id]["included_exons"]:

			right_coord_list.append(exon[1])

		final_right_coord = min(right_coord_list)

		final_left_exon = [final_left_coord, putative_ri_events[event_id]["junction"][0]]
		final_right_exon = [putative_ri_events[event_id]["junction"][1], final_right_coord]

		excluded_exons = [final_left_exon, final_right_exon]

		included_exons = [[final_left_coord, final_right_coord]]

		if new_id not in standard_event_dict:

			standard_event_dict[new_id] = {
					"event_type": event_type,
					"included_exons": copy.deepcopy(included_exons),
					"excluded_exons": copy.deepcopy(excluded_exons),
					"chrom": chrom[:],
					"strand": strand[:],
					"included_count": 0,
					"excluded_count": 0,
					"included_form_transcripts": copy.deepcopy(included_form_transcripts),
					"excluded_form_transcripts": copy.deepcopy(excluded_form_transcripts)
				}



def af_al_events(standard_event_dict, first_donor_acceptor_left_dict, first_donor_acceptor_right_dict, standard_transcript_dict, standard_donor_acceptor_indexed_transcript_dict, standard_junction_indexed_transcript_dict):

	AF_counter = 0
	AL_counter = 0

	directions = ["left", "right"]

	for direction in directions:

		for half_junction in eval("first_donor_acceptor_" + direction + "_dict"):

			putative_forms = {}

			chrom = half_junction.split("_")[0]
			strand = half_junction.split("_")[-1]

			if strand == "+" and direction == "left":

				event_type = "AF"
				case = 1

			elif strand == "+":

				event_type = "AL"
				case = 2

			elif strand == "-" and direction == "left":

				event_type = "AL"
				case = 3

			elif strand == "-":

				event_type = "AF"
				case = 4

			else:

				sys.exit("Non-standard strand type in af_al_events.  Exiting . . . ")

			#####Extract all forms/transcripts where the junction is the first inner half junction

			for transcript in eval("first_donor_acceptor_" + direction + "_dict")[half_junction]:

				if direction == "left":

					first_junction_index = 0
					first_exon_index = 0
					second_exon_index = 1

				elif direction == "right":

					first_junction_index = -1
					first_exon_index = -2
					second_exon_index = -1

				first_junction = standard_transcript_dict[transcript]["junctions"][first_junction_index]


				exon_pair = copy.deepcopy([standard_transcript_dict[transcript]["exons"][first_exon_index], standard_transcript_dict[transcript]["exons"][second_exon_index]])

				junction_key = chrom + "_" + "_".join(map(str, first_junction)) + "_" + strand

				if junction_key not in putative_forms:

					putative_forms[junction_key] = {"transcripts": [], "exon_pairs": []}

				if transcript not in putative_forms[junction_key]["transcripts"]:

					putative_forms[junction_key]["transcripts"].append(transcript)

				if exon_pair not in putative_forms[junction_key]["exon_pairs"]:

					putative_forms[junction_key]["exon_pairs"].append(exon_pair)

			#####Extract all forms where the junction need not be the first inner half junction. This allows for scenarios where multiple leading exons can be absent in one of the forms

			if half_junction in standard_donor_acceptor_indexed_transcript_dict:

				for transcript in standard_donor_acceptor_indexed_transcript_dict[half_junction]:

					if (direction == "left" and (standard_transcript_dict[transcript]["flat_junctions"].index(int(half_junction.split("_")[1])) % 2) != 0) or (direction == "right" and (standard_transcript_dict[transcript]["flat_junctions"].index(int(half_junction.split("_")[1])) % 2) == 0):  ###The index must be odd (would be even if 1-based), otherwise it's role (either as donor or acceptor) is not the same

						if direction == "left":

							query_junction_index = (standard_transcript_dict[transcript]["flat_junctions"].index(int(half_junction.split("_")[1])) - 1)/2

						elif direction == "right":

							query_junction_index = (standard_transcript_dict[transcript]["flat_junctions"].index(int(half_junction.split("_")[1])))/2



						query_junction = standard_transcript_dict[transcript]["junctions"][query_junction_index]

						query_junction_key = chrom + "_" + "_".join(map(str, query_junction)) + "_" + strand

						if query_junction_key not in putative_forms: ###If it is already in, we don't want to add any new transcripts b/c these would be transcripts in which the junction is not the first junction

							transcript_list = copy.deepcopy(standard_junction_indexed_transcript_dict[query_junction_key])

							putative_forms[query_junction_key] = {"transcripts": transcript_list, "exon_pairs": []}

							for sub_transcript in transcript_list: ###sub transcript to distinguish from outer loop transcript

								sub_transcript_query_junction_index = standard_transcript_dict[sub_transcript]["junctions"].index(query_junction)

								exon_pair = [standard_transcript_dict[sub_transcript]["exons"][sub_transcript_query_junction_index], standard_transcript_dict[sub_transcript]["exons"][sub_transcript_query_junction_index + 1]]

								if exon_pair not in putative_forms[query_junction_key]["exon_pairs"]:

									putative_forms[query_junction_key]["exon_pairs"].append(exon_pair)

			#####Collapse exon pair lists

			for junction_key in putative_forms:

				left_exon_lengths = []

				left_exons = []

				for exon_pair in putative_forms[junction_key]["exon_pairs"]:
	 
					if exon_pair[0] not in left_exons:

						left_exons.append(copy.deepcopy(exon_pair[0]))

				for exon in left_exons:

					left_exon_lengths.append(splice_lib.calc_length_exon_list([exon]))

				final_left_exon = left_exons[left_exon_lengths.index(min(left_exon_lengths))]

				right_exon_lengths = []

				right_exons = []

				for exon_pair in putative_forms[junction_key]["exon_pairs"]:
	 
					if exon_pair[1] not in right_exons:

						right_exons.append(copy.deepcopy(exon_pair[1]))

				for exon in right_exons:

					right_exon_lengths.append(splice_lib.calc_length_exon_list([exon]))

				final_right_exon = right_exons[right_exon_lengths.index(min(right_exon_lengths))]

				putative_forms[junction_key]["final_exons"] = [final_left_exon, final_right_exon]

			###Infer pairwise events

			putative_forms_keys = copy.deepcopy(putative_forms.keys())


			for index, junction_key in enumerate(putative_forms_keys):
				       
				if index < len(putative_forms_keys) - 1:

					for sub_index in range(index + 1, len(putative_forms_keys)):

						form_one_exons = copy.deepcopy(putative_forms[junction_key]["final_exons"])

						form_two_exons = copy.deepcopy(putative_forms[putative_forms_keys[sub_index]]["final_exons"])

						if direction == "left":

							constitutive_exon_index = 1

						elif direction == "right":

							constitutive_exon_index = 0

						if splice_lib.calc_length_exon_list([form_one_exons[constitutive_exon_index]]) < splice_lib.calc_length_exon_list([form_two_exons[constitutive_exon_index]]):

							form_two_exons[constitutive_exon_index] = copy.deepcopy(form_one_exons[constitutive_exon_index]) ###this is not the alternative exon here

						elif splice_lib.calc_length_exon_list([form_one_exons[constitutive_exon_index]]) > splice_lib.calc_length_exon_list([form_two_exons[constitutive_exon_index]]):

							form_one_exons[constitutive_exon_index] = copy.deepcopy(form_two_exons[constitutive_exon_index])
						###distal form is included form, (greater inter exon distance)


						if (int(form_one_exons[1][0]) - int(form_one_exons[0][1])) > (int(form_two_exons[1][0]) - int(form_two_exons[0][1])):

							included_exons = copy.deepcopy(form_one_exons)
							excluded_exons = copy.deepcopy(form_two_exons)
							included_form_transcripts = copy.deepcopy(putative_forms[junction_key]["transcripts"])
							excluded_form_transcripts = copy.deepcopy(putative_forms[putative_forms_keys[sub_index]]["transcripts"])

						elif (int(form_one_exons[1][0]) - int(form_one_exons[0][1])) < (int(form_two_exons[1][0]) - int(form_two_exons[0][1])):

							excluded_exons = copy.deepcopy(form_one_exons)
							included_exons = copy.deepcopy(form_two_exons)
							excluded_form_transcripts = copy.deepcopy(putative_forms[junction_key]["transcripts"])
							included_form_transcripts = copy.deepcopy(putative_forms[putative_forms_keys[sub_index]]["transcripts"])

						else:

							print form_one_exons
							print form_two_exons
							print direction
							print putative_forms

							sys.exit("Included and excluded forms seem to have the same inter-exon distance in af_al_events.  Exiting  . . . ")


						add_event_to_dict = True

						if direction == "left":

							test_index = 0
							exon_test_index = 1

						elif direction == "right":

							test_index = -1
							exon_test_index = 0

						for excluded_form_transcript in excluded_form_transcripts:

							if excluded_exons[test_index][exon_test_index] != standard_transcript_dict[excluded_form_transcript]["exons"][test_index][exon_test_index]:

								add_event_to_dict = False
								break


						if add_event_to_dict:

							if event_type == "AF":

								AF_counter += 1

								event_id = "AF." + str(AF_counter)

							elif event_type == "AL":

								AL_counter += 1

								event_id = "AL." + str(AL_counter)

							else:

								sys.exit("Unsupported event type in af_al_events. Exiting . . . ")

							if event_id not in standard_event_dict:

									standard_event_dict[event_id] = {
											"event_type": event_type,
											"included_exons": copy.deepcopy(included_exons),
											"excluded_exons": copy.deepcopy(excluded_exons),
											"chrom": chrom[:],
											"strand": strand[:],
											"included_count": 0,
											"excluded_count": 0,
											"included_form_transcripts": copy.deepcopy(included_form_transcripts),
											"excluded_form_transcripts": copy.deepcopy(excluded_form_transcripts)
										}
							else:
								sys.exit("Event id already found in dict - problem with naming strategy in alt_events")
				       


def mx_events(standard_transcript_dict, standard_event_dict):

	MX_counter = 0
	event_type = "MX"

	flanking_exon_dict = {}  ##In this dict transcripts (and exons) will be indexed by the donor and acceptor (e.g. chrom_donor_acceptor_strand) of the flanking exons of exon triplets

	for transcript in standard_transcript_dict:

		for index, exon in enumerate(standard_transcript_dict[transcript]["exons"]):

			if index <= (len(standard_transcript_dict[transcript]["exons"]) - 3):

				exon_triplet = standard_transcript_dict[transcript]["exons"][index:(index + 3)]

				chrom = standard_transcript_dict[transcript]["chrom"][:]
				strand = standard_transcript_dict[transcript]["strand"][:]

				donor_acceptor_index = chrom + "_" + str(exon_triplet[0][1]) + "_" + str(exon_triplet[2][0]) + "_" + strand

				putative_alternative_exon_key = "_".join(map(str, exon_triplet[1]))

				if donor_acceptor_index not in flanking_exon_dict:

					flanking_exon_dict[donor_acceptor_index] = {"chrom": chrom, "strand": strand, "left_flanking_exons": [], "right_flanking_exons": [], "putative_alternative_exons": {}}

				if exon_triplet[0] not in flanking_exon_dict[donor_acceptor_index]["left_flanking_exons"]:

					flanking_exon_dict[donor_acceptor_index]["left_flanking_exons"].append(copy.deepcopy(exon_triplet[0]))

				if exon_triplet[2] not in flanking_exon_dict[donor_acceptor_index]["right_flanking_exons"]:

					flanking_exon_dict[donor_acceptor_index]["right_flanking_exons"].append(copy.deepcopy(exon_triplet[2]))

				if putative_alternative_exon_key not in flanking_exon_dict[donor_acceptor_index]["putative_alternative_exons"]:

					flanking_exon_dict[donor_acceptor_index]["putative_alternative_exons"][putative_alternative_exon_key] = []

				if transcript not in flanking_exon_dict[donor_acceptor_index]["putative_alternative_exons"][putative_alternative_exon_key]:
					flanking_exon_dict[donor_acceptor_index]["putative_alternative_exons"][putative_alternative_exon_key].append(transcript)

	for key in flanking_exon_dict:

		chrom = flanking_exon_dict[key]["chrom"][:]
		strand = flanking_exon_dict[key]["strand"][:]

		if len(flanking_exon_dict[key]["putative_alternative_exons"]) > 1:

			left_flanking_exons = copy.deepcopy(flanking_exon_dict[key]["left_flanking_exons"])
			left_flanking_exon_lengths = []

			for exon in left_flanking_exons:

				left_flanking_exon_lengths.append(splice_lib.calc_length_exon_list([exon]))

			final_left_flanking_exon = left_flanking_exons[left_flanking_exon_lengths.index(min(left_flanking_exon_lengths))]

			right_flanking_exons = copy.deepcopy(flanking_exon_dict[key]["right_flanking_exons"])
			right_flanking_exon_lengths = []

			for exon in right_flanking_exons:

				right_flanking_exon_lengths.append(splice_lib.calc_length_exon_list([exon]))

			final_right_flanking_exon = right_flanking_exons[right_flanking_exon_lengths.index(min(right_flanking_exon_lengths))]


			alternative_exon_keys = flanking_exon_dict[key]["putative_alternative_exons"].keys()

			for index, sub_key in enumerate(alternative_exon_keys):

				if index < len(alternative_exon_keys) - 1:

					for sub_index in range(index + 1, len(alternative_exon_keys)):

						form_one_exon = map(int, sub_key.split("_"))

						form_two_exon = map(int, alternative_exon_keys[sub_index].split("_"))

						real_event = False

						if form_one_exon[0] != form_two_exon[0] and form_one_exon[1] != form_two_exon[1]:

							if form_one_exon[0] < form_two_exon[0] and form_one_exon[1] < form_two_exon[1]:

								real_event = True
								MX_counter += 1
								event_id = "MX." + str(MX_counter)

								if strand == "+":

									included_form = "form_one"
									excluded_form = "form_two"

								elif strand == "-":

									included_form = "form_two"
									excluded_form = "form_one"


							elif form_one_exon[0] > form_two_exon[0] and form_one_exon[1] > form_two_exon[1]:

								real_event = True
								MX_counter += 1
								event_id = "MX." + str(MX_counter)

								if strand == "+":

									included_form = "form_two"
									excluded_form = "form_one"

								elif strand == "-":

									included_form = "form_one"
									excluded_form = "form_two"

							if real_event:

								included_exons = [copy.deepcopy(final_left_flanking_exon), copy.deepcopy(eval(included_form + "_exon")), copy.deepcopy(final_right_flanking_exon)]
								excluded_exons = [copy.deepcopy(final_left_flanking_exon), copy.deepcopy(eval(excluded_form + "_exon")), copy.deepcopy(final_right_flanking_exon)]

								if included_form == "form_one":
									included_form_transcripts = copy.deepcopy(flanking_exon_dict[key]["putative_alternative_exons"][sub_key])
									excluded_form_transcripts = copy.deepcopy(flanking_exon_dict[key]["putative_alternative_exons"][alternative_exon_keys[sub_index]])

								elif included_form == "form_two":
									excluded_form_transcripts = copy.deepcopy(flanking_exon_dict[key]["putative_alternative_exons"][sub_key])
									included_form_transcripts = copy.deepcopy(flanking_exon_dict[key]["putative_alternative_exons"][alternative_exon_keys[sub_index]])

								if event_id not in standard_event_dict:

									standard_event_dict[event_id] = {
											"event_type": event_type,
											"included_exons": copy.deepcopy(included_exons),
											"excluded_exons": copy.deepcopy(excluded_exons),
											"chrom": chrom[:],
											"strand": strand[:],
											"included_count": 0,
											"excluded_count": 0,
											"included_form_transcripts": copy.deepcopy(included_form_transcripts),
											"excluded_form_transcripts": copy.deepcopy(excluded_form_transcripts)
										}






def generate_exon_pair_dict(standard_transcript_dict):

	exon_pair_dict = {"inner_left_missing": {}, "inner_right_missing": {}}  ###For the following exon pair [[100,200],[400,500]] with chr1, strand +, inner left missing would be chr1_100_400_500_+

	for transcript in standard_transcript_dict:

		chrom = standard_transcript_dict[transcript]["chrom"][:]
		strand = standard_transcript_dict[transcript]["strand"][:]

		def add_inner_left_missing_to_dict(inner_left_missing_key, exon_pair):

			if inner_left_missing_key not in exon_pair_dict["inner_left_missing"]:

				exon_pair_dict["inner_left_missing"][inner_left_missing_key] = {"exon_pairs": [], "transcripts": []}

			exon_pair_dict["inner_left_missing"][inner_left_missing_key]["exon_pairs"].append(exon_pair)
			exon_pair_dict["inner_left_missing"][inner_left_missing_key]["transcripts"].append(transcript)

		def add_inner_right_missing_to_dict(inner_right_missing_key, exon_pair):

			if inner_right_missing_key not in exon_pair_dict["inner_right_missing"]:

				exon_pair_dict["inner_right_missing"][inner_right_missing_key] = {"exon_pairs": [], "transcripts": []}

			exon_pair_dict["inner_right_missing"][inner_right_missing_key]["transcripts"].append(transcript)
			exon_pair_dict["inner_right_missing"][inner_right_missing_key]["exon_pairs"].append(exon_pair)


		for index, exon in enumerate(standard_transcript_dict[transcript]["exons"]):

			if index != (len(standard_transcript_dict[transcript]["exons"]) - 1): ## Avoids list index out of range for final exon

				exon_pair = [exon, standard_transcript_dict[transcript]["exons"][index + 1]]
				inner_left_missing_key = chrom + "_" + str(exon_pair[0][0]) + "_" + str(exon_pair[1][0]) + "_" + str(exon_pair[1][1]) + "_" + strand
				inner_right_missing_key = chrom + "_" + str(exon_pair[0][0]) + "_" + str(exon_pair[0][1]) + "_" + str(exon_pair[1][1]) + "_" + strand

				if index == 0:

					add_inner_right_missing_to_dict(inner_right_missing_key, exon_pair)  ###only this one bc alternative would be AF/AL

				elif index == (len(standard_transcript_dict[transcript]["exons"]) - 2):

					add_inner_left_missing_to_dict(inner_left_missing_key, exon_pair)    ###Only this one bc alternative would be AF/AL

				else:

					add_inner_left_missing_to_dict(inner_left_missing_key, exon_pair)
					add_inner_right_missing_to_dict(inner_right_missing_key, exon_pair)


	return exon_pair_dict


def output_ioe(outdir, standard_event_dict, standard_transcript_dict):

	event_centric_ioe_file = open(outdir + "/splice_lib_events.ioe", 'w')

	event_centric_ioe_file.write("seqname" + "\t" + "gene_id" + "\t" + "event_id" + "\t" + "included_transcripts" + "\t" + "total_transcripts" + "\n")

	for event in standard_event_dict:

		chrom = standard_event_dict[event]["chrom"][:]
		gene = standard_transcript_dict[standard_event_dict[event]["included_form_transcripts"][0]]["gene"][:]
		included_transcripts = copy.deepcopy(standard_event_dict[event]["included_form_transcripts"])
		included_transcripts.sort()
		total_transcripts = copy.deepcopy(standard_event_dict[event]["included_form_transcripts"]) + copy.deepcopy(standard_event_dict[event]["excluded_form_transcripts"])
		total_transcripts.sort()

		event_centric_ioe_file.write(chrom + "\t" + gene + "\t" + gene + ";" + event + "\t" + ",".join(included_transcripts) + "\t" + ",".join(total_transcripts) + "\n")

	event_centric_ioe_file.close()


def main():

	start_time = time.time()

	print "{0}: {1} seconds elapsed. Starting infer_pairwise_events.  Now parsing arguments.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	parser = argparse.ArgumentParser()
	parser.add_argument("--transcript_gtf", type = str, help = "Full transcript gtf file", required = True)
	parser.add_argument("--outdir", type = str, help = "Path to output directory", required = True)
	parser.add_argument("--min_exon_length", type = int, help = "Minimum allowable exon length in input gtf.  Transcripts with shorter exons will be filtered. (default: 3)", default = 3)
	parser.add_argument("--min_intron_length", type = int, help = "Minimum allowable intron length in input gtf. Transcripts with shorter introns will be filtered. (default: 20)", default = 20)
	parser.add_argument("--max_exon_length", type = int, help = "Maximum allowable exon length in input gtf.  Transcripts with longer exons will be filtered (default = 35000)", default = 35000)
	parser.add_argument("--max_intron_length", type = int, help = "Maximum allowable intron length in input gtf.  Transcripts with longer introns will be fitlered (default = 1000000)", default = 1000000)
	parser.add_argument("--dump_pkl_file", action = "store_true", help = "If set, program will dump pickle file of event dict.")

	args = parser.parse_args()

	transcript_gtf = args.transcript_gtf
	outdir = args.outdir
	min_exon_length = args.min_exon_length
	min_intron_length = args.min_intron_length
	max_exon_length = args.max_exon_length
	max_intron_length = args.max_intron_length
	dump_pkl_file = args.dump_pkl_file

	subprocess.call("mkdir -p " + outdir, shell = True)

	print "{0}: {1} seconds elapsed. Arguments parsed.  Now importing transcriptome gtf.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
	
	standard_transcript_dict = splice_lib.generate_standard_transcript_dict(transcript_gtf)

	print "{0}: {1} seconds elapsed. Transcriptome gtf imported.  Now sorting transcript exons.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	splice_lib.sort_transcript_dict_exons(standard_transcript_dict)

	print "{0}: {1} seconds elapsed. Sorted transcript exons.  Now adding transcript junctions to transcript dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	splice_lib.add_junctions_to_transcript_dict(standard_transcript_dict)

	print "{0}: {1} seconds elapsed. Added transcript junctions.  Now filtering transcriptome to remove transcripts with short introns/exons.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	filter_standard_transcript_dict(standard_transcript_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length)

	print "{0}: {1} seconds elapsed. Filtered transcriptome.  Now creating junction-indexed transcript dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	standard_junction_indexed_transcript_dict = splice_lib.index_transcripts_by_junctions(standard_transcript_dict)

	print "{0}: {1} seconds elapsed. Created junction-indexed transcript dict.  Now creating assorted exon subsect dicts for event inference.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	donor_acceptor_indexed_set = splice_lib.index_transcripts_by_donor_acceptor(standard_transcript_dict)

	standard_donor_acceptor_indexed_transcript_dict = donor_acceptor_indexed_set[0]

	first_donor_acceptor_left_dict = donor_acceptor_indexed_set[1]

	first_donor_acceptor_right_dict = donor_acceptor_indexed_set[2]

	exon_pair_dict = generate_exon_pair_dict(standard_transcript_dict)

	print "{0}: {1} seconds elapsed. Created assorted exon subset dicts.  Now identifying SE and MSE events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	standard_event_dict = cassette_events(standard_transcript_dict, standard_junction_indexed_transcript_dict, standard_donor_acceptor_indexed_transcript_dict)

	print "{0}: {1} seconds elapsed. SE and MSE event identification complete.  Now identifying A3 and A5 events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	alt_events(exon_pair_dict, standard_event_dict)

	print "{0}: {1} seconds elapsed. A3 and A5 event identification complete.  Now identifying AF and AL events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	af_al_events(standard_event_dict, first_donor_acceptor_left_dict, first_donor_acceptor_right_dict, standard_transcript_dict, standard_donor_acceptor_indexed_transcript_dict, standard_junction_indexed_transcript_dict)

	print "{0}: {1} seconds elapsed. AF and AL event identification complete.  Now identifying RI events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	write_intron_bedfile(standard_junction_indexed_transcript_dict, outdir)

	write_exon_bedfile(standard_transcript_dict, outdir)

	call_bedtools_intersect(outdir)

	retained_intron_events(outdir, standard_transcript_dict, standard_event_dict)

	print "{0}: {1} seconds elapsed. RI event identification complete.  Now identifying MX events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
	
	mx_events(standard_transcript_dict, standard_event_dict)

	print "{0}: {1} seconds elapsed. MX event identification complete.  Now cleaning/filtering events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	splice_lib.collapse_redundant_junction_events(standard_event_dict)

	standard_event_dict = splice_lib.rename_events(standard_event_dict)

	junction_indexed_event_dict = splice_lib.generate_junction_indexed_event_dict(standard_event_dict)

	splice_lib.filter_overlapping_se_alt_donacc_events(standard_event_dict, junction_indexed_event_dict)

	filter_event_dict(standard_event_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length)
	
	print "{0}: {1} seconds elapsed. Event cleaning/filtering complete.  Now renaming events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	standard_event_dict = splice_lib.rename_events(standard_event_dict)

	print "{0}: {1} seconds elapsed. Events renamed.  Now outputting event gtf.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	splice_lib.output_event_gtf(standard_event_dict, outdir)

	print "{0}: {1} seconds elapsed. Event gtf output.  Now outputting MISO-style gff3.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	splice_lib.output_miso_event_gff3(standard_event_dict, outdir, name="splice_lib_events")

	if dump_pkl_file:
		print "{0}: {1} seconds elapsed. MISO_style event gff3 output.  Now dumping pickle file of event dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
		pickle.dump(standard_event_dict, open(outdir + "/splice_lib_events.pkl", 'wb' ), -1 )
		print "{0}: {1} seconds elapsed. Pickle file dump complete.  Now outputting IOE file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	else:
		print "{0}: {1} seconds elapsed. MISO_style event gff3 output.  Now outputting IOE file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	output_ioe(outdir, standard_event_dict, standard_transcript_dict)

	print "{0}: {1} seconds elapsed. IOE file output.  Tabulating event type counts.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

	event_type_counts = splice_lib.assess_event_types(standard_event_dict)

	print "{0}: {1} seconds elapsed. Event type counts tabulated.  Found {2} total events with {3} MSE, {4} SE, {5} A3, {6} A5, {7} AF, {8} AL, {9} RI events and {10} MX events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)), str(event_type_counts["total"]), str(event_type_counts["MSE"]), str(event_type_counts["SE"]), str(event_type_counts["A3"]), str(event_type_counts["A5"]), str(event_type_counts["AF"]), str(event_type_counts["AL"]), str(event_type_counts["RI"]), str(event_type_counts["MX"]))

	print "{0}: {1} seconds elapsed. infer_pairwise_events complete.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
	print "Ceci n'est pas un algorithme bioinformatique."


if __name__ == 	'__main__':

	main()

