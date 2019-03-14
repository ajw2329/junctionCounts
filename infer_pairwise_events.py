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

def filter_event_dict(standard_event_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length, min_AP_AT_dist = 100):
 
    ''' input: standard event dict, constraints on exon, intron length
        output: modifies standard event dict (deletes items that fail the filtration criteria)

        Designed to deal with suspiciously short/long exons (e.g. 1 nt, 1000000 exons created by stringtie created in the former case by unclear circumstances, in the latter case by running stringtie on ribosome-depleted whole cell RNA-seq data)
    '''
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

        elif standard_event_dict[event]["event_type"] in ["AT", "AP"]:

            if splice_lib.calc_length_exon_list(standard_event_dict[event]["included_exons"]) - splice_lib.calc_length_exon_list(standard_event_dict[event]["excluded_exons"]) < min_AP_AT_dist:

                events_to_delete.append(event)

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

    filtered_event_log.close()

    del events_to_delete



def node_dir(standard_transcript_dict):
    '''
        returns node-indexed transcript dictionary (i.e. dictionary of transcripts indexed by individual exon boundaries)
    '''

    node_index_tx = {}

    def add_key(key, value):

        if key not in node_index_tx:

            node_index_tx[key] = [value]

        else:

            if value not in node_index_tx[key]:

                node_index_tx[key].append(value)


    for transcript in standard_transcript_dict:

        chrom = standard_transcript_dict[transcript]["chrom"][:]
        strand = standard_transcript_dict[transcript]["strand"][:]

        for exon in standard_transcript_dict[transcript]["exons"]:

            if strand == "+":

                upstream = 0
                downstream = 1

            elif strand == "-":

                upstream = 1
                downstream = 0

            add_key("_".join([chrom, str(exon[upstream]), strand, "upstream"]), transcript)
            add_key("_".join([chrom, str(exon[downstream]), strand, "downstream"]), transcript)

    return node_index_tx


def compare_transcripts(tx1, tx2, tx1_exons, tx2_exons, strand, chrom):

    events = {}

    if len(tx1_exons) < len(tx2_exons):

        temp_var = tx1
        tx1 = tx2
        tx2 = temp_var

        temp_exons = copy.deepcopy(tx1_exons)
        tx1_exons = copy.deepcopy(tx2_exons)
        tx2_exons = copy.deepcopy(temp_exons)

    else:

        tx1_exons = copy.deepcopy(tx1_exons)
        tx2_exons = copy.deepcopy(tx2_exons)

    def convert_to_nodes(exons):

        nodes = []
        node_exon_dict = {}

        for exon in exons:

            node_exon_dict[str(exon[0]) + "_left"] = copy.deepcopy(exon)
            node_exon_dict[str(exon[1]) + "_right"] = copy.deepcopy(exon)

            nodes.append(str(exon[0]) + "_left")
            nodes.append(str(exon[1]) + "_right")

        nodes = ["universal_left"] + nodes + ["universal_right"]
        node_exon_dict["universal_left"] = ["universal_left"]
        node_exon_dict["universal_right"] = ["universal_right"]

        return nodes, node_exon_dict

    tx1_nodes, tx1_node_exon_dict = convert_to_nodes(tx1_exons)
    tx2_nodes, tx2_node_exon_dict = convert_to_nodes(tx2_exons)


    if len(tx1_nodes) > len(tx2_nodes):

        tx2_nodes[len(tx2_nodes):len(tx1_nodes)] = (len(tx1_nodes) - len(tx2_nodes))*["filler"]

    elif len(tx1_nodes) < len(tx2_nodes):

        tx1_nodes[len(tx1_nodes):len(tx2_nodes)] = (len(tx2_nodes) - len(tx1_nodes))*["filler"]    

    tx1_index = 0
    tx2_index = 0

    def scan_for_events(tx1_index, tx2_index):

        pre_common = []
        post_common = []
        tx1_unique = []
        tx2_unique = []

        tx2_adjustment = tx1_index - tx2_index

        for i in range(tx1_index, len(tx1_nodes)):

            if i - tx2_adjustment >= len(tx2_nodes):

                return

            if tx2_nodes[i - tx2_adjustment] == tx1_nodes[i] and len(tx1_unique) == 0 and len(tx2_unique) == 0:

                pre_common.append(tx1_nodes[i])

            else:

                if tx1_nodes[i] in tx2_nodes:


                    post_common.append(tx1_nodes[i])

                    tx2_unique += tx2_nodes[tx2_nodes.index(pre_common[-1]) + 1: tx2_nodes.index(tx1_nodes[i])]

                    pre_common = [pre_common[-1]]
                    post_common = [post_common[0]]


                    form_1_nodes = pre_common + tx1_unique + post_common
                    form_2_nodes = pre_common + tx2_unique + post_common

                    form_1_exons = []
                    form_2_exons = []

                    #form_1_exons[0] = tx1_node_exon_dict[node]

                    for node in form_1_nodes:

                        if tx1_node_exon_dict[node] not in form_1_exons:

                            form_1_exons.append(copy.deepcopy(tx1_node_exon_dict[node]))

                    for node in form_2_nodes:

                        if tx2_node_exon_dict[node] not in form_2_exons:

                            form_2_exons.append(copy.deepcopy(tx2_node_exon_dict[node]))


                    #classify event type as well as included/excluded form
                    classified = classify_event(pre_common, post_common, tx1_unique, tx2_unique, form_1_exons, form_2_exons, strand, tx1, tx2)

                    key_string = chrom + ":" + strand + ":" + pre_common[-1] + "|" + "&".join(classified["included_unique"]) + "|" + "&".join(classified["excluded_unique"]) + "|" + post_common[0]


                    events[key_string] = {"included_exons": copy.deepcopy(classified["included_exons"]), "excluded_exons": copy.deepcopy(classified["excluded_exons"]), "event_type": copy.deepcopy(classified["event_type"]), "included_form_transcripts": [copy.deepcopy(classified["included_form_transcript"])], "excluded_form_transcripts": [copy.deepcopy(classified["excluded_form_transcript"])], "strand": strand, "chrom": chrom, "included_count": 0, "excluded_count": 0, "sources": [], "pre_common": copy.deepcopy(pre_common), "post_common": copy.deepcopy(post_common)}


                    if tx1_nodes.index(post_common[-1]) < len(tx1_nodes) and tx2_nodes.index(post_common[-1]) < len(tx2_nodes):
                        scan_for_events(tx1_nodes.index(post_common[-1]), tx2_nodes.index(post_common[-1]))


                    break

                else:

                    tx1_unique.append(tx1_nodes[i])

    scan_for_events(tx1_index, tx2_index)


    return events


def call_compare_transcripts(node_index_tx, standard_transcript_dict):

    standard_event_dict = {}
    compared_transcripts = {}

    for node in node_index_tx:

        for i, transcript in enumerate(node_index_tx[node]):

            if i < len(node_index_tx[node]) - 1:

                for j in range(i+1, len(node_index_tx[node])):

                    other_transcript = node_index_tx[node][j][:]

                    transcript_pair = [transcript, other_transcript]
                    transcript_pair.sort()

                    if "_".join(transcript_pair) not in compared_transcripts:

                        comparison_event_dict = compare_transcripts(transcript, other_transcript, standard_transcript_dict[transcript]["exons"], standard_transcript_dict[other_transcript]["exons"], standard_transcript_dict[transcript]["strand"], standard_transcript_dict[transcript]["chrom"])

                        compared_transcripts["_".join(transcript_pair)] = "NA"

                        for key in comparison_event_dict:

                            if key not in standard_event_dict:

                                standard_event_dict[key] = copy.deepcopy(comparison_event_dict[key])

                            else:

                                ###Combine transcript lists
                                ###Shorten outer exons if possible

                                if "universal_left" in comparison_event_dict[key]["pre_common"]:

                                    min_right = min(comparison_event_dict[key]["included_exons"][-1][-1], comparison_event_dict[key]["excluded_exons"][-1][-1], standard_event_dict[key]["included_exons"][-1][-1], standard_event_dict[key]["excluded_exons"][-1][-1])
                                    standard_event_dict[key]["included_exons"][-1][-1] = min_right
                                    standard_event_dict[key]["excluded_exons"][-1][-1] = min_right

                                elif "universal_right" in comparison_event_dict[key]["post_common"]:

                                    max_left = max(comparison_event_dict[key]["included_exons"][0][0], comparison_event_dict[key]["excluded_exons"][0][0], standard_event_dict[key]["included_exons"][0][0], standard_event_dict[key]["excluded_exons"][0][0])
                                    standard_event_dict[key]["included_exons"][0][0] = max_left
                                    standard_event_dict[key]["excluded_exons"][0][0] = max_left

                                else:

                                    max_left = max(comparison_event_dict[key]["included_exons"][0][0], comparison_event_dict[key]["excluded_exons"][0][0], standard_event_dict[key]["included_exons"][0][0], standard_event_dict[key]["excluded_exons"][0][0])
                                    standard_event_dict[key]["included_exons"][0][0] = max_left
                                    standard_event_dict[key]["excluded_exons"][0][0] = max_left

                                    min_right = min(comparison_event_dict[key]["included_exons"][-1][-1], comparison_event_dict[key]["excluded_exons"][-1][-1], standard_event_dict[key]["included_exons"][-1][-1], standard_event_dict[key]["excluded_exons"][-1][-1])
                                    standard_event_dict[key]["included_exons"][-1][-1] = min_right
                                    standard_event_dict[key]["excluded_exons"][-1][-1] = min_right

                                standard_event_dict[key]["included_form_transcripts"] = list(set(standard_event_dict[key]["included_form_transcripts"] + comparison_event_dict[key]["included_form_transcripts"]))

                                standard_event_dict[key]["excluded_form_transcripts"] = list(set(standard_event_dict[key]["excluded_form_transcripts"] + comparison_event_dict[key]["excluded_form_transcripts"]))


    return standard_event_dict



def classify_event(pre_common, post_common, tx1_unique, tx2_unique, form_1_exons, form_2_exons, strand, tx1, tx2):

    event_type = "WHOOPS"

    if "universal_left" in pre_common:

        form_1_exons = copy.deepcopy(form_1_exons)[1:]
        form_2_exons = copy.deepcopy(form_2_exons)[1:]

        min_right = min(form_1_exons[-1][-1], form_2_exons[-1][-1])

        form_1_exons[-1][-1] = min_right
        form_2_exons[-1][-1] = min_right

    elif "universal_right" in post_common:

        form_1_exons = copy.deepcopy(form_1_exons)[0:-1]
        form_2_exons = copy.deepcopy(form_2_exons)[0:-1]

        max_left = max(form_1_exons[0][0], form_2_exons[0][0])

        form_1_exons[0][0] = max_left
        form_2_exons[0][0] = max_left

    else:

        if form_1_exons[0][0] != form_2_exons[0][0]:

            max_left = max(form_1_exons[0][0], form_2_exons[0][0])

            form_1_exons[0][0] = max_left
            form_2_exons[0][0] = max_left

        if form_1_exons[-1][-1] != form_2_exons[-1][-1]:

            min_right = min(form_1_exons[-1][-1], form_2_exons[-1][-1])

            form_1_exons[-1][-1] = min_right
            form_2_exons[-1][-1] = min_right

    if "universal_left" not in pre_common:

        if "universal_right" not in post_common:

            if not (len(form_1_exons) == len(form_2_exons) and len(form_1_exons) == 3 and len(tx1_unique) == 2 and len(tx2_unique) == 2):

                if splice_lib.calc_length_exon_list(form_1_exons) > splice_lib.calc_length_exon_list(form_2_exons):

                    included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2

                elif splice_lib.calc_length_exon_list(form_1_exons) < splice_lib.calc_length_exon_list(form_2_exons):

                    included_exons = form_2_exons; excluded_exons = form_1_exons; included_unique = tx2_unique; excluded_unique = tx1_unique; included_form_transcript = tx2; excluded_form_transcript = tx1

                elif len(form_1_exons) > len(form_2_exons):

                    included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2

                elif len(form_1_exons) < len(form_2_exons):

                    included_exons = form_2_exons; excluded_exons = form_1_exons; included_unique = tx2_unique; excluded_unique = tx1_unique; included_form_transcript = tx2; excluded_form_transcript = tx1

                else: ##In this case, wherein both the lengths and exon counts of both forms are equal, I arbitrarily assign included form to form 1

                    included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2

                if len(included_exons) == 3 and len(excluded_exons) == 2 and len(included_unique) == 2 and len(excluded_unique) == 0:

                    event_type = "SE"

                elif len(included_exons) > 3 and len(excluded_exons) == 2 and len(included_unique) > 2 and (len(included_unique) % 2) == 0 and len(excluded_unique) == 0:

                    event_type = "MS"

                elif len(included_exons) == 2 and len(excluded_exons) == 2 and len(included_unique) == 1 and len(excluded_unique) == 1:

                    if int(included_unique[0].split("_")[0]) == included_exons[0][1]:

                        if strand == "+":
                            event_type = "A5"
                        elif strand == "-":
                            event_type = "A3"
                    elif int(included_unique[0].split("_")[0]) == included_exons[1][0]:
                        if strand == "+":
                            event_type = "A3"
                        elif strand == "-":
                            event_type = "A5"

                elif len(included_exons) == 1 and len(excluded_exons) == 2:

                    event_type = "RI"

                elif len(included_exons) == 1 and len(excluded_exons) > 2:

                    event_type = "MR"

                else:
                    event_type = "CO" ## for complex

            else:
                event_type = "MX"

                if strand == "+":
                    if int(tx1_unique[0].split("_")[0]) < int(tx2_unique[0].split("_")[0]):
                        included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2
                    else:
                        included_exons = form_2_exons; excluded_exons = form_1_exons; included_unique = tx2_unique; excluded_unique = tx1_unique; included_form_transcript = tx2; excluded_form_transcript = tx1

                elif strand == "-":
                    if int(tx1_unique[-1].split("_")[0]) > int(tx2_unique[-1].split("_")[0]):
                        included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2
                    else:
                        included_exons = form_2_exons; excluded_exons = form_1_exons; included_unique = tx2_unique; excluded_unique = tx1_unique; included_form_transcript = tx2; excluded_form_transcript = tx1

        else:

            if form_1_exons[-1][-1] > form_2_exons[-1][-1]:
                included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2
            else:
                included_exons = form_2_exons; excluded_exons = form_1_exons; included_unique = tx2_unique; excluded_unique = tx1_unique; included_form_transcript = tx2; excluded_form_transcript = tx1

            if len(included_exons) == 2 and len(excluded_exons) == 2 and len(included_unique) == 2 and len(excluded_unique) == 2:
                if strand == "+":
                    event_type = "AL"
                elif strand == "-":
                    event_type = "AF"

            elif (len(included_exons) > 2 and len(excluded_exons) >= 2 and len(included_unique) > 2 and len(excluded_unique) >= 2) or (len(included_exons) >= 2 and len(excluded_exons) > 2 and len(included_unique) >= 2 and len(excluded_unique) > 2):
                if strand == "+":
                    event_type = "ML"
                elif strand == "-":
                    event_type = "MF"

            elif len(included_exons) >= 2 and len(included_unique) >= 2 and len(excluded_exons) == 1 and len(excluded_unique) == 0:
                if strand == "+":
                    event_type = "UL"
                elif strand == "-":
                    event_type = "UF"

            elif len(included_exons) == 1 and len(excluded_exons) == 1 and len(included_unique) == 1 and len(excluded_unique) == 1:
                if strand == "+":
                    event_type = "AP"
                elif strand == "-":
                    event_type = "AT"

            else:
                if strand == "+":
                    event_type = "CL"
                elif strand == "-":
                    event_type = "CF"

    else:

        if form_1_exons[0][0] < form_2_exons[0][0]:
            included_exons = form_1_exons; excluded_exons = form_2_exons; included_unique = tx1_unique; excluded_unique = tx2_unique; included_form_transcript = tx1; excluded_form_transcript = tx2
        else:
            included_exons = form_2_exons; excluded_exons = form_1_exons; included_unique = tx2_unique; excluded_unique = tx1_unique; included_form_transcript = tx2; excluded_form_transcript = tx1

        if len(included_exons) == 2 and len(excluded_exons) == 2 and len(included_unique) == 2 and len(excluded_unique) == 2:
            if strand == "+":
                event_type = "AF"
            elif strand == "-":
                event_type = "AL"

        elif (len(included_exons) > 2 and len(excluded_exons) >= 2 and len(included_unique) > 2 and len(excluded_unique) >= 2 and (included_exons[0] == excluded_exons[0] or included_exons[-1] == excluded_exons[-1])) or (len(included_exons) >= 2 and len(excluded_exons) > 2 and len(included_unique) >= 2 and len(excluded_unique) > 2 and (included_exons[0] == excluded_exons[0] or included_exons[-1] == excluded_exons[-1])):
            if strand == "+":
                event_type = "MF"
            elif strand == "-":
                event_type = "ML"

        elif len(included_exons) >= 2 and len(included_unique) >= 2 and len(excluded_exons) == 1 and len(excluded_unique) == 0:
            if strand == "+":
                event_type = "UF"
            elif strand == "-":
                event_type = "UL"

        elif len(included_exons) == 1 and len(excluded_exons) == 1 and len(included_unique) == 1 and len(excluded_unique) == 1:
            if strand == "+":
                event_type = "AT"
            elif strand == "-":
                event_type = "AP"

        else:
            if strand == "+":
                event_type = "CF"
            elif strand == "-":
                event_type = "CL"

    outdict = {"included_exons": included_exons, "excluded_exons": excluded_exons, "event_type": event_type, "included_form_transcript": included_form_transcript, "excluded_form_transcript": excluded_form_transcript, "included_unique": included_unique, "excluded_unique": excluded_unique, "strand": strand}


    #for i in ["args", ",".join(map(str,pre_common)), ",".join(map(str, post_common)), ",".join(map(str, tx1_unique)), ",".join(map(str, tx2_unique)), form_1_exons, form_2_exons, strand, tx1, tx2]:

    #    print >> sys.stderr, i

    #print >> sys.stderr, outdict

    if event_type == "WHOOPS":
         print outdict

    return outdict



def write_intron_bedfile(standard_junction_indexed_transcript_dict, outdir):
    '''
        input: junction-indexed transcript dict (dict that has all junctions (which allows definition of introns); values are transcripts to which they belong)
        output: writes a bed file containing intron coordinates to the output directory

        Used for running bedtools intersect to define retained intron events.  An intron that completely overlaps an exon satisfies the definition.
    '''

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

    intron_bedfile.close()


def write_exon_bedfile(standard_transcript_dict, outdir):
    '''
        input: standard transcript dict, output directory.
        output: bedfile of all transcript exon coordinates along with transcripts to which they belong
        
        Used to create input for bedtools intersect, which in turn allows definition of 
    '''

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



def call_bedtools_intersect(outdir, bedtools_path):
    '''
        input: output directory and path to bedtools executable (in arguments, default is set to "bedtools" if path is not specified)
        output: bedtools intersect generates a bed file containing all introns completely overlapping an exon (along with both sets of respective coordinates).  This file is used later to define RI events.
    '''

    bedtools_command = bedtools_path + " intersect -a " + outdir + "/RI_introns.bed -b " + outdir + "/transcript_exon_bedfile.bed -wa -wb -s -f 1.00 > " + outdir + "/RI_introns_exon_overlaps.bed"


    try:
        subprocess.check_output(bedtools_command, shell = True, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print e.output
        sys.exit("Bedtools intersect failed in 'call_bedtools_intersect'.  RI inference failed. Exiting . . . ")

    #subprocess.call(bedtools_command, shell = True)


def retained_intron_events(outdir, standard_transcript_dict, standard_event_dict):
    '''
        input: output directory, standard transcript dict, standard event dict
        output: Modifies standard event dict to include new RI events defined by bedtools intersect of introns and exons
    '''

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

    for event_id in putative_ri_events:


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


        key = chrom + ":" + strand + ":" + str(included_exons[0][0]) + "_left" + "|" + "&".join([]) + "|" + "&".join([str(excluded_exons[0][1]) + "_right", str(excluded_exons[1][0]) + "_left"]) + "|" + str(included_exons[-1][-1]) + "_right"


        if key not in standard_event_dict:

            standard_event_dict[key] = {
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

            max_left = max(included_exons[0][0], excluded_exons[0][0], standard_event_dict[key]["included_exons"][0][0], standard_event_dict[key]["excluded_exons"][0][0])
            standard_event_dict[key]["included_exons"][0][0] = max_left
            standard_event_dict[key]["excluded_exons"][0][0] = max_left

            min_right = min(included_exons[-1][-1], excluded_exons[-1][-1], standard_event_dict[key]["included_exons"][-1][-1], standard_event_dict[key]["excluded_exons"][-1][-1])
            standard_event_dict[key]["included_exons"][-1][-1] = min_right
            standard_event_dict[key]["excluded_exons"][-1][-1] = min_right

        standard_event_dict[key]["included_form_transcripts"] = list(set(standard_event_dict[key]["included_form_transcripts"] + included_form_transcripts))

        standard_event_dict[key]["excluded_form_transcripts"] = list(set(standard_event_dict[key]["excluded_form_transcripts"] + excluded_form_transcripts))


def separate_jc_friendly_events(standard_event_dict):
    '''
        input: standard event dict
        output: two dictionaries, "friendly" and "unfriendly" containing respectively events that can be quantified by junctioncounts in its current form and events that cannot

        junctioncounts relies exclusively on junction reads, with the one exception of RI in which bedtools intersect is used directly on the reads in order to quantify reads overlapping exon-intron junctions. Consequently, AT and AP events (for example) cannot be quantified (yet).
    '''

    friendly = {}
    unfriendly = {}

    for event in standard_event_dict:

        if (len(standard_event_dict[event]["included_junctions"]) == 0 and len(standard_event_dict[event]["included_ei_junctions"]) == 0) or (len(standard_event_dict[event]["excluded_junctions"]) == 0 and len(standard_event_dict[event]["excluded_ei_junctions"]) == 0):

            unfriendly[event] = copy.deepcopy(standard_event_dict[event])

        else:

            friendly[event] = copy.deepcopy(standard_event_dict[event])

    return friendly, unfriendly


def output_ioe(outdir, standard_event_dict, standard_transcript_dict, name = "splice_lib_events"):
    '''
        input: output directory, standard event dict, standard transcript dict, name for output file (not including extension)
        output: IOE files as defined by SUPPA
    '''

    event_centric_ioe_file = open(outdir + "/" + name + ".ioe", 'w')

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


def main(args, transcript_dict = None):

    start_time = time.time()

    print "{0}: {1} seconds elapsed. Starting infer_pairwise_events.  Now parsing arguments.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    parser = argparse.ArgumentParser()
    parser.add_argument("--transcript_gtf", type = str, help = "Full transcript gtf file.  Not required, but if not provided a transcript dict must be passed as a parameter to the main function.")
    parser.add_argument("--outdir", type = str, help = "Path to output directory", required = True)
    parser.add_argument("--min_exon_length", type = int, help = "Minimum allowable exon length in input gtf.  Transcripts with shorter exons will be filtered. (default: 3)", default = 3)
    parser.add_argument("--min_intron_length", type = int, help = "Minimum allowable intron length in input gtf. Transcripts with shorter introns will be filtered. (default: 20)", default = 20)
    parser.add_argument("--max_exon_length", type = int, help = "Maximum allowable exon length in input gtf.  Transcripts with longer exons will be filtered (default = 35000)", default = 35000)
    parser.add_argument("--max_intron_length", type = int, help = "Maximum allowable intron length in input gtf.  Transcripts with longer introns will be fitlered (default = 1000000)", default = 1000000)
    parser.add_argument("--dump_pkl_file", action = "store_true", help = "If set, program will dump pickle file of event dict.")
    parser.add_argument("--bedtools_path", type = str, help = "Path to bedtools executable (default = 'bedtools')", default = "bedtools")
    parser.add_argument("--suppress_output", action = "store_true", help = "If set, GTF, GFF3, and IOE files will not be written.")
    parser.add_argument("--min_AP_AT_dist", type = int, help = "Specifies minimum distance between alternative polyadenylation sites, TSS in order for AP, AT events to be retained. Default is 100", default = 100)

    args = parser.parse_args(args)

    transcript_gtf = args.transcript_gtf
    outdir = args.outdir
    min_exon_length = args.min_exon_length
    min_intron_length = args.min_intron_length
    max_exon_length = args.max_exon_length
    max_intron_length = args.max_intron_length
    dump_pkl_file = args.dump_pkl_file
    bedtools_path = args.bedtools_path
    suppress_output = args.suppress_output
    min_AP_AT_dist = args.min_AP_AT_dist

    if transcript_gtf is None and transcript_dict is None:

        sys.exit("Neither transcript gtf nor transcript dict provided.  Exiting . . . ")

    elif not (transcript_gtf is None or transcript_dict is None):

        sys.exit("Both transcript gtf and transcript dict provided.  Please provided only one.  Exiting  . . . ")

    subprocess.call("mkdir -p " + outdir, shell = True)

    if transcript_gtf is not None:

        print "{0}: {1} seconds elapsed. Arguments parsed.  Now importing transcriptome gtf.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
        
        standard_transcript_dict = splice_lib.generate_standard_transcript_dict(transcript_gtf)

    elif transcript_dict is not None:

        print "{0}: {1} seconds elapsed. Arguments parsed.  Now copying transcript dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

        standard_transcript_dict = copy.deepcopy(transcript_dict)

    print "{0}: {1} seconds elapsed. Transcriptome imported.  Now sorting transcript exons.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    splice_lib.sort_transcript_dict_exons(standard_transcript_dict)

    print "{0}: {1} seconds elapsed. Sorted transcript exons.  Now adding transcript junctions to transcript dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    splice_lib.add_junctions_to_transcript_dict(standard_transcript_dict)

    print "{0}: {1} seconds elapsed. Added transcript junctions.  Now filtering transcriptome to remove transcripts with short introns/exons.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    filter_standard_transcript_dict(standard_transcript_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length)

    print "{0}: {1} seconds elapsed. Filtered transcriptome.  Now creating junction-indexed transcript dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    standard_junction_indexed_transcript_dict = splice_lib.index_transcripts_by_junctions(standard_transcript_dict)

    print "{0}: {1} seconds elapsed. Created junction-indexed transcript dict.  Creating node dict to associate groups of transcripts.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))


    node_index_tx = node_dir(standard_transcript_dict)


    print "{0}: {1} seconds elapsed. Created node dict.  Now comparing pairs of transcripts to create standard event dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    standard_event_dict = call_compare_transcripts(node_index_tx, standard_transcript_dict)

    print "{0}: {1} seconds elapsed. Created standard event dict.  Now searching for RI events that may lack common nodes.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    write_intron_bedfile(standard_junction_indexed_transcript_dict, outdir)

    write_exon_bedfile(standard_transcript_dict, outdir)

    call_bedtools_intersect(outdir, bedtools_path)

    retained_intron_events(outdir, standard_transcript_dict, standard_event_dict)

    print "{0}: {1} seconds elapsed. RI event identification complete.  Now cleaning/filtering events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))    

    splice_lib.complete_event_dict(standard_event_dict, suppress_unique_edges = True, suppress_eij = False, no_ends = True, inform_using_ri_events = True)
    splice_lib.collapse_redundant_junction_events(standard_event_dict, outdir)

    #splice_lib.complete_event_dict(standard_event_dict)

    standard_event_dict = splice_lib.rename_events(standard_event_dict)

    junction_indexed_event_dict = splice_lib.generate_junction_indexed_event_dict(standard_event_dict)

    splice_lib.filter_overlapping_se_alt_donacc_events(standard_event_dict, junction_indexed_event_dict)

    filter_event_dict(standard_event_dict, outdir, min_exon_length, min_intron_length, max_exon_length, max_intron_length, min_AP_AT_dist = min_AP_AT_dist)
    
    print "{0}: {1} seconds elapsed. Event cleaning/filtering complete.  Now renaming events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    standard_event_dict = splice_lib.rename_events(standard_event_dict)


    print "{0}: {1} seconds elapsed. Event renaming complete.  Now separating events into junctioncounts friendly/unfriendly categories.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    friendly, unfriendly = separate_jc_friendly_events(standard_event_dict)


    if not suppress_output:
        print "{0}: {1} seconds elapsed. Events separated.  Now outputting event gtfs and bedfiles.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

        splice_lib.output_event_gtf(friendly, outdir)
        splice_lib.output_event_gtf(unfriendly, outdir, name = "junctioncounts_unfriendly_events")
        splice_lib.output_event_bedfile(friendly, outdir)
        splice_lib.output_event_bedfile(unfriendly, outdir, name = "junctioncounts_unfriendly_events")


        print "{0}: {1} seconds elapsed. Event gtf output.  Now outputting MISO-style gff3s.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

        splice_lib.output_miso_event_gff3(friendly, outdir)
        splice_lib.output_miso_event_gff3(unfriendly, outdir, name="junctioncounts_unfriendly_events")

        if dump_pkl_file:
            print "{0}: {1} seconds elapsed. MISO_style event gff3 output.  Now dumping pickle file of junctioncounts-friendly event dict.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
            pickle.dump(friendly, open(outdir + "/splice_lib_events.pkl", 'wb' ), -1 )
            print "{0}: {1} seconds elapsed. Pickle file dump complete.  Now outputting IOE file.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

        else:
            print "{0}: {1} seconds elapsed. MISO_style event gff3 output.  Now outputting IOE files.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

        output_ioe(outdir, friendly, standard_transcript_dict)
        output_ioe(outdir, unfriendly, standard_transcript_dict, name = "junctioncounts_unfriendly_events")

        print "{0}: {1} seconds elapsed. IOE file output.  Tabulating event type counts.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    else:

        print "{0}: {1} seconds elapsed. Events separated.  Tabulating event type counts.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))

    event_type_counts = splice_lib.assess_event_types(standard_event_dict)

    print "{0}: {1} seconds elapsed. Event type counts tabulated.  Found {2} total events with {3} MS, {4} SE, {5} A3, {6} A5, {7} AF, {8} AL, {9} MF, {10} ML, {11} CF, {12} CL, {13} UF, {14} UL, {15} AT, {16} AP, {17} RI, {18} MX, {19} CO, {20} AB and {21} MR events.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)), str(event_type_counts["total"]), str(event_type_counts["MS"]), str(event_type_counts["SE"]), str(event_type_counts["A3"]), str(event_type_counts["A5"]), str(event_type_counts["AF"]), str(event_type_counts["AL"]), str(event_type_counts["MF"]), str(event_type_counts["ML"]), str(event_type_counts["CF"]), str(event_type_counts["CL"]), str(event_type_counts["UF"]), str(event_type_counts["UL"]), str(event_type_counts["AT"]), str(event_type_counts["AP"]),  str(event_type_counts["RI"]), str(event_type_counts["MX"]),  str(event_type_counts["CO"]), str(event_type_counts["AB"]), str(event_type_counts["MR"]))

    print "{0}: {1} seconds elapsed. infer_pairwise_events complete.".format(str(datetime.now().replace(microsecond = 0)), str(round(time.time() - start_time, 1)))
    print "Ceci n'est pas un algorithme bioinformatique."

    return friendly


if __name__ ==     '__main__':

    main(sys.argv[1:])

    #outdict = compare_transcripts("transcript1", "transcript2", [[100,200],[300,400],[500,600],[700,800],[900,1000], [1100, 1200], [1300, 1400], [1500, 1600]], [[20,50],[300,400],[500,600],[900,1000]], "-", "chr1")
    
    #for i in outdict:

    #    print i
    #    print outdict[i]
    #    print ""

