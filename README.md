# junctionCounts

junctionCounts (https://doi.org/10.1093/nargab/lqae093) is a tool intended to identify and quantify binary alternative splicing events from RNA-seq data. junctionCounts currently consists of two primary utilities: `infer_pairwise_events.py`, which identifies binary alternative splicing events (e.g. skipped exon events), and `junctionCounts.py`, which quantifies those events (i.e. provides a percent-spliced-in, or PSI value) from RNA-seq data (specifically from a BAM file). These utilities can be used independent of one another, and events identified independently by `infer_pairwise_events.py` and the transcripts with which they are associated (as marked in the .ioe files - a format introduced by the authors of SUPPA) may be useful to a variety of researchers, regardless of their interest in performing short read quantification. While `junctionCounts.py` could in principle be run without the use of `infer_pairwise_events.py`, it may be easiest to use both tools together for quantification due to the specific reqirements of the input GTF for the quantification. For predicting the functional outcomes of alternative splicing events, including NMD, NSD and coding-to-noncoding switch events, check out junctionCounts' partner utilities at https://github.com/ajw2329/cds_insertion.

![Screenshot](./figures/junctionCounts_flowchart.png)

# Dependencies

`python3`  
`ncls` (https://github.com/hunt-genes/ncls)  
`numpy`  
`pysam`
`splice_lib` (github.com/ajw2329/splice_lib/ - install using `python setup.py install`)  
`bedtools`  
`R (≥ 3.5.0)` `dplyr` `DEXSeq` `pacman` `optparse` 

Note that `ncls`, `numpy`, and `pysam` are all python packages that can be installed using `pip install package_name`.

## Benchmarking

Panels A-K: We evaluated junctionCounts on its performance relative to four established splicing analysis tools: MAJIQ v2.5.6.dev1+g8423f68b (https://doi.org/10.1038/s41467-023-36585-y), rMATS-turbo v4.3.0 (https://doi.org/10.1038/s41596-023-00944-2), splAdder v3.0.4 (https://doi.org/10.1093/bioinformatics/btw076) and Whippet v1.6.2 (https://doi.org/10.1016/j.molcel.2018.08.018). We generated four paired-end simulated datasets from real RNA-seq data using polyester v1.36.0 (https://doi.org/10.1093/bioinformatics/btv272) with the arguments: read.length = 100, fragment.length.min = 100, fragment.length.max = 500, fragment.length.mean = 180, fragment.length.sd = 40 and simulated.sequencing.error = TRUE. Three of the simulated datasets: 25M, 50M and 75M, were used to evaluate performance at different library depths (25, 50 and 75 million reads per library respectively). These datasets were modeled on mouse cerebellum and liver RNA-seq data in triplicate for each cell type from Vaquero-Garcia et al. (2016) (https://doi.org/10.7554/eLife.11752), who replicated the experiments in Zhang et al. (2014) (https://doi.org/10.1073/pnas.1408886111). The fourth simulated dataset was modeled on human RNA-seq data from HeLa cells treated with either spliceostatin A or DMSO (https://doi.org/10.1038/s41467-017-02007-z) at 50 million reads per library with triplicates for each condition. This dataset introduced a larger pool of potential AS events to test compared to the murine simulated datasets, and importantly increased the total number of MR events tested.

Panels L-N: To test junctionCounts' NMD event predictions, we procured RT-PCR data for well-characterized NMD targets and accompanying RNA-seq data from HEK-293 cells upon UPF1 knockdown versus non-targeting siRNA (https://doi.org/10.15252/embj.2021108898). We then ran junctionCounts and its partner utilities with default settings on the RNA-seq data to predict NMD switch events within the NMD targets verified with RT-PCR in the original study. Out of the 13 predicted NMD switch events associated with the NMD targets, 11 had the expected dPSI directionality (Panel L). We then used the same approach on a dataset in which treatment with emetine, a translation elongation inhibitor, was reported to increase the abundance of NMD substrates in HEK-293 T cells (https://doi.org/10.1093/nar/gkw1109). cdsInsertion and findSwitchEvents identified 636 potential NMD switch events with significant changes in splicing (|dPSI| ≥ 0.1 and Q-value ≤ 0.05) upon emetine treatment, out of which 472 (74.2%) exhibited the expected dPSI directionality (Panels M, N).

![Screenshot](./figures/benchmarking_results.png)  

# infer_pairwise_events.py

## Description

`infer_pairwise_events.py` takes as input a transcriptome GTF file and performs a pairwise comparison of all overlapping transcripts to identify the set of (one or more) minimal binary alternative splicing (or other forms of co/post-transcriptional processing) events that differentiate each pair. `infer_pairwise_events.py` considers an alternative event to be a collection of non-shared (i.e. unique to one transcript in the considered pair) exon boundaries that is bounded either a) by a set of shared outer exon boundaries, or b) by transcript termini. This definition includes all canonical binary alternative splicing events (of which the author is aware such as skipped exons, alternative donors/acceptors, alternative first/last exons) while also including non-canonical alternative splicing events. A schematic illustrating the events in question is below. Note that certain event "types" (specifically CO, CF, and CL) are generic terms that stand for 'complex', 'complex first', and 'complex last' respectively. These terms do not relate to modern concepts of alternative splicing complexity (i.e. non-binary events) but rather are buckets for events that both a) do not fit into the other, more specific categories, and b) do not involve the transcript termini (CO), involve different 5'-termini (CF), or involve different 3'-termini (CL). This also applies to event types MF, ML, MS, MR (multiple first, last, multiple skipped exon, and multiple retained intron respectively) - these events are highly analogous to AF, AL and SE events, but are not restricted to one pair of alternative exons in the case of AF, AL (for example, an MF or ML event could have 3 exons unique to the included form and two unique to the excluded form) or a single alternative exon in the case of SE. 

![Screenshot](./figures/event_models.png)

The meaning of each abbreviation (bearing in mind the above discussion) is as follows:   
SE - skipped exon   
MS - multiple skipped exon   
MX - mutually exclusive exons   
A5 - alternative 5'-splice site   
A3 - alternative 3'-splice site   
RI - retained intron   
MR - multiple retained intron   
AF - alternative first exon   
AL - alternative last exon   
MF - multiple first exon   
ML - multiple last exon   
CF - complex first   
CL - complex last   

## Basic usage

```
python /path/to/infer_pairwise_events.py --transcript_gtf /path/to/transcriptome.gtf --outdir /path/to/output/dir/
```

### Input file description

#### transcriptome GTF

| seqname | source    | feature    | start | end   | score | strand | frame | attribute                                                                                                                         |
|---------|-----------|------------|-------|-------|-------|--------|-------|-----------------------------------------------------------------------------------------------------------------------------------|
| chr1    | StringTie | transcript | 11869 | 14409 | 1000  | +      | .     | gene_id "MSTRG.1"; transcript_id "ENST00000456328.2_1"; gene_name "DDX11L1"; ref_gene_id "ENSG00000223972.5_2";                   |
| chr1    | StringTie | exon       | 11869 | 12227 | 1000  | +      | .     | gene_id "MSTRG.1"; transcript_id "ENST00000456328.2_1"; exon_number "1"; gene_name "DDX11L1"; ref_gene_id "ENSG00000223972.5_2";  |
| chr1    | StringTie | exon       | 12613 | 12721 | 1000  | +      | .     | gene_id "MSTRG.1"; transcript_id "ENST00000456328.2_1"; exon_number "2"; gene_name "DDX11L1"; ref_gene_id "ENSG00000223972.5_2";  |
| chr1    | StringTie | exon       | 13221 | 14409 | 1000  | +      | .     | gene_id "MSTRG.1"; transcript_id "ENST00000456328.2_1"; exon_number "3"; gene_name "DDX11L1"; ref_gene_id "ENSG00000223972.5_2";  |

### Output file descriptions

`splice_lib_events.gtf`

| seqname | source           | feature | start     | end       | score | strand | frame | attribute                                                  |
|---------|------------------|---------|-----------|-----------|-------|--------|-------|------------------------------------------------------------|
| chr3    | splice_lib_event | exon    | 52321872  | 52321892  | .     | +      | .     | gene_id "SE.0006314"; transcript_id "SE.0006314_included"; |
| chr3    | splice_lib_event | exon    | 52324320  | 52324735  | .     | +      | .     | gene_id "SE.0006314"; transcript_id "SE.0006314_included"; |
| chr3    | splice_lib_event | exon    | 52324976  | 52325127  | .     | +      | .     | gene_id "SE.0006314"; transcript_id "SE.0006314_included"; |
| chr3    | splice_lib_event | exon    | 52321872  | 52321892  | .     | +      | .     | gene_id "SE.0006314"; transcript_id "SE.0006314_excluded"; |
| chr3    | splice_lib_event | exon    | 52324976  | 52325127  | .     | +      | .     | gene_id "SE.0006314"; transcript_id "SE.0006314_excluded"; |

`splice_lib_events.bed`

| seqname | start | end | event_id | score | strand |
|------|-----------|-----------|------------|------|---|
| chr3 | 52321872  | 52325127  | SE.0006314 | 1000 | + |
| chr3 | 42605015  | 42610565  | SE.0003301 | 1000 | - |
| chr9 | 128419930 | 128469482 | MS.0000307 | 1000 | - |

`splice_lib_events.ioe`

| seqname | gene_id     | event_id               | included_transcripts              | total_transcripts                                |
|---------|-------------|------------------------|-----------------------------------|--------------------------------------------------|
| chr10   | CREM        | CREM;MF.0011380        | ENST00000354759.7_1               | ENST00000354759.7_1,ENST00000490511.1_1          |
| chr13   | EBPL        | EBPL;MF.0011381        | ENST00000378284.6_1               | ENST00000378284.6_1,MSTRG.7758.6                 |
| chr10   | FRMD4A      | FRMD4A;MF.0011382      | ENST00000357447.6_2               | ENST00000357447.6_2,MSTRG.3067.1                 |
| chr9    | MSTRG.27860 | MSTRG.27860;MF.0011383 | MSTRG.27860.2                     | ENST00000382293.7_1,MSTRG.27860.2                |
| chr19   | GATAD2A     | GATAD2A;MF.0011384     | ENST00000360315.7_1,MSTRG.14812.1 | ENST00000360315.7_1,MSTRG.14812.1,MSTRG.14812.12 |
| chr16   | ZFP1        | ZFP1;MF.0011385        | ENST00000332307.4_2               | ENST00000332307.4_2,ENST00000393430.6_1          |
| chr17   | METTL23     | METTL23;MF.0011386     | ENST00000586738.5_1               | ENST00000586200.1_1,ENST00000586738.5_1          |
| chr2    | MSTRG.16852 | MSTRG.16852;MF.0011387 | MSTRG.16852.6                     | ENST00000422956.6_1,MSTRG.16852.6                |
| chr5    | MATR3       | MATR3;MF.0011388       | ENST00000394800.6_2               | ENST00000394800.6_2,MSTRG.23433.20               |

`splice_lib_events.gff3`

| seqname | source           | feature    | start    | end      | score | strand | frame | attribute                                                |
|---------|------------------|------------|----------|----------|-------|--------|-------|----------------------------------------------------------|
| chr3    | splice_lib_event | gene       | 52321872 | 52325127 | .     | +      | .     | ID=SE.0006314;Name=SE.0006314                            |
| chr3    | splice_lib_event | transcript | 52321872 | 52325127 | .     | +      | .     | ID=SE.0006314_included;Parent=SE.0006314                 |
| chr3    | splice_lib_event | exon       | 52321872 | 52321892 | .     | +      | .     | ID=exon:SE.0006314_included:1;Parent=SE.0006314_included |
| chr3    | splice_lib_event | exon       | 52324320 | 52324735 | .     | +      | .     | ID=exon:SE.0006314_included:2;Parent=SE.0006314_included |
| chr3    | splice_lib_event | exon       | 52324976 | 52325127 | .     | +      | .     | ID=exon:SE.0006314_included:3;Parent=SE.0006314_included |
| chr3    | splice_lib_event | transcript | 52321872 | 52325127 | .     | +      | .     | ID=SE.0006314_excluded;Parent=SE.0006314                 |
| chr3    | splice_lib_event | exon       | 52321872 | 52321892 | .     | +      | .     | ID=exon:SE.0006314_excluded:1;Parent=SE.0006314_excluded |
| chr3    | splice_lib_event | exon       | 52324976 | 52325127 | .     | +      | .     | ID=exon:SE.0006314_excluded:2;Parent=SE.0006314_excluded |

`junctioncounts_unfriendly_events.gtf`

`junctioncounts_unfriendly_events.bed`

`junctioncounts_unfriendly_events.ioe`

`junctioncounts_unfriendly_events.gff3`

## Help statement

```
usage: infer_pairwise_events.py [-h] [--transcript_gtf TRANSCRIPT_GTF]
                                --outdir OUTDIR
                                [--min_exon_length MIN_EXON_LENGTH]
                                [--min_intron_length MIN_INTRON_LENGTH]
                                [--max_exon_length MAX_EXON_LENGTH]
                                [--max_intron_length MAX_INTRON_LENGTH]
                                [--dump_pkl_file]
                                [--bedtools_path BEDTOOLS_PATH]
                                [--suppress_output]
                                [--min_AP_AT_dist MIN_AP_AT_DIST]

optional arguments:
  -h, --help            show this help message and exit
  --transcript_gtf TRANSCRIPT_GTF
                        Full transcript gtf file. Not required, but if not
                        provided a transcript dict must be passed as a
                        parameter to the main function.
  --outdir OUTDIR       Path to output directory
  --min_exon_length MIN_EXON_LENGTH
                        Minimum allowable exon length in input gtf.
                        Transcripts with shorter exons will be filtered.
                        (default: 3)
  --min_intron_length MIN_INTRON_LENGTH
                        Minimum allowable intron length in input gtf.
                        Transcripts with shorter introns will be filtered.
                        (default: 20)
  --max_exon_length MAX_EXON_LENGTH
                        Maximum allowable exon length in input gtf.
                        Transcripts with longer exons will be filtered
                        (default = 35000)
  --max_intron_length MAX_INTRON_LENGTH
                        Maximum allowable intron length in input gtf.
                        Transcripts with longer introns will be fitlered
                        (default = 1000000)
  --dump_pkl_file       If set, program will dump pickle file of event dict.
  --bedtools_path BEDTOOLS_PATH
                        Path to bedtools executable (default = 'bedtools')
  --suppress_output     If set, GTF, GFF3, and IOE files will not be written.
  --min_AP_AT_dist MIN_AP_AT_DIST
                        Specifies minimum distance between alternative
                        polyadenylation sites, TSS in order for AP, AT events
                        to be retained. Default is 100

```

# junctionCounts.py

## Description

## Basic usage

```
python /path/to/junctionCounts.py --event_gtf /path/to/splice_lib_events.gtf --bam /path/to/input_file.bam --outdir /path/to/output/dir/ --sample_name sample_name
```

### Output file descriptions

`h9esc_cyto_jd004_count_psi_outfile.tsv`

| sample_name      | event_id   | event_type | min_ijc | min_ejc | avg_psi | max_gene_frac | all_ijc     | all_ejc   | avg_ijc       | avg_ejc | span_psi | min_psi | ijc_min_psi | ejc_min_psi | max_psi | ijc_max_psi | ejc_max_psi | mid_psi | bootstrap_num |
|------------------|------------|------------|---------|---------|---------|---------------|-------------|-----------|---------------|---------|----------|---------|-------------|-------------|---------|-------------|-------------|---------|---------------|
| h9esc_cyto_jd004 | RI.0024780 | RI         | 2       | 0       | 1.0     | 0.08          | 25,2        | 0         | 13.5          | 0.0     | 0.0      | 1.0     | 25          | 0           | 1.0     | 25          | 0           | 1.0     | NA            |
| h9esc_cyto_jd004 | RI.0022510 | RI         | 5       | 0       | 1.0     | 0.2           | 5,5         | 0         | 5.0           | 0.0     | 0.0      | 1.0     | 5           | 0           | 1.0     | 5           | 0           | 1.0     | NA            |
| h9esc_cyto_jd004 | CO.0032123 | CO         | 1       | 0       | 1.0     | 0.0058        | 7,6,1,30,1  | 0,0,0     | 9.0           | 0.0     | 0.0      | 1.0     | 7           | 0           | 1.0     | 7           | 0           | 1.0     | NA            |
| h9esc_cyto_jd004 | AF.0004999 | AF         | 3       | 0       | 1.0     | 1.0           | 3           | 0         | 3.0           | 0.0     | 0.0      | 1.0     | 3           | 0           | 1.0     | 3           | 0           | 1.0     | NA            |
| h9esc_cyto_jd004 | CF.0019530 | CF         | 1       | 43      | 0.0227  | 1.0           | 1           | 43        | 1.0           | 43.0    | 0.0      | 0.0227  | 1           | 43          | 0.0227  | 1           | 43          | 0.0227  | NA            |
| h9esc_cyto_jd004 | CF.0069352 | CF         | 2       | 1       | 0.7579  | 0.0053        | 3,2,6       | 1         | 3.66666666667 | 1.0     | 0.1905   | 0.6667  | 2           | 1           | 0.8571  | 6           | 1           | 0.7619  | NA            |
| h9esc_cyto_jd004 | ML.0016599 | ML         | 0       | 6       | 0.5895  | 0.0028        | 2010,2114,0 | 545,322,6 | 1374.66666667 | 291.0   | 0.9972   | 0.0     | 0           | 545         | 0.9972  | 2114        | 6           | 0.4986  | NA            |
| h9esc_cyto_jd004 | ML.0022641 | ML         | 28      | 3       | 0.915   | 0.4828        | 38,28       | 3         | 33.0          | 3.0     | 0.0236   | 0.9032  | 28          | 3           | 0.9268  | 38          | 3           | 0.915   | NA            |
| h9esc_cyto_jd004 | CL.0039596 | CL         | 9       | 4       | 0.7457  | 0.0336        | 9,21        | 5,4       | 15.0          | 4.5     | 0.1971   | 0.6429  | 9           | 5           | 0.84    | 21          | 4           | 0.7414  | NA            |


#### Field descriptions:

`sample_name` : sample name provided by user via --sample_name  
`event_id` : unique name for event  
`event_type` : type of event as discussed above  
`min_ijc` : (integer) minimum read count of all junctions in the included form  
`min_ejc` : (integer) minimum read count of all junctions in the excluded (skipped) form  
`avg_psi` : (float) average of all PSI values computed all pairs of included, excluded junction counts  
`max_gene_frac` : (float) `max(min_ijc, min_sjc)/max_gene_jc` where `max_gene_jc` is the maximum junction count for any junction in the event's gene. Only computed if `--calc_gene_frac` is passed along with an ioe file  
`all_ijc` : (comma-sep list of integers) list of junction counts for all included form junctions  
`all_ejc` : (comma-sep list of integers) list of junction counts for all excluded form junctions  
`avg_ijc` : (float) mean junction count for included form junctions  
`avg_ejc` : (float) mean junction count for excluded form junctions  
`span_psi` : (float) `max_psi - min_psi`, i.e. the difference between the largest and smallest PSI values computed using any pair of included, excluded form junction counts. Large spans may indicate unreliability of the PSI estimate for a number of reasons such as partial junction overlap with another event combined with an abundance of multiply counted reads.  
`min_psi` : (float) minimum of all PSI values computed all pairs of included, excluded junction counts  
`ijc_min_psi` : (float) included form junction count corresponding to the `min_psi` calculation  
`ejc_min_psi` : (float) excluded form junction count corresponding to the `min_psi` calculation  
`max_psi` : (float) maximum of all PSI values computed all pairs of included, excluded junction counts  
`ijc_max_psi` : (float) included form junction count corresponding to the `max_psi` calculation  
`ejc_max_psi` : (float) excluded form junction count corresponding to the `max_psi` calculation  
`mid_psi` : (float) midpoint between `min_psi` and `max_psi`  

## Help statement

```
usage: junctionCounts.py [-h] [--event_gtf EVENT_GTF] --bam BAM [--se]
                         [--forward_read FORWARD_READ] --outdir OUTDIR
                         --sample_name SAMPLE_NAME [--dump_pkl_dict]
                         [--gzipped] [--event_ioe EVENT_IOE]
                         [--calc_gene_frac] [--suppress_output]
                         [--enable_edge_use] [--turn_off_no_ends]
                         [--suppress_eij_use] [--disable_ri_extrapolation]
                         [--n_bootstraps N_BOOTSTRAPS]

optional arguments:
  -h, --help            show this help message and exit
  --event_gtf EVENT_GTF
                        GTF describing pairwise events
  --bam BAM             BAM read file for counting junction reads
  --se                  BAM file is single-ended. Default assumes paired-end
                        reads.
  --forward_read FORWARD_READ
                        Specify read that matches the 'forward' strand.
                        Options are 'R1','R2' or 'unstranded' if the library
                        is unstranded. Unstranded use is not currently
                        recommended. default = 'R2'
  --outdir OUTDIR       Path for output files
  --sample_name SAMPLE_NAME
                        Sample name
  --dump_pkl_dict       Dumps a pkl file of the splice event dict with all
                        quantifications
  --gzipped             Output files will be gzipped if set
  --event_ioe EVENT_IOE
                        Event ioe file - used to recover event-gene
                        association
  --calc_gene_frac      Requires IOE file to be passed to --event_ioe. If set,
                        a maximal event fraction of gene expression will be
                        estimated by max(min_excluded,
                        min_included)/max(gene_junctions)
  --suppress_output     Suppresses output files if set
  --enable_edge_use     Use individual exon edges that are unique to form in
                        quantification
  --turn_off_no_ends    Disable the exclusion of transcript-termini from
                        isoform-specific exon edge quantification (default is
                        to exclude ends)
  --suppress_eij_use    Don't use exon-overlapping exon-intron junctions for
                        quantification except for RI events.
  --disable_ri_extrapolation
                        Disable the use of RI/MR included forms to inform
                        quantification of other events that contain these
                        exons
  --n_bootstraps N_BOOTSTRAPS
                        Number of bootstraps. default = 0
```
# DEXSeq_comparison.R

## Description

`DEXSeq_comparison.R` statistically tests changes in PSI between two conditions with at least 2 replicates per condition.

## Basic usage

```
Rscript /path/to/DEXSeq_comparison.R -e /path/to/infer_pairwise_events -p /path/to/jc_psi_files -s /path/to/sampleTable.csv 
```

## Example sample_table.csv

| sample | condition | file                        |
|--------|-----------|-----------------------------|
| esc_1  | esc       | esc_1_count_psi_outfile.tsv |
| esc_2  | esc       | esc_2_count_psi_outfile.tsv |
| esc_3  | esc       | esc_3_count_psi_outfile.tsv |
| npc_1  | npc       | npc_1_count_psi_outfile.tsv |
| npc_2  | npc       | npc_2_count_psi_outfile.tsv |
| npc_3  | npc       | npc_3_count_psi_outfile.tsv |

## Help statement

```
usage: DEXSeq_comparison.py [-h] [-e, --events INFER_PAIRWISE_EVENTS_DIR]
                                [-p, --psi_dir PSI_TSV_DIR]
                                [-s, --sample_table SAMPLE_TABLE]
                                [--min_jc MINIMUM_JC]
                                [--min_psi MINIMUM_PSI]

optional arguments:
  -h, --help            show this help message and exit
  --min_jc MINIMUM_JC   Minimum number of total mean junction read suport across
                        conditions to report an event: i.e. for conditions A and B,
                        mean(ijc_A) + mean(ejc_A) + mean(ijc_B) + mean(ejc_B) >= min_jc.
`                       (Default: 15)
  --min_psi MINIMUM_PSI Minimum max(mean(PSI)) across conditions to report an
                        event. (Default: 0.1)
  --ri_span RI_SPAN_PSI Maximum max(mean(span_PSI)) across conditions to report
                        an RI/MR event. (Default: 0.03)
```

### Output file descriptions

`npc_esc_dpsi.csv`

| gene | event_id | event_type | chr  | start    | end      | strand | esc_mean_ijc | esc_mean_ejc | esc_mean_psi | npc_mean_ijc | npc_mean_ejc | npc_mean_psi | dpsi    | event_qval | sig |
|------|----------|------------|------|----------|----------|--------|--------------|--------------|--------------|--------------|--------------|--------------|---------|------------|-----|
| ADD1 | A5.00001 | A5         | chr4 | 2904764  | 2907844  | +      | 254.75       | 307.5        | 0.45125      | 217.5        | 255          | 0.45645      | 0.0052  | 1          | 0   |
| BRD2 | AF.00001 | AF         | chr6 | 32968660 | 32972141 | +      | 48.25        | 282          | 0.1431       | 118          | 81           | 0.54875      | 0.40565 | 0          | 1   |


#### Field descriptions:

`gene` : gene symbol  
`event_id` : unique name for event  
`event_type` : type of event  
`chr` : chromosome  
`start` : start coordinate of event                                                                                                                                                                 
`end` : end coordinate of event  
`strand` : strand  
`conditionA_mean_ijc` : (float) mean read count of all junctions in the included form across replicates in conditionA   
`conditionA_mean_ejc` : (float) mean read count of all junctions in the excluded (skipped) form across replicates in conditionA    
`conditionA_mean_psi` : (float) mean of PSI values computed for all pairs of included, excluded junction counts across replicates in condtionA    
`conditionB_mean_ijc` : (float) mean read count of all junctions in the included form across replicates in conditionB   
`conditionB_mean_ejc` : (float) mean read count of all junctions in the excluded (skipped) form across replicates in conditionB    
`conditionB_mean_psi` : (float) mean of PSI values computed for all pairs of included, excluded junction counts across replicates in condtionB  
`dpsi` : (float) difference of mean(PSI) between conditions  
`event_qval` : (float) Q-value describing the positive false discovery rate of the difference of mean(PSI) between conditions.  
`sig` : 1 = event is significantly more included (dpsi ≥ 0.1, event_qval ≤ 0.05) in condition B relative to condition A.   
         -1 = event is significantly less included (dpsi ≤ -0.1, event_qval ≤ 0.05) in condition B relative to condition A.  
          0 = event is not significantly different between conditions (|dpsi| < 0.1, event_qval > 0.05).  
