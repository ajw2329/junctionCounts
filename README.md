# junctionCounts

junctinCounts is a tool intended to identify and quantify binary alternative splicing events from RNA-seq data. junctionCounts currently consists of two primary utilities: `infer_pairwise_events.py`, which identifies binary alternative splicing events (e.g. skipped exon events), and `junctionCounts.py`, which quantifies those events (i.e. provides a percent-spliced-in, or PSI value) from Illumina RNA-seq data (specifically from a BAM file). These utilities can be used independent of one another, and events identified independently by `infer_pairwise_events.py` and the transcripts with which they are associated (as marked in the .ioe files - a format introduced by the authors of SUPPA) may be useful to a variety of researchers, regardless of their interest in performing short read quantification. While `junctionCounts.py` could in principle be run without the use of `infer_pairwise_events.py`, it may be easiest to use both tools together for quantification due to the specific reqirements of the input GTF for the quantification.

# Dependencies

`python2.7`  
`ncls` (https://github.com/hunt-genes/ncls)  
`numpy`  
`pysam`  

Note that `ncls`, `numpy`, and `pysam` are all python packages that can be installed using `pip install package_name`.

# infer_pairwise_events.py

## Description

`infer_pairwise_events.py` takes as input a transcriptome GTF file and performs a pairwise comparison of all overlapping transcripts to identify the set of (one or more) minimal binary alternative splicing (or other forms of co/post-transcriptional processing) events that differentiate each pair. `infer_pairwise_events.py` considers an alternative event to be a collection of non-shared (i.e. unique to one transcript in the considered pair) exon boundaries that is bounded either a) by a set of shared outer exon boundaries, or b) by transcript termini. This definition includes all canonical binary alternative splicing events (of which the author is aware such as skipped exons, alternative donors/acceptors, alternative first/last exons) while also including non-canonical alternative splicing events. A schematic illustrating the events in question is below. Note that certain event "types" (specifically CO, CF, and CL) are generic terms that stand for 'complex', 'complex first', and 'complex last' respectively. These terms do not relate to modern concepts of alternative splicing complexity (i.e. non-binary events) but rather are buckets for events that both a) do not fit into the other, more specific categories, and b) do not involve the transcript termini (CO), involve different 5'-termini (CF), or involve different 3'-termini (CL). This also applies to event types MF, ML, MS, MR (multiple first, last, multiple skipped exon, and multiple retained intron respectively) - these events are highly analogous to AF, AL and SE events, but are not restricted to one pair of alternative exons in the case of AF, AL (for example, an MF or ML event could have 3 exons unique to the included form and two unique to the excluded form) or a single alternative exon in the case of SE. UF and UL events are another pair of unconventional events. They are situations in which the included form transcript's penterminal exon is the excluded form transcript's terminal exon, and the included form contains one additional (up or downstream respectively) exon. These events can be found in established annotations but seem unusual and are possibly artifactual.

![Screenshot](./figures/event_schematic.png)

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
AT - alternative transcription start site (tandem UTR)   
AL - alternative polyadenylation (tandem UTR)   
UF - unique first   
UL - unique last   

## Basic usage

```
python2.7 /path/to/infer_pairwise_events.py --transcript_gtf /path/to/transcriptome.gtf --outdir /path/to/output/dir/
```

### Input file description

#### transcriptome GTF

    | chr1 | ENSEMBL | gene       | 134901 | 139379 | . | - | . | gene_id "ENSG00000237683.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; level 3; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                                                                                                                                                                                                                                    |
| chr1 | ENSEMBL | transcript | 134901 | 139379 | . | - | . | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                             |
| chr1 | ENSEMBL | exon       | 137621 | 139379 | . | - | . | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; exon_number 1; exon_id "ENSE00002221580.1"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19"; |
| chr1 | ENSEMBL | CDS        | 138533 | 139309 | . | - | 0 | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; exon_number 1; exon_id "ENSE00002221580.1"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19"; |
| chr1 | ENSEMBL | stop_codon | 138530 | 138532 | . | - | 0 | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; exon_number 1; exon_id "ENSE00002221580.1"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19"; |
| chr1 | ENSEMBL | exon       | 134901 | 135802 | . | - | . | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; exon_number 2; exon_id "ENSE00002314092.1"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19"; |
| chr1 | ENSEMBL | UTR        | 137621 | 138532 | . | - | . | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                             |
| chr1 | ENSEMBL | UTR        | 139310 | 139379 | . | - | . | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                             |
| chr1 | ENSEMBL | UTR        | 134901 | 135802 | . | - | . | gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; level 3; protein_id "ENSP00000473460.1"; tag "basic"; tag "appris_principal"; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                             |
| chr1 | ENSEMBL | gene       | 157784 | 157887 | . | - | . | gene_id "ENSG00000222623.1"; gene_type "snRNA"; gene_status "KNOWN"; gene_name "RNU6-1100P"; level 3; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                                                                                                                                                                                                                                             |
| chr1 | ENSEMBL | transcript | 157784 | 157887 | . | - | . | gene_id "ENSG00000222623.1"; transcript_id "ENST00000410691.1"; gene_type "snRNA"; gene_status "KNOWN"; gene_name "RNU6-1100P"; transcript_type "snRNA"; transcript_status "KNOWN"; transcript_name "RNU6-1100P-201"; level 3; tag "basic"; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                                                                                                       |
| chr1 | ENSEMBL | exon       | 157784 | 157887 | . | - | . | gene_id "ENSG00000222623.1"; transcript_id "ENST00000410691.1"; gene_type "snRNA"; gene_status "KNOWN"; gene_name "RNU6-1100P"; transcript_type "snRNA"; transcript_status "KNOWN"; transcript_name "RNU6-1100P-201"; exon_number 1; exon_id "ENSE00001807945.1"; level 3; tag "basic"; remap_status "automatic_gene"; remap_substituted_missing_target "V19";                                                                           |                                                                     |

### Output file descriptions

`splice_lib_events.gtf`

| chr3 | splice_lib_event | exon | 52321872  | 52321892  | . | + | . | gene_id "SE.0006314"; transcript_id "SE.0006314_included"; |
| chr3 | splice_lib_event | exon | 52324320  | 52324735  | . | + | . | gene_id "SE.0006314"; transcript_id "SE.0006314_included"; |
| chr3 | splice_lib_event | exon | 52324976  | 52325127  | . | + | . | gene_id "SE.0006314"; transcript_id "SE.0006314_included"; |
| chr3 | splice_lib_event | exon | 52321872  | 52321892  | . | + | . | gene_id "SE.0006314"; transcript_id "SE.0006314_excluded"; |
| chr3 | splice_lib_event | exon | 52324976  | 52325127  | . | + | . | gene_id "SE.0006314"; transcript_id "SE.0006314_excluded"; |
| chr3 | splice_lib_event | exon | 42605015  | 42605178  | . | - | . | gene_id "SE.0003301"; transcript_id "SE.0003301_included"; |
| chr3 | splice_lib_event | exon | 42605257  | 42605472  | . | - | . | gene_id "SE.0003301"; transcript_id "SE.0003301_included"; |
| chr3 | splice_lib_event | exon | 42610357  | 42610565  | . | - | . | gene_id "SE.0003301"; transcript_id "SE.0003301_included"; |
| chr3 | splice_lib_event | exon | 42605015  | 42605178  | . | - | . | gene_id "SE.0003301"; transcript_id "SE.0003301_excluded"; |
| chr3 | splice_lib_event | exon | 42610357  | 42610565  | . | - | . | gene_id "SE.0003301"; transcript_id "SE.0003301_excluded"; |
| chr9 | splice_lib_event | exon | 128419930 | 128420078 | . | - | . | gene_id "MS.0000307"; transcript_id "MS.0000307_included"; |
| chr9 | splice_lib_event | exon | 128432097 | 128432186 | . | - | . | gene_id "MS.0000307"; transcript_id "MS.0000307_included"; |
| chr9 | splice_lib_event | exon | 128434595 | 128434922 | . | - | . | gene_id "MS.0000307"; transcript_id "MS.0000307_included"; |
| chr9 | splice_lib_event | exon | 128469250 | 128469482 | . | - | . | gene_id "MS.0000307"; transcript_id "MS.0000307_included"; |
| chr9 | splice_lib_event | exon | 128419930 | 128420078 | . | - | . | gene_id "MS.0000307"; transcript_id "MS.0000307_excluded"; |
| chr9 | splice_lib_event | exon | 128469250 | 128469482 | . | - | . | gene_id "MS.0000307"; transcript_id "MS.0000307_excluded"; |

`splice_lib_events.bed`

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

| chr3 | splice_lib_event | gene       | 52321872  | 52325127  | . | + | . | ID=SE.0006314;Name=SE.0006314                            |
| chr3 | splice_lib_event | transcript | 52321872  | 52325127  | . | + | . | ID=SE.0006314_included;Parent=SE.0006314                 |
| chr3 | splice_lib_event | exon       | 52321872  | 52321892  | . | + | . | ID=exon:SE.0006314_included:1;Parent=SE.0006314_included |
| chr3 | splice_lib_event | exon       | 52324320  | 52324735  | . | + | . | ID=exon:SE.0006314_included:2;Parent=SE.0006314_included |
| chr3 | splice_lib_event | exon       | 52324976  | 52325127  | . | + | . | ID=exon:SE.0006314_included:3;Parent=SE.0006314_included |
| chr3 | splice_lib_event | transcript | 52321872  | 52325127  | . | + | . | ID=SE.0006314_excluded;Parent=SE.0006314                 |
| chr3 | splice_lib_event | exon       | 52321872  | 52321892  | . | + | . | ID=exon:SE.0006314_excluded:1;Parent=SE.0006314_excluded |
| chr3 | splice_lib_event | exon       | 52324976  | 52325127  | . | + | . | ID=exon:SE.0006314_excluded:2;Parent=SE.0006314_excluded |
| chr3 | splice_lib_event | gene       | 42605015  | 42610565  | . | - | . | ID=SE.0003301;Name=SE.0003301                            |
| chr3 | splice_lib_event | transcript | 42605015  | 42610565  | . | - | . | ID=SE.0003301_included;Parent=SE.0003301                 |
| chr3 | splice_lib_event | exon       | 42605015  | 42605178  | . | - | . | ID=exon:SE.0003301_included:1;Parent=SE.0003301_included |
| chr3 | splice_lib_event | exon       | 42605257  | 42605472  | . | - | . | ID=exon:SE.0003301_included:2;Parent=SE.0003301_included |
| chr3 | splice_lib_event | exon       | 42610357  | 42610565  | . | - | . | ID=exon:SE.0003301_included:3;Parent=SE.0003301_included |
| chr3 | splice_lib_event | transcript | 42605015  | 42610565  | . | - | . | ID=SE.0003301_excluded;Parent=SE.0003301                 |
| chr3 | splice_lib_event | exon       | 42605015  | 42605178  | . | - | . | ID=exon:SE.0003301_excluded:1;Parent=SE.0003301_excluded |
| chr3 | splice_lib_event | exon       | 42610357  | 42610565  | . | - | . | ID=exon:SE.0003301_excluded:2;Parent=SE.0003301_excluded |
| chr9 | splice_lib_event | gene       | 128419930 | 128469482 | . | - | . | ID=MS.0000307;Name=MS.0000307                            |
| chr9 | splice_lib_event | transcript | 128419930 | 128469482 | . | - | . | ID=MS.0000307_included;Parent=MS.0000307                 |
| chr9 | splice_lib_event | exon       | 128419930 | 128420078 | . | - | . | ID=exon:MS.0000307_included:1;Parent=MS.0000307_included |
| chr9 | splice_lib_event | exon       | 128432097 | 128432186 | . | - | . | ID=exon:MS.0000307_included:2;Parent=MS.0000307_included |
| chr9 | splice_lib_event | exon       | 128434595 | 128434922 | . | - | . | ID=exon:MS.0000307_included:3;Parent=MS.0000307_included |
| chr9 | splice_lib_event | exon       | 128469250 | 128469482 | . | - | . | ID=exon:MS.0000307_included:4;Parent=MS.0000307_included |
| chr9 | splice_lib_event | transcript | 128419930 | 128469482 | . | - | . | ID=MS.0000307_excluded;Parent=MS.0000307                 |
| chr9 | splice_lib_event | exon       | 128419930 | 128420078 | . | - | . | ID=exon:MS.0000307_excluded:1;Parent=MS.0000307_excluded |
| chr9 | splice_lib_event | exon       | 128469250 | 128469482 | . | - | . | ID=exon:MS.0000307_excluded:2;Parent=MS.0000307_excluded |

`junctioncounts_unfriendly_events.gtf`

`junctioncounts_unfriendly_events.bed`

`junctioncounts_unfriendly_events.ioe`

`junctioncounts_unfriendly_events.gff3`

## Help statement

usage: infer_pairwise_events.py [-h] [--transcript_gtf TRANSCRIPT_GTF]  
                                --outdir OUTDIR  
                                [--min_exon_length MIN_EXON_LENGTH]  
                                [--min_intron_length MIN_INTRON_LENGTH]  
                                [--max_exon_length MAX_EXON_LENGTH]  
                                [--max_intron_length MAX_INTRON_LENGTH]  
                                [--dump_pkl_file]  
                                [--bedtools_path BEDTOOLS_PATH]  
                                [--suppress_output]  

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

# junctionCounts.py

## Description

## Basic usage

```
python2.7 /path/to/junctionCounts.py --event_gtf /path/to/splice_lib_events.gtf --bam /path/to/input_file.bam --outdir /path/to/output/dir/ --sample_name sample_name
```

### Output file descriptions

`SRR536342_count_psi_outfile.tsv`

| sample_name | event_id   | event_type | min_IJC | min_SJC | IncFormLen | SkipFormLen | PSI             | max_gene_frac | all_IJC                    | all_SJC | span_PSI | PSI_lo | PSI_hi | bootstrap_num |
|-------------|------------|------------|---------|---------|------------|-------------|-----------------|---------------|----------------------------|---------|----------|--------|--------|---------------|
| SRR536342   | ML.0005073 | ML         | 0       | 4       | 1          | 1           | 0.941558441558  | 0.0286        | 124,81,1,17,0,87,105,95,70 | 4       | 0.9688   | 0      | 0.9688 | NA            |
| SRR536342   | RI.0002381 | RI         | 130     | 678     | 1          | 1           | 0.164510166359  | 1             | 130137                     | 678     | 0.0072   | 0.1609 | 0.1681 | NA            |
| SRR536342   | CF.0001650 | CF         | 24      | 0       | 1          | 1           | 1               | 0.1043        | 31,185,24                  | 0       | 0        | 1      | 1      | NA            |
| SRR536342   | SE.0009399 | SE         | 75      | 142     | 1          | 1           | 0.351598173516  | 0.4897        | 75,79                      | 142     | 0.0118   | 0.3456 | 0.3575 | NA            |
| SRR536342   | MF.0012942 | MF         | 0       | 36      | 1          | 1           | 0.0886075949367 | 1             | 7,0                        | 36      | 0.1628   | 0      | 0.1628 | NA            |
| SRR536342   | MX.0000766 | MX         | 0       | 0       | 1          | 1           | NA              | 0             | 0,0                        | 0,0     | NA       | NA     | NA     | NA            |
| SRR536342   | CL.0003508 | CL         | 0       | 36      | 1          | 1           | 0.389830508475  | 0.12          | 46,0                       | 36      | 0.561    | 0      | 0.561  | NA            |
| SRR536342   | CO.0001190 | CO         | 0       | 0       | 1          | 1           | 1               | 0             | 0,0,2                      | 0,0     | NA       | NA     | NA     | NA            |
| SRR536342   | CL.0002855 | CL         | 177     | 36      | 1          | 1           | 0.830985915493  | 0.594         | 177                        | 36      | 0        | 0.831  | 0.831  | NA            |
| SRR536342   | CL.0002857 | CL         | 0       | 1       | 1          | 1           | 0               | 0.0135        | 0                          | 1       | 0        | 0      | 0      | NA            |
| SRR536342   | ML.0009767 | ML         | 16      | 39      | 1          | 1           | 0.285714285714  | 0.9512        | 16                         | 41,39   | 0.0102   | 0.2807 | 0.2909 | NA            |
| SRR536342   | CL.0009753 | CL         | 14      | 0       | 1          | 1           | 1               | 0.4667        | 17,18,14                   | 0       | 0        | 1      | 1      | NA            |
| SRR536342   | SE.0005423 | SE         | 8       | 6       | 1          | 1           | 0.720930232558  | 0.0606        | 8,23                       | 6       | 0.2217   | 0.5714 | 0.7931 | NA            |
| SRR536342   | AF.0003423 | AF         | 12      | 60      | 1          | 1           | 0.166666666667  | 0.625         | 12                         | 60      | 0        | 0.1667 | 0.1667 | NA            |
| SRR536342   | SE.0005931 | SE         | 0       | 0       | 1          | 1           | NA              | 0             | 0,0                        | 0       | NA       | NA     | NA     | NA            |
| SRR536342   | SE.0005394 | SE         | 0       | 0       | 1          | 1           | 1               | 0             | 275,0                      | 0       | NA       | NA     | NA     | NA            |
| SRR536342   | MF.0006829 | MF         | 0       | 11      | 1          | 1           | 0.12            | 0.046         | 0,3                        | 11      | 0.2143   | 0      | 0.2143 | NA            |


### Input file descriptions

## Benchmarking

### RT-PCR

![Screenshot](./figures/combined_jc_vs_rtpcr_dpsi.png)
![Screenshot](./figures/combined_jc_vs_rtpcr_psi.png)

## Help statement
