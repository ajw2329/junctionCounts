source ~/.bash_profile

num_threads=6

export star_exe=/public/home/anjowall/STAR/bin/Linux_x86_64/STAR
export samtools_exe=/public/home/anjowall/samtools-1.8/bin/samtools

###These dirs must exist
export topdir=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/
export human_fastq_dir="$topdir"/fastq/
export human_genome_dir=/public/groups/sanfordlab/people/anjowall/genomes/GRCm38/STAR_genome_GRCm38.p6/
export human_genome_fasta=/public/groups/sanfordlab/people/anjowall/genomes/GRCm38/STAR_genome_GRCm38.p6/GRCm38.primary_assembly.genome.fa

###These dirs/files can be created
export human_bam_dir_pass_one="$topdir"/bam/

export human_fastq_prefix_list=$(ls "$human_fastq_dir"/*fastq.gz | awk -F\/ '{ print $NF}' | awk -F"_[12].fastq" '{ print $1 }' | sort | uniq)


mkdir -p "$human_bam_dir_pass_one"


#####STAR indexing
##Human

if [ ! -e "$human_genome_dir"/SA ]
        then
                echo "genome suffix array file not found . . . running STAR genome generate . . . "
                $star_exe --runMode genomeGenerate --genomeDir "$human_genome_dir" --genomeFastaFiles "$human_genome_fasta" --runThreadN "$num_threads" --genomeSAindexNbases 14 --limitGenomeGenerateRAM=140000000000 
else
        echo "Genome suffix array file found - assuming genome is indexed - attempting to proceed . . . "
fi


$star_exe --genomeDir $human_genome_dir --genomeLoad Remove
$star_exe --genomeDir $human_genome_dir --genomeLoad LoadAndExit


#####STAR mapping
###Human
echo "Attempting to map with STAR"

for prefix in ${human_fastq_prefix_list[*]}
        do
                if [ ! -e "$human_bam_dir_pass_one"/"$prefix""_pass_1_""Aligned.sortedByCoord.out.bam" ]
                        then
                                echo "Mapping sample  . . . "
                                $star_exe --runMode alignReads --readFilesIn "$human_fastq_dir""$prefix"_1.fastq.gz "$human_fastq_dir""$prefix"_2.fastq.gz --readFilesCommand zcat --genomeDir "$human_genome_dir" --outSAMtype BAM SortedByCoordinate --runThreadN "$num_threads" --alignEndsType Local --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --chimScoreMin 1 --outSAMattributes NH HI AS nM --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix "$human_bam_dir_pass_one"/"$prefix""_pass_1_" --outFilterMultimapNmax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM=60688817693
                        echo "Sample already mapped . . .  proceeding . . . "
                fi
        done


$star_exe --genomeDir $human_genome_dir --genomeLoad Remove
