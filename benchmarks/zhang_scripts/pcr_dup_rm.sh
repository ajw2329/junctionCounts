source ~/.bash_profile

num_threads=$1

export bam_dir=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/bam/

export star_exe=/public/home/anjowall/STAR/bin/Linux_x86_64/STAR
export samtools_exe=/public/home/anjowall/samtools-1.8/samtools


###These are just used for naming
export fastq_dir=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014//fastq/
export fastq_prefix_list=$(ls "$fastq_dir"/*fastq.gz | awk -F\/ '{ print $NF}' | awk -F"_[12].fastq" '{ print $1 }' | sort | uniq)
###

bam_dir_rmdup=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/bam_rmdup/

mkdir -p "$bam_dir_rmdup"


run_STAR_rmdup(){
	prefix="$1"
	in_bam_dir="$2"
	out_bam_dir="$3"
	marked_dup_file="$in_bam_dir"/"$prefix""_""Processed.out.bam"
	removed_dup_file="$out_bam_dir"/"$prefix""_""dup_removed.bam"
	removed_dup_file_index="$out_bam_dir"/"$prefix""_""dup_removed.bam.bai"

	if [ ! -e "$marked_dup_file" ]
		then
			echo "Dup marked bam file not find - running STAR for dup removal "
			$star_exe --limitBAMsortRAM 10000000000 --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType UniqueIdentical --inputBAMfile "$in_bam_dir"/"$prefix""_pass_1_""Aligned.sortedByCoord.out.bam" --outFileNamePrefix "$in_bam_dir"/"$prefix""_" --outSAMtype BAM SortedByCoordinate

	else
		echo "Dup marked file found. Continuing . . . "
	fi

	if [ ! -e "$removed_dup_file" ]
		then
			echo "Dup removed file not found.  Running samtools for dup removal."
			$samtools_exe view -hb -F 0x400 "$marked_dup_file" > "$removed_dup_file"
	else
		echo "Dup removed file found. Continuing . . . "
	fi
	if [ ! -e "$removed_dup_file_index" ]
		then
		echo "Dup removed file index not found. Running samtools index."
		$samtools_exe index "$removed_dup_file"
	else
		echo "Dup removed file index found. Continuing . . ."
	fi
}


i=""
(
for prefix in ${fastq_prefix_list[*]}
	do
		((i=i%num_threads)); ((i++==0)) && wait
		run_STAR_rmdup "$prefix" "$bam_dir" "$bam_dir_rmdup" &
done
wait
) 


