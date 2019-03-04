
bam_dir=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/bam/
ref_path=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/gencode.v29lift37.basic.annotation.gtf
stringtie_out=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/shen_2014/stringtie/

for i in $bam_dir/*.bam
do
stringtie "$i" -G $ref_path -o $stringtie_out/$(echo "$i" | awk -F\/ '{ print $NF}').gtf & done

wait; stringtie --merge -G $ref_path -o $stringtie_out/all_merged.gtf $stringtie_out/*.gtf

