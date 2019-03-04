
n_threads=10
topdir=/public/groups/sanfordlab/people/anjowall/projects/rtpcr_validated/zhang_2014/
fastq_dir=$topdir/fastq/
kallisto_dir=$topdir/kallisto/
kallisto_idx=$topdir/gencode.vM20.basic.annotation.idx

export fastq_prefix_list=$(ls "$fastq_dir"/*fastq.gz | awk -F\/ '{ print $NF}' | awk -F"_[12].fastq" '{ print $1 }' | sort | uniq)

#kallisto index -i $kallisto_idx $topdir/gencode.vM20.basic.annotation.fa

mkdir -p $kallisto_dir

for prefix in ${fastq_prefix_list[*]}
do
subdir=$kallisto_dir/$prefix
mkdir -p $subdir
kallisto quant -i $kallisto_idx -o $subdir --rf-stranded -t $n_threads "$fastq_dir""$prefix"_1.fastq.gz "$fastq_dir""$prefix"_2.fastq.gz
done


