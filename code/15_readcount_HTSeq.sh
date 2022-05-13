#!/bin/bash -l

#SBATCH -A uppmax2022-2-5
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J htseq_11
#SBATCH --mail-type=ALL
#SBATCH --mail-user nanxing.liu.3019@student.uu.se

module load bioinfo-tools
module load samtools
module load htseq

#index: have used samtool to sort at step 11 before , so use the sorted directly to index
#add index in sorted.file, so the indexed.file also in the sorted_path

sorted_path=/proj/genomeanalysis2022/nobackup/work/nanxi/analyses/11_merge_bam/sort

for file in $sorted_path/*.bam
do
	samtools index $file
done


#htseq: use indexed.file, braker.file(gff)

annotation_braker_path=/proj/genomeanalysis2022/nobackup/work/nanxi/analyses/12_annotation_braker
htseq_path=/proj/genomeanalysis2022/nobackup/work/nanxi/analyses/15_readcount_HTSeq

for file in $sorted_path/*.bam
do
	basename_file=$(basename "$file" .bam)

	htseq-count -f bam -t CDS -r pos -i ID --stranded=no \
	$file $annotation_braker_path/augustus.hints.gff3 > $htseq_path/count_htseq_$basename_file.txt
done


#  $(basename "$file" .bam)： extract the basename of $file(*.bam), before the .bam 
