#!/bin/bash -l

#SBATCH -A uppmax2022-2-5
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J trinity_11
#SBATCH --mail-type=ALL
#SBATCH --mail-user nanxing.liu.3019@student.uu.se

module load bioinfo-tools
module load samtools 
module load trinity/2.4.0

# use merged.bam after mapping(star)
Trinity --genome_guided_bam /proj/genomeanalysis2022/nobackup/work/nanxi/analyses/11_merge_bam/merge/map_sort_merge.bam \
 --max_memory 100G --genome_guided_max_intron 10000 --CPU 8 \
--output /proj/genomeanalysis2022/nobackup/work/nanxi/analyses/13_RNA_assembly_trinity/trinity_output


# use sorted.bam after mapping(star)

//cd /proj/genomeanalysis2022/nobackup/work/nanxi/analyses/11_merge_bam/sort

// for file in *_sort.bam:
// do
//	Trinity --genome_guided_bam $file ----max_memory 100G --genome_guided_max_intron 10000 --CPU 8 \
//	--output /proj/genomeanalysis2022/nobackup/work/nanxi/analyses/13_RNA_assembly_trinity/${file}_trinity_out
//done

