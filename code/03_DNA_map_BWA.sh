#  script in:   /home/nanxi/Genome_anaysis/code/ 03_DNA_map_BWA.sh: 



#!/bin/bash -l
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2                       # book 2 notes
#SBATCH -t 03:00:00
#SBATCH -J DNA_BWA_11
#SBATCH --mail-type=ALL
#SBATCH --mail-user nanxing.liu.3019@student.uu.se

module load bioinfo-tools
module load bwa/0.7.17

bwa index /home/nanxi/Genome_anaysis/analyses/01_DNA_assembly/DNA_pacbio_assembly_11.contigs.fasta

bwa mem -t 2 \                      #  [ -t nThreads ]   use 4 notes
/home/nanxi/Genome_anaysis/analyses/01_DNA_assembly/DNA_pacbio_assembly_11.contigs.fasta \
/home/nanxi/Genome_anaysis/data/raw_data/4_Tean_Teh_2017/illumina_data/SRR6058604_scaffold_11.1P.fastq.gz \
/home/nanxi/Genome_anaysis/data/raw_data/4_Tean_Teh_2017/illumina_data/SRR6058604_scaffold_11.2P.fastq.gz \
> /home/nanxi/Genome_anaysis/analyses/03_DNA_map/DNA_bwa_align.sam
