#   script in /home/nanxi/Genome_anaysis/code/02_DNA_short_QC.sh:



#!/bin/bash -l
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J Dna_shore_QC_11
#SBATCH --mail-type=ALL
#SBATCH --mail-user nanxing.liu.3019@student.uu.se

module load bioinfo-tools FastQC

for x in /home/nanxi/Genome_anaysis/data/raw_data/4_Tean_Teh_2017/illumina_data/*_scaffold_11.*.fastq.gz
do
    fastqc $x -o /home/nanxi/Genome_anaysis/analyses/02_DNA_short_QC;
done
