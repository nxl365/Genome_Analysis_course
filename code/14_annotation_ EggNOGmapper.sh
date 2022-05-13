#!/bin/bash -l

#SBATCH -A uppmax2022-2-5
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -J eggnog_11
#SBATCH --mail-type=ALL
#SBATCH --mail-user nanxing.liu.3019@student.uu.se

module load bioinfo-tools
module load eggNOG-mapper

emapper.py -i /home/nanxi/Genome_anaysis/analyses/09_DNA_masking/DNA_pilon_out.fasta.masked --itype genome \
-o /proj/genomeanalysis2022/nobackup/work/nanxi/analyses/14_annotation_eggnogmapper/eggnogmapper_out \
--cpu 8
