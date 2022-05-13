#!/bin/bash -l

#SBATCH -A uppmax2022-2-5
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J braker_11
#SBATCH --mail-type=ALL
#SBATCH --mail-user nanxing.liu.3019@student.uu.se

#BRAKER needs the following modules loaded to be able run:
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

#GeneMark key needs to be copied to the home directory
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

#AUGUSTUS needs writing access to the configuration file:
source $AUGUSTUS_CONFIG_COPY
chmod a+w -R /home/nanxi/Genome_anaysis/code/12_annotation_braker/augustus_config/species/



#BRAKER needs the following uppercase flags of the paths to the softwares
#genome is masked, use the --softmasking

braker.pl \
--AUGUSTUS_CONFIG_PATH=/home/nanxi/Genome_anaysis/code/12_annotation_braker/augustus_config \
--AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin \
--AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts/ \
--GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy \
--species=Durio_zibethinus \
--useexisting \
--genome=/home/nanxi/Genome_anaysis/analyses/09_DNA_masking/DNA_pilon_out.fasta.masked \
--bam=11_merge_bam/merge/map_sort_merge.bam \
--workingdir=/proj/genomeanalysis2022/nobackup/work/nanxi/analyses/12_annotation_braker \
--softmasking \
--cores=8 \
--gff3
