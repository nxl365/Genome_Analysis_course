## 1)	download the count results of HTSeq from uppmax to local:

scp -r nanxi@rackham.UPPMAX.uu.se:/proj/genomeanalysis2022/nobackup/work/nanxi/analyses/15_readcount_HTSeq .


2)R studio

# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.14")

browseVignettes("DESeq2") #we can get a intro link from this

install.packages("tidyverse")

library(DESeq2)
library(tidyverse)


# read htseq count file
# search the SRR number with NCBI, get the simple name of these rna

count_htseq_SRR6040092 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040092_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: leaf
count_htseq_SRR6040093 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040093_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: root
count_htseq_SRR6040096 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040096_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: stem

count_htseq_SRR6040094 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040094_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: aril 2
count_htseq_SRR6040095 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040095_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: aril 1
count_htseq_SRR6040097 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040097_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: aril 3

count_htseq_SRR6156066 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6156066_scaffold_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Monthong: aril 2
count_htseq_SRR6156067 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6156067_scaffold_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Monthong: aril 3
count_htseq_SRR6156069 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6156069_scaffold_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Monthong: aril 1
