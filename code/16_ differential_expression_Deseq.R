# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")

browseVignettes("DESeq2") #we can get a intro link from this

install.packages("tidyverse")
install.packages('pheatmap')

# load library

library(DESeq2)
library(pheatmap)
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(genefilter)


# 1. read htseq count files
# and search the SRR number with NCBI, get the simple name of these rna
count_htseq_SRR6040092 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040092_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: leaf
count_htseq_SRR6040093 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040093_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: root
count_htseq_SRR6040096 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040096_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: stem

count_htseq_SRR6040094 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040094_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: aril 2
count_htseq_SRR6040095 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040095_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: aril 1
count_htseq_SRR6040097 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6040097_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Musang King: aril 3

count_htseq_SRR6156066 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6156066_scaffold_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Monthong: aril 2
count_htseq_SRR6156067 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6156067_scaffold_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Monthong: aril 3
count_htseq_SRR6156069 <- read.delim("15_readcount_HTSeq/count_htseq_SRR6156069_scaffold_11_Aligned.out.bam_sort.txt", header=FALSE) #RNA-seq of Durio zibethinus Monthong: aril 1

# 2. merge all the htseq data of files to one table, by shared column names"V1"(genes)
count_rna <- merge(count_htseq_SRR6040092,count_htseq_SRR6040093,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6040096,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6040094,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6040095,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6040097,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6156066,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6156067,by="V1")
count_rna <- merge(count_rna,count_htseq_SRR6156069,by="V1")

# 3. name columns
names(count_rna)<-c("gene","MK_leaf","MK_root","MK_stem","MK_aril2","MK_aril1","MK_aril3","MT_aril2","MT_aril3","MT_aril1")

#4. remove the first 5 unuseful rows,keep all columns
count_rna <- count_rna[-c(1:5),]

#5. creat two dataframes: `count_rna`& `metadata` 

  # make a new table which drop the first column (gene) of `Count_rna`
for_col_names <- count_rna[,-1] 

  # create the conditions(non-fruit & fruit) for deseq, PCA
condition <- c("Non-fruit","Non-fruit","Non-fruit","MK-fruit","MK-fruit","MK-fruit","MT-fruit","MT-fruit","MT-fruit")

  # create a metadata (data frame) with row name as `count_rna`'s column name(except "gene"), and add conditions to the frame
metadata <- data.frame(row.names=colnames(for_col_names), condition) 

  #transfer the count_rna to a data frame, with row name as gene name, column name as columnnames of count_rna
count_rna_f <- data.frame(row.names=count_rna[,1], count_rna[2], 
                          count_rna[3], count_rna[4], 
                          count_rna[5], count_rna[6], 
                          count_rna[7],count_rna[8],count_rna[9], count_rna[10])


# 6. do Deseq2
 # build dds matrix
dds_table <- DESeqDataSetFromMatrix(countData=count_rna_f, colData=metadata, design= ~ condition)
 # run analysis
dds <- DESeq(dds_table)

# 7. do PCA (Principal component analysis)
   # count_rna has done normalization from Htseq
  
rld<-rlog(dds)  # Transform data (regularized-logarithm transformation)
colData(dds)
plotPCA(rld,intgroup=c("condition","sizeFactor"))

# 8. Volcano plot (MA)
res <- results(dds)
plotMA(res)

# 9. Ordering for the clusters of heatmap
head(res)
summary(res)
   # order the dds result by padj(adjusted p-value),smaller P is more significant
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
head(resOrdered)
   # delete NA values
resOrdered=na.omit(resOrdered)  
   # output result into file
write.csv(resOrdered,"resOrdered.csv")

# 10. heatmap
order_select<- head(order(res$padj), 10)  #for top 10 smallest p-val genes (most significiant)
nt<-normTransform(dds)
log2.counts<-assay(nt)[order_select,]
dataf<-as.data.frame(colData(dds))

pheatmap(log2.counts,cluster_rows = TRUE,show_rownames = TRUE,cluster_cols = TRUE,annotation_col = dataf)


## another way of heatmap
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 10)

heatmap.2( assay(rld)[topVarGenes, ], scale = "row",
           trace="none", dendrogram="column",
           col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), margins = c(6,10))
