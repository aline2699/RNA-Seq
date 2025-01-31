library(DESeq2)
library(ggplot2)
library(pheatmap)


#----------------------------------------------
# Step 5: Exploratory data analysis
#----------------------------------------------

# Read the precleaned counts-table into R
cleaned_counts <- read.csv("cleaned_counts.csv", header = TRUE, sep=";", row.names = 1)
# Remove NA's from the counts-table
cleaned_counts <- na.omit(cleaned_counts)
# Make a dataframe that has the column names of the counts-table as rownames and assignes a group to each of the samples
condition_frame <- data.frame(group = rep(c("HER2", "NonTNBC", "Normal", "TNBC"), each=3), row.names = colnames(cleaned_counts))

# Create the DESeqDataSet object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cleaned_counts, colData = condition_frame, design = ~group)
dds <- DESeq2::DESeq(dds)
# Remove the dependency that the variance has on the mean
dds_vst <- DESeq2::vst(dds, blind=TRUE)

# PLOT PCA plot

# Make a PCA plot to see how the gene expression profiles of the samples cluster
DESeq2::plotPCA(dds_vst, intgroup = "group")
# factor extra


library("pheatmap")
# Find the indices of the top 20 most highly expressed genes in the dds object after the data has been normalized
# 1. counts: extracts a matrix of the normalized counts from the dds data object
# 2: rowMeans: calculates the mean of the normalized counts across samples for each gene. This gives a measure of the average expression level for each gene across all samples
# 3: order[1:20]: order the normalized counts in a decreasing order, then pick the first 10
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:10]

# Make a normalization transfomation on the data object dds
ntd <- normTransform(dds)

# PLOT heatmat

# Plot the heatmap of the count matrix, when the data is normalization transformed
# 1: assay: accesses the data of ntd
# 2: [select,]: we take the 20 most expressed genes of the ntd set
# 3: cluster_rows: whether or not the tree thing should be shown on the left side
# 4: show_rownames: whether or not the gene names should be shown on the right side
# 5: cluster_cols: whether or not the tree thing should be shown on top
# 6: annotation_col: dataframe that contains the sample names that will be displayed on the bottom
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=condition_frame)
# Plot the heatmap of the count matrix, when the the dependency of the variance on the mean has been removed from the data
pheatmap(assay(dds_vst)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=condition_frame)




library("RColorBrewer")
# Create a matrix that contains the distances between the samples by comparing the values of the samples
# 1. assay: Extract the variance-stabilized values of the samples from the DESeqTransform object dds_vst.
# 2. t: transpose the matrix-like entries, so that the samples will be the rows
# 3. dist: computes the Euclidean distance (measure of the straight-line) between samples in the expression space of genes
sampleDists <- dist(t(assay(dds_vst)))
# Because dist doesn't output a matrix, but just a bunch of distance values, the sampleDists object has to be transformed into a matrix
sampleDistMatrix <- as.matrix(sampleDists)
# Give names to the rows and columns, so the heatmap will have the indication of which group is represented
rownames(sampleDistMatrix) <- dds_vst$group
colnames(sampleDistMatrix) <- dds_vst$group

# Choose colorpalette
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# PLOT heatmap sample vs sample
# Plot the heatmap of the sample-to-sample distances
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



#----------------------------------------------
# Step 6: Differential expression analysis
#----------------------------------------------

# Normal vs HER2
Results_Normal_vs_HER2 <- DESeq2::results(dds, contrast = c("group","HER2","Normal"))
N_vs_H2_sign <- cleaned_counts[which(Results_Normal_vs_HER2$padj<0.05),]
N_vs_H2_Hupreg <- cleaned_counts[which(Results_Normal_vs_HER2$padj<0.05 & Results_Normal_vs_HER2$log2FoldChange>0),]
N_vs_H2_Hdownreg <- cleaned_counts[which(Results_Normal_vs_HER2$padj<0.05 & Results_Normal_vs_HER2$log2FoldChange<0),]


# Normal vs TNBC
Results_Normal_vs_TNBC <- DESeq2::results(dds, contrast = c("group","TNBC","Normal"))
N_vs_T_sign <- cleaned_counts[which(Results_Normal_vs_TNBC$padj<0.05),]
N_vs_T_Tupreg <- cleaned_counts[which(Results_Normal_vs_TNBC$padj<0.05 & Results_Normal_vs_TNBC$log2FoldChange>0),]
N_vs_T_Tdownreg <- cleaned_counts[which(Results_Normal_vs_TNBC$padj<0.05 & Results_Normal_vs_TNBC$log2FoldChange<0),]

# Normal vs NonTNBC
Results_Normal_vs_NonTNBC <- DESeq2::results(dds, contrast = c("group","NonTNBC","Normal"))
N_vs_NT_sign <- cleaned_counts[which(Results_Normal_vs_NonTNBC$padj<0.05),]
N_vs_NT_NTupreg <- cleaned_counts[which(Results_Normal_vs_NonTNBC$padj<0.05 & Results_Normal_vs_NonTNBC$log2FoldChange>0),]
N_vs_NT_NTdownreg <- cleaned_counts[which(Results_Normal_vs_NonTNBC$padj<0.05 & Results_Normal_vs_NonTNBC$log2FoldChange<0),]

# TNBC vs NonTNBC
Results_TNBC_vs_NonTNBC <- DESeq2::results(dds, contrast = c("group","NonTNBC","TNBC"))
T_vs_NT_sign <- cleaned_counts[which(Results_TNBC_vs_NonTNBC$padj<0.05),]
T_vs_NT_NTupreg <- cleaned_counts[which(Results_TNBC_vs_NonTNBC$padj<0.05 & Results_TNBC_vs_NonTNBC$log2FoldChange>0),]
T_vs_NT_NTdownreg <- cleaned_counts[which(Results_TNBC_vs_NonTNBC$padj<0.05 & Results_TNBC_vs_NonTNBC$log2FoldChange<0),]

a <- Results_TNBC_vs_NonTNBC[which(Results_TNBC_vs_NonTNBC$padj<0.05 & Results_TNBC_vs_NonTNBC$log2FoldChange<0),]

# TNBC vs HER2
Results_TNBC_vs_HER2 <- DESeq2::results(dds, contrast = c("group","HER2","TNBC"))
T_vs_H_sign <- cleaned_counts[which(Results_TNBC_vs_HER2$padj<0.05),]
T_vs_H_Hupreg <- cleaned_counts[which(Results_TNBC_vs_HER2$padj<0.05 & Results_TNBC_vs_HER2$log2FoldChange>0),]
T_vs_H_Hdownreg <- cleaned_counts[which(Results_TNBC_vs_HER2$padj<0.05 & Results_TNBC_vs_HER2$log2FoldChange<0),]

# HER2 vs NonTNBC
Results_HER2_vs_NonTNBC <- DESeq2::results(dds, contrast = c("group","NonTNBC","HER2"))
H_vs_NT_sign <- cleaned_counts[which(Results_HER2_vs_NonTNBC$padj<0.05),]
H_vs_NT_NTupreg <- cleaned_counts[which(Results_HER2_vs_NonTNBC$padj<0.05 & Results_HER2_vs_NonTNBC$log2FoldChange>0),]
H_vs_NT_NTdownreg <- cleaned_counts[which(Results_HER2_vs_NonTNBC$padj<0.05 & Results_HER2_vs_NonTNBC$log2FoldChange<0),]

# Remove effect of between-sample differences in sequencing depth
normalized_dds <- DESeq2::counts(dds)


#----------------------------------------------
# Step 7: Overrepresentation analysis
#----------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

T_vs_NT <- clusterProfiler::enrichGO(gene = row.names(T_vs_NT_sign), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

# PLOT dotplot
library(enrichplot)
barplot(T_vs_NT, showCategory=10) 

library(dplyr)
mutate(T_vs_NT, qscore = -log(p.adjust, base=10))%>% 
  barplot(x="qscore")

dotplot(T_vs_NT, showCategory=10) + ggtitle("dotplot for Overrepresentation Analysis")


# Which genes/geneclusters are overrepresented in TNBC group when compared to the nonTNBC?
T_overrep_vs_NT <- clusterProfiler::enrichGO(gene = row.names(T_vs_NT_NTdownreg), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

library(enrichplot)
barplot(T_overrep_vs_NT, showCategory=10) 

library(dplyr)
mutate(T_overrep_vs_NT, qscore = -log(p.adjust, base=10))%>% 
  barplot(x="qscore")

dotplot(T_overrep_vs_NT, showCategory=10) + ggtitle("Biological Pathways: \nTNBC upregulated compared to NonTNBC")

# Which genes/geneclusters are overrepresented in NonTNBC group when compared to the TNBC?
NT_overrep_vs_T <- clusterProfiler::enrichGO(gene = row.names(T_vs_NT_NTupreg), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

dotplot(NT_overrep_vs_T, showCategory=10) + ggtitle("Biological Pathways: \nNonTNBC upregulated compared to TNBC")


# Which genes/geneclusters are overrepresented in HER2 group when compared to the TNBC?
H_overrep_vs_T <- clusterProfiler::enrichGO(gene = row.names(T_vs_H_Hupreg), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

dotplot(H_overrep_vs_T, showCategory=10) + ggtitle("Biological Pathways: \nHER2 upregulated compared to TNBC")

# Which genes/geneclusters are overrepresented in TNBC group when compared to the HER2?
T_overrep_vs_H <- clusterProfiler::enrichGO(gene = row.names(T_vs_H_Hdownreg), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

dotplot(T_overrep_vs_H, showCategory=10) + ggtitle("Biological Pathways: \nTNBC upregulated compared to HER2")




# Which genes/geneclusters are overrepresented in NonTNBC group when compared to the HER2?
NT_overrep_vs_H <- clusterProfiler::enrichGO(gene = row.names(H_vs_NT_NTupreg), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

dotplot(NT_overrep_vs_H, showCategory=10) + ggtitle("Biological Pathways: \nNonTNBC upregulated compared to HER2")

# Which genes/geneclusters are overrepresented in HER2 group when compared to the NonTNBC?
H_overrep_vs_NT <- clusterProfiler::enrichGO(gene = row.names(H_vs_NT_NTdownreg), universe = row.names(cleaned_counts), OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "ENSEMBL")

dotplot(H_overrep_vs_NT, showCategory=10) + ggtitle("Biological Pathways: \nHER2 upregulated compared to NonTNBC")




# --------------------------------------------------------------------------------
# Plots for different expression levels for specific genes

# Normalized counts table
a <- DESeq2::counts(dds)

# define names on x axis
c <- c("HER2", "NonTNBC", "Normal", "TNBC")
library(lattice)
library(ggplot2)

# plots for the genes that were picked
RAB21 <- a["ENSG00000080371",]
averagesRAB21 <- c(sum(RAB21[1:3])/3, sum(RAB21[4:6])/3, sum(RAB21[7:9])/3, sum(RAB21[10:12])/3)
dotplot(averagesRAB21~c, type="b", ylab="Average reads of RAB21 per group", xlab="Groups", main="Number of RAB21 reads per group", pch=19, cex=2, lwd=1, lty=1, col="darkblue")

TMEM219 <- a["ENSG00000149932",]
averagesTMEM219 <- c(sum(TMEM219[1:3])/3, sum(TMEM219[4:6])/3, sum(TMEM219[7:9])/3, sum(TMEM219[10:12])/3)
dotplot(averagesTMEM219~c, type="b", ylab="Average reads of TMEM219 per group", xlab="Groups", main="Number of TMEM219 reads per group", pch=19, cex=2, lwd=1, lty=1, col="darkblue")

RACK1 <- a["ENSG00000204628",]
averagesRACK1 <- c(sum(RACK1[1:3])/3, sum(RACK1[4:6])/3, sum(RACK1[7:9])/3, sum(RACK1[10:12])/3)
dotplot(averagesRACK1~c, type="b", ylab="Average reads of RACK1 per group", xlab="Groups", main="Number of RACK1 reads per group", pch=19, cex=2, lwd=1, lty=1, col="darkblue")

B2M <- a["ENSG00000166710",]
averagesB2M <- c(sum(B2M[1:3])/3, sum(B2M[4:6])/3, sum(B2M[7:9])/3, sum(B2M[10:12])/3)
dotplot(averagesB2M~c, type="b", ylab="Average reads of B2M per group", xlab="Groups", main="Number of B2M reads per group", pch=19, cex=2, lwd=1, lty=1, col="darkblue")

