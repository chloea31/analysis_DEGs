#===============================================================================

#-------------------------------------------------------------------------------
# Description: Analysis of differential gene expression between control and 
# disease from cancer dataset
#-------------------------------------------------------------------------------
# Data frame containing the read counts for each genes (obtained with 
# feature-count?)
# 
# Steps of the script:
# 1) Data Acquisition, formating and filtering
# 2) Exploratory analysis and visualisation with DESeq2 data transformation
# 3) Differential Expression with DESeq2
# 4) Final DE gene list extraction

#===============================================================================


###############################
### Import packages
###############################


library("BiocManager")
library("dplyr")
library("gridBase")
library("gridExtra")

library("DESeq2")
library(scales) # needed for oob parameter
library(viridis)

# libraries for plotting results
library("ggplot2")
library("NMF")
library("ggbeeswarm")
library("genefilter")
library("pheatmap")
library("beeswarm")


#============== Data Acquisition ==================

# import tables of read counts
## Reading a pre-formatting count matrix for all samples

setwd() # to complete

expression <- read.table("expression_counts.txt", 
                         header = TRUE, 
                         row.names = 1, 
                         check.names = FALSE)

View(expression)
dim(expression)
names(expression)
colnames(expression)
row.names(expression)
class(expression)

# transpose matrix
expression_T = t(expression)
View(expression_T)
dim(expression_T)
class(expression_T)
names(expression_T)
colnames(expression_T)

# re-order columns in the transpose matrix
col_order <- c("control_1", "control_2", "control_3", "disease_1", "disease_2", "disease_3")
ReadCounts <- expression_T[, col_order]
View(ReadCounts)
dim(ReadCounts)
class(ReadCounts)
names(ReadCounts)
colnames(ReadCounts)

#================= Data filtering =======================================================

# Dataset filtering: gene with very low counts across all libraries 
# provide little evidence for differential expression. Step required by the
# EdgeR package, DESeq2 has its own internal filtering process, done for 
# all the analysis for comparison)
table(rowSums(ReadCounts)==0) #count the number of genes expressed in none of the samples

# A good threshold can be chosen by identifying the CPM (Count Per 
# million calculated by the EdgeR package) 
# that corresponds to a count of 10
cpm(10, mean(colSums(ReadCounts))) #tiens en compte la taille des librairies en terme de reads mappes
# discard low counts entries (count must be sufficient in at least two 
# samples)

cutoff <- 0.8000456
#keep <- rowSums(cpm(ReadCounts)>cutoff) >= 2
keep <- rowSums(cpm(ReadCounts)>10) >= 2
summary(keep)
ReadCounts <- ReadCounts[keep,] # on garde uniquement les genes respectant la condition
# permet d'enlever les genes qui sont trop faiblement exprimes
dim(ReadCounts)
class(ReadCounts)
names(ReadCounts)
colnames(ReadCounts)

boxplot(ReadCounts,
        main = "Representation of the dispersion of the data from the counting matrix for the conditions: control and disease",
        xlab = "Samples",
        ylab = "Countings")

#-----------------------------------------------------------------------
#============   Differential Expression with DESeq2     ================
#-----------------------------------------------------------------------

# Warning : The DGE analysis has to be performed on the raw read counts (untransformed, not normalized
# for sequencing depth)

#============== Differential Expression Gene Analysis ==================

# In our data set, we want to identify the genes differentially between the different groups
# TODO : chose a reference group (first-level factor)


SampleInfo  <- data.frame(group = c("control", 
                                    "control", 
                                    "control", 
                                    "disease", 
                                    "disease", 
                                    "disease"),
                          row.names = names(ReadCounts))

ref_group <- "control"
colData(DESeq.ds)$group <- relevel(colData(DESeq.ds)$group, reference_group)

# in order to detect difference between control and siDDX5_17
DESeq.ds <- DESeqDataSetFromMatrix (countData = ReadCounts ,
                                    colData = SampleInfo ,
                                    design = ~ group) 

ref_group <- "control"
colData(DESeq.ds)$group <- relevel(colData(DESeq.ds)$group, ref_group)

# The DESeq function will performed three successive steps of analysis :
# sequencing depth normalization between the samples
# gene-wise dispersion estimation across all samples
# fits a negative  binomial  GLM  and applies  Wald  statistics  to each  gene
# estimation de la dispersion

DESeq.ds <- DESeq(DESeq.ds)
summary(DESeq.ds)
#plotDispEsts(DESeq.ds)

# Building the results table
DESeqResults <- results(DESeq.ds, contrast=c("group", "control", "disease"), pAdjustMethod = "BH")

summary(DESeqResults)
write.table(DESeqResults, "DESeq2_results.txt", quote=FALSE, row.names=TRUE, col.names=TRUE)

DESeqResDF <- as.data.frame(DESeqResults)

# Set a boolean column for significance
DESeqResDF$significant <- ifelse(DESeqResDF$padj < .05, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(DESeqResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Let's add some more detail
ggplot(DESeqResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)

# Read counts of the top gene and plot, mainly useful that there was no mistake when setting the constrasts (which groups to compare)
topGene <- rownames(DESeqResDF)[which.min(DESeqResults$padj)]
plotCounts(DESeq.ds, gene = topGene, intgroup="group")

# sort by adjusted p-value
DESeqResDF  <- DESeqResDF[order(DESeqResDF$padj), ]

DESeqResDF$significant
DESeqResDF$significant <- "NO"

DESeqResDF$significant[DESeqResDF$log2FoldChange > 0.6 & DESeqResDF$padj < 0.05] <- "UP"
DESeqResDF$significant[DESeqResDF$log2FoldChange < -0.6 & DESeqResDF$padj < 0.05] <- "DOWN"
p <- ggplot(data=DESeqResDF, aes(x=log2FoldChange, y=-log10(padj), col=significant)) + geom_point() + theme_minimal()
p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p1 <- p + scale_color_manual(values=c("blue", "black", "red"))
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p1 <- p + scale_colour_manual(values = mycolors)

DESeqResDF$absLogFC <- abs(DESeqResDF$log2FoldChange)
keep <- DESeqResDF$padj < 0.05
DESeqResDF.sens <- DESeqResDF[keep,]
DESeqResDF.sens <- DESeqResDF.sens[order(DESeqResDF.sens$padj), ]


# identify  genes  with  the  desired  adjusted p-value cut -off

#DGEgenes  <- rownames(subset(DESeqResults.sorted, padj < alpharisk))
DGEgenes  <- rownames(subset(DESeqResDF.sens, padj < .05))

cutoff <- 0.0001
#topresults.DESeq2 <- DESeqResults.sorted[which(DESeqResults.sorted$padj < cutoff), ]
topresults.DESeq2 <- DESeqResDF.sens[which(DESeqResDF.sens$padj < cutoff), ]

write.table(topresults.DESeq2, "topgenes.DESeq2.txt", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Extract genes significantly up-regulated or down-regulated, put a condition on the fold change
UP_topresults <- topresults.DESeq2[topresults.DESeq2$log2FoldChange>0,]
UP_topgenes <- rownames(UP_topresults)
write.table(UP_topresults, "up_regulated.txt", quote=FALSE, row.names = TRUE, col.names = TRUE)
# Extract genes significantly up-regulated in group1 vs group2
Down_topresults <- topresults.DESeq2[topresults.DESeq2$log2FoldChange<0,]
Down_topgenes <- rownames(Down_topresults)
write.table(Down_topresults, "down_regulated.txt", quote=FALSE, row.names = TRUE, col.names = TRUE)

liste <- DESeqResDF.sens[1:10,  c(2,5,6,7)]
write.table(liste, "results.tsv", quote=FALSE, row.names = TRUE, col.names = TRUE)

hist(topresults.DESeq2$log2FoldChange,
     main = "Histogram representing the dispersion of the log2 fold change of DESeq2 genes")


#----------------------------------------------------------------------------------
#==  Exploratory analysis and visualisation with DESeq2 data transformation    ====
#----------------------------------------------------------------------------------

# investigate  library  sizes
colSums(counts(DESeq.ds))

# normalisation des donnees en faisant le logarithme des donnees sur la taille des librairies

# calculate  the  size  factors (for library size normalization)
DESeq.ds <- estimateSizeFactors(DESeq.ds) #facteur correctif qui tient compte de la taille de la librairie
DESeq.ds@colData$sizeFactor

# make a data.frame  with meta-data  where  row.names  should  match  the  individual sample  names
# group provides informations about replicates, ind identifies uniquely each library, sex and dev allow to test for different effects
# in practive here we will only test for differences between groups, but it is useful to know that it is possible to test for sex only or dev stage only, or for the effect of both
# TODO group = () : to complete with /data/share/MADT/TP_hbadouin/2020-2021/count_matrix/samples_decription.lst
# Replicates must have the same group identifier, so that the package can use them to estimate dispersion


#### rlog transformations for exploratory analysis and visualization
# PCA and clustering should be done on normalized and preferably transformed read counts, so
# that the high variability of read counts does not occlude potentially informative data trends

# Here I used the regularized log-transformation of the DESeq2 package that reduce variance 
# of low read counts by using the dispestion-mean trend seen for the entire data set as a reference
# permettent de recuperer les valeurs normalises par DESeq
rld <- rlog(DESeq.ds, blind = FALSE)
rlog.norm.counts <- assay(rld)

# Read count correlations
# we compute 1 minus the correlation of read counts to obtain a distance
# cor permet de calculer la correlation entre les librairies via le coef de correlation de Pearson
# on calcule la distance car cela permet de faire du clastering d=0 "identiques"  d=1 "aucune similarite"
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))

# we use the distance to perform a hierarchical clustering and plot the resulting tree
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
     main = "rlog  transformed  read  counts\ndistance: Pearson  correlation")

# TODO : 
# -check that replicates are clustered in the tree
# -see if some branches are unexpectedly long 

# PCA : the PCA will inform on the structure of the signal. It helps to verify that the overall experiment produced a signal.
P <- plotPCA(rld, intgroup = "group")
P + theme_bw() + ggtitle("Rlog transformed counts")

#P <- plotPCA(rld, intgroup = "sex")
#P + theme_bw() + ggtitle("Rlog transformed counts")

# PC1 distingue les etapes du dev et PC2 separe les hermaphordites du reste
#P <- plotPCA(rld, intgroup = "dev")
#P + theme_bw() + ggtitle("Rlog transformed counts")

write.table(rlog.norm.counts, "rlog.norm.counts.tsv", quote=FALSE, row.names = TRUE, col.names = TRUE)
