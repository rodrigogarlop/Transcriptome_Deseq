library("DESeq2")
df <- read.table("input_RSEM_freq_1_samples_2_groups_2.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE) #We'll now reload the new matrix to load gene names as first column
meta <- read.table("metadata_obese.tsv",comment.char = "", sep = "\t", header = TRUE,skip = 0,quote="",fill = FALSE)
dds <- DESeqDataSetFromMatrix(countData=df, colData=meta, design=~Group, tidy = TRUE) #create deseq2 object where input data is df and metadata is meta. Using column Group in the metadata
dds <- DESeq(dds) #now carry out deseq differential analysis
#This carries out the following procedures
#estimateSizeFactors: This calculates the relative library depth of each sample 
#estimateDispersions: estimates the dispersion of counts for each gene 
#nbinomWaldTest: calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs
OvM <- results(dds, contrast=c("Group","O","M")) #Change to compare O vs M (3220 items)
pdf("Differential_abundances_q-values.pdf")
with(subset(OvM, padj>=0.05 ), plot(log2FoldChange, padj, pch=20, main="Fold change associated to groups by q-value", xlab="Fold Change (log2 scale)",ylab="q-value", ylim=c(0,1),xlim=c(-10,10), xaxt='n',yaxt='n')) #Create a volcano plot with all pvalues
mtext("Left: Obese with Metabolic Syndrome, Right: Only obese")
axis(1, las=2, at=(seq(-10,10,1)),labels=c(rev(2^(1:10)),2^(0:10)))
axis(2,las=2, at=seq(0,1,.05))
with(subset(OvM, padj<.05 ), points(log2FoldChange, padj, pch=20, col="turquoise")) #paint points with <0.05 qvalue coral
with(subset(OvM, padj<.01 ), points(log2FoldChange, padj, pch=20, col="coral")) #paint points with <0.01
dev.off()
pdf("Differential_abundances_p-values.pdf")
with(subset(OvM, pvalue>=.05 ), plot(log2FoldChange, pvalue, pch=20, main="Fold change associated to groups by p-value", xlab="Fold Change (log2 scale)",ylab="p-value", ylim=c(0,1),xlim=c(-10,10), xaxt='n',yaxt='n'))
mtext("Left: Obese with Metabolic Syndrome, Right: Only obese")
axis(1, las=2, at=(seq(-10,10,1)),labels=c(rev(2^(1:10)),2^(0:10)))
axis(2,las=2, at=seq(0,1,.05))#Create a volcano plot with all pvalues
with(subset(OvM, pvalue<.05 ), points(log2FoldChange, pvalue, pch=20, col="turquoise")) #paint points with <0.05 pvalue coral
with(subset(OvM, pvalue<.01 ), points(log2FoldChange, pvalue, pch=20, col="coral")) #paint points with <0.01 pval
dev.off()
#Get those with qval <0.05
with_q_val_OvM <- OvM[which(!is.na(OvM$padj)),] #select those that actually have a calculated q.value (571)
with_q_val_OvM <- with_q_val_OvM[with_q_val_OvM$padj <0.05,] #only keep those with q-value < 0.05 (26)
with_q_val_OvM <- with_q_val_OvM[order(with_q_val_OvM$log2FoldChange),] #sort by lfc
#Get those with pval <0.05
with_p_val_OvM <- OvM[which(!is.na(OvM$pvalue)),] #select those that actually have a calculated q.value (571)
with_p_val_OvM <- with_p_val_OvM[with_p_val_OvM$pvalue <0.05,] #only keep those with q-value < 0.05 (26)
with_p_val_OvM <- with_p_val_OvM[order(with_p_val_OvM$log2FoldChange),] #sort by lfc

# Use the raw table (which includes the healthy discarded samples) to evaluate the actual original counts
dforg <- read.table("input_matrix_for_deseq.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #We can load the original table again, to retain those
# dim(OvM_raw)
# [1] 26  6
library("viridis")
# First those with pval <0.05
OvM_raw <- dforg[rownames(with_p_val_OvM),]
pdf("Differential_abundances_pval_lt_0.5.pdf")
heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",labRow="",cexCol=0.8)
dev.off()
# Again, but with those with qval <0.05
OvM_raw <- dforg[rownames(with_q_val_OvM),]
pdf("Differential_abundances_qval_lt_0.5.pdf")
heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",labRow="",cexCol=0.8)
dev.off()

# Repeat with normalized values (extracted from the DSEq object)
dfnorm <- counts(dds, normalized=TRUE)
OvM_raw <- dfnorm[rownames(with_p_val_OvM),]
pdf("Norm_Differential_abundances_pval_lt_0.5.pdf")
heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",labRow="",cexCol=0.8)
dev.off()
# Again, but with those with qval <0.05
OvM_raw <- dfnorm[rownames(with_q_val_OvM),]
pdf("Norm_Differential_abundances_qval_lt_0.5.pdf")
heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",cexCol=0.8)
dev.off()
