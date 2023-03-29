# Started on 2019-12-12 by Rodrigo García López for the 8a lab (iBT-UNAM, Mexico)
# This is based on the 00_Commands_13k_secreted_prots.R from the same project
# These commands are for analyzing proteins predicted from Obese, Obese with Metabolic Syndrome and Control samples from mexican children
#IN BASH
cd /home/rod/Desktop/Secretome/01_Data_diff_abdnc_30k
sed 's/TRINITY_//' /home/rod/Desktop/Secretome/01_Data_diff_abdnc_30k/30024_RSEM.isoforms.results.gene.counts.matrix.secretome >filtered.tsv
R #Initialize R
setwd("/home/rod/Desktop/Secretome/01_Data_diff_abdnc_30k/") #we are working locally with R version 3.6.1 (2019-07-05) -- "Action of the Toes"
list.files() #we have two files, both are preducted protien matrices without normalization
# [1] "30024_RSEM.isoforms.results.gene.counts.matrix.secretome"
# [2] "filtered.tsv"                                            
# [3] "metadata_trinity_2.txt"                                               
df <- read.table("filtered.tsv",comment.char = "", sep = "\t", header = TRUE, row.names = 1,skip = 0,quote="",fill = FALSE) #We'll start by reading the fpkm matrix (as it has not been normalized yet)
# dim(df)
# 30024     8
map <- read.table("metadata_trinity_2.txt", sep="\t",header=F, skip=0, comment.char='',quote="",fill=F) #load actual names and groups
# map
#             V1      V2      V3  V4
# 1 healthy_rep1  015_B4  NW.015  NW
# 2 healthy_rep2  164_B1  NW.164  NW
# 3   obese_rep1  024_B2   O.024   O
# 4   obese_rep2  074_B2   O.074   O
# 5   obese_rep3 090_25M   O.090   O
# 6     oms_rep1  146_B1 OMS.146 OMS
# 7     oms_rep2  153_B3 OMS.153 OMS
# 8     oms_rep3  258_B3 OMS.258 OMS

df <- df[as.character(factor(map[,1]))] #sort samples as they appear in the mapping file
names(df) <- as.character(factor(map[,3])) # update the names to the actual sample names (in 2nd col of map file)
# names(df)
# [1] "NW.015"  "NW.164"  "O.024"   "O.074"   "O.090"   "OMS.146" "OMS.153"
# [8] "OMS.258"
df <- as.matrix(df[order(rowSums(df),decreasing=T),]) # sort input by most abundant first and convert to matrix type object #we are no longer using rounding of the numbers
df <- df[which(rowSums(df)>0),] # Filter empty rows
# dim(df)
# [1] 30022     8 #Two items were empty in the table
write.table(as.matrix(df),"no_empty_rows_matrix.tsv", sep="\t", quote=FALSE, col.names=NA)
df <- read.table("no_empty_rows_matrix.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #We'll now reload the new matrix to load gene names as first column
# dim(df) #Just checking the outfile was ok
# [1] 30022     8
col_s <- names(df) #create a vector for coloring the samples by group
#col_s #oms stands for obese metabolic syndrome
# [1] "healthy_rep1" "healthy_rep2" "obese_rep1"   "obese_rep2"   "obese_rep3"  
# [6] "oms_rep1"     "oms_rep2"     "oms_rep3" 
col_s[grep("^NW\\.",col_s)] <- "black" #define group colors here
col_s[grep("^O\\.",col_s)] <- "darkorange"
col_s[grep("^OMS.",col_s)] <- "darkorchid2"
# # col_s
# [1] "black"       "black"       "darkorange"  "darkorange"  "darkorange" 
# [6] "darkorchid2" "darkorchid2" "darkorchid2"
pdf("01_Protein_reads_per_sample_raw_counts.pdf")
matplot(t(df), type="l",col=1:4,lty=2, main="All proteins per sample totals",xaxt='n') # a plot with all observations overlaped and shown per sample
axis(1, at=1:8,names(df),las=2)
barplot(colSums(df),col=col_s,xaxt='n') #bars for total items per samples
axis(1, at=seq(.7,9.1,1.2),labels=names(df),las=2)
dev.off()

#UPDATE 2019-12-03: We decided to leave the healthy patients out of the equation, thus we'll be working with 6 items only

### Remove the healthy samples
#From Bash
cut -f 1,4- no_empty_rows_matrix.tsv >filtered_obese.tsv 
grep -v "NW" metadata_trinity_2.txt >metadata_obese.tsv
# In R
df <- read.table("filtered_obese.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #We'll now reload the new matrix to load gene names as first column
# dim(df)
# [1] 29727     6
df <- df[which(rowSums(df)>0),] # Filter empty rows
# dim(df)
# 12951     6 # 259 items were only found in healthy samples
write.table(as.matrix(df),"obese_input_matrix_for_deseq.tsv", sep="\t", quote=FALSE, col.names=NA) 

### Apply stringent filters to the input matrix (sparsity reduction)
df <- read.table("obese_input_matrix_for_deseq.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #Load the input matrix (filtered has no "TRINITY" prefixes and is the result of retaining only the obese samples)
# dim(df) #just checking the outfile
# [1] 29727     6
library(Matrix) #load the matrix package
nnzero(df)/length(unlist(df)) #Get the current proportion of non-zeroes in the table
# [1] 0.2447494
(length(unlist(df))-nnzero(df))/nnzero(df) # get
# [1] 2.716568
#Set cutoffs
valid_count <- 1 # Items below this minimum are not considered
min_freq <- 2 # Set the minimum frequency sum of a feature across the whole table
min_samples <- 3  # Set the minimum number of samples expected
in_x_per_group <-2 # Only considered "in a group" if at least found in this number of samples within the group
min_groups <- 2 #Minimum groups that have this
output_name <- "input_RSEM_freq_1_samples_3_groups_2.tsv"
group_sizes <- names(df) #create a vector for detecting groups
# group_sizes[grep("healthy",group_sizes)] <- "H" #define group colors here
group_sizes[grep("O\\.",group_sizes)] <- "O"
group_sizes[grep("OMS\\.",group_sizes)] <- "M"
group_sizes <- table(group_sizes)
groups <- length(group_sizes)
df <- df[order(rowSums(df),decreasing=T),] # sort input by most abundant first
# Filter 1: min frequency
# Prior to any filter, make those with less than valid_count equal to 0 (so that they are not considered)
df[df<valid_count]=0
# length(which(rowSums(df)>=2))
# [1] 29691 # Most have at least 2
df <- df[rowSums(df)>=min_freq,] # Whole table filter features with less than min_freq total
# dim(df)
# [1] 29691     6
bin_df <- df; bin_df[bin_df>0] = 1 #Copy and create a binary table (absense/presence; use all columns
# Filter 2: min_samples
# > length(which(rowSums(bin_df)>=1))
# [1] 29691
# > length(which(rowSums(bin_df)>=2))
# [1] 7887
# > length(which(rowSums(bin_df)>=3))
# [1] 2283
# length(which(rowSums(bin_df)>=4))
# [1] 722
# length(which(rowSums(bin_df)>=5))
# [1] 127
# length(which(rowSums(bin_df)>=6))
# [1] 12

df <- df[rowSums(bin_df)>=min_samples,] #remove items only seen in less than min_samples in the binary matrix
bin_df <- df; bin_df[bin_df>0] = 1 #Recalculate the binary matrix
# dim(df)
# [1] 2283    6
init <- 0
compare <- data.frame(matrix(NA,nrow=nrow(df),ncol=groups)) # Create an empty matrix with the same number of otus with columns
for(i in 1:groups){ # For each group, sum the columns in the binary table to compare them with the min_samples cutoff
	end <- init+group_sizes[i] #set the advancing window for columns
	compare[i] <- as.numeric(rowSums(bin_df[,(init+1):end])>=in_x_per_group) # Create a binary vector with absence/presence per group, only count them as valid if they have at least in_x_per_group
	init <- init+group_sizes[i] #reset the starting column for next iteration
}
df <- df[rowSums(compare)>=min_groups,] # Filter the original set based on rows (features) where in at least n groups passes the filter
# dim(df)
# [1] 683    6
nnzero(df)/length(unlist(df)) #Get the current proportion of non-zeroes in the table
# [1] 0.7005857 #for 3 samples
(length(unlist(df))-nnzero(df))/nnzero(df) # zero to non-zero ratio
# [1] 0.4273772 #for 3 samples
write.table(round(df),output_name, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # and export it as a tsv file (using integers as deseq requires them)

#Repeat with 4 groups
output_name <- "input_RSEM_freq_1_samples_4_groups_2.tsv"
df <- read.table("obese_input_matrix_for_deseq.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #Load the input matrix (filtered has no "TRINITY" prefixes and is the result of retaining only the obese samples)
#Set cutoffs
valid_count <- 1 # Items below this minimum are not considered
min_freq <- 2 # Set the minimum frequency sum of a feature across the whole table
min_samples <- 4  # Set the minimum number of samples expected
in_x_per_group <-2 # Only considered "in a group" if at least found in this number of samples within the group
min_groups <- 2 #Minimum groups that have this
group_sizes <- names(df) #create a vector for detecting groups
# group_sizes[grep("healthy",group_sizes)] <- "H" #define group colors here
group_sizes[grep("O\\.",group_sizes)] <- "O"
group_sizes[grep("OMS\\.",group_sizes)] <- "M"
group_sizes <- table(group_sizes)
groups <- length(group_sizes)
df <- df[order(rowSums(df),decreasing=T),] # sort input by most abundant first
# Filter 1: min frequency
# Prior to any filter, make those with less than valid_count equal to 0 (so that they are not considered)
df[df<valid_count]=0
# length(which(rowSums(df)>=2))
# [1] 29691 # Most have at least 2
df <- df[rowSums(df)>=min_freq,] # Whole table filter features with less than min_freq total
bin_df <- df; bin_df[bin_df>0] = 1 
df <- df[rowSums(bin_df)>=min_samples,] #remove items only seen in less than min_samples in the binary matrix
bin_df <- df; bin_df[bin_df>0] = 1 #Recalculate the binary matrix
# dim(df)
# [1] 722   6
init <- 0
compare <- data.frame(matrix(NA,nrow=nrow(df),ncol=groups)) # Create an empty matrix with the same number of otus with columns
for(i in 1:groups){ # For each group, sum the columns in the binary table to compare them with the min_samples cutoff
	end <- init+group_sizes[i] #set the advancing window for columns
	compare[i] <- as.numeric(rowSums(bin_df[,(init+1):end])>=in_x_per_group) # Create a binary vector with absence/presence per group, only count them as valid if they have at least in_x_per_group
	init <- init+group_sizes[i] #reset the starting column for next iteration
}
df <- df[rowSums(compare)>=min_groups,] # Filter the original set based on rows (features) where in at least n groups passes the filter
# dim(df)
# [1] 683    6
nnzero(df)/length(unlist(df)) #Get the current proportion of non-zeroes in the table
# [1] 0.7005857 #for 3 samples
(length(unlist(df))-nnzero(df))/nnzero(df) # zero to non-zero ratio
# [1] 0.4273772 #for 3 samples
write.table(round(df),output_name, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # and export it as a tsv file (using integers as deseq requires them)

# Selecting 3 or 4 samples results in the same 683 proteins
# This was tested by executing the following command in bash
# diff input_RSEM_freq_1_samples_4_groups_2.tsv input_RSEM_freq_1_samples_3_groups_2.tsv

col_s <- names(df) #create a vector for coloring the samples by group
col_s[grep("^O\\.",col_s)] <- "darkorange"
col_s[grep("^OMS.",col_s)] <- "darkorchid2"
pdf("01_Protein_reads_per_sample_filtered_counts.pdf")
matplot(t(df), type="l",col=1:4,lty=2, main="All proteins per sample totals",xaxt='n') # a plot with all observations overlaped and shown per sample
axis(1, at=1:6,names(df),las=2)
barplot(colSums(df),col=col_s,xaxt='n') #bars for total items per samples
axis(1, at=seq(.7,6.7,1.2),labels=names(df),las=2)
dev.off()

### deseq with filtered set
library("DESeq2")
df <- read.table("input_RSEM_freq_1_samples_3_groups_2.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE) #We'll now reload the new matrix to load gene names as first column
# dim(df) # Just checking the outfile
# [1] 683   7
meta <- read.table("Deseq_metadata_obese.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F)
dds <- DESeqDataSetFromMatrix(countData=df, colData=meta, design=~Group, tidy = TRUE) #create deseq2 object where input data is df and metadata is meta. Using column Group in the metadata
dds <- DESeq(dds) #now carry out deseq differential analysis
#This carries out the following procedures
#estimateSizeFactors: This calculates the relative library depth of each sample 
#estimateDispersions: estimates the dispersion of counts for each gene 
#nbinomWaldTest: calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs
OvM <- results(dds, contrast=c("Group","O","OMS")) #Change to compare O vs M (1219 items with 3)
#Get those with qval <0.05
with_q_val_OvM <- OvM[which(!is.na(OvM$padj)),] #select those that actually have a calculated q.value (571)
with_q_val_OvM <- with_q_val_OvM[with_q_val_OvM$padj <0.05,] #only keep those with q-value < 0.05 (26)
with_q_val_OvM <- with_q_val_OvM[order(with_q_val_OvM$log2FoldChange),] #sort by lfc
#Get those with pval <0.05
with_p_val_OvM <- OvM[which(!is.na(OvM$pvalue)),] #select those that actually have a calculated q.value (571)
with_p_val_OvM <- with_p_val_OvM[with_p_val_OvM$pvalue <0.05,] #only keep those with q-value < 0.05 (26)
with_p_val_OvM <- with_p_val_OvM[order(with_p_val_OvM$log2FoldChange),] #sort by lfc
# Volcano plots
pdf("Differential_abundances_q-values.pdf")
with(subset(OvM, padj>=0.05 ), plot(-log2FoldChange, -padj, pch=20, main="Fold change associated to groups by q-value", xlab="Fold Change (log2 scale)",ylab="q-value", ylim=c(-1,0),xlim=c(-8,8), xaxt='n',yaxt='n')) #Create a volcano plot with all pvalues
mtext("Left: Only obese; Right: Obese with Metabolic Syndrome")
axis(1, las=2, at=(seq(-10,10,1)),labels=c(rev(2^(1:10)),2^(0:10)))
axis(2,las=2, at=seq(-1,-0,.05),labels=format(seq(-1,-0,.05)*-1,digits=2))
with(subset(OvM, padj<.05 ), points(-log2FoldChange, -padj, pch=20, col="turquoise")) #paint points with <0.05 qvalue coral1
with(subset(OvM, padj<.01 ), points(-log2FoldChange, -padj, pch=20, col="coral1")) #paint points with <0.01
dev.off()
pdf("Differential_abundances_p-values.pdf")
with(subset(OvM, pvalue>=.05 ), plot(-log2FoldChange, -pvalue, pch=20, main="Fold change associated to groups by p-value", xlab="Fold Change (log2 scale)",ylab="p-value", ylim=c(-1,0),xlim=c(-7,7), xaxt='n',yaxt='n'))
mtext("Left: Only obese; Right: Obese with Metabolic Syndrome")
axis(1, las=2, at=(seq(-10,10,1)),labels=c(rev(2^(1:10)),2^(0:10)))
axis(2,las=2, at=seq(-1,-0,.05),labels=format(seq(-1,-0,.05)*-1,digits=2))
with(subset(OvM, pvalue<.05 ), points(-log2FoldChange, -pvalue, pch=20, col="turquoise3")) #paint points with <0.05 pvalue coral1
with(subset(OvM, pvalue<.01 ), points(-log2FoldChange, -pvalue, pch=20, col="firebrick1")) #paint points with <0.01 pval
legend("bottomright", pch=20, title="p-values",col=c("turquoise3","coral1"),legend=c("<0.05:31","<0.01:  5"))
dev.off()

# # # with(subset(OvM, padj<.001 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
# # # #paint points with < .001 pvalue red
# # # vsdata <- vst(dds, blind=FALSE) #use varianceStabilizingTransformation on the data
# # # vst_expr <- assay(vsdata) #Now extract the observations
# # # pca <- prcomp(vst_expr) #create principal component analysis
# # # pvar_expl <- round(((pca$sdev ^ 2) / sum(pca$sdev ^ 2)) * 100, 2)
# # # plotPCA(vsdata, intgroup="Group")

# Use the raw table (which includes the healthy discarded samples) to evaluate the actual original counts
dforg <- read.table("no_empty_rows_matrix.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #We can load the original table again, to retain those
# dim(OvM_raw)
# [1] 26  6
library("viridis")
# First those with pval <0.05
OvM_raw <- dforg[rownames(with_p_val_OvM),]
pdf("Differential_abundances_pval_lt_0.5.pdf")
heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs OMS", xlab="", ylab="Differentially abundant features",cexCol=0.8)
dev.off()
# # Again, but with those with qval <0.05
# OvM_raw <- dforg[rownames(with_q_val_OvM),]
# pdf("Differential_abundances_qval_lt_0.5.pdf")
# heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",labRow="",cexCol=0.8)
# dev.off()

# Repeat with normalized values (extracted from the DSEq object)
dfnorm <- counts(dds, normalized=TRUE)
OvM_std <- dfnorm[rownames(with_p_val_OvM),]
pdf("Norm_Differential_abundances_pval_lt_0.5.pdf")
heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",cexCol=0.8,cexRow=0.5)
dev.off()
# # Again, but with those with qval <0.05
# OvM_std <- dfnorm[rownames(with_q_val_OvM),]
# pdf("Norm_Differential_abundances_qval_lt_0.5.pdf")
# heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",cexCol=0.8)
# dev.off()

### get the data for plotting the differential abundances
write.table(OvM_raw,"Differential_features_with_pval_lt_0.05-raw_table.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # export the item observations in the raw table
write.table(OvM_std,"Differential_features_with_pval_lt_0.05-std_table.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # export the item observations in the standardized table (No normal samples present)
dseq_diff <- as.data.frame(with_p_val_OvM)
write.table(dseq_diff,"Differential_features_with_pval_lt_0.05.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # and export it as a tsv file

### Plot the fold change
colors <- c(rep("firebrick1",length(which(dseq_diff[,2]<0))),rep("darkorange1",length(which(dseq_diff[,2]>0))))
pdf("log2FoldChange_extra.pdf")
barplot(abs(dseq_diff[,2]),horiz=TRUE, xlim=c(2,8), xpd = FALSE, col=colors,border=colors, main="Log2 Fold Change of differentially abundance features", xlab="log2 Fold Change", ylab="")
axis(2,las=1,at=seq(0.7, 37,1.2), labels=rownames(dseq_diff),cex=0.5)
abline(v=seq(1,10,.5), lty= 2, col="darkgray")
legend("right",pch=15, col=c("darkorange1", "firebrick1"), legend=c("Obese > OMS","OMS > Obese"), bg="white",pt.cex=2, title="Associated Group")
dev.off()

### extract the fold change and get the raw values (2 exp value)
matrix <- cbind("Log2_Fold_Change"=abs(dseq_diff[,2]),"Raw_Fold_Change"=2^abs(dseq_diff[,2]))
row.names(matrix) <- row.names(dseq_diff)
write.table(matrix,"Differential_features_with_pval_lt_0.05_l2fc.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # and export it as a tsv file

### Create a new subset to normalize the 683 items in used for O and OMS DESeq to include the NW samples
df <- read.table("input_RSEM_freq_1_samples_3_groups_2.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1)
# dim(df)
# [1] 683   6
df_all <- read.table("no_empty_rows_matrix.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1)
# dim(df_all)
# [1] 30022     8
subset <- df_all[rownames(df),] # extract the corresponding 683 items in the original table (so that we get the NW samples)
subset[subset<1]=0
write.table(round(subset),"Raw_subset_683_prots_all_samples.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # and export it as a tsv file
# 138 out of 683 do not have at least 1 observation in the NW samples
# Only 85 items had a sum >100 in the NW samples (both sets considered)
library("DESeq2")
# [1] 683   7
df <- read.table("Raw_subset_683_prots_all_samples.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE) #We'll now reload the new matrix to load gene names as first column
meta2 <- read.table("Deseq_metadata.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F)
dds2 <- DESeqDataSetFromMatrix(countData=df, colData=meta2, design=~Group, tidy = TRUE) #create deseq2 object where input data is df and metadata is meta. Using column Group in the metadata
dds2 <- DESeq(dds2) #now carry out deseq differential analysis
#This carries out the following procedures
#estimateSizeFactors: This calculates the relative library depth of each sample 
#estimateDispersions: estimates the dispersion of counts for each gene 
#nbinomWaldTest: calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs
OvM2 <- results(dds2, contrast=c("Group","O","OMS")) #Change to compare O vs M (1219 items with 3)
#Get those with qval <0.05
with_q_val_OvM2 <- OvM2[which(!is.na(OvM2$padj)),] #select those that actually have a calculated q.value (571)
with_q_val_OvM2 <- with_q_val_OvM2[with_q_val_OvM2$padj <0.05,] #only keep those with q-value < 0.05 (26)
with_q_val_OvM2 <- with_q_val_OvM2[order(with_q_val_OvM2$log2FoldChange),] #sort by lfc
#Get those with pval <0.05
with_p_val_OvM2 <- OvM2[which(!is.na(OvM2$pvalue)),] #select those that actually have a calculated q.value (571)
with_p_val_OvM2 <- with_p_val_OvM2[with_p_val_OvM2$pvalue <0.05,] #only keep those with q-value < 0.05 (26)
with_p_val_OvM2 <- with_p_val_OvM2[order(with_p_val_OvM2$log2FoldChange),] #sort by lfc
# Volcano plots
pdf("Differential_abundances_q-values_with_NW.pdf")
with(subset(OvM2, padj>=0.05 ), plot(-log2FoldChange, -padj, pch=20, main="Fold change associated to groups by q-value", xlab="Fold Change (log2 scale)",ylab="q-value", ylim=c(-1,0),xlim=c(-10,10), xaxt='n',yaxt='n')) #Create a volcano plot with all pvalues
mtext("Left: Only obese; Right: Obese with Metabolic Syndrome")
axis(1, las=2, at=(seq(-10,10,1)),labels=c(rev(2^(1:10)),2^(0:10)))
axis(2,las=2, at=seq(-1,-0,.05),labels=format(seq(-1,-0,.05)*-1,digits=2))
with(subset(OvM2, padj<.05 ), points(-log2FoldChange, -padj, pch=20, col="turquoise")) #paint points with <0.05 qvalue coral1
with(subset(OvM2, padj<.01 ), points(-log2FoldChange, -padj, pch=20, col="coral1")) #paint points with <0.01
dev.off()
pdf("Differential_abundances_p-values_with_NW.pdf")
with(subset(OvM2, pvalue>=.05 ), plot(-log2FoldChange, -pvalue, pch=20, main="Fold change associated to groups by p-value", xlab="Fold Change (log2 scale)",ylab="p-value", ylim=c(-1,0),xlim=c(-8,8), xaxt='n',yaxt='n'))
mtext("Left: Only obese; Right: Obese with Metabolic Syndrome")
axis(1, las=2, at=(seq(-10,10,1)),labels=c(rev(2^(1:10)),2^(0:10)))
axis(2,las=2, at=seq(-1,-0,.05),labels=format(seq(-1,-0,.05)*-1,digits=2))
with(subset(OvM2, pvalue<.05 ), points(-log2FoldChange, -pvalue, pch=20, col="turquoise3")) #paint points with <0.05 pvalue coral1
with(subset(OvM2, pvalue<.01 ), points(-log2FoldChange, -pvalue, pch=20, col="firebrick1")) #paint points with <0.01 pval
legend("bottomright", pch=20, title="p-values",col=c("turquoise3","coral1"),legend=c("<0.05:62","<0.01:  14"))
dev.off()

# # # with(subset(OvM, padj<.001 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
# # # #paint points with < .001 pvalue red
# # # vsdata <- vst(dds, blind=FALSE) #use varianceStabilizingTransformation on the data
# # # vst_expr <- assay(vsdata) #Now extract the observations
# # # pca <- prcomp(vst_expr) #create principal component analysis
# # # pvar_expl <- round(((pca$sdev ^ 2) / sum(pca$sdev ^ 2)) * 100, 2)
# # # plotPCA(vsdata, intgroup="Group")

# Use the raw table (which includes the healthy discarded samples) to evaluate the actual original counts
dforg <- read.table("no_empty_rows_matrix.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1) #We can load the original table again, to retain those
library("viridis")
# First those with pval <0.05
OvM_raw2 <- dforg[rownames(with_p_val_OvM2),]
pdf("Differential_abundances_pval_lt_0.5_with_NW.pdf")
heatmap(as.matrix(OvM_raw2),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs OMS", xlab="", ylab="Differentially abundant features",cexCol=0.8,cexRow=0.5)
dev.off()
# # Again, but with those with qval <0.05
# OvM_raw <- dforg[rownames(with_q_val_OvM),]
# pdf("Differential_abundances_qval_lt_0.5.pdf")
# heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",labRow="",cexCol=0.8)
# dev.off()

# Repeat with normalized values (extracted from the DSEq object)
dfnorm2 <- counts(dds2, normalized=TRUE)
OvM_std2 <- dfnorm2[rownames(with_p_val_OvM2),]
pdf("Norm_Differential_abundances_pval_lt_0.5_with_NW.pdf")
heatmap(as.matrix(OvM_std2),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",cexCol=0.8,cexRow=0.5)
dev.off()
# # Again, but with those with qval <0.05
# OvM_std <- dfnorm[rownames(with_q_val_OvM),]
# pdf("Norm_Differential_abundances_qval_lt_0.5.pdf")
# heatmap(as.matrix(OvM_raw),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",cexCol=0.8)
# dev.off()

write.table(OvM_raw2,"Differential_features_with_pval_lt_0.05-raw_table_with_NW.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # export the item observations in the raw table
write.table(OvM_std2,"Differential_features_with_pval_lt_0.05-std_table_with_NW.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # export the item observations in the standardized table (No normal samples present)
dseq_diff2 <- as.data.frame(with_p_val_OvM2)
write.table(dseq_diff2,"Differential_features_with_pval_lt_0.05.tsv_with_NW", sep="\t", quote=FALSE, row.names=T,col.names =NA) # and export it as a tsv file
nrow(na.omit(dseq_diff2[rownames(dseq_diff),]))
write.table(na.omit(dseq_diff2[rownames(dseq_diff),]),"Differential_features_in_both_DESeq_O_OMS_and_DESeq_All_groups.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # and export it as a tsv file. Values are from the DESeq_All set
# [1] 11 # 11 markers were present in both the DESeq_O_OMS and the DESeq_All_groups sets


# Let's now try something different: create a heatmap plot with the 31 markers in the first DESeq analysis using the standardized tables of the DESeq_All trial
# nrow(dseq_diff)
# [1] 31
dfnorm2 <- counts(dds2, normalized=TRUE)
dfnorm2 <- counts(dds2, normalized=TRUE)
# dim(dfnorm2)
# [1] 683   8
DESeq_O_OMS_in_DESeq_All_STD <- dfnorm2[rownames(dseq_diff),]
pdf("DESeq_O_OMS_in_DESeq_All_STD.pdf")
heatmap(as.matrix(DESeq_O_OMS_in_DESeq_All_STD),col=plasma(50),scale="row",Colv = NA, Rowv = NA, main="Differentially abundant transcripts O vs M", xlab="", ylab="Differentially abundant features",cexCol=0.8,cexRow=0.5)
dev.off()
#Finally, output this data
write.table(DESeq_O_OMS_in_DESeq_All_STD,"DESeq_O_OMS_in_DESeq_All_STD_per_sample.tsv", sep="\t", quote=FALSE, row.names=T,col.names =NA) # and export it as a tsv file



### Other tests
# Use the normalized values to evaluate distances
# dist(vst_expr)
samples <- colnames(dfnorm2)
groups <- NULL
groups[grep("^NW\\.",samples)] <- "chartreuse3"
groups[grep("^O\\.",samples)] <- "darkorange1"
groups[grep("^OMS\\.",samples)] <- "firebrick1"
colors <- levels(factor(groups))

library("vegan")
bray <- vegdist(t(dfnorm2),method="bray",upper=TRUE,diag=TRUE) # calculate similarity matrix
PCOA <- cmdscale(bray, eig = T, k=ncol(dfnorm2)-1)
PCOA.axes <- round(PCOA$eig*100/sum(PCOA$eig),1) #And determine how much variation each set of linear combinations can explain
adonis_grp <- adonis(formula = bray ~ as.factor(groups), permutations = 10000, method = t(df)) #get some statistics
grp_R2 <- round(adonis_grp$aov.tab[1,5],4) #retain R2
grp_pval <- round(adonis_grp$aov.tab[1,6],4) #and pvalue
pdf("DESeq_Normalized_Secretome_PCOA_filtered_items.pdf")
# par(mfrow=c(1,3))
# Plot coordinates 1,2
ordiplot(PCOA,type="n", main="", choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%")) #Create empty plot with % explained by axes
mtext(paste("Dimensions 1 and 2 - Adonis R2:",grp_R2,"p =",grp_pval)) #add the adonis test results
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=groups, choices=c(1,2)) #add the samples
ordihull(PCOA, groups=groups, draw="polygon", col=colors, border=colors, label=F, choices=c(1,2)) #Paint polygons
# Plot coordinates 1,3
ordiplot(PCOA,type="n", main="", choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",grp_R2,"p =",grp_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=groups, choices=c(1,3))
ordihull(PCOA, groups=groups, draw="polygon", col=colors, border=colors, label=F, choices=c(1,3))
# Plot coordinates 2,3
ordiplot(PCOA,type="n", main="", choices=c(2,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",grp_R2,"p =",grp_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=groups, choices=c(2,3))
ordihull(PCOA, groups=groups, draw="polygon", col=colors, border=colors, label=F, choices=c(2,3))
dev.off()
#Now plot some overlying metadata
meta <- read.table("metadata_obese.tsv",comment.char = "", sep = "\t", header = TRUE,skip = 0,quote="",fill = FALSE)
unique <- apply(dfnorm>0,2,sum) # Get a vector of unique features per sample
NMDS <- metaMDS(bray,distance=bray,k=5,trymax=10000,autotransform=FALSE,wascores = FALSE)
pdf("PCOA_filtered_items.pdf")
# par(mfrow=c(1,3))
# Plot coordinates 1,2
ordiplot(NMDS,type="n", main="", choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%")) #Create empty plot with % explained by axes
mtext(paste("Dimensions 1 and 2 - Adonis R2:",grp_R2,"p =",grp_pval)) #add the adonis test results
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=groups, choices=c(1,2)) #add the samples
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE) #Now plot the metadata
# Plot coordinates 1,3
ordiplot(PCOA,type="n", main="", choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",grp_R2,"p =",grp_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=groups, choices=c(1,3))
ordihull(PCOA, groups=col_s, draw="polygon", col=colors, border=colors, label=F, choices=c(1,3))
# Plot coordinates 2,3
ordiplot(PCOA,type="n", main="", choices=c(2,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",grp_R2,"p =",grp_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=groups, choices=c(2,3))
ordihull(PCOA, groups=col_s, draw="polygon", col=colors, border=colors, label=F, choices=c(2,3))
dev.off()

### Correlations
# We will be using file DESeq_O_OMS_in_DESeq_All_STD_per_sample.tsv to correlate its items to the different clinical variables in a modified version of Clinical_data.tsv (named Clinical_data_secretome_samples_only.tsv) contructed for the secretome and only using our 8 samples with the same tags and order.
setwd("/home/rod/Desktop/Resp_IBt/Resp_IBt_Secretome_local_2019_12-13/00_All_tables/Metadata/")
meta <- read.table("Clinical_data_secretome_samples_only.tsv",comment.char = "", sep = "\t", header = TRUE,skip = 0,quote="",fill = FALSE,row.names=1,check.names=F)
df <- as.matrix(read.table("DESeq_O_OMS_in_DESeq_All_STD_per_sample.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1)) #This data have been normalized with DESeq
df <- t(apply(df,1, function(x) x/colSums(df))) # change into proportions summing up to 1
# dim(df)
# [1] 31  8
# meta
#          Org_name Sample     Name Group  Age Sex Weight  Size BMI Waist
# NW.015         15     15 Norm-015  Norm  9.1   f     34 142.0  56  58.5
# NW.164       164B    164 Norm-164  Norm  7.3   m     26 128.6  58  56.3
# O.024          24     24    O-024    Ob 10.1   m     43 137.6  95  80.4
# O.074          74     74    O-074    Ob  9.1   m     44 138.5  96  87.0
# O.090   VE15.090R     90    O-090    Ob  9.4   m     45 143.2  95  84.4
# OMS.146       146    146  OMS-146   OMS  9.0   m     59 146.1  98  91.9
# OMS.153      153I    153  OMS-153   OMS  7.9   m     52 136.4  99  90.2
# OMS.258      258B    258  OMS-258   OMS  9.3   f     44 131.2  98  89.4
#         Waist..percentile.     TA Systolic.BP..percentile.
# NW.015                  25  95/70                       22
# NW.164                  25  98/65                       41
# O.024                  >75 100/68                       44
# O.074                  >90  99/62                       37
# O.090                  >90  91/64                       11
# OMS.146                >90 124/80                       96
# OMS.153                >90  98/62                       33
# OMS.258                >90  88/60                       13
#         Diastolic.BP..percentile. Glucose..mg.dl. Triglycerides..mg.dl.
# NW.015                         79              81                    79
# NW.164                         68              92                    35
# O.024                          74              90                    55
# O.074                          52              90                    56
# O.090                          56              96                    72
# OMS.146                        93              81                   135
# OMS.153                        53              89                   152
# OMS.258                        53              79                   126
#         HDL..mg.dl.
# NW.015           50
# NW.164           61
# O.024            52
# O.074            59
# O.090            40
# OMS.146          47
# OMS.153          33
# OMS.258          42

# pdf("Test.pdf")
# pairs(meta)
# pairs(df)
# dev.off()

# library("GGally")
# library(tidyr)
library("ggplot2")

for (i in c(5,7,8,9,10,13:18)){
    clin <- names(meta)[i]
    clin_out <- gsub(" ", "_",clin)
    clin_out <- gsub("/", "_",clin_out)
    clin_out <- gsub("\\(", "_",clin_out)
    clin_out <- gsub("\\)", "_",clin_out)
    correlation <- apply(df, 1, function(x) cor(x,meta[,i],method="spearman")) # Calculate the correlation of the age of all differentially abundant proteins vs all numerical data
    corrp <- apply(df, 1, function(x) cor.test(x,meta[,i],method="spearman")$p.value) #and get a pvalue for the correlation test
    out <- cbind(correlation,corrp,df)
    out <- out[order(out[,2],decreasing=F),] #sort by pvalue
    write.table(out,paste("Spearman_cor_prots_vs",clin_out,"31_diff_abd.tsv",sep="_"), sep="\t", quote=FALSE, col.names=NA)
    plot_list = list()
    for (j in 1:nrow(df)){
        x <- out[j,3:ncol(out)]
        y <- meta[,i]
        R <- round(cor(x,y,method="spearman"),4)
        model <- lm(y ~ x) #get 1st polynomial regression
        b <- trunc(model$coefficients[1]*100)/100
        sl <- trunc(model$coefficients[2]*100)/100
        pval <- round(cor.test(x,y,method="spearman")$p.value,4)
        p <- ggplot(data.frame(x,y),aes(x, y)) + geom_point() + theme_bw() + geom_smooth(method=lm) + labs(title =  substitute(paste(italic(R)^2 == r2,   "   p-value" == pvalue,"  ", italic(y) == inter, " + ",slope, italic(x)),list(r2 = R, pvalue=pval, slope=sl, inter=b))) + labs(x = paste(rownames(df)[j], "relative abundance"),y=clin)
        plot_list[[j]] = p
#         plot(x,y, main="1st deg. regression fit to data", xlab=paste(rownames(df)[j], "relative abundance"), ylab=clin,pch=20)
#         mtext(paste("Spearman R2:",R," p-val:",pval))
#         lines(x,predict(model),lwd=1.5, col="blue")
    }
    pdf(paste("Spearman_cor_prots_vs",clin_out,"31_diff_abd.pdf",sep="_"))
    for (j in 1:nrow(df)) {
        print(plot_list[[j]])
    }
    dev.off()
}
# apply(df, 1, function(x) plot(x,meta[,i],xlab=names(x)))#)
# dev.off()

### Now repeat with the LEfSe results
setwd("/home/rod/Desktop/Resp_IBt/Resp_IBt_Secretome_local_2019_12-13/00_All_tables/LEfSe_vs_clinical")
meta <- read.table("Clinical_data.tsv",comment.char = "", sep = "\t", header = TRUE,skip = 0,quote="",fill = FALSE,row.names=1,check.names=F)
df <- as.matrix(read.table("name_otu_table_sorted_L2.tsv",comment.char = "", sep = "\t", header = TRUE, skip = 0,quote="",fill = FALSE,row.names=1)) #This data have been normalized with DESeq
# # # df <- t(apply(df,1, function(x) x/colSums(df))) # change into proportions summing up to 1 # They were already in relative numbers
library("ggplot2")
K="L2"
for (i in c(5,7,8,9,10,13:17)){
    clin <- names(meta)[i]
    clin_out <- gsub(" ", "_",clin)
    clin_out <- gsub("/", "_",clin_out)
    clin_out <- gsub("\\(", "_",clin_out)
    clin_out <- gsub("\\)", "_",clin_out)
    correlation <- apply(df, 1, function(x) cor(x,meta[,i],method="spearman")) # Calculate the correlation of the age of all differentially abundant proteins vs all numerical data
    corrp <- apply(df, 1, function(x) cor.test(x,meta[,i],method="spearman")$p.value) #and get a pvalue for the correlation test
    out <- cbind(correlation,corrp,df)
    out <- out[order(out[,2],decreasing=F),] #sort by pvalue
    write.table(out,paste(K,"Spearman_cor_prots_vs",clin_out,"7_diff_abd.tsv",sep="_"), sep="\t", quote=FALSE, col.names=NA)
    plot_list = list()
    for (j in 1:nrow(df)){
        x <- out[j,3:ncol(out)]
        y <- meta[,i]
        R <- round(cor(x,y,method="spearman"),4)
        model <- lm(y ~ x) #get 1st polynomial regression
        b <- trunc(model$coefficients[1]*100)/100
        sl <- trunc(model$coefficients[2]*100)/100
        pval <- round(cor.test(x,y,method="spearman")$p.value,4)
        p <- ggplot(data.frame(x,y),aes(x, y)) + geom_point() + theme_bw() + geom_smooth(method=lm) + labs(title =  substitute(paste(italic(R)^2 == r2,   "   p-value" == pvalue,"  ", italic(y) == inter, " + ",slope, italic(x)),list(r2 = R, pvalue=pval, slope=sl, inter=b))) + labs(x = paste(rownames(df)[j], "relative abundance"),y=clin)
        plot_list[[j]] = p
#         plot(x,y, main="1st deg. regression fit to data", xlab=paste(rownames(df)[j], "relative abundance"), ylab=clin,pch=20)
#         mtext(paste("Spearman R2:",R," p-val:",pval))
#         lines(x,predict(model),lwd=1.5, col="blue")
    }
    pdf(paste(K,"Spearman_cor_prots_vs",clin_out,"7_diff_abd.pdf",sep="_"))
    for (j in 1:nrow(df)) {
        print(plot_list[[j]])
    }
    dev.off()
}
i=18 # LDL is a special case
# OMS.064 and OMS.446 are missing and are removed before processing correlations (rows 21 and 27).
df <- df[,colnames(df[,-21])]
df <- df[,colnames(df[,-26])]
meta <- meta[rownames(meta[-21,],),]
meta <- meta[rownames(meta[-26,],),]
clin <- names(meta)[i]
clin_out <- gsub(" ", "_",clin)
clin_out <- gsub("/", "_",clin_out)
clin_out <- gsub("\\(", "_",clin_out)
clin_out <- gsub("\\)", "_",clin_out)
correlation <- apply(df, 1, function(x) cor(x,meta[,i],method="spearman")) # Calculate the correlation of the age of all differentially abundant proteins vs all numerical data
corrp <- apply(df, 1, function(x) cor.test(x,meta[,i],method="spearman")$p.value) #and get a pvalue for the correlation test
out <- cbind(correlation,corrp,df)
out <- out[order(out[,2],decreasing=F),] #sort by pvalue
write.table(out,paste(K,"Spearman_cor_prots_vs",clin_out,"diff_abd.tsv",sep="_"), sep="\t", quote=FALSE, col.names=NA)
plot_list = list()
for (j in 1:nrow(df)){
    x <- out[j,3:ncol(out)]
    y <- meta[,i]
    R <- round(cor(x,y,method="spearman"),4)
    model <- lm(y ~ x) #get 1st polynomial regression
    b <- trunc(model$coefficients[1]*100)/100
    sl <- trunc(model$coefficients[2]*100)/100
    pval <- round(cor.test(x,y,method="spearman")$p.value,4)
    p <- ggplot(data.frame(x,y),aes(x, y)) + geom_point() + theme_bw() + geom_smooth(method=lm) + labs(title =  substitute(paste(italic(R)^2 == r2,   "   p-value" == pvalue,"  ", italic(y) == inter, " + ",slope, italic(x)),list(r2 = R, pvalue=pval, slope=sl, inter=b))) + labs(x = paste(rownames(df)[j], "relative abundance"),y=clin)
    plot_list[[j]] = p
    # plot(x,y, main="1st deg. regression fit to data", xlab=paste(rownames(df)[j], "relative abundance"), ylab=clin,pch=20)
    # mtext(paste("Spearman R2:",R," p-val:",pval))
    # lines(x,predict(model),lwd=1.5, col="blue")
}
pdf(paste(K,"Spearman_cor_prots_vs",clin_out,"7_diff_abd.pdf",sep="_"))
    for (j in 1:nrow(df)) {
    print(plot_list[[j]])
}
dev.off()
