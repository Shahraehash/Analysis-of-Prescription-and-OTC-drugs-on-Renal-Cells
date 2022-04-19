########################################################################################################################################################################
library(DESeq2)
library(clusterProfiler)
library(org.Rn.eg.db)
########################################################################################################################################################################
#Cyclosporine Data Analysis
#Load Cyclosporine Data
Cyclosporine_data <- read.csv("Cyclosporine_RNA_Seq_Counts.csv", sep = ",", row.names = NULL, header = TRUE)
Cyclosporine_data_2 <- Cyclosporine_data[,1:7]
Cyclosporine_metadata <- read.csv(file = "Cyclosporine_RNA_Seq_Counts_MetaData.csv")



#DESeq On Cyclosporine
dds_Cyclosporine <- DESeqDataSetFromMatrix(countData = Cyclosporine_data_2, colData = Cyclosporine_metadata, design =~CsA_Healthy, tidy = TRUE)
dds_Cyclosporine <- DESeq(dds_Cyclosporine)
results_dds_Cyclosporine <- results(dds_Cyclosporine)
results_dds_Cyclosporine <- results_dds_Cyclosporine[order(results_dds_Cyclosporine$log2FoldChange),]

#Convert from S4 to Df
write.csv(results_dds_Cyclosporine, 'results_Cyclosporine.csv')
dataframe_results_Cyclosporine <- read.csv('results_Cyclosporine.csv', sep = ",", row.names = NULL, header = TRUE)

#filter based on log2FoldChange value
threshold <- 0.1
dataframe_results_Cyclosporine_least_log_fold_change <- dataframe_results_Cyclosporine[which(dataframe_results_Cyclosporine$log2FoldChange > (-1*threshold)),]
dataframe_results_Cyclosporine_least_log_fold_change <- dataframe_results_Cyclosporine_least_log_fold_change[which(dataframe_results_Cyclosporine_least_log_fold_change$log2FoldChange < threshold),]

#GSEA
Cyclosporine_PA <- dataframe_results_Cyclosporine_least_log_fold_change
gene_name_PA <- Cyclosporine_PA$X

Cyclosporine_original_gene_list <- Cyclosporine_PA$log2FoldChange
names(Cyclosporine_original_gene_list) <- Cyclosporine_PA$X
Cyclosporine_gene_list <- na.omit(Cyclosporine_original_gene_list)
Cyclosporine_gene_list = sort(Cyclosporine_gene_list, decreasing = TRUE)

Cyclosporine_gse <- gseGO(geneList=Cyclosporine_gene_list, 
                          ont ="ALL", 
                          keyType = "GENENAME", 
                          nPerm = 10000, 
                          minGSSize = 3, 
                          maxGSSize = 800, 
                          pvalueCutoff = 1, 
                          verbose = TRUE, 
                          OrgDb = org.Hs.eg.db, 
                          pAdjustMethod = "none")
require(DOSE)
dotplot(Cyclosporine_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
########################################################################################################################################################################
#Acetaminophen Data Analysis
#Load Acetaminophen Data
Acetaminophen_data <- read.csv("Acetaminophen_RNA_Seq_Counts.csv", sep = ",", row.names = NULL, header = TRUE)
Acetaminophen_data_2 <- Acetaminophen_data[,1:7]
Acetaminophen_metadata <- read.csv(file = "Acetaminophen_RNA_Seq_Counts_MetaData.csv")

#DESeq on Acetaminophen
dds_Acetaminophen <- DESeqDataSetFromMatrix(countData = Acetaminophen_data_2, colData = Acetaminophen_metadata, design =~Control.Treatment, tidy = TRUE)
dds_Acetaminophen <- DESeq(dds_Acetaminophen)
results_dds_Acetaminophen <- results(dds_Acetaminophen)
results_dds_Acetaminophen <- results_dds_Acetaminophen[order(results_dds_Acetaminophen$log2FoldChange),]

#convert S4 to Df
write.csv(results_dds_Acetaminophen, 'results_Acetaminophen.csv')
dataframe_results_Acetaminophen <- read.csv('results_Acetaminophen.csv', sep = ",", row.names = NULL, header = TRUE)

#filter based on log2FoldChange value 
dataframe_results_Acetaminophen_least_log_fold_change <- dataframe_results_Acetaminophen[which(dataframe_results_Acetaminophen$log2FoldChange > (-1*threshold)),]
dataframe_results_Acetaminophen_least_log_fold_change <- dataframe_results_Acetaminophen_least_log_fold_change[which(dataframe_results_Acetaminophen_least_log_fold_change$log2FoldChange < threshold),]

#GSEA
Acetaminophen_PA <- dataframe_results_Acetaminophen_least_log_fold_change
Acetaminophen_gene_name_PA <- Acetaminophen_PA$X

Acetaminophen_original_gene_list <- Acetaminophen_PA$log2FoldChange
names(Acetaminophen_original_gene_list) <- Acetaminophen_PA$X
Acetaminophen_gene_list <- na.omit(Acetaminophen_original_gene_list)
Acetaminophen_gene_list = sort(Acetaminophen_gene_list, decreasing = TRUE)


Acetaminophen_gse <- gseGO(Acetaminophen_gene_list, 
                           ont ="ALL", 
                           keyType = "ENSEMBL", 
                           nPerm = 10000, 
                           minGSSize = 3, 
                           maxGSSize = 800, 
                           pvalueCutoff = 1, 
                           verbose = TRUE, 
                           OrgDb = org.Rn.eg.db, 
                           pAdjustMethod = "none")
require(DOSE)
dotplot(Acetaminophen_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


