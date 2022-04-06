########################################
# R script for running differntial expression analysis using DESeq2
# Usage: Rscript DE.R RNA_expr_matrix_file Sample_mate_info_file output_folder
#######################################
# install packages
require('tidyverse',quietly=T, warn.conflicts=F) || install.packages('tidyverse', repo='http://cran.rstudio.com');
require('RCurl',quietly=T, warn.conflicts=F) || install.packages('RCurl', repo='http://cran.rstudio.com');
require('hexbin',quietly=T, warn.conflicts=F) || install.packages('hexbin', repo='http://cran.rstudio.com');
require('pheatmap',quietly=T, warn.conflicts=F) || install.packages('pheatmap', repo='http://cran.rstudio.com');
require('RColorBrewer',quietly=T, warn.conflicts=F) || install.packages('RColorBrewer', repo='http://cran.rstudio.com');
require('hwriter',quietly=T, warn.conflicts=F) || install.packages('hwriter', repo='http://cran.rstudio.com');

if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
require('vsn',quietly=T, warn.conflicts=F) || BiocManager::install('vsn');
require('DESeq2',quietly=T, warn.conflicts=F) || BiocManager::install('DESeq2');
require('ReportingTools',quietly=T, warn.conflicts=F) || BiocManager::install('ReportingTools');
require('BiocParallel',quietly=T, warn.conflicts=F) || BiocManager::install('BiocParallel');
require('limma',quietly=T, warn.conflicts=F) || BiocManager::install('limma');
require('EnhancedVolcano',quietly=T, warn.conflicts=F) || BiocManager::install('EnhancedVolcano');

args<-commandArgs(TRUE)

expr_file=args[1]  # TPM data file
expr_file_cts=args[2]
mate_file=args[3]
prefix=args[4]

# expr_file="PPMI_CaseAll_CtrlNonCarrier_allgene_train70_test30_tpm_Train.txt"  # either filename or stdin
# mate_file="PPMI_CaseAll_CtrlNonCarrier_allgene_train70_test30_cov_Train.txt"
#
# prefix="PPMI_CaseAll_CtrlNonCarrier_allgene_train70_test30"

output_dir = paste0("../results/", prefix,"/DE")

# Create folder if the directory doesn't exist
dir.create(output_dir, recursive =T, showWarnings = FALSE)

##==================================
# 1. load data
##==================================
print("[INFO] 1.1 Loading quantification data...")
cts = read.delim(paste0("../run_inout/",expr_file_cts), row.names = 1, check.names = FALSE)
# read meta table
print("[INFO] 1.2 Loading meta data table...")
covarianceTable = read.delim(paste0("../run_inout/",mate_file), stringsAsFactors = F, row.names = 1, check.names = FALSE)
covarianceTable$CONDITION = factor(covarianceTable$case_control_other_latest, levels=c("Control","Case"))
# covarianceTable$CONDITION = factor(covarianceTable$Mutation, levels=c("NonCarrier","Carrier"))
covarianceTable$age = covarianceTable$age_at_baseline
covarianceTable$sex = factor(covarianceTable$sex)
covarianceTable$Plate = factor(covarianceTable$Plate)

print("Origanl covariates dim...")
dim(covarianceTable);
print("PRS filtering...(skipped)")
# covarianceTable = covarianceTable[!is.na(covarianceTable$PRS),]
dim(covarianceTable);
table(covarianceTable$CONDITION)

# filter and re-order
cts = cts[,rownames(covarianceTable)]
all(rownames(covarianceTable) %in% colnames(cts))
dim(cts); cts = cts[, rownames(covarianceTable)]; dim(cts)
all(rownames(covarianceTable) == colnames(cts))

## Filtration
print("[INFO] 1.3. Data filtration...")
# keep genes with read count > 5 in >10% samples
keep = (rowMeans(cts>5) > 0.1)
dim(cts); cts <- cts[keep,]; dim(cts);

## filter out genes with 0 in >90% samples
# notAllZero <- (rowMeans(cts > 0) > 0.1)
# dim(cts); cts <- cts[notAllZero,]; dim(cts);

# keep genes with variances > 1
dim(cts);
cts = cts[apply(cts, 1, var, na.rm=TRUE) > 1,];
dim(cts);


##======================================
# 2. load data to DEseq
##======================================
# Note: With no arguments to results, the results will be for the last variable in the design formula,
# and if this is a factor, the comparison will be the last level of this variable over the first level.
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
print("[INFO] 3. Load data to DEseq...")

## For ctrl-Carrier vs ctrl-NonCarrier:Dod are all 0, remove it from design
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                colData = covarianceTable,
                                design = ~ DoD + age + sex + Plate + CONDITION)

##======================================
# 3. QA of the data [optional]
##======================================
# Note: This part is not necessary for DEseq, but important for data QA

##--------------------------------------
## 3.1: compare different variance stabilization methods
##--------------------------------------
print("[INFO] 3.1 Data QA...vsd")
ntd <- normTransform(dds) # log2(x+1)
vsd <- vst(dds, blind=F) # Note: blind to the design, equal to design = ~ 1
vsd_T <- vst(dds, blind=T)

# using limma to remove covariates, it returns adjusted values in log2 scale
library(limma)
vsd_adjusted_log2 <- removeBatchEffect(assay(vsd), batch=vsd$age, batch2=vsd$sex)

##https://support.bioconductor.org/p/83286/
##https://support.bioconductor.org/p/134901/
# design <- model.matrix(~TREAT + A + B + C, data=covariates_df)
# treatment.design <- design[,1:2]
# batch.design <- design[,-(1:2)]
# corrected_vst_df <- removeBatchEffect(t(vst_df),design=treatment.design,covariates=batch.design)

design0 = model.matrix(~CONDITION + sex + DoD + age + Plate, data= covarianceTable)
condition.design =design0[,2, drop=F]
batch.design= design0[,-2]
vsd_removeBE = removeBatchEffect(assay(vsd),design=condition.design,covariates=batch.design)

pdf(file.path(output_dir,"data_diagnosis.pdf"))
msd <- meanSdPlot(counts(dds), ranks = FALSE); msd$gg + ggtitle("no transformation") +theme(plot.title = element_text(hjust = 0.5))
msd <- meanSdPlot(assay(ntd), ranks = FALSE); msd$gg + ggtitle("log2(x+1) transform")+ theme(plot.title = element_text(hjust = 0.5))
msd <- meanSdPlot(assay(vsd), ranks = FALSE); msd$gg + ggtitle("VST")+theme(plot.title = element_text(hjust = 0.5))
msd <- meanSdPlot(vsd_adjusted_log2, ranks = FALSE); msd$gg + ggtitle("vsd_adjusted_log2")+theme(plot.title = element_text(hjust = 0.5))
msd <- meanSdPlot(vsd_removeBE, ranks = FALSE); msd$gg + ggtitle("vsd_removeBE")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

##--------------------------------------
## 3.2: save normalized reads count
##--------------------------------------

## save the raw reads count
write.table(counts(dds), file.path(output_dir,"raw.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)

## save the variance-stabilized data
write.table(assay(vsd), file.path(output_dir,"vst.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)
write.table(assay(vsd_T), file.path(output_dir,"vst.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)

## save the variance-stabilized data with covariates adjusted
write.table(vsd_adjusted_log2, file.path(output_dir,"vsd_adjusted_log2.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)
write.table(vsd_removeBE, file.path(output_dir,"vsd_removeBE.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)

##--------------------------------------
## 3.3: clustering of samples
##--------------------------------------
print("Data clustering...")
sampleDists <- dist(t(assay(vsd)))

pdf(file.path(output_dir,"clustering.tree.pdf"))
## heatmap
par(cex=0.5, mar=c(5, 8, 4, 1))
sampleDistMatrix <- as.matrix( sampleDists )

colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(8, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         clustering_method = 'ward.D', fontsize = 8,
         col = colors,main = "Heatmap - Counts-vsd")
## tree
plot(hclust(sampleDists,method = "ward.D"), xlab='', main="Cluster Dendrogram - vsd")

#########################
sampleDists <-  dist(t(assay(vsd_T)))
## heatmap
par(cex=0.5, mar=c(5, 8, 4, 1))
sampleDistMatrix <- as.matrix( sampleDists )

colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(8, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         clustering_method = 'ward.D', fontsize = 8,
         col = colors,main = "Heatmap - vsd_T")
## tree
plot(hclust(sampleDists,method = "ward.D"), xlab='', main="Cluster Dendrogram - vsd_T")

## PCA on raw
se=SummarizedExperiment(counts(dds), colData=colData(dds))
pcaData <- plotPCA(DESeqTransform(se), intgroup = c("CONDITION"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION)) +
  geom_point(size =3) +
  # geom_text(nudge_x = 0.5, nudge_y=0.5, size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA based on raw count")+theme(plot.title = element_text(hjust = 0.5))
print(p)

## PCA on vsd
pcaData <- plotPCA(vsd, intgroup = c("CONDITION", "sex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggrepel)
p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION, shape=sex)) +
  geom_point(size =3) +
  # geom_text_repel(nudge_x = 0.5, nudge_y=0.5, size=2, alpha=1, max.overlaps =100) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  xlim(c(-50,40)) +
  ggtitle("PCA based on vsd (vst(dds, blind=F))")+theme(plot.title = element_text(hjust = 0.5))
print(p)

## PCA on vsd_T
pcaData <- plotPCA(vsd_T, intgroup = c("CONDITION", "sex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggrepel)
p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION, shape=sex)) +
  geom_point(size =3) +
  # geom_text_repel(nudge_x = 0.5, nudge_y=0.5, size=2, alpha=1, max.overlaps =100) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  xlim(c(-50,40)) +
  ggtitle("PCA based on vsd (vst(dds, blind=T))")+theme(plot.title = element_text(hjust = 0.5))
print(p)

## PCA on vsd_removeBE
vsd_removeBE_se=SummarizedExperiment(vsd_removeBE, colData=colData(vsd))
pcaData <- plotPCA(DESeqTransform(vsd_removeBE_se), intgroup = c("CONDITION","sex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION,shape=sex)) +
  geom_point(size =3) +
  # geom_text(nudge_x = 0.5, nudge_y=0.5, size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA based on vsd_removeBE")+theme(plot.title = element_text(hjust = 0.5))
print(p)

dev.off()

##==============================
# 4. Run DESeq2
##==============================
print("[INFO] 4.1.  Run DESeq2...")
## MulticoreParam() is not supported on Windows, use SnowParam()
register(MulticoreParam(8))

## It takes a while to finish......
# dds_subset <- DESeq(dds,fitType = "glmGamPoi", parallel=TRUE, BPPARAM=MulticoreParam(4))
dds_subset <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(8))
resultsNames(dds_subset)

res <- results(dds_subset, contrast=c("CONDITION","Case","Control"),
               alpha = 0.05, parallel=TRUE, BPPARAM=MulticoreParam(8))

# decimal value of Fold-change
res$FoldChange <- 2**res$log2FoldChange

genes_annotation = read.table("../../data/RNAseq/genes_annotation.tsv", header = F, stringsAsFactors = F,
                              col.names = c("chr","start","end","geneID","score","strand","geneSymbol","geneType"));
res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
res$geneDescription <- genes_annotation$chr[match(row.names(res), genes_annotation$geneID)]
# add additional columns in the output
output_additonal_columns="mi"
if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group,
  baseMeanPerLvl <- sapply( levels(dds_subset$CONDITION), function(lvl) rowMeans( counts(dds_subset,normalized=FALSE)[,dds_subset$CONDITION == lvl] ) )
  colnames(baseMeanPerLvl) = paste0("baseMean_raw.", colnames(baseMeanPerLvl))
  res = cbind(res, baseMeanPerLvl)
}
if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of normalized expression values for each group,
  baseMeanPerLvl <- sapply( levels(dds_subset$CONDITION), function(lvl) rowMeans( counts(dds_subset,normalized=TRUE)[,dds_subset$CONDITION == lvl] ) )
  colnames(baseMeanPerLvl) = paste0("baseMean_norm.", colnames(baseMeanPerLvl))
  res = cbind(res, baseMeanPerLvl)
}
if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample
  individual <- counts(dds_subset,normalized=FALSE)
  colnames(individual) = paste0("ind_raw.", colnames(individual))
  res = cbind(res, individual)
}
if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
  individual <- counts(dds_subset,normalized=TRUE)
  colnames(individual) = paste0("ind_norm.", colnames(individual))
  res = cbind(res, individual)
}
print("[INFO] 4.4 Save DESeq2 result...")
res <- res[order(res$padj),]
print("--- Res data")
dim(res)

write.table(as.data.frame(res),
            file=gzfile(file.path(output_dir, "DEresult.all.xls.gz")),
            sep="\t", quote =F, na="", row.names=T, col.names = NA)

write.table(as.data.frame(subset(res, padj<=0.05 & abs(log2FoldChange)>=1)),
            file=file.path(output_dir, "DEresult.padj05_log2FCgt1.xls"),
            sep="\t", quote =F, na="", row.names=T, col.names = NA)

write.table(as.data.frame(subset(res, padj<=0.05)),
            file=file.path(output_dir, "DEresult.padj05.xls"),
            sep="\t", quote =F, na="", row.names=T, col.names = NA)
res0 = res
res=DESeqResults(res)


##=================================
##      LFC
##=================================
# print("[INFO] 4.2 DESeq2...Run shrinkage...")
# resLFC <- lfcShrink(dds_subset, coef="CONDITION_Case_vs_Control", type="apeglm",svalue=TRUE,
#                     parallel=TRUE, BPPARAM=MulticoreParam(8))
# resLFC$FoldChange <- 2**resLFC$log2FoldChange
#
# resLFC$symbol <- genes_annotation$geneSymbol[match(row.names(resLFC), genes_annotation$geneID)]
# resLFC$geneType <- genes_annotation$geneType[match(row.names(resLFC), genes_annotation$geneID)]
# resLFC$geneDescription <- genes_annotation$chr[match(row.names(resLFC), genes_annotation$geneID)]
# if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group,
#   resLFC = cbind(resLFC, baseMeanPerLvl)
# }
# if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of normalized expression values for each group,
#   resLFC = cbind(resLFC, baseMeanPerLvl)
# }
# if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample
#   resLFC = cbind(resLFC, individual)
# }
# if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
#   resLFC = cbind(resLFC, individual)
# }
#
# print("--- resLFC data")
# resLFC <- resLFC[order(resLFC$padj),]
# # head(resLFC)
# dim(resLFC)
# write.table(as.data.frame(resLFC),
#             file=gzfile(file.path(output_dir, "DEresult.all.resLFC.xls.gz")),
#             sep="\t", quote =F, na="", row.names=T, col.names = NA)
#
# write.table(as.data.frame(subset(resLFC, padj<=0.05 & abs(log2FoldChange)>=1)),
#             file=file.path(output_dir, "DEresult.padj05_log2FCgt1.resLFC.xls"),
#             sep="\t", quote =F, na="", row.names=T, col.names = NA)
#
# write.table(as.data.frame(subset(resLFC, padj<=0.05)),
#             file=file.path(output_dir, "DEresult.padj05.resLFC.xls"),
#             sep="\t", quote =F, na="", row.names=T, col.names = NA)
# resLFC0=resLFC
# resLFC=DESeqResults(resLFC)

resLFC = res

print("[INFO] 4.3 DESeq2...Saving results...")
save.image(file.path(output_dir,"DESeq2_res.RData"))

##==============================
#  5. MAKING PLOTS
##==============================
print("[INFO] 5 Making plots...")
pdf(file.path(output_dir, "DEresult.padj_05_volcano.pdf"), paper = 'USr')
##------------------
# MA plot
##------------------
print("--- MA plot ")
DESeq2::plotMA(res, alpha = 0.05, colNonSig = "gray", main=prefix)
DESeq2::plotMA(resLFC, alpha = 0.05, colNonSig = "gray", main=paste0("resLFC: ",prefix))

##------------------
# vocano plot
##------------------
print("--- Volcano plot ")
## FDR
topT <- as.data.frame(res)
with(topT, plot(log2FoldChange, -log10(padj),
                pch=20, cex=0.5, main=paste0("res: ",prefix), col='gray',
                xlab=bquote(~Log[2]~fold~change),
                ylab=bquote(~-log[10]~FDR)))

if(nrow(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1))>0){
  with(subset(topT, padj<=0.05 & log2FoldChange>=1),
       points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1))
  with(subset(topT, padj<=0.05 & log2FoldChange<=-1),
       points(log2FoldChange, -log10(padj), pch=20, col="blue", cex=1))
  with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1),
       text(log2FoldChange, -log10(padj), labels=symbol, col="gray", cex=0.5, pos=1, offset=0.2))
}
abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)

###### pvalue
topT <- as.data.frame(res)
with(topT, plot(log2FoldChange, -log10(pvalue),
                pch=20, cex=0.5, main=paste0("res: ",prefix), col='gray',
                xlab=bquote(~Log[2]~fold~change),
                ylab=bquote(~-log[10]~pvalue)))

if(nrow(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1))>0){
  with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1),
       text(log2FoldChange, -log10(pvalue), labels=symbol, col="gray", cex=0.5, pos=1, offset=0.2))
  with(subset(topT, padj<=0.05 & log2FoldChange>=1),
       points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
  with(subset(topT, padj<=0.05 & log2FoldChange<=-1),
       points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex=1))
  
}
abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
pCutoff = max(res$pvalue[res$padj<=0.05],na.rm =T)
abline(h=-log10(pCutoff), col="black", lty=4, lwd=2.0)


####### FDR
topT_LFC <- as.data.frame(resLFC)
with(topT_LFC, plot(log2FoldChange, -log10(padj),
                pch=20, cex=0.5, main=paste0("resLFC: ",prefix), col='gray',
                xlab=bquote(~Log[2]~fold~change),
                ylab=bquote(~-log[10]~FDR)))

if(nrow(subset(topT_LFC, padj<=0.05 & abs(log2FoldChange)>=1))>0){
  with(subset(topT_LFC, padj<=0.05 & log2FoldChange>=1),
       points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1))
  with(subset(topT_LFC, padj<=0.05 & log2FoldChange<=-1),
       points(log2FoldChange, -log10(padj), pch=20, col="blue", cex=1))
  with(subset(topT_LFC, padj<=0.05 & abs(log2FoldChange)>=1),
       text(log2FoldChange, -log10(padj), labels=symbol, col="gray", cex=0.5, pos=1, offset=0.2))
}
abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)

###### pvalue
topT_LFC <- as.data.frame(resLFC)
with(topT_LFC, plot(log2FoldChange, -log10(pvalue),
                pch=20, cex=0.5, main=paste0("resLFC: ",prefix), col='gray',
                xlab=bquote(~Log[2]~fold~change),
                ylab=bquote(~-log[10]~pvalue)))

if(nrow(subset(topT_LFC, padj<=0.05 & abs(log2FoldChange)>=1))>0){
  with(subset(topT_LFC, padj<=0.05 & log2FoldChange>=1),
       points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
  with(subset(topT_LFC, padj<=0.05 & log2FoldChange<=-1),
       points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex=1))
  with(subset(topT_LFC, padj<=0.05 & abs(log2FoldChange)>=1),
       text(log2FoldChange, -log10(pvalue), labels=symbol, col="gray", cex=0.5, pos=1, offset=0.2))
}
abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
pCutoff = max(resLFC$pvalue[resLFC$padj<=0.05],na.rm =T)
abline(h=-log10(pCutoff), col="black", lty=4, lwd=2.0)

## enhanced Vocano plot
library(EnhancedVolcano)
EnhancedVolcano(res, lab = res$symbol, x = 'log2FoldChange', y = 'padj', drawConnectors = TRUE)
EnhancedVolcano(res, lab = res$symbol, x = 'log2FoldChange', y = 'pvalue',
                    pCutoff = max(res$pvalue[res$padj<=0.05],na.rm =T),
                    xlim = range(res$log2FoldChange, na.rm = T)*1.15,
                    ylim = max(-log10(res$pvalue))*1.15)

EnhancedVolcano(resLFC, lab = resLFC$symbol, x = 'log2FoldChange', y = 'padj', drawConnectors = TRUE)
EnhancedVolcano(resLFC, lab = resLFC$symbol, x = 'log2FoldChange', y = 'pvalue',
                    pCutoff = max(resLFC$pvalue[res$padj<=0.05],na.rm =T),
                    xlim = range(resLFC$log2FoldChange, na.rm = T)*1.15,
                    ylim = max(-log10(resLFC$pvalue))*1.15)
dev.off()

##------------------
# heatmap plot
##------------------
print("--- Heatmap plot ")
heatmap_flag = 0
pdf(file.path(output_dir, paste0("DEresult.padj_05.heatmap.pdf")), paper = 'USr')
if(nrow(subset(res, abs(log2FoldChange)>=1 & padj<0.05))>2){
  message("--- >>heatmap for top 20 DE genes")
  topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)

  topT=topT[order(topT$padj),] # sort by padj in increasing order
  topT=rbind(head(subset(topT, log2FoldChange<0),10),head(subset(topT, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(topT),]
  rownames(topDE) = topT$symbol
  annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
  rownames(annotation_row) = topT$symbol
  annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    CONDITION = c(Control = "green", Case = "red")
  )

  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN

  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 8,
           main =paste0("heatmap for top 20 DE genes in ", prefix),
           #fontsize_col = 5,
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none",
           clustering_method = 'ward.D',
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,cutree_cols=2,
           cluster_cols = TRUE,
           clustering_distance_cols = "correlation")

  heatmap_flag = 1;
}
print("Plot Heatmap resLFC")
if(nrow(subset(resLFC, abs(log2FoldChange)>=1 & padj<0.05))>2){
  ## ==============================
  message("# heatmap for top 10 DE genes")
  ## ==============================
  topT=subset(resLFC, abs(log2FoldChange)>=1 & padj<0.05)
  topT=topT[order(topT$padj),] # sort by padj in increasing order
  topT=rbind(head(subset(topT, log2FoldChange<0),10),head(subset(topT, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(topT),]; rownames(topDE) = topT$symbol
  annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
  rownames(annotation_row) = topT$symbol
  annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    CONDITION = c(Control = "green", Case = "red")
  )

  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN

  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           main =paste0("heatmap for top 10 DE genes in resLFC: ", prefix),
           fontsize = 8,
           #fontsize_col = 5,
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none",
           clustering_method = 'ward.D',
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,cutree_cols=2,
           cluster_cols = TRUE,
           clustering_distance_cols = "correlation")

  heatmap_flag = 1;

}

## Regardless the FoldChange
if(nrow(subset(resLFC, padj<0.05))>2){
  ## ==============================
  message("# heatmap for top 10 DE genes")
  ## ==============================
  topT=subset(resLFC, padj<0.05)
  topT=topT[order(topT$padj),] # sort by padj in increasing order
  topT=rbind(head(subset(topT, log2FoldChange<0),10),head(subset(topT, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(topT),]; rownames(topDE) = topT$symbol
  annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
  rownames(annotation_row) = topT$symbol
  annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    CONDITION = c(Control = "green", Case = "red")
  )

  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN

  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           main =paste0("heatmap for top 10 DE genes in resLFC: ", prefix),
           fontsize = 8,
           #fontsize_col = 5,
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none",
           clustering_method = 'ward.D',
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,cutree_cols=2,
           cluster_cols = TRUE,
           clustering_distance_cols = "correlation")

  heatmap_flag = 1;

}

if(heatmap_flag == 0){
  print("Heatmap flag==0")
  plot.new()
 text(.5,.5, "No Sig DEGs",font=4, cex=2)
}

dev.off()

## tapply(as.numeric(covariates_df$TREAT),covariates_df$B, table)

## SNCA:ENSG00000145335.15 padj=0.38
## MARCH8:ENSG00000165406.15 padj=0.48
## RANBP10:ENSG00000141084.10 padj=0.45
## VDR ENSG00000111424.10
## CFLAR ENSG00000003402.19
## COMT ENSG00000093010.13


## Check
## smallest padj:
# gene_ids = c("ENSG00000111424.10","ENSG00000003402.19","ENSG00000093010.13")
# for (gene_i in gene_ids){
#   # gene_i = "ENSG00000141084.10"
#   gene = as.numeric(as.vector(tpm[rownames(tpm)==gene_i,]))
#   condition = covarianceTable$CONDITION
#   df = data.frame(gene,condition)
#   df$gene = log((df$gene + 1), 2)
#   pdf(file.path(output_dir, paste0(gene_i ,"_boxplot.pdf")), paper = 'USr')
#   boxplot(gene~condition,
#           main = gene_i,
#           data=df,
#           xlab="Condition",
#           ylab="Expression log2(tpm+1)"
#   )
#   dev.off()
#   print(resLFC[gene_i,"padj"])
# }
