########################################
# R script for running differntial expression analysis using DESeq2
#######################################

# install packages
library("tidyverse")
library('pheatmap')
library('RColorBrewer')
library('hwriter')

library('vsn')
library('DESeq2')
library('BiocParallel')
library('limma')
library('EnhancedVolcano')

args<-commandArgs(TRUE)

prefix=args[1]
expr_file_cts= args[2]

# prefix="PDBF_all_samples"
# expr_file_cts="merge.exon.allSamples.PDBF.circRNArawcount.wide.txt"

mate_file=paste0(prefix,"_cov_PCs.tsv")
output_dir = paste0("../../results/", prefix,"_circRNAexon/DE")
print(output_dir)

# Create folder if the directory doesn't exist
dir.create(output_dir, recursive =T, showWarnings = FALSE)

##==================================
# 1. load data
##==================================
print("[INFO] 1.1 Loading quantification data...") 
cts0 = read.delim(expr_file_cts,row.names = 1, check.names = FALSE)

# read meta table
print("[INFO] 1.2 Loading meta data table...")
covarianceTable = read.delim(paste0(mate_file), stringsAsFactors = F, row.names = 1, check.names = FALSE)
covarianceTable$CONDITION = factor(covarianceTable$case_control_other_latest, levels=c("Control","Case"))
covarianceTable$age = covarianceTable$age_at_visit
covarianceTable$sex = factor(covarianceTable$sex)
covarianceTable$Plate = factor(covarianceTable$Plate)

print("Origanl covariates dim...")
dim(covarianceTable);
print("PRS filtering...(skipped)")
# covarianceTable = covarianceTable[!is.na(covarianceTable$PRS),]
dim(covarianceTable);
table(covarianceTable$CONDITION)

cts = cts0
samples = intersect(rownames(covarianceTable),colnames(cts))

# filter and re-order
cts = cts[, samples]
covarianceTable = covarianceTable[samples,]
dim(covarianceTable)

cts[cts<0] =  0

all(rownames(covarianceTable) %in% colnames(cts))
dim(cts); cts = cts[, rownames(covarianceTable)]; dim(cts) #1030286     983
all(rownames(covarianceTable) == colnames(cts))

## Filtration
print("[INFO] 1.3. Data filtration...")
# keep genes with read count > 5 in >10% samples
# keep = (rowMeans(cts>5) > 0.1) # 1171  983
# keep = (rowMeans(cts>2) > 0.1) # 3182  983
# keep = (rowMeans(cts>1) > 0.1) # 5800  983
# keep = (rowMeans(cts>0) > 0.1) # 16153  983
# keep = (rowMeans(cts>0) > 0.05) # 26765  983

## Rosan, exon, basline
# keep = (rowMeans(cts>5) > 0.1) # 1152  983
keep = (rowMeans(cts>2) > 0.1) # 3126  983
# keep = (rowMeans(cts>1) > 0.1) # 5707  983
# keep = (rowMeans(cts>0) > 0.1) # 15881  983
# keep = (rowMeans(cts>0) > 0.05) # 26336  983

sum(keep)
dim(cts); cts <- cts[keep,]; dim(cts);

# keep genes with variances > 1
# dim(cts);
# cts = cts[apply(cts, 1, var, na.rm=TRUE) > 1,];
# dim(cts);


##======================================
# 2. load data to DEseq
##======================================
# Note: With no arguments to results, the results will be for the last variable in the design formula,
# and if this is a factor, the comparison will be the last level of this variable over the first level.
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
print("[INFO] 3. Load data to DEseq...")

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                colData = covarianceTable,
                                design = ~age+sex+Plate+RIN+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+CONDITION)

##======================================
# 3. QA of the data [optional]
##======================================
# Note: This part is not necessary for DEseq, but important for data QA

##--------------------------------------
## 3.1: compare different variance stabilization methods
##--------------------------------------
print("[INFO] 3.1 Data QA...vsd")
ntd <- normTransform(dds) # log2(x+1)
vsd <- varianceStabilizingTransformation(dds, blind=F) # Note: blind to the design, equal to design = ~ 1
vsd_T <- varianceStabilizingTransformation(dds, blind=T)

# using limma to remove covariates, it returns adjusted values in log2 scale
library(limma)
vsd_adjusted_log2 <- removeBatchEffect(assay(vsd), batch=vsd$age, batch2=vsd$sex)

##https://support.bioconductor.org/p/83286/
##https://support.bioconductor.org/p/134901/
# design <- model.matrix(~TREAT + A + B + C, data=covariates_df)
# treatment.design <- design[,1:2]


# batch.design <- design[,-(1:2)]
# corrected_vst_df <- removeBatchEffect(t(vst_df),design=treatment.design,covariates=batch.design)

design0 = model.matrix(~CONDITION + sex + age + Plate, data= covarianceTable)
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
write.table(assay(vsd_T), file.path(output_dir,"vst_blindT.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)

## save the variance-stabilized data with covariates adjusted
write.table(vsd_adjusted_log2, file.path(output_dir,"vsd_adjusted_log2.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)
write.table(vsd_removeBE, file.path(output_dir,"vsd_removeSexBE.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)


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
       text(log2FoldChange, -log10(padj), col="gray", cex=0.5, pos=1, offset=0.2))
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
       text(log2FoldChange, -log10(pvalue), col="gray", cex=0.5, pos=1, offset=0.2))
  with(subset(topT, padj<=0.05 & log2FoldChange>=1),
       points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
  with(subset(topT, padj<=0.05 & log2FoldChange<=-1),
       points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex=1))
  
}
abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
pCutoff = max(res$pvalue[res$padj<=0.05],na.rm =T)
abline(h=-log10(pCutoff), col="black", lty=4, lwd=2.0)

## enhanced Vocano plot
library(EnhancedVolcano)
EnhancedVolcano(res, x = 'log2FoldChange', y = 'padj', lab=FALSE,drawConnectors = TRUE)
EnhancedVolcano(res, x = 'log2FoldChange', y = 'pvalue',lab=FALSE,
                    pCutoff = max(res$pvalue[res$padj<=0.05],na.rm =T),
                    xlim = range(res$log2FoldChange, na.rm = T)*1.15,
                    ylim = max(-log10(res$pvalue))*1.15)
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
  annotation_row = dplyr::select(as.data.frame(topT), log2FoldChange) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown)
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

if(heatmap_flag == 0){
  print("Heatmap flag==0")
  plot.new()
 text(.5,.5, "No Sig DEGs",font=4, cex=2)
}

dev.off()
