###########################################
# Rscript to do the qulity control(QC) processes
# Usage: Rscript QC.R RNA_expr_matrix_file Sample_mate_info_file output_folder
###########################################
if (!require(ape)) install.packages('ape'); library(ape)
if (!require(reshape2)) install.packages('reshape2'); library(reshape2)
require('RColorBrewer',quietly=T, warn.conflicts=F) || install.packages('RColorBrewer', repo='http://cran.rstudio.com');
library(ggplot2)

args<-commandArgs(TRUE)

expr_file=args[1]  # either filename or stdin
mate_file=args[2]
prefix=args[3]

output_folder = paste0("../results/", prefix,"/QC")
expr_file_cts = paste0(prefix,"_read_Train.txt")

dir.create(output_folder, recursive = T, showWarnings = FALSE)

##==================================
# 1. load data
##===================================
message("[INFO] 1. Loading expression data...")
tpm = read.delim(paste0("../run_inout/",expr_file), row.names = 1, check.names = FALSE)
cts = read.delim(paste0("../run_inout/",expr_file_cts), row.names = 1, check.names = FALSE)
message("#Expression data dim:"); dim(tpm);

covarianceTable = read.delim(paste0("../run_inout/",mate_file), stringsAsFactors = F, row.names = 1, check.names = FALSE)
covarianceTable$CONDITION = covarianceTable$case_control_other_latest

## re-order
if(all(rownames(covarianceTable) %in% colnames(tpm))){
  tpm = tpm[, rownames(covarianceTable)]
}else{
  stop("Error,Check the sample IDs")
}

if(all(rownames(covarianceTable) %in% colnames(tpm))){
  tpm = tpm[, rownames(covarianceTable)]
}else{
  stop("[Error] Check the sample IDs")
}

if(all(rownames(covarianceTable) != colnames(tpm))) {
  stop("[Error] sample ID orders in Expr and Meta file do not match.")
}

##===================================
# 2. Filtration
##===================================
message("[INFO] 2. Data filtration...(based on read counts data)")
#keep genes with 5 reads in >90% samples
message("-- >> genes with >5 reads in >90% samples")
keep = (rowMeans(cts>5) > 0.1)
dim(cts); cts <- cts[keep,]; dim(cts);

message("-- >> Variance > 1")
cts = cts[apply(cts, 1, var, na.rm=TRUE) > 1,]; dim(cts);

message("-- >> TPM data filtration...")
dim(tpm); tpm = tpm[rownames(cts), colnames(cts)]; dim(tpm);


##===================================
# 3. RLE plot
##===================================
message("[INFO] 3. Generating RLE plot...")
# RLE: For each gene and each sample, ratios are calculated between the expression of a gene and the median expression of this gene across all samples of the experiment.
# For each sample, these relative expression values are displayed as a box plot.
# Since it is assumed that in most experiments only relatively few genes are differentially expressed, the boxes should be similar in range and be centered close to 0.
# Two effects may characterize arrays with lower quality: 1) the spread is greater than that of other arrays from this experiment, and 2) the box is not centered near 0.

logtpm = log10(tpm + 1e-4)  # so row value of 0 will be -2 in the transformed value
rle_m = logtpm - apply(logtpm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.

rle = melt(cbind(ID = rownames(rle_m), rle_m), variable.name = "Sample", value.name = "tpm_val", id = "ID")
bymedian <- with(rle, reorder(Sample, tpm_val, IQR))  # sort by IQR

#########
bymedian_lvl = levels(bymedian)
outputfile = paste0(output_folder, "/3.1.QC.RLE_multi_color.pdf")
pdf(outputfile, width = 100, height = 20)

for(i in c('CONDITION','sex','Mutation','Plate')){
  x_i = covarianceTable[,i]
  myLabels <- c(sort(unique(x_i)))
  myColors <- c(rainbow(length(unique(myLabels))))
  condition_ordered = x_i[match(bymedian_lvl, rownames(covarianceTable))]
  boxColors <- myColors[match(condition_ordered, myLabels, nomatch = 1)]

  op = par(mar = c(7, 3, 3, 1))
  boxplot(tpm_val ~ bymedian, data = rle, outline = F,
          las = 2, boxwex = 1, col = boxColors,
          cex.axis = 0.2, main = paste0("Relative Log Expression (color: ",i,")"),
          xlab = "", ylab = "RLE", frame = F, axes=FALSE)

  ## set the lable ticks color
  axis(2,cex.axis = 1, font =20)
  Map(axis, side=1,  at=1:length(boxColors), col.axis=boxColors, labels=bymedian_lvl, lwd=0, las=2, cex.axis = 0.8,font=1)
  axis(1,at=1:length(boxColors),labels=FALSE)

  legend("bottomleft",
         c(paste0("--",i,"--"), sort(unique(x_i))),
         text.col = c("black",myColors),
         bty = 'n', cex = .5)
  # ggplot(rle, aes(x=bymedian, y=tpm)) +
  #   geom_boxplot(fill=boxColors, alpha=0.9)+
  #   theme(axis.text.y = element_text(color=boxColors,size = 1))
}
dev.off()

##===================================
# 3.1. RLE plot on right 100 samples
##===================================
message("[INFO] 3.1 Generating RLE plot...Zoom into the 100 on right side")
## Zoom into the 100 on right side
r_n = 100
bymedian_lvl_R = bymedian_lvl[c((length(bymedian_lvl)-r_n+1):length(bymedian_lvl))]
rle_R = rle_m[,colnames(rle_m) %in% bymedian_lvl_R]

rle_R = melt(cbind(ID = rownames(rle_R), rle_R), variable.name = "Sample", value.name = "tpm_val", id = "ID")
bymedian_R <- with(rle_R, reorder(Sample, tpm_val, IQR))  # sort by IQR

outputfile = paste0(output_folder, "/3.2.QC.RLE_multi_color_",r_n,".pdf")
pdf(outputfile, width = 8, height = 8)

for(i in c('CONDITION','sex','Mutation','Plate')){
  x_i = covarianceTable[,i]
  myLabels <- c(sort(unique(x_i)))
  myColors <- c(rainbow(length(unique(myLabels))))
  condition_ordered = x_i[match(bymedian_lvl_R, rownames(covarianceTable))]
  boxColors <- myColors[match(condition_ordered, myLabels, nomatch = 1)]
  op = par(mar = c(7, 3, 3, 1))
  boxplot(tpm_val ~ bymedian_R, data = rle_R, outline = F,
          las = 2, boxwex = 1, col = boxColors,
          main = paste0("Relative Log Expression (100 samples on most right side, color: ",i,")"),
          xlab = "", ylab = "RLE", frame = F,axes=FALSE)

  ## set the lable ticks color
  axis(2,cex.axis = 1, font =14)
  Map(axis, side=1,  at=1:length(boxColors), col.axis=boxColors, labels=bymedian_lvl_R, lwd=0, las=2,cex.axis = 0.5,font=1)
  axis(1,at=1:length(boxColors),labels=FALSE)

  legend("bottomleft",
        c(paste0("--",i,"--"), sort(unique(x_i))),
         text.col = c("black",myColors),
         bty = 'n', cex = .5)

  # ggplot(rle, aes(x=bymedian, y=tpm)) +
  #   geom_boxplot(fill=boxColors, alpha=0.9)+
  #   theme(axis.text.y = element_text(color=boxColors,size = 1))
}
dev.off()

# The other graphical representation (NUSE) represents normalized standard error (SE) estimates from the PLM fit.
# The SE estimates are normalized such that for each probe set, the median standard error across all arrays is equal to 1.
# A box plot of NUSE values is drawn for each array.
# On the NUSE plot, arrays with lower quality will have boxes that are centered higher and/or have a larger spread than the other good quality arrays from the same experiment.
# Typically, boxes centered above 1.1 represent arrays that have quality problems which are often detected in several of the other QC measures presented in this chapter.


##===================================
# 4. Clustering
##===================================
message("[INFO] 4. Clustering...")
sampleDists = 1 - cor(tpm, method = 'spearman')

outputfile = paste0(output_folder, "/4.1.QC.clustring_hc.pdf")
pdf(outputfile, width = 100, height = 20)
message("generating clustering plot...")

hc = hclust(as.dist(sampleDists), method = "complete")
plot(hc, cex = 0.7, xlab = '', main = "Cluster Dendrogram (dis = 1 - Spearman_rank_correlation, linkage = complete)")
dev.off()

outputfile = paste0(output_folder, "/4.2.QC.clustring_multi_color_fan.pdf")
pdf(outputfile, width = 7, height = 7)
subject = hc$labels
tree = as.phylo(hc)

for(i in c('CONDITION','sex','Mutation','Plate')){
  x_i = covarianceTable[,i]

  ## Update: fix edge.color bug. See https://stackoverflow.com/a/22102420/951718
  myLabels <- c('node', sort(unique(x_i)))
  myColors <- c("black", rainbow(length(unique(myLabels)) - 1))
  ## match colors and labels (nomatch == node => select idx 1)
  ## (myLabels are reordered by edge ordering
  edgeColors <- myColors[match(x_i[tree$edge[, 2]], myLabels, nomatch = 1)]
  tipColors <- myColors[match(x_i, myLabels, nomatch = 1)]

  par(mar = c(1, 1, 1, 1))
  plot(tree, type = "fan",
       cex = .3, lab4ut = 'axial', underscore = T,
       tip.color = tipColors,
       edge.color = edgeColors,
       main = paste0("Clustering of samples based on Spearman correlation (color: ",i,")")
       )
  legend("bottomleft",
        c(paste0("--",i,"--"), sort(unique(x_i))),
         text.col = c(myColors),
         bty = 'n', cex = .5)
}
dev.off()

outputfile = paste0(output_folder, "/4.3.QC.clustring_multi_color_regular.pdf")
pdf(outputfile, width = 100, height = 20)
subject = hc$labels
tree = as.phylo(hc)

for(i in c('CONDITION','sex','Mutation','Plate')){
  x_i = covarianceTable[,i]

  ## Update: fix edge.color bug. See https://stackoverflow.com/a/22102420/951718
  myLabels <- c('node', sort(unique(x_i)))
  myColors <- c("black", rainbow(length(unique(myLabels)) - 1))
  ## match colors and labels (nomatch == node => select idx 1)
  ## (myLabels are reordered by edge ordering
  edgeColors <- myColors[match(x_i[tree$edge[, 2]], myLabels, nomatch = 1)]
  tipColors <- myColors[match(x_i, myLabels, nomatch = 1)]

  par(mar = c(1, 1, 1, 1))
  plot(tree, direction = "downwards",
       cex = .3, lab4ut = 'axial', underscore = T,
       tip.color = tipColors,
       edge.color = edgeColors,
       main = paste0("Clustering of samples based on Spearman correlation (color: ",i,")")
       )
  legend("bottomleft",
        c(paste0("--",i,"--"), sort(unique(x_i))),
         text.col = c(myColors),
         bty = 'n', cex = .5)
}
dev.off()

##===================================
# 5. D-statistic
##====================================
message("[INFO] 5. D-statistic...")
outputfile = paste0(output_folder, "/5.QC.D-statistic.pdf")
pdf(outputfile, width = 7, height = 7)
message("generating D-statistic plot...")

# D-statistic
par(op)
D = apply(1 - sampleDists, 1, median)
hist(D, breaks = 100, ylab = "Number of samples", xlab = "D-statistic", main = "Histogram of D-statistic")
cutoffD = quantile(D, probs = 0.05) # 5% quantitle
if (sum(D < cutoffD)) legend("topleft", paste(names(sort(D[which(D < cutoffD)])), round(sort(D[which(D < cutoffD)]), 2)), bty = 'n', cex = .5)

dev.off()

D_lt85 = names(D[D<0.85])
# Chart

##===================================
# 6. Venn plot: RLE_r100 & D<0.85
##====================================
message("[INFO] 6. Venn: RLE_r100 & D<0.85...")
library('eulerr')
set1 <- bymedian_lvl_R
set2 <- D_lt85
myCol <- c(rainbow(2))
outputfile = paste0(output_folder, "/6.QC.Venn.pdf")
pdf(outputfile, width = 7, height = 7)

plot(euler(list("RLE_100"=set1,'D_lt85'=set2), shape = "ellipse"), quantities = TRUE)
dev.off()

##===================================
# 7. Check gender-match
##====================================
message("[INFO] 7. Generating gender-match plot...")

# Gender-specific expression
# chrY-unique regions (chrY - pseudoautosomal regions)
# ref: http://useast.ensembl.org/info/genome/genebuild/assembly.html
# grep -w chrY gencode.v19.annotation.bed12 | intersectBed -a - -b chrY.PAR.bed -v -wa | cut -f4 | sed 's/__/\t/g' | cut -f2 | sort -u > chrY-unique.geneID.txt
#chrY=read.table('/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/chrY-unique.geneID.txt')
#chrY=logtpm[chrY$V1,]
#chrY=apply(chrY, 2, mean)
#chrX=logtpm[chrX,]
chrX = 'ENSG00000229807'  # XIST
# chrY = 'ENSG00000129824'  # RPS4Y1
chrY = 'ENSG00000067048'  # DDX3Y

sex = covarianceTable$sex
if(any(grepl(chrX, rownames(logtpm))) && any(grepl(chrY, rownames(logtpm))))
{
    outputfile = paste0(output_folder, "/7.QC.gender.pdf")
    pdf(outputfile, width = 7, height = 7)
    ind = c(grep(chrX, rownames(logtpm))[1], grep(chrY, rownames(logtpm))[1])

    d = as.data.frame(t(logtpm[ind,])); colnames(d) = c("chrX", "chrY"); d$SEX = sex;
    with(d,
         plot(chrX, chrY,
              xlab = "Expression of XIST",
              ylab = "Expression of DDX3Y",
              col = 'white',
              bg = ifelse(SEX == "Female", 'red', 'blue'),
              pch = 21, bty = "n",
              main = "Gender-specific expression"
         )
    )
    #if(nrow(subset(d, chrX>0 & chrX<1.5 & chrY<0.2))>0) text(subset(d, chrX>0 & chrX<1.5 & chrY<0.2), rownames(subset(d, chrX>0 & chrX<1.5 & chrY<0.2)),pos=2, cex=0.5)
    mean2sd_chrX = mean(d$chrX[sex == "Female"]) - sd(d$chrX[sex == "Female"])
    mean2sd_chrY = mean(d$chrY[sex == "Male"]) - sd(d$chrY[sex == "Male"])
    if (nrow(subset(d, (SEX == "Female" & chrX < mean2sd_chrX) | (SEX == "Male" & chrY < mean2sd_chrY))) > 0) {
      text(subset(d, (SEX == "Female" & chrX < mean2sd_chrX) | (SEX == "Male" & chrY < mean2sd_chrY)),
           gsub("-BLM0T1", "", rownames(subset(d, (SEX == "Female" & chrX < mean2sd_chrX) | (SEX == "Male" & chrY < mean2sd_chrY)))),
           pos = 2, cex = 0.5)
    }
    legend('bottomleft', pch = 21, c("Female", "Male"), col = 'white', pt.bg = c("red", "blue"), bty = 'n', cex = 1)

    dev.off()
}

##===================================
# 8. Plot Case/Control in plate
##====================================
message("[INFO] 8. Plot Case/Control in plate plot...")
outputfile = paste0(output_folder, "/8.QC.plate_carrier_non.pdf")
pdf(outputfile, width = 7, height = 7)
ggplot(covarianceTable, aes(x=Plate, fill=CONDITION)) + geom_histogram(stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##===================================
# 9. preliminary: PCA
##====================================
message("[INFO] 9. PCA plot...")
pca_res = prcomp(t(log(tpm+1,2)))
pdf( paste0(output_folder,"/9.QC.PCA_plot.pdf"), 8, 6)
screeplot(pca_res, type = "l", npcs = 10, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"), col=c("red"), lty=5, cex=0.6)

cumpro <- cumsum(pca_res$sdev^2 / sum(pca_res$sdev^2))
plot(cumpro[0:20], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 10, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC10"), col=c("blue"), lty=5, cex=0.6)

## PCA plot
library("factoextra")
fviz_pca_ind(pca_res, geom.ind = "point", pointshape = 21,
             pointsize = 2,
             fill.ind = covarianceTable$CONDITION,
             col.ind = "black",
             palette = "jco",
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Group") +
  ggtitle("PCA-plot of PPMI using log2(TPM+1)") +
  theme(plot.title = element_text(hjust = 0.5))

####
percentVar = round(100 * pca_res$sdev^2/ sum(pca_res$sdev^2))
pcdata = as.data.frame(pca_res$x[,1:2])
pcdata = cbind(pcdata,covarianceTable$CONDITION,covarianceTable$sex,covarianceTable$race)
mycolor <- c(brewer.pal(9,"Set1"),brewer.pal(8, "Set2"),brewer.pal(12,"Set3")[c(-2,-12)],brewer.pal(12,"Paired"))

for (i in c(3,4,5)){
  pcdata_clean = pcdata[,c(1:2,i)]
  title = colnames(pcdata_clean)[3]
  print(title)
  colnames(pcdata_clean)[3] = "Category"
  p <- ggplot(data = pcdata_clean, aes(x = PC1, y = PC2, color = Category)) + geom_point(size=1,shape=16)
  p = p + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
  p <- p + theme_bw() + scale_color_manual(values = mycolor[1:20])
  p <- p + theme(axis.title.x = element_text(size = 10, color = "black", angle = 0))
  p <- p + theme(axis.title.y = element_text(size = 10, color = "black", angle = 90))
  p <- p + theme(axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))
  p <- p + theme(axis.text.y = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))
  p = p + ggtitle("PCA-plot of PPMI using log2(TPM+1)") + theme(plot.title = element_text(hjust = 0.5))

  p1 <- p + stat_ellipse(level = 0.95, linetype = 4)
  # p2 <- p + geom_text(aes(x = PC1, y = PC2, label= File_ID))

  print(p1)
  # print(p2)
}
dev.off()

##===================================
# 10. preliminary: tSNE plot
##====================================
message("[INFO] 10. tSNE plot...")
library(Rtsne)
tsne <- Rtsne(t(tpm), dims = 2, perplexity=30, verbose=TRUE, max_iter = 1000)
pdf (paste0(output_folder,"/10.QC.tSNE.pdf"), 5, 5)
col_sex = rainbow(length(unique(covarianceTable$sex)))
names(col_sex) = unique(covarianceTable$sex)
plot(tsne$Y, col=col_sex[covarianceTable$sex], pch = 16, cex = 0.5, main = "tSNE - Sex")

col_condition = rainbow(length(unique(covarianceTable$CONDITION)))
names(col_condition) = unique(covarianceTable$CONDITION)
plot(tsne$Y, col=col_condition[covarianceTable$CONDITION], pch = 16, cex = 0.5,main = "tSNE - Condition")
dev.off()

message("QC done")
##End
