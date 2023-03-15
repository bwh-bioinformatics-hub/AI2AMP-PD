# https://support.bioconductor.org/p/49864/,https://www.biostars.org/p/128931/, qvalue, FDR
# https://github.com/YuLab-SMU/clusterProfiler/issues/124
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(KEGG.db)
library(ReactomePA)
library(DOSE)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(DESeq2)
library(enrichplot)
library('fgsea')
library(httr)

set_config(config(ssl_verifypeer = 0L))

prefix="PDBF_CaseA_CtrlA_PCs"
output_dir = paste0("../results/", prefix,"/Enrich_profiler")
# Create folder if the directory doesn't exist
file.exists(output_dir) || dir.create(output_dir, recursive = T)

input_DE_result=  paste0("../results/",prefix,"/DE/DEresult.all.xls.gz")
DE_result0=read.delim(input_DE_result, stringsAsFactors=F, row.names = 1, header=T, check.names =F)
DE_result0$gseaRANK = -log10(DE_result0$pvalue) * sign(DE_result0$log2FoldChange)

PC0 = DE_result0[DE_result0$geneType=="protein_coding",]
PC_ranked0 = PC0[order(PC0$gseaRANK,decreasing = TRUE),]
PC_ranked0$logpvalue = -log10(PC0$pvalue)

ranked_gene_list0 = PC_ranked0$gseaRANK
names(ranked_gene_list0) = as.character(PC_ranked0$symbol)
ranked_gene_list0 = sort(ranked_gene_list0, decreasing = TRUE)

DE_result=subset(DE_result0, padj<0.05) # padj<0.05
PC = DE_result[DE_result$geneType=="protein_coding",]
gene_list = PC$symbol


top_n = 30
## FDR correction
## BP
go_enrich <- enrichGO(gene = gene_list,
                      OrgDb = 'org.Hs.eg.db',
                      keyType="SYMBOL",
                      pAdjustMethod = "fdr",
                      readable = F,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

X_BP = go_enrich@result
write.table(X_BP,paste0(output_dir,"/enrich_of_genes_BP_FDR.txt"))
X_BP=X_BP[X_BP$p.adjust<0.05,]
X_BP <- X_BP[order(X_BP$p.adjust),]
X_BP <- X_BP[0:min(nrow(X_BP),top_n),]

## CC
go_enrich <- enrichGO(gene = gene_list,
                      OrgDb = 'org.Hs.eg.db',
                      keyType="SYMBOL",
                      pAdjustMethod = "fdr",
                      readable = F,
                      ont = "CC",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

X_CC = go_enrich@result
write.table(X_CC,paste0(output_dir,"/enrich_of_genes_CC_FDR.txt"))
X_CC=X_CC[X_CC$p.adjust<0.05,]
X_CC <- X_CC[order(X_CC$p.adjust),]
X_CC <- X_CC[0:min(nrow(X_CC),top_n),]


## MF
go_enrich <- enrichGO(gene = gene_list,
                      OrgDb = 'org.Hs.eg.db',
                      keyType="SYMBOL",
                      pAdjustMethod = "fdr",
                      readable = F,
                      ont = "MF",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

X_MF = go_enrich@result
write.table(X_MF,paste0(output_dir,"/enrich_of_genes_MF_FDR.txt"))
X_MF=X_MF[X_MF$p.adjust<0.05,]
X_MF <- X_MF[order(X_MF$p.adjust),]
X_MF <- X_MF[0:min(nrow(X_MF),top_n),]

## bonferroni correction
go_enrich <- enrichGO(gene = gene_list,
                       OrgDb = 'org.Hs.eg.db',
                       keyType="SYMBOL",
                       pAdjustMethod = "bonferroni",
                       readable = F,
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

 X_BP_ben = go_enrich@result
 write.table(X_BP_ben,paste0(output_dir,"/enrich_of_genes_BP_bonferroni.txt"))
 X_BP_ben=X_BP_ben[X_BP_ben$p.adjust<0.05,]
 X_BP_ben <- X_BP_ben[order(X_BP_ben$p.adjust),]
 X_BP_ben <- X_BP_ben[0:min(nrow(X_BP_ben),top_n),]
 
 ## bonferroni correction
 go_enrich <- enrichGO(gene = gene_list,
                       OrgDb = 'org.Hs.eg.db',
                       keyType="SYMBOL",
                       pAdjustMethod = "bonferroni",
                       readable = F,
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

 X_CC_ben = go_enrich@result
 write.table(X_CC_ben,paste0(output_dir,"/enrich_of_genes_CC_bonferroni.txt"))
 X_CC_ben=X_CC_ben[X_CC_ben$p.adjust<0.05,]
 X_CC_ben <- X_CC_ben[order(X_CC_ben$p.adjust),]
 X_CC_ben <- X_CC_ben[0:min(nrow(X_CC_ben),top_n),]

eg = bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)
## save 
write.table(eg,paste0(output_dir,"/gene_id_symbol_map.txt"))
## KEGG
kegg_enrich <- enrichKEGG(gene = eg$ENTREZID,
                          organism = 'human',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

X_kegg = kegg_enrich@result
write.table(X_kegg,paste0(output_dir,"/enrich_of_genes_KEGG_FDR.txt"))
X_kegg=X_kegg[X_kegg$p.adjust<0.05,]
X_kegg <- X_kegg[order(X_kegg$p.adjust),]
X_kegg <- X_kegg[0:min(nrow(X_kegg),top_n),]


source("wikiPathway_hu.R")
## WikiPathway
WP_enrich <- enrichWP_hu(gene = eg$ENTREZID,
                          organism = "Homo sapiens",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
X_WP = WP_enrich@result
write.table(X_WP,paste0(output_dir,"/enrich_of_genes_WP_FDR.txt"))
X_WP=X_WP[X_WP$p.adjust<0.05,]
X_WP <- X_WP[order(X_WP$p.adjust),]
X_WP <- X_WP[0:min(nrow(X_WP),top_n),]

####################
library(ggplot2)
# Basic barplot
top_n = 30
pdf(paste0(output_dir,"/Enrich.pdf"))
p<-ggplot(data=X_BP, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-60,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="GO - BP")
# Horizontal bar plot
p + coord_flip()

p<-ggplot(data=X_CC, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-50,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="GO - CC")
# Horizontal bar plot
p + coord_flip()

p<-ggplot(data=X_MF, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-50,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="GO - MF")
# Horizontal bar plot
p + coord_flip()

p<-ggplot(data=X_kegg, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-50,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="KEGG")
# Horizontal bar plot
p + coord_flip()

p<-ggplot(data=X_WP, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-50,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="WikiPathway")
# Horizontal bar plot
p + coord_flip()

p<-ggplot(data=X_BP_ben, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-60,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="GO - Biological processes")
# Horizontal bar plot
p + coord_flip()

p<-ggplot(data=X_CC_ben, aes(x=reorder(Description, -log(p.adjust, base=10)), y= -log(p.adjust, base=10),fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_gradient(low = "red", high = "blue",limits=c(1e-50,0.055))+
  theme(axis.text.x= element_text(colour="black", size=10))+
  theme(axis.text.y= element_text(colour="black",size=10)) +
  labs(y="-log10(p.adjust)",x="GO - Cellular components")
# Horizontal bar plot
p + coord_flip()

dev.off()



