library(stats)
library(ROCR)
library(pscl)

##For repeatability of samples
set.seed(0)

## min max scaling/normalization
norm_minmax <- function(x){
  (x- min(x)) /(max(x)-min(x))
}

args<-commandArgs(TRUE)

expr_data_train = args[1]  # either filename or stdin
meta_data_train = args[2]
expr_data_test = args[3]
meta_data_test = args[4]
prefix = args[5]
feature_selection = args[6]

cat("======================\n","Train set:\n",expr_data_train,"\n",meta_data_train,"\n")
cat("======================\n","Test set:\n",expr_data_test,"\n",meta_data_test,"\n")
cat("======================\n","Feature selection:\n",feature_selection,"\n")
cat("======================\n","Run model Random forest\n")

cv.fold = 5
cv.time = 5 
scale = "seq"
cv.step = 1

###############
# expr_data_train = "CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest_tpm_Train.txt"
# meta_data_train = "CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest_cov_Train.txt"
# expr_data_test = "CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest_tpm_Test.txt"
# meta_data_test = "CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest_cov_Test.txt"
# prefix ="CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest"
# feature_selection ="lasso"

output_dir = paste0("../results/", prefix,"/MultiG/RF_",feature_selection)
# Create folder if the directory doesn't exist
file.exists(output_dir) || dir.create(output_dir, recursive = T)

DE_gene_file=paste0("../results/",prefix,"/DE/DEresult.all.xls.gz")
DE_gene_data = read.delim(DE_gene_file, stringsAsFactors=F, row.names = 1, header=T, check.names =F)
DE_gene_data = DE_gene_data[(DE_gene_data$padj < 0.1) & (!is.na(DE_gene_data$padj)),]

feature_ls = rownames(DE_gene_data)
cat("Total feature number:", length(feature_ls),"\n")

# data
expr_data_train = read.table(paste0("../run_inout/",expr_data_train),sep="\t",head=T,row.names=1,check.names=F)
expr_data_test = read.table(paste0("../run_inout/",expr_data_test),sep="\t",head=T,row.names=1,check.names=F)
meta_data_train = read.table(paste0("../run_inout/",meta_data_train),sep="\t",head=T,row.names=1,check.names=F)
meta_data_test = read.table(paste0("../run_inout/",meta_data_test),sep="\t",head=T,row.names=1,check.names=F)

cat("Train sample size:", dim(expr_data_train),"\n")
cat("Test sample size:", dim(expr_data_test),"\n")

expr_data_train = expr_data_train[,rownames(meta_data_train)]
train_data = as.data.frame(t(expr_data_train))[,feature_ls]
train_data$Y = ifelse(meta_data_train$case_control_other_latest =='Case', 1, 0)

expr_data_test = expr_data_test[,rownames(meta_data_test)]
test_data = as.data.frame(t(expr_data_test))[,feature_ls]
test_data$Y = ifelse(meta_data_test$case_control_other_latest =='Case', 1, 0)


# package
library(randomForest)
# function
source("functions_lib.R")

# data
train.x <- train_data[, !names(train_data) %in% c("Y")]
train.y <- as.factor(train_data$Y)
train.l <- levels(train.y)
levels(train.y) <- 0:1

test.x <- test_data[, !names(test_data) %in% c("Y")]
test.y <- as.factor(test_data$Y)
test.l <- levels(test.y)
levels(test.y) <- 0:1

pdf.dir <- paste0(output_dir, "/randomForest.pdf")
pdf(pdf.dir, width = 7 , height = 7)

source("functions_lib.R")
#Feature selection cross-validation
if(feature_selection=="mrmr"){
  marker.p <- mRMR_features(train_data, output_dir, cv.fold = cv.fold,cv.time=cv.time, scale = scale, step = cv.step,marker.num=0)
}else if(feature_selection=="lasso"){
  marker.p <- LASSO_features(train_data, output_dir)
}else{
  marker.p = feature_ls
}
cat("Selected feature number:", length(marker.p),"\n")
write.table(marker.p,paste0(output_dir, "/RandomForest_Features.txt"), sep = "\t", quote = F, col.names = NA)

# train model
train.rf <- randomForest(train.x[, marker.p], train.y, importance = T)

#########################################
##train boxplot
train.p <- predict(train.rf, type = "prob")
boxplot(train.p[, 2] ~ train.y, col = c("#F8766D","#00BFC4"),boxwex=0.6, ylab = "Probability", names = train.l,main="Train Set")
#########################################
###predict result
pr.dir <- paste0(output_dir, "/randomForest_train_probability.txt")

train.pp=array(,c(nrow(train.p),3))

for(i in 1:nrow(train.p)){
        train.pp[i,1]=train.p[i,2]
        if(train.pp[i,1]>=0.5){
                train.pp[i,2]=1
        }
        else{
                train.pp[i,2]=0
        }
        train.pp[i,3]=as.numeric(as.character(train.y[i]))
}
rownames(train.pp)=rownames(train.p)
colnames(train.pp)=c("predict_probability","predict_group","true_group")

write.table(train.pp, pr.dir, sep = "\t", quote = F, col.names = NA)

# train ROC
plot_roc(train.y, train.p[, 2],main="Train set")

# test predict
test.p <- predict(train.rf, test.x, type = "prob")
pr.dir <- paste0(output_dir, "/randomForest_test_probability.txt")

test.pp=array(,c(nrow(test.p),3))
for(i in 1:nrow(test.p)){
        test.pp[i,1]=test.p[i,2]
        if(test.pp[i,1]>=0.5){
                test.pp[i,2]=1
        }
        else{
                test.pp[i,2]=0
        }
        test.pp[i,3]=as.numeric(as.character(test.y[i]))
}
rownames(test.pp)=rownames(test.p)
colnames(test.pp)=c("predict_probability","predict_group","true_group")
write.table(test.pp, pr.dir, sep = "\t", quote = F, col.names = NA)

# predict plot
p.col <- ifelse(is.na(test.y), "#00BA38", ifelse(test.y==0,"#F8766D","#00BFC4"))

plot(rank(test.p[, 2]), test.p[, 2], col = p.col, pch = 16, xlab = "", ylab = "Probability", main = "Test set")
txt <- train.l
if (length(test.l) > 2) {
  txt <- c(txt, "the rest")
}
legend("bottomright", txt, col = c("#F8766D","#00BFC4"), pch = 16)
abline(h = 0.5)

# test ROC
plot_roc(test.y, test.p[, 2],main="Test set")
dev.off()

#####################################importance
imp=train.rf$importance
pdf(paste0(output_dir,"/RF_importance.pdf"),10,4)

par(mar=c(5,25,2,2)+0.01)
barplot(imp[,3],space=0.5,width=0.5,horiz=T,yaxt="n",xlab="Mean Decrease Accuracy",col="skyblue",border="skyblue")
axis(side=2,at=seq(0.5,nrow(imp[,3,drop=F])*0.75,0.75),labels=rownames(imp),las=1,cex.axis=0.8)

dev.off()

