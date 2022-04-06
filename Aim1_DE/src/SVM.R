library(stats)
library(ROCR)
library(pscl)

library(kernlab)      # SVM methodology
library(e1071)

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
cat("======================\n","Run model SVM\n")

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

output_dir = paste0("../results/", prefix,"/MultiG/SVM_",feature_selection)
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

cat("Train sample:", dim(expr_data_train),"\n")
cat("Test sample:", dim(expr_data_test),"\n")

expr_data_train = expr_data_train[,rownames(meta_data_train)]
train_data = as.data.frame(t(expr_data_train))[,feature_ls]
train_data$Y = ifelse(meta_data_train$case_control_other_latest =='Case', 1, 0)

expr_data_test = expr_data_test[,rownames(meta_data_test)]
test_data = as.data.frame(t(expr_data_test))[,feature_ls]
test_data$Y = ifelse(meta_data_test$case_control_other_latest =='Case', 1, 0)

cat("Train size:", dim(train_data),"\n")
cat("Test size:", dim(test_data),"\n")

pdf.dir <- paste0(output_dir, "/SVM.pdf")
pdf(pdf.dir, width = 7, height = 7)

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
write.table(marker.p,paste0(output_dir, "/SVM_Features.txt"), sep = "\t", quote = F, col.names = NA)

train_dataset = train_data[,c(marker.p,"Y")]
## Do normalization on dataset
train_dataset <- log((train_dataset+1), 2)
train_dataset <- as.data.frame(lapply(train_dataset, norm_minmax))

test_dataset = test_data[,c(marker.p,"Y")]
test_dataset <- log((test_dataset+1), 2)
test_dataset <- as.data.frame(lapply(test_dataset, norm_minmax))

cat("Final train size:", dim(train_dataset),"\n")
cat("Final test size:", dim(test_dataset),"\n")

##############################
## using the linear kernel
## Define the model
model <- svm(Y~., data = train_dataset, kernel = "linear", scale = FALSE,cost = 0.1,gamma=0.1)

# Plot Results
fitted.results <- predict(model,newdata=train_dataset,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != train_dataset$Y)
print(paste('Trains accuracy:',1-misClasificError))

fitted.results <- predict(model,newdata=test_dataset,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_dataset$Y)
print(paste0('Test accuracy:',1-misClasificError))

# pR2(model)
# step(model, test="LRT")

pred_train <- predict(model, newdata=train_dataset, type="response")
pr_train <- prediction(pred_train, train_dataset$Y)
prf_train <- performance(pr_train, measure = "tpr", x.measure = "fpr")

pred <- predict(model, newdata=test_dataset, type="response")
pr <- prediction(pred, test_dataset$Y)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")

plot(prf,col="blue", main="SVM: linear kernel")
plot(prf_train,add = TRUE, col="red")
lines(0:1, 0:1, type = "l", lty = 2)

auc_train <- performance(pr_train, measure = "auc")
auc_train <- auc_train@y.values[[1]]

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]

cat("AUC-Train:",auc_train,"\n","AUC-Test:",auc,"\n")

txt <- c(paste0("AUC-Test=", round(auc*100, 2),"%"), paste0("AUC-Train=", round(auc_train*100, 2),"%"))
legend("bottomright", txt, col= c("blue","red"),lwd = 1, lty = c(1,1))
#######################
## using the linear kernel
## Define the model
model <- svm(Y~., data = train_dataset, kernel = "radial", scale = FALSE,cost = 0.1,gamma=0.1)
# Plot Results0.
fitted.results <- predict(model,newdata=train_dataset,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != train_dataset$Y)
print(paste('Trains accuracy:',1-misClasificError))

fitted.results <- predict(model,newdata=test_dataset,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_dataset$Y)
print(paste0('Test accuracy:',1-misClasificError))

# pR2(model)
# step(model, test="LRT")

pred_train <- predict(model, newdata=train_dataset, type="response")
pr_train <- prediction(pred_train, train_dataset$Y)
prf_train <- performance(pr_train, measure = "tpr", x.measure = "fpr")

pred <- predict(model, newdata=test_dataset, type="response")
pr <- prediction(pred, test_dataset$Y)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")

plot(prf,col="blue", main="SVM: non-linear kernel")
plot(prf_train,add = TRUE, col="red")
lines(0:1, 0:1, type = "l", lty = 2)

auc_train <- performance(pr_train, measure = "auc")
auc_train <- auc_train@y.values[[1]]

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]

cat("AUC-Train:",auc_train,"\n","AUC-Test:",auc,"\n")

txt <- c(paste0("AUC-Test=", round(auc*100, 2),"%"), paste0("AUC-Train=", round(auc_train*100, 2),"%"))
legend("bottomright", txt, col= c("blue","red"),pch=16)

dev.off()