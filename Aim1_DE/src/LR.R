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

output_dir = paste0("../results/", prefix,"/MultiG/LR_",feature_selection)
# Create folder if the directory doesn't exist
file.exists(output_dir) || dir.create(output_dir, recursive = T)

DE_gene_file=paste0("../results/",prefix,"/DE/DEresult.padj05.resLFC.xls")
DE_gene_data = read.delim(DE_gene_file, stringsAsFactors=F, row.names = 1, header=T, check.names =F)
feature_ls = rownames(DE_gene_data)

# data
expr_data_train = read.table(paste0("../run_inout/",expr_data_train),sep="\t",head=T,row.names=1,check.names=F)
expr_data_test = read.table(paste0("../run_inout/",expr_data_test),sep="\t",head=T,row.names=1,check.names=F)
meta_data_train = read.table(paste0("../run_inout/",meta_data_train),sep="\t",head=T,row.names=1,check.names=F)
meta_data_test = read.table(paste0("../run_inout/",meta_data_test),sep="\t",head=T,row.names=1,check.names=F)

train_data = as.data.frame(t(expr_data_train))[,feature_ls]
train_data$Y = ifelse(meta_data_train$case_control_other_latest =='Case', 1, 0)

test_data = as.data.frame(t(expr_data_test))[,feature_ls]
test_data$Y = ifelse(meta_data_test$case_control_other_latest =='Case', 1, 0)


pdf.dir <- paste0(output_dir, "/LogisticRgress.pdf")
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

train_dataset = train_data[,c(marker.p,"Y")]
## Do normalization on dataset
train_dataset <- log((train_dataset+1), 2)
train_dataset <- as.data.frame(lapply(train_dataset, norm_minmax))

test_dataset = test_data[,c(marker.p,"Y")]
test_dataset <- log((test_dataset+1), 2)
test_dataset <- as.data.frame(lapply(test_dataset, norm_minmax))

## Define the model
model <- glm(Y ~.,family=binomial(link='logit'),data=train_dataset)
summary(model)

fitted.results <- predict(model,newdata=test_dataset,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_dataset$Y)
print(paste('Accuracy',1-misClasificError))

pR2(model)

step(model, test="LRT")

pred <- predict(model, newdata=test_dataset, type="response")
pr <- prediction(pred, test_dataset$Y)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
lines(0:1, 0:1, type = "l", lty = 2)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
print(auc)
txt <- c(paste("AUC=", round(auc*100, 4),"%"))
legend("bottomright", txt)

dev.off()