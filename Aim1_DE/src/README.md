### Source code folder
>* 1  00_process_participant.py: Merge duplicated data of each participant

>* 2  01_filtration.py: Do filtration based on RNAseq sample information for DE analysis.

>* 3  10_covariate_table.py: Build covriates table for filtered samples as DEseq2 input.

>* 4  20_extract_quant.py: Build expression tables for selected samples as DEseq2 input.

>* 5  30_split_train_test.py: Prepare the Training and Testing dataset

>* 6  QC.R, DE.R, Enrich.R, LR.R, RF.R, SVM.R can be run separately OR using run them as a pipeline using amppd_analysis.sh

>* 7  lib.R, functions_lib.R are supporting functions used in other R scripts.

### Example of amppd_analysis.sh
```shell
prefix="CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest"

expr_file_train=$prefix"_tpm_Train.txt"  # either filename or stdin
mate_file_train=$prefix"_cov_Train.txt"

expr_file_test=$prefix"_tpm_Test.txt"  # either filename or stdin
mate_file_test=$prefix"_cov_Test.txt"

###
#expr_file_train="PPMI_CaseAll_CtrlNonCarrier_allgene_DEanalysis_tpm_Train.txt"  # either filename or stdin
#mate_file_train="PPMI_CaseAll_CtrlNonCarrier_allgene_DEanalysis_cov_Train.txt"
#prefix="PPMI_CaseAll_CtrlNonCarrier_allgene_DEanalysis"

feature_select="lasso"  # mrmr OR lasso

###QC
Rscript QC.R $expr_file_train $mate_file_train $prefix
#
###DE analysis
Rscript DE.R $expr_file_train $mate_file_train $prefix
#
###Enrichment abalysis
Rscript Enrich.R $prefix

Rscript LR.R $expr_file_train $mate_file_train $expr_file_test $mate_file_test $prefix $feature_select
Rscript RF.R $expr_file_train $mate_file_train $expr_file_test $mate_file_test $prefix $feature_select
Rscript SVM.R $expr_file_train $mate_file_train $expr_file_test $mate_file_test $prefix $feature_select

```

> When running 30_split_train_test.py, the names of  training/testing files have been defined and unified. 
> * expr_file_train = $prefix"_tpm_Train.txt" 
> * mate_file_train = $prefix"_cov_Train.txt"
> * expr_file_test = $prefix"_tpm_Test.txt" 
> * mate_file_test = $prefix"_cov_Test.txt"
> 
> So only a prefix string is needed in the amppd_analysis.sh. the program can automaticly find the input files 




