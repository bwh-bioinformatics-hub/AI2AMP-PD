#!/usr/bin/env bash
###########################################
# Shell script to run the amppd pipeline
###########################################

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
#Rscript QC.R $expr_file_train $mate_file_train $prefix
#
###DE analysis
#Rscript DE.R $expr_file_train $mate_file_train $prefix
#
###Enrichment abalysis
#Rscript Enrich.R $prefix

Rscript LR.R $expr_file_train $mate_file_train $expr_file_test $mate_file_test $prefix $feature_select
Rscript RF.R $expr_file_train $mate_file_train $expr_file_test $mate_file_test $prefix $feature_select
#Rscript SVM.R $expr_file_train $mate_file_train $expr_file_test $mate_file_test $prefix $feature_select
