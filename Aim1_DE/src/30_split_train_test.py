import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split


output_folder = "../run_inout"

################################################################################################
# expr_read_file = "../PPMI/run_inout/PPMI_CaseAll_CtrlNonCarrier_expr_matrix_reads_all.txt"
# expr_tpm_file = "../PPMI/run_inout/PPMI_CaseAll_CtrlNonCarrier_expr_matrix_tpm_all.txt"
# meta_file = "../PPMI/run_inout/PPMI_CaseAll_CtrlNonCarrier_covariates.tsv"
# prefix = "PPMI_CaseAll_CtrlNonCarrier_allgene_DEanalysis"
#
# expr_read_df = pd.read_csv(expr_read_file,sep="\t", index_col=0, header=0)
# expr_tpm_df = pd.read_csv(expr_tpm_file,sep="\t", index_col=0, header=0)
# meta_df = pd.read_csv(meta_file,sep="\t", index_col=0, header=0)
#
# meta_df = meta_df.sample(frac= 1.0)
#
# meta_df_train,meta_df_test = train_test_split(meta_df,test_size=0,shuffle=True,random_state=0)
# # meta_df_train,meta_df_test = meta_df, meta_df
#
# print(meta_df_train["case_control_other_latest"].value_counts())
# print(meta_df_test["case_control_other_latest"].value_counts())
#
# expr_tpm_df_train = expr_tpm_df.loc[:,meta_df_train.index]
# expr_tpm_df_test = expr_tpm_df.loc[:,meta_df_test.index]
#
# expr_read_df_train = expr_read_df.loc[:,meta_df_train.index]
# expr_read_df_test = expr_read_df.loc[:,meta_df_test.index]
#
# meta_df_train.to_csv(output_folder + "/" + prefix+"_cov_Train.txt",sep="\t")
# meta_df_test.to_csv(output_folder + "/" + prefix+"_cov_Test.txt",sep="\t")
#
# expr_tpm_df_train.to_csv(output_folder + "/" + prefix+"_tpm_Train.txt",sep="\t")
# expr_tpm_df_test.to_csv(output_folder + "/" + prefix+"_tpm_Test.txt",sep="\t")
#
# expr_read_df_train.to_csv(output_folder + "/" + prefix+"_read_Train.txt",sep="\t")
# expr_read_df_test.to_csv(output_folder + "/" + prefix+"_read_Test.txt",sep="\t")

################################################################################################
expr_read_file_train = "../PPMI/run_inout/PPMI_CaseAll_CtrlNonCarrier_expr_matrix_reads_all.txt"
expr_tpm_file_train = "../PPMI/run_inout/PPMI_CaseAll_CtrlNonCarrier_expr_matrix_tpm_all.txt"
meta_file_train = "../PPMI/run_inout/PPMI_CaseAll_CtrlNonCarrier_covariates.tsv"

expr_read_file_test = "../PDBP/run_inout/PDBP_CaseAll_CtrlNonCarrier_expr_matrix_reads_all.txt"
expr_tpm_file_test = "../PDBP/run_inout/PDBP_CaseAll_CtrlNonCarrier_expr_matrix_tpm_all.txt"
meta_file_test= "../PDBP/run_inout/PDBP_CaseAll_CtrlNonCarrier_covariates.tsv"

prefix = "CaseAll_CtrlNonCarrier_allgene_PPMITrain_PDBPTest"

####
expr_read_df_train = pd.read_csv(expr_read_file_train,sep="\t", index_col=0, header=0)
expr_tpm_df_train = pd.read_csv(expr_tpm_file_train,sep="\t", index_col=0, header=0)
meta_df_train = pd.read_csv(meta_file_train,sep="\t", index_col=0, header=0)

meta_df_train = meta_df_train.sample(frac= 1.0)

####
expr_read_df_test = pd.read_csv(expr_read_file_test,sep="\t", index_col=0, header=0)
expr_tpm_df_test = pd.read_csv(expr_tpm_file_test,sep="\t", index_col=0, header=0)
meta_df_test = pd.read_csv(meta_file_test,sep="\t", index_col=0, header=0)

meta_df_test = meta_df_test.sample(frac= 1.0)

#####
print(meta_df_train["case_control_other_latest"].value_counts())
print(meta_df_test["case_control_other_latest"].value_counts())

#####
expr_tpm_df_train = expr_tpm_df_train.loc[:,meta_df_train.index]
expr_tpm_df_test = expr_tpm_df_test.loc[:,meta_df_test.index]

expr_read_df_train = expr_read_df_train.loc[:,meta_df_train.index]
expr_read_df_test = expr_read_df_test.loc[:,meta_df_test.index]

meta_df_train.to_csv(output_folder + "/" + prefix+"_cov_Train.txt",sep="\t")
meta_df_test.to_csv(output_folder + "/" + prefix+"_cov_Test.txt",sep="\t")

expr_tpm_df_train.to_csv(output_folder + "/" + prefix+"_tpm_Train.txt",sep="\t")
expr_tpm_df_test.to_csv(output_folder + "/" + prefix+"_tpm_Test.txt",sep="\t")

expr_read_df_train.to_csv(output_folder + "/" + prefix+"_read_Train.txt",sep="\t")
expr_read_df_test.to_csv(output_folder + "/" + prefix+"_read_Test.txt",sep="\t")
