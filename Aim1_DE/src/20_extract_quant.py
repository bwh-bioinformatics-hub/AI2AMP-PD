import pandas as pd

## load gene quantification data
gene_read_df = pd.read_csv("../../data/RNAseq/gene_expr_matrix_reads_row_genes.txt", index_col=0, header=0, sep="\t")
print(gene_read_df.shape)

gene_tpm_df = pd.read_csv("../../data/RNAseq/gene_expr_matrix_tpm_row_genes.txt", index_col=0, header=0, sep="\t")
print(gene_tpm_df.shape)

annotation_df = pd.read_csv("../../data/RNAseq/genes_annotation.tsv", index_col=None, header=0, sep="\t")
protein_coding_gene_df = annotation_df.loc[annotation_df['gene_type'] == "protein_coding", :]
protein_coding_gene_id = list(protein_coding_gene_df["gene_id"])

######################
Study = "PPMI"
output_folder = "../" + Study

full_cov_table_df = pd.read_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlA_cov.tsv", index_col=None, header=0, sep="\t")
CaseA_CtrlA_samples = list(full_cov_table_df['sample_id'])

## Ctrl_CN
condition = (full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"].isin(["Carrier","NonCarrier"]))
full_cov_table_ctrlCN_df = full_cov_table_df.loc[condition,]
full_cov_table_ctrlCN_df.to_csv(output_folder+"/run_inout/"+Study+"_CtrlCN_cov.tsv", sep="\t", index=False)
ctrlCN_samples = list(full_cov_table_ctrlCN_df['sample_id'])

## CaseC vs CtrlN
condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="Carrier")) | \
            ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
full_cov_table_CaseC_CtrlN_df = full_cov_table_df.loc[condition,]
full_cov_table_CaseC_CtrlN_df.to_csv(output_folder+"/run_inout/"+Study+"_CaseC_CtrlN_cov.tsv", sep="\t", index=False)
CaseC_CtrlN_samples = list(full_cov_table_CaseC_CtrlN_df['sample_id'])

## CaseN vs CtrlN
condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="NonCarrier")) | \
            ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
full_cov_table_CaseN_CtrlN_df = full_cov_table_df.loc[condition,]
full_cov_table_CaseN_CtrlN_df.to_csv(output_folder+"/run_inout/"+Study+"_CaseN_CtrlN_cov.tsv", sep="\t", index=False)
CaseN_CtrlN_samples = list(full_cov_table_CaseN_CtrlN_df['sample_id'])

## CaseA vs CtrlN
condition = (full_cov_table_df["case_control_other_latest"]=="Case")  | \
            ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
full_cov_table_CaseA_CtrlN_df = full_cov_table_df.loc[condition,]
full_cov_table_CaseA_CtrlN_df.to_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlN_cov.tsv", sep="\t", index=False)
CaseA_CtrlN_samples = list(full_cov_table_CaseA_CtrlN_df['sample_id'])

## Case_rmNA vs CtrlN
condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]!="NoValue") ) | \
            ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
full_cov_table_Case_rmNA_CtrlN_df = full_cov_table_df.loc[condition,]
full_cov_table_Case_rmNA_CtrlN_df.to_csv(output_folder+"/run_inout/"+Study+"_Case_rmNA_CtrlN_cov.tsv", sep="\t", index=False)
Case_rmNA_CtrlN_samples = list(full_cov_table_Case_rmNA_CtrlN_df['sample_id'])



#######################
## 1.1 All gene
# Case All - Ctrl All
RNAseq = gene_read_df.loc[:, CaseA_CtrlA_samples]
print("CA_CA:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlA_expr_matrix_reads_all.txt", sep="\t")

## CaseC_CtrlN
RNAseq = gene_read_df.loc[:, CaseC_CtrlN_samples]
print("CC_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseC_CtrlN_expr_matrix_reads_all.txt", sep="\t")

## CaseN_CtrlN
RNAseq = gene_read_df.loc[:, CaseN_CtrlN_samples]
print("CN_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseN_CtrlN_expr_matrix_reads_all.txt", sep="\t")

# Ctrl_samples
RNAseq = gene_read_df.loc[:, ctrlCN_samples]
print("Ctrl_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CtrlCN_expr_matrix_reads_all.txt", sep="\t")

## CaseA_CtrlN
RNAseq = gene_read_df.loc[:, CaseA_CtrlN_samples]
print("CA_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlN_expr_matrix_reads_all.txt", sep="\t")

## Case_rmNA_CtrlN
RNAseq = gene_read_df.loc[:, Case_rmNA_CtrlN_samples]
print("C_rmNA_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Case_rmNA_CtrlN_expr_matrix_reads_all.txt", sep="\t")


## TPM
#######################
## 2.1 All gene
# Case All - Ctrl All
RNAseq = gene_tpm_df.loc[:, CaseA_CtrlA_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlA_expr_matrix_tpm_all.txt", sep="\t")

## CaseC_CtrlN
RNAseq = gene_tpm_df.loc[:, CaseC_CtrlN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseC_CtrlN_expr_matrix_tpm_all.txt", sep="\t")

## CaseN_CtrlN
RNAseq = gene_tpm_df.loc[:, CaseN_CtrlN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseN_CtrlN_expr_matrix_tpm_all.txt", sep="\t")

# Ctrl_samples
RNAseq = gene_tpm_df.loc[:, ctrlCN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CtrlCN_expr_matrix_tpm_all.txt", sep="\t")

## CaseA_CtrlN
RNAseq = gene_tpm_df.loc[:, CaseA_CtrlN_samples]
print("CA_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlN_expr_matrix_tpm_all.txt", sep="\t")

## Case_rmNA_CtrlN
RNAseq = gene_tpm_df.loc[:, Case_rmNA_CtrlN_samples]
print("C_rmNA_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Case_rmNA_CtrlN_expr_matrix_tpm_all.txt", sep="\t")

