import pandas as pd

## load gene quantification data
# gene_read_df = pd.read_csv("../../data/RNAseq/gene_expr_matrix_reads_row_genes.txt", index_col=0, header=0, sep="\t")
# print(gene_read_df.shape)
#
# gene_tpm_df = pd.read_csv("../../data/RNAseq/gene_expr_matrix_tpm_row_genes.txt", index_col=0, header=0, sep="\t")
# print(gene_tpm_df.shape)

######################
Study = "PDBF"
output_folder = "../" + Study

case_control_pd = pd.read_csv(output_folder+"/run_inout/"+Study+"_CaseAll_CtrlNonCarrier_covariates.tsv", index_col=None, header=0, sep="\t")
samples = list(case_control_pd['sample_id'])

control_carrier_non_pd = pd.read_csv(output_folder+"/run_inout/"+Study+"_Ctrl_Carriers_Non_covariates.tsv", index_col=None, header=0, sep="\t")
control_samples = list(control_carrier_non_pd['sample_id'])

CC_CN_pd = pd.read_csv(output_folder+"/run_inout/"+Study+"_CC_CN_covariates.tsv", index_col=None, header=0, sep="\t")
CC_CN_samples = list(CC_CN_pd['sample_id'])

CA_CA_pd = pd.read_csv(output_folder+"/run_inout/"+Study+"_Case_Ctrl_all_covariates.tsv", index_col=None, header=0, sep="\t")
CA_CA_samples = list(CA_CA_pd['sample_id'])

annotation_df = pd.read_csv("../../data/RNAseq/genes_annotation.tsv", index_col=None, header=0, sep="\t")
protein_coding_gene_df = annotation_df.loc[annotation_df['gene_type'] == "protein_coding", :]
protein_coding_gene_id = list(protein_coding_gene_df["gene_id"])

Stop
#######################
## 1.1 All gene
## Case - Control
RNAseq = gene_read_df.loc[:, samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseAll_CtrlNonCarrier_expr_matrix_reads_all.txt", sep="\t")

## Carrier - NonCarrier
RNAseq = gene_read_df.loc[:, control_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Ctrl_Carriers_Non_expr_matrix_reads_all.txt", sep="\t")

## CC - CN
RNAseq = gene_read_df.loc[:, CC_CN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CC_CN_expr_matrix_reads_all.txt", sep="\t")

## CA - CA
RNAseq = gene_read_df.loc[:, CA_CA_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Case_Ctrl_all_expr_matrix_reads_all.txt", sep="\t")

################
## 1.2 protein coding genes
## Case - Control
RNAseq = gene_read_df.loc[gene_read_df.index.isin(protein_coding_gene_id), samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseAll_CtrlNonCarrier_expr_matrix_reads_pc.txt", sep="\t")

## Carrier - NonCarrier
RNAseq = gene_read_df.loc[gene_read_df.index.isin(protein_coding_gene_id), control_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Ctrl_Carriers_Non_expr_matrix_reads_pc.txt", sep="\t")

## CC - CN
RNAseq = gene_read_df.loc[gene_read_df.index.isin(protein_coding_gene_id), CC_CN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CC_CN_expr_matrix_reads_pc.txt", sep="\t")

## CA - CA
RNAseq = gene_read_df.loc[gene_read_df.index.isin(protein_coding_gene_id), CA_CA_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Case_Ctrl_all_expr_matrix_reads_pc.txt", sep="\t")


########################
## 2.1 TPM: All genes
## Case - control
RNAseq = gene_tpm_df.loc[:, samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseAll_CtrlNonCarrier_expr_matrix_tpm_all.txt", sep="\t")

## Carrier - NonCarrier
RNAseq = gene_tpm_df.loc[:, control_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Ctrl_Carriers_Non_expr_matrix_tpm_all.txt", sep="\t")

## CC - CN
RNAseq = gene_tpm_df.loc[:, CC_CN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CC_CN_expr_matrix_tpm_all.txt", sep="\t")

## CA - CA
RNAseq = gene_tpm_df.loc[:, CA_CA_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Case_Ctrl_all_expr_matrix_tpm_all.txt", sep="\t")


################
## 2.2 TPM: Protein coding genes
## Case - control
RNAseq = gene_tpm_df.loc[gene_tpm_df.index.isin(protein_coding_gene_id), samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CaseAll_CtrlNonCarrier_expr_matrix_tpm_pc.txt", sep="\t")

## Carrier - NonCarrier
RNAseq = gene_tpm_df.loc[gene_tpm_df.index.isin(protein_coding_gene_id), control_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Ctrl_Carriers_Non_expr_matrix_tpm_pc.txt", sep="\t")

## CC - CN
RNAseq = gene_tpm_df.loc[gene_tpm_df.index.isin(protein_coding_gene_id), CC_CN_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_CC_CN_expr_matrix_tpm_pc.txt", sep="\t")

## CA - CA
RNAseq = gene_tpm_df.loc[gene_tpm_df.index.isin(protein_coding_gene_id), CA_CA_samples]
print(RNAseq.shape)
RNAseq.to_csv(output_folder+"/run_inout/"+Study+"_Case_Ctrl_all_expr_matrix_tpm_pc.txt", sep="\t")


## save data in Salmon format
import gzip
f_out = gzip.open(output_folder+"/run_inout/"+Study+"_Case_Ctrl_all_expr_matrix_Salmon.txt.gz", "wt")
with gzip.open("../../data/RNAseq/aggregated.transcripts.tsv.gz",'rt') as f:
    header = next(f)
    f_out.write(header)
    for line in f:
        line_ls = line.split("\t")
        sample_id = line_ls[1]
        if sample_id in CA_CA_samples:
            f_out.write(line)
f_out.close()