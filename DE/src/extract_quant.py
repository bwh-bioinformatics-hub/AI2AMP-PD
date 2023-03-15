import pandas as pd

## load gene quantification data
gene_read_df = pd.read_csv("../../../../../data/RNAseq/gene_expr_matrix_reads_row_genes.txt", index_col=0, header=0, sep="\t")
print(gene_read_df.shape)

gene_tpm_df = pd.read_csv("../../../../../data/RNAseq/gene_expr_matrix_tpm_row_genes.txt", index_col=0, header=0, sep="\t")
print(gene_tpm_df.shape)

annotation_df = pd.read_csv("../../../../../data/RNAseq/genes_annotation.tsv", index_col=None, header=0, sep="\t")
protein_coding_gene_df = annotation_df.loc[annotation_df['gene_type'] == "protein_coding", :]
protein_coding_gene_id = list(protein_coding_gene_df["gene_id"])

AJ_info = pd.read_csv("../../../../../data/variants/output_ancestry_genetic_ancestry_all_pca_plus.csv", index_col=0, header=0)
AJ_ls = AJ_info.loc[AJ_info["InfAJ"]==1,"IID"].tolist()

stop

######################
Study = "PPMI"
output_folder = "../run_inout/"
full_cov_table_df = pd.read_csv("../" + Study+"/run_inout/"+Study+"_CaseA_CtrlA_with_PCs.tsv", index_col=None, header=0, sep="\t")
# full_cov_table_df = full_cov_table_df.loc[~full_cov_table_df["participant_id"].isin(AJ_ls),:]

print("=====================================")
# full_cov_table_df.to_csv(output_folder+Study+"_CaseA_CtrlA_cov_with_PCs.tsv", sep="\t", index=False)

# CaseA_CtrlA_samples = list(full_cov_table_df['sample_id'])
# RNAseq = gene_read_df.loc[:, CaseA_CtrlA_samples]
# print("CA_CA:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseA_CtrlA_read_with_PCs.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseA_CtrlA_samples]
# print("CA_CA:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseA_CtrlA_tpm_with_PCs.tsv", sep="\t")
#
# print(full_cov_table_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))

# print("=====================================")
# ## Ctrl_CN
# condition = (full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"].isin(["Carrier","NonCarrier"]))
# full_cov_table_ctrlCN_df = full_cov_table_df.loc[condition,]
# full_cov_table_ctrlCN_df.to_csv(output_folder+Study+"_CtrlCN_cov.tsv", sep="\t", index=False)
# ctrlCN_samples = list(full_cov_table_ctrlCN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, ctrlCN_samples]
# print("Ctrl_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlCN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, ctrlCN_samples]
# print("CtrlCN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlCN_tpm.tsv", sep="\t")
#
# print(full_cov_table_ctrlCN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# print("=====================================")
# ## Case_CN
# condition = (full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"].isin(["Carrier","NonCarrier"]))
# full_cov_table_caseCN_df = full_cov_table_df.loc[condition,]
# full_cov_table_caseCN_df.to_csv(output_folder+Study+"_CaseCN_cov.tsv", sep="\t", index=False)
# caseCN_samples = list(full_cov_table_caseCN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, caseCN_samples]
# print("Case_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseCN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, caseCN_samples]
# print("CaseCN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseCN_tpm.tsv", sep="\t")
#
# print(full_cov_table_caseCN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##--------------------------------------------
# ## vs CtrlN
# print("=====================================")
# ## CaseC vs CtrlN
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="Carrier")) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# full_cov_table_CaseC_CtrlN_df = full_cov_table_df.loc[condition,]
# full_cov_table_CaseC_CtrlN_df.to_csv(output_folder+Study+"_CaseC_CtrlN_cov.tsv", sep="\t", index=False)
# CaseC_CtrlN_samples = list(full_cov_table_CaseC_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseC_CtrlN_samples]
# print("CC_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseC_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseC_CtrlN_samples]
# print("CC_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseC_CtrlN_tpm.tsv", sep="\t")
#
# print(full_cov_table_CaseC_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))

print("=====================================")
## CaseN vs CtrlN
condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="NonCarrier")) | \
            ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
full_cov_table_CaseN_CtrlN_df = full_cov_table_df.loc[condition,]
full_cov_table_CaseN_CtrlN_df.to_csv(output_folder+Study+"_CaseN_CtrlN_cov_with_PCs.tsv", sep="\t", index=False)
CaseN_CtrlN_samples = list(full_cov_table_CaseN_CtrlN_df['sample_id'])

RNAseq = gene_read_df.loc[:, CaseN_CtrlN_samples]
print("CN_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+Study+"_CaseN_CtrlN_read_with_PCs.tsv", sep="\t")

RNAseq = gene_tpm_df.loc[:, CaseN_CtrlN_samples]
print("CN_CN:",RNAseq.shape)
RNAseq.to_csv(output_folder+Study+"_CaseN_CtrlN_tpm_with_PCs.tsv", sep="\t")

print(full_cov_table_CaseN_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# print("=====================================")
# ## CaseA vs CtrlN
# condition = (full_cov_table_df["case_control_other_latest"]=="Case")  | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# full_cov_table_CaseA_CtrlN_df = full_cov_table_df.loc[condition,]
# full_cov_table_CaseA_CtrlN_df.to_csv(output_folder+Study+"_CaseA_CtrlN_cov.tsv", sep="\t", index=False)
# CaseA_CtrlN_samples = list(full_cov_table_CaseA_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseA_CtrlN_samples]
# print("CA_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseA_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseA_CtrlN_samples]
# print("CA_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseA_CtrlN_tpm.tsv", sep="\t")
#
# print(full_cov_table_CaseA_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# print("=====================================")
# ## Case_rmNA vs CtrlN
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]!="NoValue") ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# full_cov_table_Case_rmNA_CtrlN_df = full_cov_table_df.loc[condition,]
# full_cov_table_Case_rmNA_CtrlN_df.to_csv(output_folder+Study+"_Case_rmNA_CtrlN_cov.tsv", sep="\t", index=False)
# Case_rmNA_CtrlN_samples = list(full_cov_table_Case_rmNA_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, Case_rmNA_CtrlN_samples]
# print("C_rmNA_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_Case_rmNA_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, Case_rmNA_CtrlN_samples]
# print("C_rmNA_CN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_Case_rmNA_CtrlN_tpm.tsv", sep="\t")
#
# print(full_cov_table_Case_rmNA_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##--------------------------------------------
# ## vs CtrlC
# print("=====================================")
# ## CaseC vs CtrlC
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="Carrier")) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="Carrier"))
# full_cov_table_CaseC_CtrlC_df = full_cov_table_df.loc[condition,]
# full_cov_table_CaseC_CtrlC_df.to_csv(output_folder+Study+"_CaseC_CtrlC_cov.tsv", sep="\t", index=False)
# CaseC_CtrlC_samples = list(full_cov_table_CaseC_CtrlC_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseC_CtrlC_samples]
# print("CC_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseC_CtrlC_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseC_CtrlC_samples]
# print("CC_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseC_CtrlC_tpm.tsv", sep="\t")
#
# print(full_cov_table_CaseC_CtrlC_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
#
# print("=====================================")
# ## CaseN vs CtrlC
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="NonCarrier")) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="Carrier"))
# full_cov_table_CaseN_CtrlC_df = full_cov_table_df.loc[condition,]
# full_cov_table_CaseN_CtrlC_df.to_csv(output_folder+Study+"_CaseN_CtrlC_cov.tsv", sep="\t", index=False)
# CaseN_CtrlC_samples = list(full_cov_table_CaseN_CtrlC_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseN_CtrlC_samples]
# print("CN_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseN_CtrlC_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseN_CtrlC_samples]
# print("CN_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseN_CtrlC_tpm.tsv", sep="\t")
#
# print(full_cov_table_CaseN_CtrlC_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# print("=====================================")
# ## CaseA vs CtrlC
# condition = (full_cov_table_df["case_control_other_latest"]=="Case")  | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="Carrier"))
# full_cov_table_CaseA_CtrlC_df = full_cov_table_df.loc[condition,]
# full_cov_table_CaseA_CtrlC_df.to_csv(output_folder+Study+"_CaseA_CtrlC_cov.tsv", sep="\t", index=False)
# CaseA_CtrlC_samples = list(full_cov_table_CaseA_CtrlC_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseA_CtrlC_samples]
# print("CA_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseA_CtrlC_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseA_CtrlC_samples]
# print("CA_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseA_CtrlC_tpm.tsv", sep="\t")
#
# print(full_cov_table_CaseA_CtrlC_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
#
# print("=====================================")
# ## Case_rmNA vs CtrlC
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]!="NoValue") ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="Carrier"))
# full_cov_table_Case_rmNA_CtrlC_df = full_cov_table_df.loc[condition,]
# full_cov_table_Case_rmNA_CtrlC_df.to_csv(output_folder+Study+"_Case_rmNA_CtrlC_cov.tsv", sep="\t", index=False)
# Case_rmNA_CtrlC_samples = list(full_cov_table_Case_rmNA_CtrlC_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, Case_rmNA_CtrlC_samples]
# print("C_rmNA_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_Case_rmNA_CtrlC_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, Case_rmNA_CtrlC_samples]
# print("C_rmNA_CC:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_Case_rmNA_CtrlC_tpm.tsv", sep="\t")
#
# print(full_cov_table_Case_rmNA_CtrlC_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##------------------------------------------------
# ##Drug
# ## CaseCDrugY_CtrlN
# print("===================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") &
#              (full_cov_table_df["Mutation"]=="Carrier") &
#              (full_cov_table_df["treat"]==1) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseCDrugY_CtrlN_df = full_cov_table_df.loc[condition,]
# CaseCDrugY_CtrlN_df.to_csv(output_folder+Study+"_CaseCDrugY_CtrlN_cov.tsv", sep="\t", index=False)
# CaseCDrugY_CtrlN_samples = list(CaseCDrugY_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseCDrugY_CtrlN_samples]
# print("CaseCDrugY_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseCDrugY_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseCDrugY_CtrlN_samples]
# print("CaseCDrugY_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseCDrugY_CtrlN_tpm.tsv", sep="\t")
# print(CaseCDrugY_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CaseCDrugN_CtrlN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") &
#              (full_cov_table_df["Mutation"]=="Carrier") &
#              (full_cov_table_df["treat"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseCDrugN_CtrlN_df = full_cov_table_df.loc[condition,]
# CaseCDrugN_CtrlN_df.to_csv(output_folder+Study+"_CaseCDrugN_CtrlN_cov.tsv", sep="\t", index=False)
# CaseCDrugN_CtrlN_samples = list(CaseCDrugN_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseCDrugN_CtrlN_samples]
# print("CaseCDrugN_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseCDrugN_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseCDrugN_CtrlN_samples]
# print("CaseCDrugN_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseCDrugN_CtrlN_tpm.tsv", sep="\t")
# print(CaseCDrugN_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##------------------------------
# ## CaseNDrugY_CtrlN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") &
#              (full_cov_table_df["Mutation"]=="NonCarrier") &
#              (full_cov_table_df["treat"]==1) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseNDrugY_CtrlN_df = full_cov_table_df.loc[condition,]
# CaseNDrugY_CtrlN_df.to_csv(output_folder+Study+"_CaseNDrugY_CtrlN_cov.tsv", sep="\t", index=False)
# CaseNDrugY_CtrlN_samples = list(CaseNDrugY_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseNDrugY_CtrlN_samples]
# print("CaseNDrugY_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseNDrugY_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseNDrugY_CtrlN_samples]
# print("CaseNDrugY_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseNDrugY_CtrlN_tpm.tsv", sep="\t")
# print(CaseNDrugY_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CaseNDrugN_CtrlN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") &
#              (full_cov_table_df["Mutation"]=="NonCarrier") &
#              (full_cov_table_df["treat"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseNDrugN_CtrlN_df = full_cov_table_df.loc[condition,]
# CaseNDrugN_CtrlN_df.to_csv(output_folder+Study+"_CaseNDrugN_CtrlN_cov.tsv", sep="\t", index=False)
# CaseNDrugN_CtrlN_samples = list(CaseNDrugN_CtrlN_df['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseNDrugN_CtrlN_samples]
# print("CaseNDrugN_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseNDrugN_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseNDrugN_CtrlN_samples]
# print("CaseNDrugN_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseNDrugN_CtrlN_tpm.tsv", sep="\t")
# print(CaseNDrugN_CtrlN_df.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##-----------------------------------------------------------------
# ## mutation
# ## CaseGBA_CaseN
# print("====================================")
# condition = (full_cov_table_df["case_control_other_latest"]=="Case") & \
#             (((full_cov_table_df["GBA"]==1) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==0) ) |
#              (full_cov_table_df["Mutation"]=="NonCarrier")
#              )
# CaseGBA_CaseN = full_cov_table_df.loc[condition,]
# CaseGBA_CaseN.to_csv(output_folder+Study+"_CaseGBA_CaseN_cov.tsv", sep="\t", index=False)
# CaseGBA_CaseN_samples = list(CaseGBA_CaseN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseGBA_CaseN_samples]
# print("CaseGBA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseGBA_CaseN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseGBA_CaseN_samples]
# print("CaseGBA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseGBA_CaseN_tpm.tsv", sep="\t")
# print(CaseGBA_CaseN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CaseSNCA_CaseN
# print("====================================")
# condition = (full_cov_table_df["case_control_other_latest"]=="Case") & \
#             (((full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==1) & (full_cov_table_df["LRRK2"]==0) ) |
#              (full_cov_table_df["Mutation"]=="NonCarrier")
#              )
# CaseSNCA_CaseN = full_cov_table_df.loc[condition,]
# CaseSNCA_CaseN.to_csv(output_folder+Study+"_CaseSNCA_CaseN_cov.tsv", sep="\t", index=False)
# CaseSNCA_CaseN_samples = list(CaseSNCA_CaseN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseSNCA_CaseN_samples]
# print("CaseSNCA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseSNCA_CaseN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseSNCA_CaseN_samples]
# print("CaseSNCA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseSNCA_CaseN_tpm.tsv", sep="\t")
# print(CaseSNCA_CaseN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CaseLRRK2_CaseN
# print("====================================")
# condition = (full_cov_table_df["case_control_other_latest"]=="Case") & \
#             (((full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==1) ) |
#              (full_cov_table_df["Mutation"]=="NonCarrier")
#              )
# CaseLRRK2_CaseN = full_cov_table_df.loc[condition,]
# CaseLRRK2_CaseN.to_csv(output_folder+Study+"_CaseLRRK2_CaseN_cov.tsv", sep="\t", index=False)
# CaseLRRK2_CaseN_samples = list(CaseLRRK2_CaseN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseLRRK2_CaseN_samples]
# print("CaseLRRK2_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseLRRK2_CaseN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseLRRK2_CaseN_samples]
# print("CaseLRRK2_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseLRRK2_CaseN_tpm.tsv", sep="\t")
# print(CaseLRRK2_CaseN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##------------------------------------------
# ## CaseGBA_CtrlN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["GBA"]==1) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseGBA_CtrlN = full_cov_table_df.loc[condition,]
# CaseGBA_CtrlN.to_csv(output_folder+Study+"_CaseGBA_CtrlN_cov.tsv", sep="\t", index=False)
# CaseGBA_CtrlN_samples = list(CaseGBA_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseGBA_CtrlN_samples]
# print("CaseGBA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseGBA_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseGBA_CtrlN_samples]
# print("CaseGBA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseGBA_CtrlN_tpm.tsv", sep="\t")
# print(CaseGBA_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CaseSNCA_CtrlN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==1) & (full_cov_table_df["LRRK2"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseSNCA_CtrlN = full_cov_table_df.loc[condition,]
# CaseSNCA_CtrlN.to_csv(output_folder+Study+"_CaseSNCA_CtrlN_cov.tsv", sep="\t", index=False)
# CaseSNCA_CtrlN_samples = list(CaseSNCA_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseSNCA_CtrlN_samples]
# print("CaseSNCA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseSNCA_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseSNCA_CtrlN_samples]
# print("CaseSNCA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseSNCA_CtrlN_tpm.tsv", sep="\t")
# print(CaseSNCA_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CaseLRRK2_CtrlN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Case") &(full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==1) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CaseLRRK2_CtrlN = full_cov_table_df.loc[condition,]
# CaseLRRK2_CtrlN.to_csv(output_folder+Study+"_CaseLRRK2_CtrlN_cov.tsv", sep="\t", index=False)
# CaseLRRK2_CtrlN_samples = list(CaseLRRK2_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CaseLRRK2_CtrlN_samples]
# print("CaseLRRK2_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseLRRK2_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CaseLRRK2_CtrlN_samples]
# print("CaseLRRK2_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CaseLRRK2_CtrlN_tpm.tsv", sep="\t")
# print(CaseLRRK2_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##------------------------------------------
# ## CtrlGBA_CtrlN
# print("====================================")
# condition = (full_cov_table_df["case_control_other_latest"]=="Control") & \
#             (((full_cov_table_df["GBA"]==1) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==0) ) |
#              (full_cov_table_df["Mutation"]=="NonCarrier")
#              )
# CtrlGBA_CtrlN = full_cov_table_df.loc[condition,]
# CtrlGBA_CtrlN.to_csv(output_folder+Study+"_CtrlGBA_CtrlN_cov.tsv", sep="\t", index=False)
# CtrlGBA_CtrlN_samples = list(CtrlGBA_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CtrlGBA_CtrlN_samples]
# print("CtrlGBA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlGBA_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CtrlGBA_CtrlN_samples]
# print("CtrlGBA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlGBA_CtrlN_tpm.tsv", sep="\t")
# print(CtrlGBA_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CtrlSNCA_CtrlN
# print("====================================")
# condition = (full_cov_table_df["case_control_other_latest"]=="Control") & \
#             (((full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==1) & (full_cov_table_df["LRRK2"]==0) ) |
#              (full_cov_table_df["Mutation"]=="NonCarrier")
#              )
# CtrlSNCA_CtrlN = full_cov_table_df.loc[condition,]
# CtrlSNCA_CtrlN.to_csv(output_folder+Study+"_CtrlSNCA_CtrlN_cov.tsv", sep="\t", index=False)
# CtrlSNCA_CtrlN_samples = list(CtrlSNCA_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CtrlSNCA_CtrlN_samples]
# print("CtrlSNCA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlSNCA_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CtrlSNCA_CtrlN_samples]
# print("CtrlSNCA_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlSNCA_CtrlN_tpm.tsv", sep="\t")
# print(CtrlSNCA_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CtrlLRRK2_CtrlN
# print("====================================")
# condition = (full_cov_table_df["case_control_other_latest"]=="Control") & \
#             (((full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==1) ) |
#              (full_cov_table_df["Mutation"]=="NonCarrier")
#              )
# CtrlLRRK2_CtrlN = full_cov_table_df.loc[condition,]
# CtrlLRRK2_CtrlN.to_csv(output_folder+Study+"_CtrlLRRK2_CtrlN_cov.tsv", sep="\t", index=False)
# CtrlLRRK2_CtrlN_samples = list(CtrlLRRK2_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CtrlLRRK2_CtrlN_samples]
# print("CtrlLRRK2_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlLRRK2_CtrlN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CtrlLRRK2_CtrlN_samples]
# print("CtrlLRRK2_CtrlN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlLRRK2_CtrlN_tpm.tsv", sep="\t")
# print(CtrlLRRK2_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ##====================================================================
# ## CtrlLRRK2_CaseN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Control") &(full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==1) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CtrlLRRK2_CaseN = full_cov_table_df.loc[condition,]
# CtrlLRRK2_CaseN.to_csv(output_folder+Study+"_CtrlLRRK2_CaseN_cov.tsv", sep="\t", index=False)
# CtrlLRRK2_CaseN_samples = list(CtrlLRRK2_CaseN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CtrlLRRK2_CaseN_samples]
# print("CtrlLRRK2_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlLRRK2_CaseN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CtrlLRRK2_CaseN_samples]
# print("CtrlLRRK2_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlLRRK2_CaseN_tpm.tsv", sep="\t")
# print(CtrlLRRK2_CaselN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CtrlBGA_CaseN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Control") &(full_cov_table_df["GBA"]==1) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CtrlLGBA_CtrlN = full_cov_table_df.loc[condition,]
# CtrlLGBA_CtrlN.to_csv(output_folder+Study+"_CtrlGBA_CaseN_cov.tsv", sep="\t", index=False)
# CtrlLGBA_CtrlN_samples = list(CtrlLGBA_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CtrlLGBA_CtrlN_samples]
# print("CtrlGBA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlGBA_CaseN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CtrlLGBA_CtrlN_samples]
# print("CtrlGBA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlGBA_CaseN_tpm.tsv", sep="\t")
# print(CtrlLGBA_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ## CtrlSNCA_CaseN
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Control") &(full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==1) & (full_cov_table_df["LRRK2"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Case") & (full_cov_table_df["Mutation"]=="NonCarrier"))
# CtrlSNCA_CtrlN = full_cov_table_df.loc[condition,]
# CtrlSNCA_CtrlN.to_csv(output_folder+Study+"_CtrlSNCA_CaseN_cov.tsv", sep="\t", index=False)
# CtrlSNCA_CtrlN_samples = list(CtrlSNCA_CtrlN['sample_id'])
#
# RNAseq = gene_read_df.loc[:, CtrlSNCA_CtrlN_samples]
# print("CtrlSNCA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlSNCA_CaseN_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, CtrlSNCA_CtrlN_samples]
# print("CtrlGBA_CaseN:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_CtrlSNCA_CaseN_tpm.tsv", sep="\t")
# print(CtrlSNCA_CtrlN.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
# ###########################################
#
# ##====================================================================
# ## LRRK2_Case_Ctrl
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==1) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Case") &   (full_cov_table_df["GBA"]==0) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==1) )
# LRRK2_Case_Ctrl = full_cov_table_df.loc[condition,]
# LRRK2_Case_Ctrl.to_csv(output_folder+Study+"_LRRK2_Case_Ctrl_cov.tsv", sep="\t", index=False)
# LRRK2_Case_Ctrl_samples = list(LRRK2_Case_Ctrl['sample_id'])
#
# RNAseq = gene_read_df.loc[:, LRRK2_Case_Ctrl_samples]
# print("LRRK2_Case_Ctrl:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_LRRK2_Case_Ctrl_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, LRRK2_Case_Ctrl_samples]
# print("LRRK2_Case_Ctrl:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_LRRK2_Case_Ctrl_tpm.tsv", sep="\t")
# print(LRRK2_Case_Ctrl.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
#
# ##====================================================================
# ## GBA_Case_Ctrl
# print("====================================")
# condition = ((full_cov_table_df["case_control_other_latest"]=="Control") & (full_cov_table_df["GBA"]==1) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==0) ) | \
#             ((full_cov_table_df["case_control_other_latest"]=="Case") &   (full_cov_table_df["GBA"]==1) & (full_cov_table_df["SNCA"]==0) & (full_cov_table_df["LRRK2"]==0) )
# GBA_Case_Ctrl = full_cov_table_df.loc[condition,]
# GBA_Case_Ctrl.to_csv(output_folder+Study+"_GBA_Case_Ctrl_cov.tsv", sep="\t", index=False)
# GBA_Case_Ctrl_samples = list(GBA_Case_Ctrl['sample_id'])
#
# RNAseq = gene_read_df.loc[:, GBA_Case_Ctrl_samples]
# print("GBA_Case_Ctrl:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_GBA_Case_Ctrl_read.tsv", sep="\t")
#
# RNAseq = gene_tpm_df.loc[:, GBA_Case_Ctrl_samples]
# print("GBA_Case_Ctrl:",RNAseq.shape)
# RNAseq.to_csv(output_folder+Study+"_GBA_Case_Ctrl_tpm.tsv", sep="\t")
# print(GBA_Case_Ctrl.loc[:,["case_control_other_latest","Mutation"]].value_counts(sort=False))
#
