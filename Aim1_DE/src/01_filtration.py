import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3_unweighted

## Load mutation records
mutation_df = pd.read_csv("../../data/participant_mutations.csv", index_col=0, header=0)
print(mutation_df.shape)

## Only keep genetics carriers
mutation_df = mutation_df.loc[mutation_df['has_known_PD_mutation_in_WGS'] == 'Yes', :]
print(mutation_df.shape)

## load cohort information
cohort = pd.read_csv("../../data/Clinical/amp-pd-participants.csv", index_col=None)
print("Cohort shape: ", cohort.shape)
print(cohort.info())

## Load all the sample information, the sample table was generated from the SQL
########################################
# # V2.5
# SELECT
# t0.participant_id AS participant_id,
# t0.visit_month AS visit_month,
# t0.sample_id AS sample_id,
# t0.BAM AS BAM,
#
# t2.diagnosis as diagnosis,
# t2.initial_diagnosis as initial_diagnosis,
# t2.most_recent_diagnosis as most_recent_diagnosis,
# t2.change_in_diagnosis as change_in_diagnosis,
# t2.age_at_diagnosis as age_at_diagnosis,
# t2.pd_diagnosis_months_after_baseline as pd_diagnosis_after_baseline,
# t2.pd_medication_initiation_months_after_baseline as pd_medication_initiation_after_baseline,
# t2.pd_medication_start_months_after_baseline as pd_medication_start_date_after_baseline,
# t2.use_of_pd_medication as use_of_pd_medication,
# t2.pd_medication_recent_use_months_after_baseline  as date_of_pd_medication_recent_use_after_baseline,
# t2.on_levodopa as on_levodopa,
# t2.on_dopamine_agonist as on_dopamine_agonist,
# t2.on_other_pd_medications as on_other_pd_medications,
#
# t3.updrs1_ment_behav_mood_score AS updrs1_ment_behav_mood_score,
# t3.updrs2_adl_score AS updrs2_adl_score,
# t3.updrs3_motor_examination_score AS updrs3_motor_examination_score,
# t3.updrs4_therapy_complications_score AS updrs4_therapy_complications_score,
# t4.upsit_total_score AS upsit_total_score,
# t5.moca_total_score AS moca_total_score,
# t6.case_control_other_latest AS case_control_other_latest,
# t6.case_control_other_at_baseline AS case_control_other_at_baseline,
# t6.diagnosis_latest AS diagnosis_latest,
# t6.diagnosis_at_baseline AS diagnosis_at_baseline,
#
# t7.Plate as Plate,
# t7.RIN_Value as RIN,
# t7.Position as Position,
#
# t8.study as Study
#
# FROM
# (
#     SELECT participant_id, sample_id, visit_month, BAM
#     FROM `amp-pd-research.2021_v2_5release_0510_transcriptomics.rna_seq_samples`
# ) t0
# LEFT JOIN
# (
#     SELECT participant_id, visit_month, diagnosis,
#     initial_diagnosis, most_recent_diagnosis, change_in_diagnosis,
#     pd_diagnosis_months_after_baseline, pd_medication_initiation_months_after_baseline, age_at_diagnosis,
#     pd_medication_start_months_after_baseline, use_of_pd_medication, pd_medication_recent_use_months_after_baseline,
#     on_levodopa, on_dopamine_agonist, on_other_pd_medications
#     FROM `amp-pd-research.2021_v2_5release_0510.PD_Medical_History`
# ) t2 USING(participant_id, visit_month)
# LEFT JOIN
# (
#     SELECT
#     DISTINCT participant_id, visit_month, updrs1_ment_behav_mood_score, updrs2_adl_score, updrs3_motor_examination_score, updrs4_therapy_complications_score
#     FROM `amp-pd-research.2021_v2_5release_0510.UPDRS`
# )t3 USING(participant_id, visit_month)
# LEFT JOIN
# (
#     SELECT DISTINCT participant_id, visit_month, upsit_total_score
#     FROM `amp-pd-research.2021_v2_5release_0510.UPSIT`
# )t4 USING(participant_id, visit_month)
# LEFT JOIN
# (
#     SELECT DISTINCT participant_id, visit_month, moca_total_score
#     FROM `amp-pd-research.2021_v2_5release_0510.MOCA`
# )t5 USING(participant_id, visit_month)
# LEFT JOIN
# (
#     SELECT participant_id, case_control_other_latest, case_control_other_at_baseline,diagnosis_latest, diagnosis_at_baseline
#     FROM `amp-pd-research.2021_v2_5release_0510.amp_pd_case_control`
# ) t6 USING(participant_id)
# left join
# (
#     select sample_id, RIN_Value, Plate, Position
#     from `amp-pd-research.2021_v2_5release_0510_transcriptomics.rna_quality_metrics`
# ) t7 USING(sample_id)
# left join
# (
#     select participant_id, study
#     from `amp-pd-research.2021_v2_5release_0510.amp_pd_participants`
# ) t8 USING(participant_id)

RNA_samples25 = pd.read_csv("../../data/sample_info_with_plate_rin_study_06082021_v25.csv", index_col=None, header=0)
RNA_samples20 = pd.read_csv("../../data/sample_info_with_plate_rin_study_06022021.csv", index_col=None, header=0)
RNA_samples = RNA_samples25
print("Samples shape: ", RNA_samples.shape)
print(RNA_samples.info()) ## 8461 samples

## Test if there are samples with date_of_pd_medication_initiation earlier than date_of_pd_diagnosis
# SELECT date_of_pd_diagnosis,date_of_pd_medication_initiation
# FROM `amp-pd-research.2020_v2release_1218.PD_Medical_History`
# where (date_of_pd_diagnosis is not null) and (date_of_pd_medication_initiation is not null) and
# (date_of_pd_medication_initiation < date_of_pd_diagnosis)
## Return: 0

## Load participant information
## the table was first generated using the SQL code below,
## then processed using the python script: process_participnt.py
# SELECT
#     DISTINCT t0.participant_id AS participant_id,
#     t1.sex AS sex,
#     t1.age_at_baseline AS age_at_baseline,
#     t1.race AS race,
#     t1.ethnicity AS ethnicity,
#     t2.age_at_diagnosis AS age_at_diagnosis,
# FROM
# (
#     SELECT DISTINCT participant_id
#     FROM `amp-pd-research.2020_v2release_1218_transcriptomics.rna_seq_samples`
# ) t0 LEFT JOIN
# (
#     SELECT DISTINCT participant_id, sex, age_at_baseline, race, ethnicity
#     FROM `amp-pd-research.2020_v2release_1218.Demographics`
# )t1 USING(participant_id) LEFT JOIN
# (
#     SELECT DISTINCT participant_id, date_of_pd_diagnosis, age_at_diagnosis
# FROM `amp-pd-research.2020_v2release_1218.PD_Medical_History`
# ) t2 USING(participant_id)

participants = pd.read_csv("../../data/participant_info_full_results_04212021_unique_id.tsv", index_col=None, header=0,
                           sep="\t")
print("Participant shape: ", participants.shape)
print(participants.info())

PRS = pd.read_csv("../../data/variants/PPB_prs_all.best", index_col=None, header=0,sep=" ")

participants = pd.merge(participants,PRS[["IID","PRS"]], how="left", left_on="participant_id", right_on="IID")
participants.drop(['IID'], axis="columns", inplace=True)
# participants.to_csv("../../data/PPB_participants_covariance_08092021.xls", sep="\t",index=False)

## Merge samples and participant, columns with the same name were renamed with a suffix "_y"
samples_info = pd.merge(RNA_samples, participants, how="left", on="participant_id", suffixes=("_x", ""))
print(samples_info.info())
samples_info.index=samples_info["sample_id"]
samples_info.drop(['sample_id',"BAM","age_at_diagnosis_x", "case_control_other_latest_x"], axis="columns", inplace=True)
# samples_info.to_csv("../../data/PPB_RNAseq_samples_covariance_08092021.xls", sep="\t")


# ## all genes
# gene_read_df = pd.read_csv("../data/RNAseq/gene_expr_matrix_reads_row_genes.txt", index_col=0, header=0, sep="\t")
# print(gene_read_df.shape)
#
# gene_tpm_df = pd.read_csv("../data/RNAseq/gene_expr_matrix_tpm_row_genes.txt", index_col=0, header=0, sep="\t")
# print(gene_tpm_df.shape)

#######################
Study="PDBF"
output_folder = "../" + Study
# SELECT participant_id FROM `amp-pd-research.2020_v2release_1218.amp_pd_participants` WHERE study="PPMI"
## Retrun: 1923
if Study == "PDBF":
    Study_cohort = cohort.loc[cohort['study'].isin(["PDBP","BioFIND"]), :]
elif Study == "PDBP":
    Study_cohort = cohort.loc[cohort['study'].isin(["PDBP"]), :]
else:
    Study_cohort = cohort.loc[cohort['study'].isin(["PPMI"]), :]

print(f'{Study} size: {Study_cohort.shape}')

Study_participant_ls = Study_cohort['participant_id']
print(len(Study_participant_ls))  ## 1923 participants

condition = (samples_info["participant_id"].isin(Study_participant_ls))
Study_RNA = samples_info.loc[condition, :]
print(f'{Study} RNAseq size: {Study_RNA.shape}')
print(len(set(Study_RNA['participant_id'].values)))  ## 1582 participants

## Study Baseline samples
condition = (Study_RNA["visit_month"].isin([0.0, 0.5]))
Study_RNA_baseline = Study_RNA.loc[condition, :]
print(f'Visit at M0: {Study_RNA_baseline.shape}')
# Study_RNA_baseline.to_csv(Study+"_RNA_baseline_info.csv")

## For the participants having no M0 RNAseq, using the RNAseq data at the earliest available time point.
participants_id_without_0month_ls = [i for i in set(Study_RNA['participant_id'].values) if
                                     i not in set(Study_RNA_baseline['participant_id'].values)]
participants_without_0month = Study_RNA.loc[
    Study_RNA['participant_id'].isin(participants_id_without_0month_ls), ['participant_id', 'visit_month']]

participants_without_0month_earliest_visit = participants_without_0month.groupby(['participant_id']).min()
print(f'Without M0，esarlist visit:\n{participants_without_0month_earliest_visit["visit_month"].value_counts()}')

participants_6month_visit = participants_without_0month.loc[
                            participants_without_0month['visit_month'] == 6.0, :]

participants_6month_visit_full = pd.merge(Study_RNA, participants_6month_visit, how='inner',
                                          left_index=True,right_index=True, suffixes=("", "_y"))
participants_6month_visit_full.drop(list(participants_6month_visit_full.filter(regex='_y$')), axis=1, inplace=True)

Study_RNA_baseline = Study_RNA_baseline.append(participants_6month_visit_full, ignore_index=False)
print(f'Baseline data shape: {Study_RNA_baseline.shape}')

## Get samples with Plate Or RIN  are NA
Plate_na_ls = list(Study_RNA_baseline.loc[Study_RNA_baseline['Plate'].isna(),].index)
RIN_na_ls = list(Study_RNA_baseline.loc[Study_RNA_baseline['RIN'].isna(),].index)
Plate_RIN_na_ls = list(
    Study_RNA_baseline.loc[(Study_RNA_baseline['RIN'].isna()) | (Study_RNA_baseline['Plate'].isna()),].index)

condition = (Study_RNA_baseline.index.isin(Plate_RIN_na_ls))
Study_RNA_baseline = Study_RNA_baseline.loc[~condition, :]
print(f'Final Baseline data shape: {Study_RNA_baseline.shape}')

## Study Cases/Controls
condition = (Study_RNA_baseline['case_control_other_latest'].isin(["Case", "Control"]))
Study_RNA_baseline_PDvsNO = Study_RNA_baseline.loc[condition, :]
print(f'Case_Control: {Study_RNA_baseline_PDvsNO.shape}')

print(Study_RNA_baseline_PDvsNO["race"].value_counts())

# fig = plt.figure()
# ax = Study_RNA_baseline_PDvsNO["race"].hist()
# # plt.draw()
# ax.set_xticklabels(ax.get_xticklabels(), rotation=-45, ha='left')
# plt.xlabel("Race")
# plt.title(f"Distribution of Race in {Study} Case/Control samples")
# plt.close()

## Keep white race only
Study_RNA_baseline_PDvsNO_white = Study_RNA_baseline_PDvsNO.loc[Study_RNA_baseline_PDvsNO['race'] == "White", :]
print(f'RNA_baseline_PDvsNO_white: {Study_RNA_baseline_PDvsNO_white.shape}')

## RIN > 5
condition = Study_RNA_baseline_PDvsNO_white["RIN"] > 5.0
Study_RNA_baseline_PDvsNO_white_RIN = Study_RNA_baseline_PDvsNO_white.loc[condition, :]
print(f'RNA_baseline_PDvsNO_white_RIN: {Study_RNA_baseline_PDvsNO_white_RIN.shape}')
###############################
Study_RNA_baseline_Case = Study_RNA_baseline_PDvsNO_white_RIN.loc[
                         Study_RNA_baseline_PDvsNO_white_RIN['case_control_other_latest'] == "Case", :]
Study_RNA_baseline_Ctrl = Study_RNA_baseline_PDvsNO_white_RIN.loc[
                            Study_RNA_baseline_PDvsNO_white_RIN['case_control_other_latest'] == "Control", :]

## keep case samples with on changes batween basline and latest
condition = Study_RNA_baseline_Case["case_control_other_at_baseline"] == "Case"
Study_RNA_baseline_Case_nochanges = Study_RNA_baseline_Case.loc[condition, :]

Study_RNA_baseline_Case_nochanges_white_final = Study_RNA_baseline_Case_nochanges
Study_RNA_baseline_Case_nochanges_white_final.loc[
    Study_RNA_baseline_Case_nochanges_white_final["participant_id"].isin(mutation_df.index),
    ['Mutation']] = "Carrier"

Study_RNA_baseline_Case_nochanges_white_final.loc[
    ~Study_RNA_baseline_Case_nochanges_white_final["participant_id"].isin(mutation_df.index),
    ['Mutation']] = "NonCarrier"

Study_RNA_baseline_Case_nochanges_white_Carrier = \
    Study_RNA_baseline_Case_nochanges_white_final.loc[
    Study_RNA_baseline_Case_nochanges_white_final['Mutation'] == "Carrier",:]

## Check
# Study_RNA_baseline_Case["case_control_other_at_baseline"].value_counts()
# Study_RNA_baseline_Case["diagnosis_at_baseline"].value_counts()
# Study_RNA_baseline_Case["diagnosis_latest"].value_counts()

######################
##Control
condition = Study_RNA_baseline_Ctrl["case_control_other_at_baseline"] == "Control"
Study_RNA_baseline_Ctrl_nochanges = Study_RNA_baseline_Ctrl.loc[condition, :]

## Check
# Study_RNA_baseline_Ctrl["case_control_other_at_baseline"].value_counts()
# Study_RNA_baseline_Ctrl["diagnosis_at_baseline"].value_counts()
# Study_RNA_baseline_Ctrl["diagnosis_latest"].value_counts()

Ctrl_with_mut = Study_RNA_baseline_Ctrl_nochanges.loc[
                   Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(mutation_df.index), :]
Ctrl_with_mut['Mutation'] = "Carrier"

family_df = pd.read_csv("../../data/Clinical/clinical/Family_History_PD.csv", index_col=0, header=0)
print("family_df shape: ", family_df.shape)
print(family_df.info())
# print(family_df.head(n=5))

family_PD_father = family_df.loc[family_df["biological_father_with_pd"] == "Yes", :]
family_PD_mother = family_df.loc[family_df["biological_mother_with_pd"] == "Yes", :]
family_PD_other = family_df.loc[family_df["other_relative_with_pd"] == "Yes", :]

Ctrl_with_father_mut = Study_RNA_baseline_Ctrl_nochanges.loc[
                          Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(family_PD_father.index), :]
Ctrl_with_mother_mut = Study_RNA_baseline_Ctrl_nochanges.loc[
                          Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(family_PD_mother.index), :]
Ctrl_with_other_mut = Study_RNA_baseline_Ctrl_nochanges.loc[
                         Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(family_PD_other.index), :]

## Venn diagram
# vd3 = venn3_unweighted([set(Ctrl_with_father_mut['participant_id']), set(Ctrl_with_mother_mut['participant_id']),
#                         set(Ctrl_with_mut['participant_id'])],
#                        set_labels=('Ctrl_Fa_PD', 'Ctrl_Mo_PD', 'Ctrl_Mut(Carrier)'))

# for text in vd3.set_labels:
#     text.set_fontsize(16);
#
# plt.title('Venn Diagram', fontsize=20, pad=30);
# plt.savefig("Venn.pdf")
# plt.close()

condition = (Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(Ctrl_with_mut['participant_id']))
Study_RNA_baseline_Ctrl_nochanges_nomut = Study_RNA_baseline_Ctrl_nochanges.loc[~condition, :]

Study_RNA_baseline_Ctrl_nochanges_nomut['Mutation'] = 'NonCarrier'
Study_RNA_baseline_Ctrl_nochanges_nomut_white_final = Study_RNA_baseline_Ctrl_nochanges_nomut

## Combine the final cases, controls-nonCarrier, save the result into a file.
CA_CN_pd = pd.concat(
    [Study_RNA_baseline_Case_nochanges_white_final, Study_RNA_baseline_Ctrl_nochanges_nomut_white_final])
CA_CN_pd.to_csv(output_folder +  "/run_inout/" + Study+"_CaseAll_CtrlNonCarrier.tsv", sep="\t", index=True)

## Combine the control carriers, non-carriers, save the result into a file.
ctrl_carriers_pd = pd.concat([Ctrl_with_mut, Study_RNA_baseline_Ctrl_nochanges_nomut_white_final])
ctrl_carriers_pd.to_csv(output_folder +  "/run_inout/" +Study+"_Ctrl_Carriers_Non.tsv", sep="\t", index=True)

## Combine the final cases, controls, save the result into a file.
CA_CA_pd = pd.concat([Study_RNA_baseline_Case_nochanges_white_final,Ctrl_with_mut, Study_RNA_baseline_Ctrl_nochanges_nomut_white_final])
CA_CA_pd.to_csv(output_folder +  "/run_inout/" +Study+"_Case_Ctrl_all.tsv", sep="\t", index=True)

## Combine the final cases-carrier, controls-nonCarrier, save the result into a file.
CC_CN_pd = pd.concat([Study_RNA_baseline_Case_nochanges_white_Carrier, Study_RNA_baseline_Ctrl_nochanges_nomut_white_final])
CC_CN_pd.to_csv(output_folder +  "/run_inout/" +Study+"_CC_CN.tsv", sep="\t", index=True)


##END
CA_CN_pd.loc[CA_CN_pd['case_control_other_latest']=='Control',['Mutation']].value_counts()
CA_CN_pd.loc[CA_CN_pd['case_control_other_latest']=='Case',['Mutation']].value_counts()
