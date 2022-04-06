import pandas as pd

## Load mutation records
mutation_df = pd.read_csv("../../data/participant_mutations.csv", index_col=0, header=0)
print(mutation_df.shape)

## genetics carriers
mutation_yes_df = mutation_df.loc[mutation_df['has_known_PD_mutation_in_WGS'] == 'Yes', :]
mutation_no_df = mutation_df.loc[mutation_df['has_known_PD_mutation_in_WGS'] == 'No', :]
mutation_NA_df = mutation_df.loc[mutation_df['has_known_PD_mutation_in_WGS'].isna(), :]

mutation_SNCA_df = mutation_df.loc[mutation_df['has_known_SNCA_mutation_in_WGS'] == 'Yes', :]
mutation_GBA_df = mutation_df.loc[mutation_df['has_known_GBA_mutation_in_WGS'] == 'Yes', :]
mutation_LRRK2_df = mutation_df.loc[mutation_df['has_known_LRRK2_mutation_in_WGS'] == 'Yes', :]

mutation_SNCA_ls = mutation_SNCA_df.index.tolist()
mutation_GBA_ls = mutation_GBA_df.index.tolist()
mutation_LRRK2_ls = mutation_LRRK2_df.index.tolist()

x = set([i for i in mutation_SNCA_df.index] + [i for i in mutation_GBA_df.index] +[i for i in mutation_LRRK2_df.index])

print(mutation_df.shape)

## load cohort information
cohort = pd.read_csv("../../data/Clinical/amp-pd-participants.csv", index_col=None)
print("Cohort shape: ", cohort.shape)

## enrollment
enrollment_df = pd.read_csv("../../data/Clinical/clinical/Enrollment.csv",sep=",", header=0, index_col=0)

## PD_Medication
medication_df = pd.read_csv("../../data/Clinical2.5/2021_v2-5release_0510/clinical/PD_Medical_History.csv",sep=",", header=0, index_col=0)
medication_df = medication_df.loc[medication_df["visit_name"]=="LOG",["use_of_pd_medication","pd_medication_start_months_after_baseline"]]

## PRS
PRS = pd.read_csv("../../data/variants/PPB_prs_CaseCtrl.best", index_col=None, header=0,sep=" ")

RNA_samples25 = pd.read_csv("../../data/sample_info_with_plate_rin_study_06082021_v25.csv", index_col=None, header=0)
RNA_samples20 = pd.read_csv("../../data/sample_info_with_plate_rin_study_06022021.csv", index_col=None, header=0)
RNA_samples = RNA_samples25
print("Samples shape: ", RNA_samples.shape)
# print(RNA_samples.info()) ## 8461 samples

participants = pd.read_csv("../../data/participant_info_full_results_04212021_unique_id.tsv", index_col=None, header=0,
                           sep="\t")

participants = pd.merge(participants,PRS[["IID","PRS"]], how="left", left_on="participant_id", right_on="IID")
participants.drop(['IID'], axis="columns", inplace=True)

participants = pd.merge(participants,medication_df, how="left", left_on="participant_id", right_index=True)
participants.loc[:,"treat"] = -1
participants.loc[participants["pd_medication_start_months_after_baseline"] <= 0,"treat"] = 1
participants.loc[participants["pd_medication_start_months_after_baseline"] > 0,"treat"] = 0

# participants_df.to_csv("../../data/PPB_participants_covariates_08092021.xls", sep="\t")

participants_with_study_arm = pd.merge(participants,enrollment_df["study_arm"],left_on="participant_id",how="left",right_index=True)
print("Participant shape: ", participants_with_study_arm.shape)

#######################################################################################################
## Merge samples and participant, columns with the same name were renamed with a suffix "_y"
samples_info = pd.merge(RNA_samples, participants_with_study_arm, how="left", on="participant_id", suffixes=("_x", ""))
samples_info.index=samples_info["sample_id"]
samples_info.drop(['sample_id',"BAM","age_at_diagnosis_x", "case_control_other_latest_x"], axis="columns", inplace=True)
# samples_info.to_csv("../../data/PPB_RNAseq_samples_covariance_08092021.xls", sep="\t")
print("samples_info shape: ", samples_info.shape)


#######################
Study="PPMI"
output_folder = "../" + Study
# SELECT participant_id FROM `amp-pd-research.2020_v2release_1218.amp_pd_participants` WHERE study="PPMI"
## Retrun: 1923
if Study == "PDBF":
    Study_cohort = cohort.loc[cohort['study'].isin(["PDBP","BioFIND"]), :]
elif Study == "PDBP":
    Study_cohort = cohort.loc[cohort['study'].isin(["PDBP"]), :]
elif Study == "PPB":
    Study_cohort = cohort.loc[cohort['study'].isin(["PDBP","BioFIND","PPMI"]), :]
else:
    Study_cohort = cohort.loc[cohort['study'].isin(["PPMI"]), :]

print(f'{Study} size: {Study_cohort.shape}')

Study_participant_ls = Study_cohort['participant_id']
print(len(Study_participant_ls))  ## 1923 participants

## Has RNA-seq data
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
participants_without_0month_id_ls = [i for i in set(Study_RNA['participant_id'].values) if
                                     i not in set(Study_RNA_baseline['participant_id'].values)]
participants_without_0month = Study_RNA.loc[
    Study_RNA['participant_id'].isin(participants_without_0month_id_ls), ['participant_id', 'visit_month']]

participants_without_0month_earliest_visit = participants_without_0month.groupby(['participant_id']).min()
print(f'Without M0，esarlist visit:\n{participants_without_0month_earliest_visit["visit_month"].value_counts()}')

participants_6month_visit = participants_without_0month.loc[
                            participants_without_0month['visit_month'] == 6.0, :]

participants_6month_visit_full = pd.merge(Study_RNA, participants_6month_visit, how='inner',
                                          left_index=True,right_index=True, suffixes=("", "_y"))
participants_6month_visit_full.drop(list(participants_6month_visit_full.filter(regex='_y$')), axis=1, inplace=True)

Study_RNA_baseline = Study_RNA_baseline.append(participants_6month_visit_full, ignore_index=False)
print(f'Baseline data shape: {Study_RNA_baseline.shape}')## 1568

## RIN > 5
condition = Study_RNA_baseline["RIN"] > 5.0
Study_RNA_baseline_RIN = Study_RNA_baseline.loc[condition, :]
print(f'Study_RNA_RIN: {len(set(Study_RNA_baseline_RIN["participant_id"].values))}') ## 1461

## Get samples with Plate Or RIN  are NA
Plate_na_ls = list(Study_RNA_baseline_RIN.loc[Study_RNA_baseline_RIN['Plate'].isna(),].index)
RIN_na_ls = list(Study_RNA_baseline_RIN.loc[Study_RNA_baseline_RIN['RIN'].isna(),].index)
Plate_RIN_na_ls = list(
    Study_RNA_baseline_RIN.loc[(Study_RNA_baseline_RIN['RIN'].isna()) | (Study_RNA_baseline_RIN['Plate'].isna()),].index)

condition = (Study_RNA_baseline_RIN.index.isin(Plate_RIN_na_ls))
Study_RNA_baseline_RIN = Study_RNA_baseline_RIN.loc[~condition, :]
print(f'Study_RNA_RIN Baseline data shape: {Study_RNA_baseline_RIN.shape}')

## Keep white race only
Study_RNA_baseline_RIN_white = Study_RNA_baseline_RIN.loc[Study_RNA_baseline_RIN['race'] == "White", :]
print(f'Study_RNA_RIN_white: {len(set(Study_RNA_baseline_RIN_white["participant_id"].values))}')

## Study Cases/Controls
condition = (Study_RNA_baseline_RIN_white['case_control_other_latest'].isin(["Case", "Control"]))
Study_RNA_baseline_PDvsNO = Study_RNA_baseline_RIN_white.loc[condition, :]
print(f'Case_Control: {Study_RNA_baseline_PDvsNO.shape}')

# if Study == "skip":
#     condition = Study_RNA_baseline_PDvsNO['study_arm'].isin(["PD","Healthy Control"])
#     Study_RNA_baseline_PDvsNO = Study_RNA_baseline_PDvsNO.loc[~condition, :]
#     print(len(set(Study_RNA_baseline_PDvsNO['participant_id'].values)))  ## 1582 participants

Study_RNA_baseline_PDvsNO.loc[:,['GBA']] = "0"
Study_RNA_baseline_PDvsNO.loc[:,['LRRK2']] = "0"
Study_RNA_baseline_PDvsNO.loc[:,['SNCA']] = "0"
Study_RNA_baseline_PDvsNO.loc[
    Study_RNA_baseline_PDvsNO["participant_id"].isin(mutation_GBA_df.index), ['GBA']] = "1"

Study_RNA_baseline_PDvsNO.loc[
    Study_RNA_baseline_PDvsNO["participant_id"].isin(mutation_LRRK2_df.index),['LRRK2']] = "1"

Study_RNA_baseline_PDvsNO.loc[
    Study_RNA_baseline_PDvsNO["participant_id"].isin(mutation_SNCA_df.index),['SNCA']] = "1"


###############################################################
Study_RNA_baseline_Case = Study_RNA_baseline_PDvsNO.loc[
                         Study_RNA_baseline_PDvsNO['case_control_other_latest'] == "Case", :]
Study_RNA_baseline_Case = Study_RNA_baseline_Case.loc[
                          Study_RNA_baseline_Case["study_arm"].isin(["PD","Genetic Cohort PD","Genetic Registry PD"]),:]

Study_RNA_baseline_Ctrl = Study_RNA_baseline_PDvsNO.loc[
                            Study_RNA_baseline_PDvsNO['case_control_other_latest'] == "Control", :]
Study_RNA_baseline_Ctrl = Study_RNA_baseline_Ctrl.loc[
                          Study_RNA_baseline_Ctrl["study_arm"].isin(["Healthy Control","Genetic Cohort Unaffected","Genetic Registry Unaffected"]),:]


## keep case samples with no changes between basline and latest
condition = Study_RNA_baseline_Case["case_control_other_at_baseline"] == "Case"
Study_RNA_baseline_Case_nochanges = Study_RNA_baseline_Case.loc[condition, :]

Study_RNA_baseline_Case_nochanges.loc[
    Study_RNA_baseline_Case_nochanges["participant_id"].isin(mutation_yes_df.index),
    ['Mutation']] = "Carrier"

Study_RNA_baseline_Case_nochanges.loc[
    Study_RNA_baseline_Case_nochanges["participant_id"].isin(mutation_no_df.index),
    ['Mutation']] = "NonCarrier"

Study_RNA_baseline_Case_nochanges.loc[
    Study_RNA_baseline_Case_nochanges["participant_id"].isin(mutation_NA_df.index),
    ['Mutation']] = "NoValue"

Study_RNA_baseline_Case_nochanges_Carrier = \
    Study_RNA_baseline_Case_nochanges.loc[
    Study_RNA_baseline_Case_nochanges['Mutation'] == "Carrier",:]

Study_RNA_baseline_Case_nochanges_NoCarrier = \
    Study_RNA_baseline_Case_nochanges.loc[
    Study_RNA_baseline_Case_nochanges['Mutation'] == "NonCarrier",:]

Study_RNA_baseline_Case_nochanges_NA = \
    Study_RNA_baseline_Case_nochanges.loc[
    Study_RNA_baseline_Case_nochanges['Mutation'] == "NoValue",:]

## Check
# Study_RNA_baseline_Case["case_control_other_at_baseline"].value_counts()
# Study_RNA_baseline_Case["diagnosis_at_baseline"].value_counts()
# Study_RNA_baseline_Case["diagnosis_latest"].value_counts()

######################
##Control
condition = Study_RNA_baseline_Ctrl["case_control_other_at_baseline"] == "Control"
Study_RNA_baseline_Ctrl_nochanges = Study_RNA_baseline_Ctrl.loc[condition, :]

Study_RNA_baseline_Ctrl_nochanges.loc[
                   Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(mutation_yes_df.index),
                   ['Mutation']] = "Carrier"
Study_RNA_baseline_Ctrl_nochanges.loc[
                   Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(mutation_no_df.index),
                   ['Mutation']] = "NonCarrier"

Study_RNA_baseline_Ctrl_nochanges.loc[
                   Study_RNA_baseline_Ctrl_nochanges["participant_id"].isin(mutation_NA_df.index),
                   ['Mutation']] = "NoValue"


Study_RNA_baseline_Ctrl_nochanges_Carrier = \
    Study_RNA_baseline_Ctrl_nochanges.loc[
    Study_RNA_baseline_Ctrl_nochanges['Mutation'] == "Carrier",:]

Study_RNA_baseline_Ctrl_nochanges_NoCarrier = \
    Study_RNA_baseline_Ctrl_nochanges.loc[
    Study_RNA_baseline_Ctrl_nochanges['Mutation'] == "NonCarrier",:]

Study_RNA_baseline_Ctrl_nochanges_NA = \
    Study_RNA_baseline_Ctrl_nochanges.loc[
    Study_RNA_baseline_Ctrl_nochanges['Mutation'] == "NoValue",:]

################################
# 0. Combine the final cases, controls, save the result into a file.
CA_CA_pd = pd.concat([Study_RNA_baseline_Case_nochanges,Study_RNA_baseline_Ctrl_nochanges])
CA_CA_pd.to_csv(output_folder +  "/run_inout/" +Study+"_CaseA_CtrlA.tsv", sep="\t", index=True)
