import matplotlib.pyplot as plt
import pandas as pd

Study="PPMI"
output_folder = "../" + Study

## One sample may have multiple diagnosis age, the AMP-PD admin suggests using the earliest one.
fout = open(output_folder+"/run_inout/participant_info_unique_id_disease_age.tsv", "w")
with open("../../data/participant_info_full_results_04212021_unique_id.tsv", "r") as f:
    header = next(f)  ## header line
    fout.write(header[:-1] + '\tdisease_year\tDoD\tmissing_flag\n')
    for line in f:
        line_ls = line[:-1].split("\t")
        age_at_baseline = line_ls[3]

        ## Format of age_at_diagnosis: 20|30|40
        age_at_diagnosis = line_ls[9]

        case_latest = line_ls[11]

        ## check the disease age, if it is 'NA', the final earliest diagnosis age will be set to 'NA'
        if age_at_diagnosis == "NA":
            earliest = "NA"
        else:
            print(age_at_diagnosis)
            ## check the disease age, if it is not 'NA', the final earliest diagnosis age will be set to the earlist Non-zero value
            age_ls = sorted(list(set([int(i) for i in age_at_diagnosis.split("|")])))
            earliest = age_ls[0]
            if earliest == 0 and len(age_ls) > 1:
                earliest = age_ls[1]

        if earliest != "NA" and age_at_baseline != "NA":
            ## if both basline age and disease age are avalibale, Duration of Disease(DoD) = baseline age - disease age
            DoD = int(age_at_baseline) - earliest
            missing_flag = DoD
        else:
            ## if neither basline age nor disease age is avalibale, DoD = 0
            if case_latest.strip('"') == "Case":
                print(line_ls[0])
                DoD = 0
                if earliest == "NA": missing_flag = "Diag_Age_NA"
                if age_at_baseline == "NA": missing_flag = "Base_Age_NA"
            else:
                ## Controls were set to 0
                DoD = 0
                missing_flag = "Control"

        fout.write(line[:-1] + "\t" + str(earliest) + "\t" + str(DoD) + "\t" + str(missing_flag) + "\n")
fout.close()

## load cohort information
cohort = pd.read_csv("../../data/Clinical/amp-pd-participants.csv", index_col=None)
print("Cohort shape: ", cohort.shape)

participants = pd.read_csv(output_folder+"/run_inout/participant_info_unique_id_disease_age.tsv", index_col=None, header=0, sep="\t")
print("Participants shape: ", participants.shape)

participants_with_cohort = pd.merge(participants, cohort[['participant_id','study']], how="left", on="participant_id")
participants_with_cohort.to_csv(output_folder+"/run_inout/participant_info_unique_id_disease_age.tsv", index=False, sep="\t")
print(participants_with_cohort.info())

columns = ['participant_id', 'disease_year', 'DoD', "missing_flag"]
participants = participants.loc[:, columns]

col_ls = ["sample_id","participant_id","Study","visit_month","case_control_other_latest",
          "upsit_total_score", "moca_total_score","Plate","RIN","sex","age_at_baseline","race",
          "PRS","GBA","LRRK2","SNCA","Mutation","DoD","study_arm","pd_medication_start_months_after_baseline","treat"]
# 0. CaseA vs ControlA
case_control_df = pd.read_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlA.tsv", index_col=None, header=0, sep="\t")
print("CA_CA:",case_control_df.shape)
full_cov_table_df = pd.merge(case_control_df, participants, how="left", on="participant_id", suffixes=("", "_y"))
full_cov_table_df.loc[(full_cov_table_df['case_control_other_latest'] == 'Control'), ["DoD"]] = 0
full_cov_table_df= full_cov_table_df.loc[~(full_cov_table_df['missing_flag'].isin(['Diag_Age_NA','Base_Age_NA'])),:]
full_cov_table_df = full_cov_table_df.loc[:,col_ls]
print(full_cov_table_df.info())
#####################################
## check
# Case_DoD = full_cov_table.loc[(full_cov_table['case_control_other_latest'] == 'Case'), ["missing_flag"]]
# Case_DoD.value_counts()
# fig = plt.figure()
# ax = Case_DoD["missing_flag"].value_counts().sort_index().plot(kind='bar')
# ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
# plt.xlabel("Duration of Disease (year)")
# plt.title("Distribution of DoD PPMI Case baseline samples")
# plt.close()

if Study=="PDBP" or Study=="PDBF" or Study=="PPB" :
    # 2 samples with very big DoD, By checking the reacord in BigQuery tables, there are multiple diagnosis ages, and I made the following changes
    full_cov_table_df.loc[full_cov_table_df['sample_id'] == "PD-PDDR399XYQ-BLM0T1", ['DoD']] = 6
    full_cov_table_df.loc[full_cov_table_df['sample_id'] == "PD-PDDR399XYQ-BLM0T1", ['missing_flag']] = 6

    full_cov_table_df.loc[full_cov_table_df['sample_id'] == "PD-PDFU000KK6-BLM0T1", ['DoD']] = 2
    full_cov_table_df.loc[full_cov_table_df['sample_id'] == "PD-PDFU000KK6-BLM0T1", ['missing_flag']] = 2

    condition = full_cov_table_df['DoD'] >= 0
    full_cov_table_df = full_cov_table_df.loc[condition, :]

# print(full_cov_table_caseA_ctrlN_df.shape)
full_cov_table_df.to_csv(output_folder+"/run_inout/"+Study+"_CaseA_CtrlA_cov.tsv", sep="\t", index=False)

print("=====================================")
print("CA_CA_cov:",full_cov_table_df.shape)
