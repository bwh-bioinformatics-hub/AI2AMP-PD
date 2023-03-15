######
# %%
import os
import sys

import pandas as pd
import numpy as np
# from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="ticks")

from sklearn.preprocessing import StandardScaler, MinMaxScaler, minmax_scale
from sklearn.model_selection import StratifiedKFold, train_test_split

from sklearn.linear_model import LogisticRegression, Lasso, LassoLars, lasso_path, lars_path, LinearRegression

from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score, average_precision_score
from sklearn.metrics import mean_squared_error,roc_curve, precision_recall_curve, confusion_matrix, balanced_accuracy_score

from matplotlib.backends.backend_pdf import PdfPages

from funcs import *

##############
# %%
## Global params
# ["SVM","LR","SGD","ABC","XGB","RF","KNN","MLP"]
# model = sys.argv[1]
scale_method = "logz"
feature_selection = "DT"
dataset = "only_PPMI"
feature_folder = 1692
###################
# %%
print("=====================")
output_dir = "../results/" + dataset + "/LASSO_pool_" + str(feature_folder) + "/" + scale_method + "/LASSO_LN"
print(output_dir)
# Create folder if the directory doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

################
# %%
## prepare dataset
print("*******************************************")
# train
train_data_file = "data/PPMIvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
train_data_df = pd.read_csv(train_data_file, sep="\t", index_col=0, header=0)

train_DEG = train_data_df.loc[:, train_data_df.columns.str.startswith("DEG")]
train_eRNA = train_data_df.loc[:, train_data_df.columns.str.startswith("eRNA")]
train_circRNA = train_data_df.loc[:, train_data_df.columns.str.startswith("circRNA")]
train_PRS = train_data_df.loc[:, train_data_df.columns.str.startswith("PRS")]
train_Clinical = train_data_df.loc[:, train_data_df.columns.str.startswith("Clin")]
DEGs = train_DEG.columns.tolist()
eRNAs = train_eRNA.columns.tolist()
circRNAs = train_circRNA.columns.tolist()
PRS_features = ["PRS"]
clinicals_features = train_Clinical.columns.tolist()  ## ["upsit_total_score","sex_int","age_at_baseline"]

# test
test_data_file = "data/PDBFvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
test_data_df = pd.read_csv(test_data_file, sep="\t", index_col=0, header=0)

####
best_feature_df = pd.read_csv(output_dir+"/best_features_LN.csv",index_col=0)
best_features = best_feature_df.index.tolist()

kept_features = best_features
print("Total feature number:", len(kept_features))

######################
# %%
train_X = train_data_df.loc[:, kept_features]
train_X = train_X.loc[~train_X.isna().any(axis=1), :]
train_y = train_data_df.loc[train_X.index.tolist(), ["Y"]]
train_i_X, val_i_X, train_i_y, val_i_y = train_test_split(train_X, train_y, test_size=0.2, shuffle=True, random_state=42)

train_index_ls = train_i_X.index.tolist()
val_index_ls = val_i_X.index.tolist()

print("Selected feature number:", len(kept_features))
train_i_X = train_i_X.loc[:, kept_features]
val_i_X = val_i_X.loc[:, kept_features]

## log transform
if scale_method.startswith("log"):
    print("log")
    train_i_X = np.log10(train_i_X + 1)
    val_i_X = np.log10(val_i_X + 1)

## MinMax
if scale_method.endswith("mm"):
    print("mm")
    train_i_X = minmax_scale(train_i_X, feature_range=(0, 1))
    val_i_X = minmax_scale(val_i_X, feature_range=(0, 1))
    # scaler = MinMaxScaler()
    # train_i_X = scaler.fit_transform(train_i_X)
    # val_i_X = scaler.transform(val_i_X)
    train_i_X = pd.DataFrame(data=train_i_X, index=train_index_ls, columns=kept_features)
    val_i_X = pd.DataFrame(data=val_i_X, index=val_index_ls, columns=kept_features)

## Standardization
if scale_method.endswith("z"):
    print("z")
    train_i_X = stats.zscore(train_i_X)
    val_i_X = stats.zscore(val_i_X)
    # scaler = StandardScaler()
    # train_i_X = scaler.fit_transform(train_i_X)
    # val_i_X = scaler.transform(val_i_X)

print("Final train size:", train_i_X.shape)
print("Final val size:", val_i_X.shape)

###################################################################
# %%
print("Test data...")
## Test Data

print("Selected feature number:", len(kept_features))
test_X = test_data_df.loc[:, kept_features]
test_X = test_X.loc[~test_X.isna().any(axis=1), :]
test_y = test_data_df.loc[test_X.index.tolist(), ["Y"]]
test_index_ls = test_X.index.tolist()
print("Final test size:", test_X.shape)

# if "Clin_sex_int" in kept_features:
#     test_X.drop(['Clin_sex_int'], axis=1, inplace=True)

## log transform
if scale_method.startswith("log"):
    print("log")
    test_X = np.log10(test_X + 1)

## Scale
if scale_method.endswith("mm"):
    print("mm")
    test_X = minmax_scale(test_X, feature_range=(0, 1))
    # test_X = scaler.transform(test_X)
    test_X = pd.DataFrame(data=test_X, index=test_index_ls, columns=kept_features)

if scale_method.endswith("z"):
    print("z")
    test_X = stats.zscore(test_X)
    # test_X = scaler.transform(test_X)

# if "Clin_sex_int" in kept_features:
#     # test_X=np.append(test_X, test_data_df.loc[test_index_ls,["Clin_sex_int"]].values, axis=1)
#     test_X.loc[:, "Clin_sex_int"] = test_data_df.loc[test_index_ls, "Clin_sex_int"].values

######
# %%
f_perform = open(output_dir+"/stepped_performances.txt","a")

train_i_X = train_i_X.loc[:, best_features]
val_i_X = val_i_X.loc[:, best_features]
test_X = test_X.loc[:, best_features]

clf = Lasso(alpha=0.028,max_iter=3000,random_state=12)
clf.fit(train_i_X.to_numpy(), train_i_y.to_numpy().ravel())
# train performance
train_pred_res = clf.predict(train_i_X.to_numpy())
train_y_pred_x = train_pred_res
train_y_pred = train_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(train_i_y, train_y_pred).ravel()
train_sensitivity=tp / (tp + fn)
train_specificity=tn / (tn + fp)
train_accuracy=accuracy_score(train_i_y, train_y_pred)
train_bal_acc=balanced_accuracy_score(train_i_y, train_y_pred)
train_precision=precision_score(train_i_y, train_y_pred)
recall_train = recall_score(train_i_y, train_y_pred)
##AUC: cv-train
train_fpr_roc, train_tpr_roc, train_thresholds_roc = roc_curve(train_i_y, train_y_pred_x)
train_roc = roc_auc_score(train_i_y, train_y_pred_x)
##PR: cv-train
train_precision_pr, train_recall_pr, train_thresholds_pr = precision_recall_curve(train_i_y, train_y_pred_x)
train_average_precision = average_precision_score(train_i_y, train_y_pred_x)

f_perform.write("\t".join([str(feature_folder),str(train_accuracy),str(train_sensitivity),str(train_specificity),str(train_precision),str(train_bal_acc),str(train_roc),str(train_average_precision)])+"\t")

print("train_roc AUC:",train_roc)
train_dy = train_i_y.to_numpy().ravel()-train_y_pred_x
train_y_df = pd.DataFrame({"trueY":train_i_y.values.ravel(),"predY": train_y_pred_x,"dY":train_dy},index=train_index_ls)
train_y_df.to_csv(output_dir+"/train_y.csv")

# val performance
val_pred_res = clf.predict(val_i_X.to_numpy())
val_y_pred_x = val_pred_res
val_y_pred = val_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(val_i_y, val_y_pred).ravel()
val_sensitivity=tp / (tp + fn)
val_specificity=tn / (tn + fp)
val_accuracy=accuracy_score(val_i_y, val_y_pred)
val_bal_acc=balanced_accuracy_score(val_i_y, val_y_pred)
val_precision=precision_score(val_i_y, val_y_pred)
recall_val = recall_score(val_i_y, val_y_pred)
##AUC: cv-val
val_fpr_roc, val_tpr_roc, val_thresholds_roc = roc_curve(val_i_y, val_y_pred_x)
val_roc = roc_auc_score(val_i_y, val_y_pred_x)
##PR: cv-val
val_precision_pr, val_recall_pr, val_thresholds_pr = precision_recall_curve(val_i_y, val_y_pred_x)
val_average_precision = average_precision_score(val_i_y, val_y_pred_x)

f_perform.write("\t".join([str(val_accuracy),str(val_sensitivity),str(val_specificity),str(val_precision),str(val_bal_acc),str(val_roc),str(val_average_precision)])+"\t")

print("val_roc AUC:",val_roc)
val_dy = val_i_y.to_numpy().ravel()-val_y_pred_x
val_y_df = pd.DataFrame({"trueY":val_i_y.values.ravel(),"predY": val_y_pred_x,"dY":val_dy},index=val_index_ls)
val_y_df.to_csv(output_dir+"/val_y.csv")

# test performance
test_pred_res = clf.predict(test_X.to_numpy())
test_y_pred_x = test_pred_res
test_y_pred = test_y_pred_x > 0.5
test_accuracy=accuracy_score(test_y, test_y_pred)
test_roc = roc_auc_score(test_y, test_y_pred_x)
print("test_roc AUC:",test_roc)
test_dy = test_y.to_numpy().ravel()-test_y_pred_x
test_y_df = pd.DataFrame({"trueY":test_y.values.ravel(),"predY": test_y_pred_x,"dY":test_dy},index=test_index_ls)
test_y_df.to_csv(output_dir+"/test_y.csv")

tn, fp, fn, tp = confusion_matrix(test_y, test_y_pred).ravel()
test_sensitivity=tp / (tp + fn)
test_specificity=tn / (tn + fp)
test_accuracy=accuracy_score(test_y, test_y_pred)
test_bal_acc=balanced_accuracy_score(test_y, test_y_pred)
test_precision=precision_score(test_y, test_y_pred)
recall_val = recall_score(test_y, test_y_pred)
##AUC: cv-val
test_fpr_roc, test_tpr_roc, test_thresholds_roc = roc_curve(test_y, test_y_pred_x)
test_roc = roc_auc_score(test_y, test_y_pred_x)
##PR: cv-val
test_precision_pr, test_recall_pr, test_thresholds_pr = precision_recall_curve(test_y, test_y_pred_x)
test_average_precision = average_precision_score(test_y, test_y_pred_x)

f_perform.write("\t".join([str(test_accuracy),str(test_sensitivity),str(test_specificity),str(test_precision),str(test_bal_acc),str(test_roc),str(test_average_precision)])+"\n")

coeff_list = clf.coef_
intercept = clf.intercept_

df1 = pd.DataFrame({"feature":best_features,"coef":coeff_list })

# scale_method = "mm"
#############################################################################################
# %%
train_data_df_new = train_data_df.loc[train_index_ls,PRS_features]
train_data_df_new = train_data_df_new.loc[~train_data_df_new.isna().any(axis=1), :]
train_index_ls = train_data_df_new.index.tolist()
train_y_df = train_y_df.loc[train_index_ls,:]
train_i_y = train_y_df["trueY"].values.ravel()

val_data_df_new = train_data_df.loc[val_index_ls,PRS_features]
val_data_df_new = val_data_df_new.loc[~val_data_df_new.isna().any(axis=1), :]
val_index_ls = val_data_df_new.index.tolist()
val_y_df = val_y_df.loc[val_index_ls,:]
val_i_y = val_y_df["trueY"].values.ravel()

test_data_df_new = test_data_df.loc[test_index_ls,PRS_features]
test_data_df_new = test_data_df_new.loc[~test_data_df_new.isna().any(axis=1), :]
test_index_ls = test_data_df_new.index.tolist()
test_y_df = test_y_df.loc[test_index_ls,:]
test_y = test_y_df["trueY"].values.ravel()

## log transform
if scale_method.startswith("log"):
    print("log")
    train_data_df_new = np.log10(train_data_df_new + 1)
    val_data_df_new = np.log10(val_data_df_new + 1)
    test_data_df_new = np.log10(test_data_df_new + 1)

## MinMax
if scale_method.endswith("mm"):
    print("mm")
    train_data_df_new = minmax_scale(train_data_df_new, feature_range=(0, 1))
    val_data_df_new = minmax_scale(val_data_df_new, feature_range=(0, 1))
    test_data_df_new = minmax_scale(test_data_df_new, feature_range=(0, 1))

    train_data_df_new = pd.DataFrame(data=train_data_df_new, index=train_index_ls, columns=PRS_features)
    val_data_df_new = pd.DataFrame(data=val_data_df_new, index=val_index_ls, columns=PRS_features)
    test_data_df_new = pd.DataFrame(data=test_data_df_new, index=test_index_ls, columns=PRS_features)

## Standardization
if scale_method.endswith("z"):
    print("z")
    train_data_df_new = stats.zscore(train_data_df_new)
    val_data_df_new = stats.zscore(val_data_df_new)
    test_data_df_new = stats.zscore(test_data_df_new)

# clf = Lasso(alpha=0.0,max_iter=3000,random_state=12)
clf = LinearRegression()
clf.fit(train_data_df_new.to_numpy(), train_y_df["dY"])

# train performance
train_pred_res_new = clf.predict(train_data_df_new.to_numpy())
train_y_pred_x = train_pred_res_new +  train_y_df["predY"].values
train_y_pred = train_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(train_i_y, train_y_pred).ravel()
train_sensitivity=tp / (tp + fn)
train_specificity=tn / (tn + fp)
train_accuracy=accuracy_score(train_i_y, train_y_pred)
train_bal_acc=balanced_accuracy_score(train_i_y, train_y_pred)
train_precision=precision_score(train_i_y, train_y_pred)
recall_train = recall_score(train_i_y, train_y_pred)
##AUC: cv-train
train_fpr_roc, train_tpr_roc, train_thresholds_roc = roc_curve(train_i_y, train_y_pred_x)
train_roc = roc_auc_score(train_i_y, train_y_pred_x)
##PR: cv-train
train_precision_pr, train_recall_pr, train_thresholds_pr = precision_recall_curve(train_i_y, train_y_pred_x)
train_average_precision = average_precision_score(train_i_y, train_y_pred_x)

f_perform.write("\t".join([str(feature_folder+1),str(train_accuracy),str(train_sensitivity),str(train_specificity),str(train_precision),str(train_bal_acc),str(train_roc),str(train_average_precision)])+"\t")


print("train_roc AUC:",train_roc)
train_dy = train_y_df["trueY"].to_numpy().ravel()-train_y_pred_x
train_y_df = pd.DataFrame({"trueY":train_y_df["trueY"].values.ravel(),"predY": train_y_pred_x,"dY":train_dy},index=train_index_ls)
train_y_df.to_csv(output_dir+"/train_y_PRS.csv")

# val performance
val_pred_res_new = clf.predict(val_data_df_new.to_numpy())
val_y_pred_x = val_pred_res_new + val_y_df["predY"].values
val_y_pred = val_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(val_i_y, val_y_pred).ravel()
val_sensitivity=tp / (tp + fn)
val_specificity=tn / (tn + fp)
val_accuracy=accuracy_score(val_i_y, val_y_pred)
val_bal_acc=balanced_accuracy_score(val_i_y, val_y_pred)
val_precision=precision_score(val_i_y, val_y_pred)
recall_val = recall_score(val_i_y, val_y_pred)
##AUC: cv-val
val_fpr_roc, val_tpr_roc, val_thresholds_roc = roc_curve(val_i_y, val_y_pred_x)
val_roc = roc_auc_score(val_i_y, val_y_pred_x)
##PR: cv-val
val_precision_pr, val_recall_pr, val_thresholds_pr = precision_recall_curve(val_i_y, val_y_pred_x)
val_average_precision = average_precision_score(val_i_y, val_y_pred_x)

f_perform.write("\t".join([str(val_accuracy),str(val_sensitivity),str(val_specificity),str(val_precision),str(val_bal_acc),str(val_roc),str(val_average_precision)])+"\t")

print("val_roc AUC:",val_roc)
val_dy = val_y_df["trueY"].to_numpy().ravel()-val_y_pred_x
val_y_df = pd.DataFrame({"trueY":val_y_df["trueY"].values.ravel(),"predY": val_y_pred_x,"dY":val_dy},index=val_index_ls)
val_y_df.to_csv(output_dir+"/val_y_PRS.csv")

# test performance
test_pred_res_new = clf.predict(test_data_df_new.to_numpy())
test_y_pred_x = test_pred_res_new + test_y_df["predY"].values
test_y_pred = test_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(test_y, test_y_pred).ravel()
test_sensitivity=tp / (tp + fn)
test_specificity=tn / (tn + fp)
test_accuracy=accuracy_score(test_y, test_y_pred)
test_bal_acc=balanced_accuracy_score(test_y, test_y_pred)
test_precision=precision_score(test_y, test_y_pred)
recall_val = recall_score(test_y, test_y_pred)
##AUC: cv-val
test_fpr_roc, test_tpr_roc, test_thresholds_roc = roc_curve(test_y, test_y_pred_x)
test_roc = roc_auc_score(test_y, test_y_pred_x)
##PR: cv-val
test_precision_pr, test_recall_pr, test_thresholds_pr = precision_recall_curve(test_y, test_y_pred_x)
test_average_precision = average_precision_score(test_y, test_y_pred_x)

f_perform.write("\t".join([str(test_accuracy),str(test_sensitivity),str(test_specificity),str(test_precision),str(test_bal_acc),str(test_roc),str(test_average_precision)])+"\n")

print("test_roc AUC:",test_roc)
test_dy = test_y_df["trueY"].to_numpy().ravel()-test_y_pred_x
test_y_df = pd.DataFrame({"trueY":test_y_df["trueY"].values.ravel(),"predY": test_y_pred_x,"dY":test_dy},index=test_index_ls)
test_y_df.to_csv(output_dir+"/test_y_PRS.csv")

coeff_list = clf.coef_
intercept = clf.intercept_

df2 = pd.DataFrame({"feature":PRS_features,"coef":coeff_list })


#############################################################################################
# %%
train_data_df_new = train_data_df.loc[train_index_ls,clinicals_features]
train_data_df_new = train_data_df_new.loc[~train_data_df_new.isna().any(axis=1), :]
train_index_ls = train_data_df_new.index.tolist()
train_y_df = train_y_df.loc[train_index_ls,:]
train_i_y = train_y_df["trueY"].values.ravel()


val_data_df_new = train_data_df.loc[val_index_ls,clinicals_features]
val_data_df_new = val_data_df_new.loc[~val_data_df_new.isna().any(axis=1), :]
val_index_ls = val_data_df_new.index.tolist()
val_y_df = val_y_df.loc[val_index_ls,:]
val_i_y = val_y_df["trueY"].values.ravel()

test_data_df_new = test_data_df.loc[test_index_ls,clinicals_features]
test_data_df_new = test_data_df_new.loc[~test_data_df_new.isna().any(axis=1), :]
test_index_ls = test_data_df_new.index.tolist()
test_y_df = test_y_df.loc[test_index_ls,:]
test_y = test_y_df["trueY"].values.ravel()

# train_data_df_new.drop(['Clin_sex_int'], axis=1, inplace=True)
# val_data_df_new.drop(['Clin_sex_int'], axis=1, inplace=True)
# test_data_df_new.drop(['Clin_sex_int'], axis=1, inplace=True)

## log transform
if scale_method.startswith("log"):
    print("log")
    train_data_df_new = np.log10(train_data_df_new + 1)
    val_data_df_new = np.log10(val_data_df_new + 1)
    test_data_df_new = np.log10(test_data_df_new + 1)

## MinMax
if scale_method.endswith("mm"):
    print("mm")
    train_data_df_new = minmax_scale(train_data_df_new, feature_range=(0, 1))
    val_data_df_new = minmax_scale(val_data_df_new, feature_range=(0, 1))
    test_data_df_new = minmax_scale(test_data_df_new, feature_range=(0, 1))
    train_data_df_new = pd.DataFrame(data=train_data_df_new, index=train_index_ls, columns=["upsit","age","sex"])
    val_data_df_new = pd.DataFrame(data=val_data_df_new, index=val_index_ls, columns=["upsit","age","sex"])
    test_data_df_new = pd.DataFrame(data=test_data_df_new, index=test_index_ls, columns=["upsit","age","sex"])

## Standardization
if scale_method.endswith("z"):
    print("z")
    train_data_df_new = stats.zscore(train_data_df_new)
    val_data_df_new = stats.zscore(val_data_df_new)
    test_data_df_new = stats.zscore(test_data_df_new)


# train_data_df_new.loc[:, "Clin_sex_int"] = train_data_df.loc[train_index_ls, "Clin_sex_int"].values
# val_data_df_new.loc[:, "Clin_sex_int"] = train_data_df.loc[val_index_ls, "Clin_sex_int"].values
# test_data_df_new.loc[:, "Clin_sex_int"] = test_data_df.loc[test_index_ls, "Clin_sex_int"].values

# clf = Lasso(alpha=0.0,max_iter=3000,random_state=12)
clf = LinearRegression()
clf.fit(train_data_df_new.to_numpy(), train_y_df["dY"])

###### train performance
train_pred_res_new = clf.predict(train_data_df_new.to_numpy())
train_y_pred_x = train_pred_res_new +  train_y_df["predY"].values
train_y_pred = train_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(train_i_y, train_y_pred).ravel()
train_sensitivity=tp / (tp + fn)
train_specificity=tn / (tn + fp)
train_accuracy=accuracy_score(train_i_y, train_y_pred)
train_bal_acc=balanced_accuracy_score(train_i_y, train_y_pred)
train_precision=precision_score(train_i_y, train_y_pred)
recall_train = recall_score(train_i_y, train_y_pred)
##AUC: cv-train
train_fpr_roc, train_tpr_roc, train_thresholds_roc = roc_curve(train_i_y, train_y_pred_x)
train_roc = roc_auc_score(train_i_y, train_y_pred_x)
##PR: cv-train
train_precision_pr, train_recall_pr, train_thresholds_pr = precision_recall_curve(train_i_y, train_y_pred_x)
train_average_precision = average_precision_score(train_i_y, train_y_pred_x)

f_perform.write("\t".join([str(feature_folder+4),str(train_accuracy),str(train_sensitivity),str(train_specificity),str(train_precision),str(train_bal_acc),str(train_roc),str(train_average_precision)])+"\t")

print("train_roc AUC:",train_roc)
train_dy = train_y_df["trueY"].to_numpy().ravel()-train_y_pred_x
train_y_df = pd.DataFrame({"trueY":train_y_df["trueY"].values.ravel(),"predY": train_y_pred_x,"dY":train_dy},index=train_index_ls)
train_y_df.to_csv(output_dir+"/train_y_PRS_clinical.csv")

# val performance
val_pred_res_new = clf.predict(val_data_df_new.to_numpy())
val_y_pred_x = val_pred_res_new + val_y_df["predY"].values
val_y_pred = val_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(val_i_y, val_y_pred).ravel()
val_sensitivity=tp / (tp + fn)
val_specificity=tn / (tn + fp)
val_accuracy=accuracy_score(val_i_y, val_y_pred)
val_bal_acc=balanced_accuracy_score(val_i_y, val_y_pred)
val_precision=precision_score(val_i_y, val_y_pred)
recall_val = recall_score(val_i_y, val_y_pred)
##AUC: cv-val
val_fpr_roc, val_tpr_roc, val_thresholds_roc = roc_curve(val_i_y, val_y_pred_x)
val_roc = roc_auc_score(val_i_y, val_y_pred_x)
##PR: cv-val
val_precision_pr, val_recall_pr, val_thresholds_pr = precision_recall_curve(val_i_y, val_y_pred_x)
val_average_precision = average_precision_score(val_i_y, val_y_pred_x)

f_perform.write("\t".join([str(val_accuracy),str(val_sensitivity),str(val_specificity),str(val_precision),str(val_bal_acc),str(val_roc),str(val_average_precision)])+"\t")

print("val_roc AUC:",val_roc)
val_dy = val_y_df["trueY"].to_numpy().ravel()-val_y_pred_x
val_y_df = pd.DataFrame({"trueY":val_y_df["trueY"].values.ravel(),"predY": val_y_pred_x,"dY":val_dy},index=val_index_ls)
val_y_df.to_csv(output_dir+"/val_y_PRS_clinical.csv")

# test performance
test_pred_res_new = clf.predict(test_data_df_new.to_numpy())
test_y_pred_x = test_pred_res_new + test_y_df["predY"].values
test_y_pred = test_y_pred_x > 0.5

tn, fp, fn, tp = confusion_matrix(test_y, test_y_pred).ravel()
test_sensitivity=tp / (tp + fn)
test_specificity=tn / (tn + fp)
test_accuracy=accuracy_score(test_y, test_y_pred)
test_bal_acc=balanced_accuracy_score(test_y, test_y_pred)
test_precision=precision_score(test_y, test_y_pred)
recall_val = recall_score(test_y, test_y_pred)
##AUC: cv-val
test_fpr_roc, test_tpr_roc, test_thresholds_roc = roc_curve(test_y, test_y_pred_x)
test_roc = roc_auc_score(test_y, test_y_pred_x)
##PR: cv-val
test_precision_pr, test_recall_pr, test_thresholds_pr = precision_recall_curve(test_y, test_y_pred_x)
test_average_precision = average_precision_score(test_y, test_y_pred_x)

f_perform.write("\t".join([str(test_accuracy),str(test_sensitivity),str(test_specificity),str(test_precision),str(test_bal_acc),str(test_roc),str(test_average_precision)])+"\n")

print("test_roc AUC:",test_roc)
test_dy = test_y_df["trueY"].to_numpy().ravel()-test_y_pred_x
test_y_df = pd.DataFrame({"trueY":test_y_df["trueY"].values.ravel(),"predY": test_y_pred_x,"dY":test_dy},index=test_index_ls)
test_y_df.to_csv(output_dir+"/test_y_PRS_clinical.csv")
f_perform.close()

coeff_list = clf.coef_
intercept = clf.intercept_

df3 = pd.DataFrame({"feature":clinicals_features,"coef":coeff_list })

df = pd.concat([df1,df2,df3])
df.to_csv(output_dir+"/coef_df.csv")
