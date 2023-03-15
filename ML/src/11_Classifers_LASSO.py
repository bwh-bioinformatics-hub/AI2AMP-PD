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
from sklearn.model_selection import StratifiedKFold,train_test_split

from sklearn.linear_model import LogisticRegression,Lasso,LassoLars,lasso_path,lars_path

from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score, average_precision_score
from sklearn.metrics import mean_squared_error

from matplotlib.backends.backend_pdf import PdfPages

from funcs import *

##############
# %%
## Global params
# ["SVM","LR","SGD","ABC","XGB","RF","KNN","MLP"]
# model = sys.argv[1]
scale_method = "logz"
feature_selection= "DT"
dataset = "only_PPMI"

################
# %%
## prepare dataset
print("*******************************************")
# train
train_data_file = "data/PPMIvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
train_data_df = pd.read_csv(train_data_file,sep="\t",index_col=0,header=0)

train_DEG = train_data_df.loc[:,train_data_df.columns.str.startswith("DEG")]
train_eRNA = train_data_df.loc[:,train_data_df.columns.str.startswith("eRNA")]
train_circRNA = train_data_df.loc[:,train_data_df.columns.str.startswith("circRNA")]
train_PRS = train_data_df.loc[:,train_data_df.columns.str.startswith("PRS")]
train_Clinical = train_data_df.loc[:,train_data_df.columns.str.startswith("Clin")]
DEGs = train_DEG.columns.tolist()
eRNAs = train_eRNA.columns.tolist()
circRNAs = train_circRNA.columns.tolist()
PRS_features = ["PRS"]
clinicals_features = train_Clinical.columns.tolist() ## ["upsit_total_score","sex_int","age_at_baseline"]

# test
test_data_file = "data/PDBFvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
test_data_df = pd.read_csv(test_data_file,sep="\t",index_col=0,header=0)


## prepare features
if feature_selection.startswith("LASSO"):
    if os.path.exists("../results/"+dataset+"/mm/DEG/lasso_feature_importances.txt"):
        feature_df = pd.read_csv("../results/"+dataset+"/mm/DEG/lasso_feature_importances.txt", sep="\t")
        DEG_features = list(feature_df.loc[feature_df["importance"] > 0, "gene"].values)

        feature_df = pd.read_csv("../results/"+dataset+"/mm/eRNA/lasso_feature_importances.txt", sep="\t")
        eRNA_features = list(feature_df.loc[feature_df["importance"] > 0, "gene"].values)
        
        feature_df = pd.read_csv("../results/"+dataset+"/mm/circRNA/lasso_feature_importances.txt", sep="\t")
        circRNA_features = list(feature_df.loc[feature_df["importance"] > 0, "gene"].values)
    else:
        raise("Error: not feature files")

elif feature_selection.startswith("DT"):
    if os.path.exists("../results/"+dataset+"/mm/DEG/dt_feature_importances.txt"):
        feature_df = pd.read_csv("../results/"+dataset+"/mm/DEG/dt_feature_importances.txt", sep="\t")
        DEG_features = list(feature_df.loc[feature_df["importance"] > 0, "gene"].values)
        print("DEG:",len(DEG_features))

        feature_df = pd.read_csv("../results/"+dataset+"/mm/eRNA/dt_feature_importances.txt", sep="\t")
        eRNA_features = list(feature_df.loc[feature_df["importance"] > 0, "gene"].values)
        print("eRNA:",len(eRNA_features))

        feature_df = pd.read_csv("../results/"+dataset+"/mm/circRNA/dt_feature_importances.txt", sep="\t")
        circRNA_features = list(feature_df.loc[feature_df["importance"] > 0, "gene"].values)
        print("circRNA:",len(circRNA_features))

    else:
        raise("Error: not feature files")
    
kept_features = circRNAs
# kept_features = list(DEG_features) + list(eRNA_features) + list(circRNA_features)
####
# best_feature_df_mRNA = pd.read_csv("../results/only_PPMI/LASSO_pool_874/logz/LASSO_LN/best_features_LN.csv",index_col=0)
# best_feature_df_ecirc = pd.read_csv("../results/only_PPMI/LASSO_pool_818/logz/LASSO_LN/best_features_LN.csv",index_col=0)
# best_features = best_feature_df_mRNA.index.tolist() + best_feature_df_ecirc.index.tolist()
# kept_features = best_features

print("Total feature number:", len(kept_features))

###################
# %%
print("=====================")
output_dir = "../results/"+dataset+"/LASSO_pool_"+str(len(kept_features))+"/"+scale_method+"/LASSO_LN"
print(output_dir)
# Create folder if the directory doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir,exist_ok=True)

######################
# %%
train_X = train_data_df.loc[:,kept_features]
train_X = train_X.loc[~train_X.isna().any(axis=1),:]
train_y = train_data_df.loc[train_X.index.tolist(),["Y"]]
train_i_X, val_i_X, train_i_Y, val_i_Y = train_test_split(train_X, train_y, test_size=0.2, shuffle=True, random_state=42)

train_index_ls =train_i_X.index.tolist()
val_index_ls = val_i_X.index.tolist()

print("Selected feature number:", len(kept_features))
train_i_X =train_i_X.loc[:,kept_features]
val_i_X =val_i_X.loc[:,kept_features]

# if "Clin_sex_int" in kept_features:
#     train_i_X.drop(['Clin_sex_int'], axis=1,inplace=True)
#     val_i_X.drop(['Clin_sex_int'], axis=1,inplace=True)

## log transform
if scale_method.startswith("log"):
    print("log")
    train_i_X = np.log10(train_i_X+1)
    val_i_X = np.log10(val_i_X+1)

## MinMax
if scale_method.endswith("mm"):
    print("mm")
    train_i_X = minmax_scale(train_i_X, feature_range=(0, 1))
    val_i_X = minmax_scale(val_i_X, feature_range=(0, 1))
    # scaler = MinMaxScaler()
    # train_i_X = scaler.fit_transform(train_i_X)
    # val_i_X = scaler.transform(val_i_X)
    train_i_X = pd.DataFrame(data=train_i_X,index=train_index_ls,columns=kept_features)
    val_i_X = pd.DataFrame(data=val_i_X,index=val_index_ls,columns=kept_features)

## Standardization
if scale_method.endswith("z"):
    print("z")
    train_i_X = stats.zscore(train_i_X)
    val_i_X = stats.zscore(val_i_X)
    # scaler = StandardScaler()
    # train_i_X = scaler.fit_transform(train_i_X)
    # val_i_X = scaler.transform(val_i_X)

# if "Clin_sex_int" in kept_features:
#     train_i_X.loc[:,"Clin_sex_int"] = train_Clinical.loc[train_index_ls,"Clin_sex_int"].values
#     val_i_X.loc[:,"Clin_sex_int"] = train_Clinical.loc[val_index_ls,"Clin_sex_int"].values

print("Final train size:", train_i_X.shape)
print("Final val size:", val_i_X.shape)

###################################################################
 # %%
print("Test data...")
## Test Data

print("Selected feature number:", len(kept_features))
test_X =test_data_df.loc[:,kept_features]
test_X = test_X.loc[~test_X.isna().any(axis=1),:]
test_y = test_data_df.loc[test_X.index.tolist(),["Y"]]
test_index_ls = test_X.index.tolist()
print("Final test size:", test_X.shape)

# if "Clin_sex_int" in kept_features:
#     test_X.drop(['Clin_sex_int'], axis=1,inplace=True)

## log transform
if scale_method.startswith("log"):
    print("log")
    test_X = np.log10(test_X+1)

## Scale
if scale_method.endswith("mm"):
    print("mm")
    test_X = minmax_scale(test_X, feature_range=(0, 1))
    # test_X = scaler.transform(test_X)
    test_X = pd.DataFrame(data=test_X,index=test_index_ls,columns=kept_features)


if scale_method.endswith("z"):
    print("z")
    test_X = stats.zscore(test_X)
    # test_X = scaler.transform(test_X)

# if "Clin_sex_int" in kept_features:
#     # test_X=np.append(test_X, test_data_df.loc[test_index_ls,["Clin_sex_int"]].values, axis=1)
#     test_X.loc[:, "Clin_sex_int"] = test_data_df.loc[test_index_ls, "Clin_sex_int"].values


######
# %%
# Compute paths
print("Computing regularization path using the lasso...")

# eps = 5e-3  # the smaller it is the longer is the path
# alpha_ls = np.linspace(0.0005,1,1000)
# alphas_lasso, coefs_lasso, _ = lasso_path(train_i_X, train_i_Y, alphas=alpha_ls)
alphas_lasso, active, coefs_lasso = lars_path(train_i_X.to_numpy(), train_i_Y.to_numpy().ravel(),method='lasso',max_iter=1000)

from itertools import cycle
fig = plt.figure()
colors = cycle(["b", "r", "g", "c", "k"])
log_alphas = np.log10(alphas_lasso)
for coef_l, c in zip(coefs_lasso, colors):
    l1 = plt.plot(log_alphas, coef_l, c=c)
plt.xlabel("Log10(alpha)")
plt.ylabel("coefficients")
plt.title("Lasso Paths")
plt.legend((l1[-1], ), ("Lasso",), loc="lower right")
plt.axis("tight")
plt.savefig(output_dir+"/coefficients_alpha.pdf")

alpha_ls_str = [str(round(i,3)) for i in alphas_lasso]
feature_coef = pd.DataFrame(data=coefs_lasso,index=kept_features,columns=alpha_ls_str)
feature_coef.to_csv(output_dir+"/feature_coefficients_alpha_LN.csv")

##################################
#%%
train_scores ={"auc":[],"accuracy":[]}
val_scores ={"auc":[],"accuracy":[]}
test_scores ={"auc":[],"accuracy":[]}

n = 320
alpha_ls = np.linspace(0.001,0.321,n,endpoint=False)
coefs_matrix = np.zeros(shape=(train_i_X.shape[1],n))
n_i = 0
for alpha_i in alpha_ls:
    print(n_i)
    clf = Lasso(alpha=alpha_i,max_iter=3000,random_state=12)
    clf.fit(train_i_X.to_numpy(), train_i_Y.to_numpy().ravel())
    # train performance
    train_pred_res = clf.predict(train_i_X.to_numpy())
    train_y_pred_x = train_pred_res
    train_y_pred = train_y_pred_x > 0.5
    train_accuracy=accuracy_score(train_i_Y, train_y_pred)
    train_roc = roc_auc_score(train_i_Y, train_y_pred_x)
    train_scores["auc"].append(train_roc)
    train_scores["accuracy"].append(train_accuracy)

    # test performance
    test_pred_res = clf.predict(test_X.to_numpy())
    test_y_pred_x = test_pred_res
    test_y_pred = test_y_pred_x > 0.5
    test_accuracy=accuracy_score(test_y, test_y_pred)
    test_roc = roc_auc_score(test_y, test_y_pred_x)
    test_scores["auc"].append(test_roc)
    test_scores["accuracy"].append(test_accuracy)

    # val performance
    val_pred_res = clf.predict(val_i_X.to_numpy())
    val_y_pred_x = val_pred_res
    val_y_pred = val_y_pred_x > 0.5
    val_accuracy=accuracy_score(val_i_Y, val_y_pred)
    val_roc = roc_auc_score(val_i_Y, val_y_pred_x)
    val_scores["auc"].append(val_roc)
    val_scores["accuracy"].append(val_accuracy)


    coefs_i = clf.coef_
    coefs_matrix[:,n_i] = coefs_i
    n_i += 1

from itertools import cycle
fig = plt.figure()
colors = cycle(["b", "r", "g", "c", "k"])
log_alphas = np.log10(alpha_ls)
for coef_l, c in zip(coefs_matrix, colors):
    l1 = plt.plot(log_alphas, coef_l, c=c)
plt.xlabel("Log10(alpha)")
plt.ylabel("coefficients")
plt.title("Lasso Paths")
plt.legend((l1[-1], ), ("Lasso",), loc="lower right")
plt.axis("tight")
plt.savefig(output_dir+"/coefficients_alpha_loop_LN.pdf")

fig = plt.figure()
l1 = plt.plot(log_alphas,train_scores["auc"], c="b",label="Train")
l2 = plt.plot(log_alphas,val_scores["auc"], c="y",label="Val")
l3 = plt.plot(log_alphas,test_scores["auc"], c="g",label="Test")
plt.xlabel("Log10(alpha)")
plt.ylabel("AUROC")
plt.title("AUC on different alphas")
plt.legend(loc="upper right")
plt.axis("tight")
plt.savefig(output_dir+"/alpha_performances_AUC_loop_LN.pdf")

fig = plt.figure()
l1 = plt.plot(log_alphas,train_scores["accuracy"], c="b",label="Train")
l2 = plt.plot(log_alphas,val_scores["accuracy"], c="y",label="Val")
l3 = plt.plot(log_alphas,test_scores["accuracy"], c="g",label="Test")
plt.xlabel("Log10(alpha)")
plt.ylabel("Accuracy")
plt.title("Accuracy on different alphas")
plt.legend(loc="upper right")
plt.axis("tight")
plt.savefig(output_dir+"/alpha_performances_accuracy_loop_LN.pdf")

alpha_ls_str = [str(round(i,3)) for i in alpha_ls]
data = {"alpha": alpha_ls_str,"log_alpha": log_alphas,
        "train_auc":train_scores["auc"],"train_accuracy":train_scores["accuracy"],
        "val_auc":val_scores["auc"],"val_accuracy":val_scores["accuracy"],
        "test_auc":test_scores["auc"],"test_accuracy":test_scores["accuracy"],}
feature_coef = pd.DataFrame(data=data)
feature_coef.to_csv(output_dir+"/perforance_alpha_loop_LN.csv")

alpha_ls_str = [str(round(i,3)) for i in alpha_ls]
feature_coef = pd.DataFrame(data=coefs_matrix,index=kept_features,columns=alpha_ls_str)
feature_coef.to_csv(output_dir+"/feature_coefficients_alpha_loop_LN.csv")


# %%
#################################
print("=====================")
from sklearn.svm import l1_min_c
output_dir = "../results/"+dataset+"/LASSO_pool_"+str(len(kept_features))+"/"+scale_method+"/LASSO_LR"
print(output_dir)
# Create folder if the directory doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir,exist_ok=True)

train_scores ={"auc":[],"accuracy":[]}
val_scores ={"auc":[],"accuracy":[]}
test_scores ={"auc":[],"accuracy":[]}

n=100
# cs = l1_min_c(train_i_X.to_numpy(), train_i_Y.to_numpy().ravel(), loss="log") * np.logspace(0, 5, n,endpoint=False)
cs = 0.005 * np.logspace(0,5,n,endpoint=False)
coefs_matrix = np.zeros(shape=(train_i_X.shape[1],n))
n_i = 0
for c_i in cs:
    print(n_i)
    clf = LogisticRegression(random_state=42,penalty="l1",solver="liblinear",C=c_i)
    clf.fit(train_i_X.to_numpy(), train_i_Y.to_numpy().ravel())

    # train performance
    train_pred_res = clf.predict(train_i_X.to_numpy())
    train_y_pred_x = train_pred_res
    train_y_pred = train_y_pred_x > 0.5
    train_accuracy=accuracy_score(train_i_Y, train_y_pred)
    train_roc = roc_auc_score(train_i_Y, train_y_pred_x)
    train_scores["auc"].append(train_roc)
    train_scores["accuracy"].append(train_accuracy)

    # test performance
    test_pred_res = clf.predict(test_X.to_numpy())
    test_y_pred_x = test_pred_res
    test_y_pred = test_y_pred_x > 0.5
    test_accuracy=accuracy_score(test_y, test_y_pred)
    test_roc = roc_auc_score(test_y, test_y_pred_x)
    test_scores["auc"].append(test_roc)
    test_scores["accuracy"].append(test_accuracy)

    # val performance
    val_pred_res = clf.predict(val_i_X.to_numpy())
    val_y_pred_x = val_pred_res
    val_y_pred = val_y_pred_x > 0.5
    val_accuracy=accuracy_score(val_i_Y, val_y_pred)
    val_roc = roc_auc_score(val_i_Y, val_y_pred_x)
    val_scores["auc"].append(val_roc)
    val_scores["accuracy"].append(val_accuracy)



    coefs_i = clf.coef_
    coefs_matrix[:,n_i] = coefs_i
    n_i += 1

fig = plt.figure()
log_cs = np.log10(cs)
colors = cycle(["b", "r", "g", "c", "k"])
for coef_l, c in zip(coefs_matrix, colors):
    l1 = plt.plot(log_cs, coef_l, c=c)
plt.xlabel("Log10(C)")
plt.ylabel("coefficients")
plt.title("Lasso Paths")
plt.legend((l1[-1], ), ("Lasso",), loc="lower left")
plt.axis("tight")
plt.savefig(output_dir+"/coefficients_alpha_loop_LR.pdf")


fig = plt.figure()
l1 = plt.plot(log_cs,train_scores["auc"], c="b",label="Train")
l2 = plt.plot(log_cs,val_scores["auc"], c="y",label="Val")
l3 = plt.plot(log_cs,test_scores["auc"], c="g",label="Test")
plt.xlabel("Log10(C)")
plt.ylabel("AUROC")
plt.title("AUC on different alphas")
plt.legend(loc="best")
plt.axis("tight")
plt.savefig(output_dir+"/alpha_performances_AUC_loop_LR.pdf")

fig = plt.figure()
l1 = plt.plot(log_cs,train_scores["accuracy"], c="b",label="Train")
l2 = plt.plot(log_cs,val_scores["accuracy"], c="y",label="Val")
l3 = plt.plot(log_cs,test_scores["accuracy"], c="g",label="Test")
plt.xlabel("Log10(C)")
plt.ylabel("Accuracy")
plt.title("Accuracy on different alphas")
plt.legend(loc="best")
plt.axis("tight")
plt.savefig(output_dir+"/alpha_performances_accuracy_loop_LR.pdf")

logcs_str = [str(round(i,3)) for i in log_cs]
data = {"logC": logcs_str,
        "train_auc":train_scores["auc"],"train_accuracy":train_scores["accuracy"],
        "val_auc":val_scores["auc"],"val_accuracy":val_scores["accuracy"],
        "test_auc":test_scores["auc"],"test_accuracy":test_scores["accuracy"],}
feature_coef = pd.DataFrame(data=data)
feature_coef.to_csv(output_dir+"/perforance_alpha_loop_LR.csv")


logcs_str = [str(round(i,3)) for i in log_cs]
feature_coef = pd.DataFrame(data=coefs_matrix,index=kept_features,columns=logcs_str)
feature_coef.to_csv(output_dir+"/feature_coefficients_alpha_loop_LR.csv")


# %%
