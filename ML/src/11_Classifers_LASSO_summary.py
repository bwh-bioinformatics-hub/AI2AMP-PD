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

def get_corr_pairs(df):
    '''Get diagonal and lower triangular pairs of correlation matrix'''
    pairs_df = {"f1":[],"f2":[],"corr_val":[]}
    cols = df.columns
    for i in range(len(cols)):
        for j in range(0, i):
            pairs_df["f1"].append(cols[i])
            pairs_df["f2"].append(cols[j])
            pairs_df["corr_val"].append(df.iloc[i, j])
    return pd.DataFrame(data=pairs_df)


## train data
train_data_file = "data/PPMIvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
train_data_df = pd.read_csv(train_data_file,sep="\t",index_col=0,header=0)

## test
test_data_file = "data/PDBFvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
test_data_df = pd.read_csv(test_data_file,sep="\t",index_col=0,header=0)


##############
# %%
## Global params
# ["SVM","LR","SGD","ABC","XGB","RF","KNN","MLP"]
# model = sys.argv[1]
scale_method = "logz"
feature_selection= "DT"
dataset = "only_PPMI"
kept_features = 1692
print("=====================")
output_dir = "../results/"+dataset+"/LASSO_pool_"+str(kept_features)+"/"+scale_method+"/LASSO_LN"
print(output_dir)

# %%
performance_coef_df = pd.read_csv(output_dir+"/perforance_alpha_loop_LN.csv",index_col=0)
best_test_auc = performance_coef_df["test_auc"].max()
best_alpha = performance_coef_df.loc[performance_coef_df["test_auc"]==best_test_auc,"alpha"].values[0]
best_logalpha = performance_coef_df.loc[performance_coef_df["test_auc"]==best_test_auc,"log_alpha"].values[0]

x = performance_coef_df.log_alpha.values
y = performance_coef_df.test_auc.values
plt.plot(x,y,label="Test_AUC",color="g")
plt.axvline(x=best_logalpha,color="red",linestyle="dashed")
s = "log$_{10}$(alpha)="+str(best_logalpha)[:7] + "\n(alpha="+str(best_alpha)+")"
plt.text(x=best_logalpha,y=0.5,s=s,horizontalalignment='center')
s = "AUC="+str(best_test_auc)[:7]
plt.text(x=best_logalpha,y=best_test_auc,s=s,horizontalalignment='center')
plt.xlabel("log10(alpha)")
plt.ylabel("AUROC")
plt.title("AUC with different alphas on Test dataset")
plt.legend()
plt.savefig(output_dir+"/best_alpha_loop_LN.pdf")

# %%
feature_coef_df = pd.read_csv(output_dir+"/feature_coefficients_alpha_loop_LN.csv",index_col=0)
best_feature_df = feature_coef_df.loc[feature_coef_df[str(best_alpha)]!=0,[str(best_alpha)]]
best_feature_df.to_csv(output_dir+"/best_features_LN.csv")
best_DEG = [i for i in best_feature_df.index if i.startswith("DEG")]
best_eRNA = [i for i in best_feature_df.index if i.startswith("eRNA")]
best_circRNA = [i for i in best_feature_df.index if i.startswith("circRNA")]

feature_count = {"total":[],"mRNA":[],"eRNA":[],"circRNA":[],"genes":[],"max_corr_train":[],"max_corr_test":[]}
for alpha_i in feature_coef_df.columns.tolist():
    sub_df = feature_coef_df.loc[feature_coef_df[alpha_i]!=0,:]
    feature_ls = sub_df.index.tolist()
    gene_str = "|".join(feature_ls)
    feature_count["genes"].append(gene_str)

    DEG_ls = [i for i in feature_ls if i.startswith("DEG")]
    eRNA_ls = [i for i in feature_ls if i.startswith("eRNA")]
    circRNA_ls = [i for i in feature_ls if i.startswith("circRNA")]
    feature_count["total"].append(len(feature_ls))
    feature_count["mRNA"].append(len(DEG_ls))
    feature_count["eRNA"].append(len(eRNA_ls))
    feature_count["circRNA"].append(len(circRNA_ls))

    ## corr
    if len(feature_ls) > 1:
        ## train
        sub_data_df = train_data_df.loc[:,feature_ls]
        sub_data_df = sub_data_df.loc[~sub_data_df.isna().any(axis=1),:]
        corr_df=sub_data_df.corr().abs()
        corr_df_long = get_corr_pairs(corr_df)
        corr_df_long.sort_values(by="corr_val",ascending=False,inplace=True)
        max_corr_train = corr_df_long["corr_val"].values[0]

        ##test
        sub_data_df = test_data_df.loc[:, feature_ls]
        sub_data_df = sub_data_df.loc[~sub_data_df.isna().any(axis=1), :]
        corr_df = sub_data_df.corr().abs()
        corr_df_long = get_corr_pairs(corr_df)
        corr_df_long.sort_values(by="corr_val", ascending=False, inplace=True)
        max_corr_test = corr_df_long["corr_val"].values[0]
    else:
        max_corr_test = 0
        max_corr_train = 0

    feature_count["max_corr_train"].append(max_corr_train)
    feature_count["max_corr_test"].append(max_corr_test)

feature_count = pd.DataFrame(data=feature_count, index=feature_coef_df.columns.tolist())
feature_count["test_auc"]=performance_coef_df["test_auc"].values
feature_count["log10alpha"]=performance_coef_df["log_alpha"].values
feature_count.to_csv(output_dir + "/feature_count_alpha_loop_LN.csv")

feature_count_log = np.log10(feature_count.loc[:,["total","mRNA","eRNA","circRNA"]]+1)
feature_count_log.loc[:,"log10alpha"] = feature_count["log10alpha"]

feature_count_log.plot.line(y=["mRNA"],x="log10alpha")
plt.axvline(x=best_logalpha,color="red",linestyle="dashed")
s = "log$_{10}$(alpha)="+str(best_logalpha)[:7] + "\n(alpha="+str(best_alpha)+")"
# plt.text(x=best_logalpha,y=2,s=s,horizontalalignment='center')
plt.xlabel("log10(alpha)")
plt.ylabel("log10(Count+1)")
plt.title("Feature counts with different alphas")
plt.savefig(output_dir+"/feature_count_alpha_loop_LN_lineplot.pdf")

feature_count_log.plot.area(y=["mRNA"],x="log10alpha")
plt.axvline(x=best_logalpha,color="red",linestyle="dashed")
s = "log$_{10}$(alpha)="+str(best_logalpha)[:7] + "\n(alpha="+str(best_alpha)+")"
# plt.text(x=best_logalpha,y=2,s=s,horizontalalignment='center')
plt.xlabel("log10(alpha)")
plt.ylabel("log10(Count+1)")
plt.title("Feature counts with different alphas")
plt.savefig(output_dir+"/feature_count_alpha_loop_LN_areaplot.pdf")

feature_count.plot.line(y=["max_corr_train","max_corr_test",],x="log10alpha")
plt.axvline(x=best_logalpha,color="red",linestyle="dashed")
s = "log$_{10}$(alpha)="+str(best_logalpha)[:7] + "\n(alpha="+str(best_alpha)+")"
# plt.text(x=best_logalpha,y=0.5,s=s,horizontalalignment='center')
plt.savefig(output_dir+"/max_corr_alpha_loop_LN.pdf")

######
# %%
## Order the features
alpha_ls = feature_coef_df.columns.tolist()[::-1]

feature_order = {}
order_x = 0
for a_i in alpha_ls:
    sub_df = feature_coef_df.loc[feature_coef_df[a_i] != 0, :]
    feature_ls = sub_df.index.tolist()
    new_feature = [i for i in feature_ls if i not in feature_order]
    if len(new_feature) > 0:
        order_x+=1
        for f_i in new_feature:
            if f_i not in feature_order:
                feature_order[f_i] = order_x

order_x+=1
for f_i in feature_coef_df.index.tolist():
    if f_i not in feature_order:
        feature_order[f_i] =order_x

feature_orde_df = pd.DataFrame.from_dict(data=feature_order,orient='index',columns=["Order"])
feature_orde_df.to_csv(output_dir+"/feature_ordering_LN.csv")

