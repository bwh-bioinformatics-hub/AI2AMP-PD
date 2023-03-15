"""
Ruifeng Hu
10-13-2022
BWH-Boston
"""
import os
import sys

import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.metrics import roc_curve, precision_recall_curve
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score, average_precision_score

from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns




folder = "../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN"
# ==============================================================================
# The PDF document of AUROC
with PdfPages(folder + "/AUROC_merged.pdf") as auroc_pdf:
    # Plot whatever you wish to plot

    ####################
    fig1, ax_roc = plt.subplots(figsize=(8,8))
    ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=1.2, color='grey', label='Chance', alpha=.8)
    ## train
    de_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/test_y.csv", index_col=0,header=0)
    de_prs_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/test_y_PRS.csv", index_col=0,header=0)
    de_prs_clinical_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/test_y_PRS_clinical.csv", index_col=0,header=0)
    #####
    dfs = {"DEGs":de_df,"DEGs_PRS":de_prs_df,"DEGs_PRS_Clinical":de_prs_clinical_df}
    colors = {"DEGs":'#f8a289',"DEGs_PRS":"#b365a6","DEGs_PRS_Clinical":'#7cae01'}
    for df_i in dfs:
        y_pred_x = dfs[df_i].loc[:, "predY"].to_numpy()
        y_pred = y_pred_x > 0.5
        y_true = dfs[df_i].loc[:, "trueY"].to_numpy()
        AUROC = roc_auc_score(np.array(y_true, dtype=np.float32), y_pred_x)
        fpr_roc, tpr_roc, thresholds_roc = roc_curve(y_true, y_pred_x)
        label_str = df_i+"-{0:.2f}".format(AUROC * 100) + "%"
        ax_roc.plot(fpr_roc, tpr_roc, color=colors[df_i], label=label_str)
    # Done with the page

    ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="AUROC - Testing")
    ax_roc.set_xlabel("False positive rate (1-specificuty)")
    ax_roc.set_ylabel("True positive rate (sensitive)")
    ax_roc.legend(loc="lower right", fontsize=8, frameon=False)
    auroc_pdf.savefig(fig1)
    plt.close(fig1)

    ####################
    fig2, ax_roc = plt.subplots(figsize=(8, 8))
    ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=1.2, color='grey', label='Chance', alpha=.8)
    ## train
    de_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/val_y.csv", index_col=0, header=0)
    de_prs_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/val_y_PRS.csv", index_col=0, header=0)
    de_prs_clinical_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/val_y_PRS_clinical.csv",index_col=0, header=0)
    #####
    dfs = {"DEGs": de_df, "DEGs_PRS": de_prs_df, "DEGs_PRS_Clinical": de_prs_clinical_df}
    colors = {"DEGs": '#f8a289', "DEGs_PRS": "#b365a6", "DEGs_PRS_Clinical": '#7cae01'}
    for df_i in dfs:
        y_pred_x = dfs[df_i].loc[:, "predY"].to_numpy()
        y_pred = y_pred_x > 0.5
        y_true = dfs[df_i].loc[:, "trueY"].to_numpy()
        AUROC = roc_auc_score(np.array(y_true, dtype=np.float32), y_pred_x)
        fpr_roc, tpr_roc, thresholds_roc = roc_curve(y_true, y_pred_x)
        label_str = df_i + "-{0:.2f}".format(AUROC * 100) + "%"
        ax_roc.plot(fpr_roc, tpr_roc, color=colors[df_i], label=label_str)
    # Done with the page

    ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="AUROC - Validation" )
    ax_roc.set_xlabel("False positive rate (1-specificuty)")
    ax_roc.set_ylabel("True positive rate (sensitive)")
    ax_roc.legend(loc="lower right", fontsize=8, frameon=False)
    auroc_pdf.savefig(fig2)
    plt.close(fig2)

    ####################
    fig3, ax_roc = plt.subplots(figsize=(8, 8))
    ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=1.2, color='grey', label='Chance', alpha=.8)
    ## train
    de_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/train_y.csv", index_col=0, header=0)
    de_prs_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/train_y_PRS.csv", index_col=0, header=0)
    de_prs_clinical_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/train_y_PRS_clinical.csv", index_col=0, header=0)
    #####
    dfs = {"DEGs": de_df, "DEGs_PRS": de_prs_df, "DEGs_PRS_Clinical": de_prs_clinical_df}
    colors = {"DEGs": '#f8a289', "DEGs_PRS": "#b365a6", "DEGs_PRS_Clinical": '#7cae01'}
    for df_i in dfs:
        y_pred_x = dfs[df_i].loc[:, "predY"].to_numpy()
        y_pred = y_pred_x > 0.5
        y_true = dfs[df_i].loc[:, "trueY"].to_numpy()
        AUROC = roc_auc_score(np.array(y_true, dtype=np.float32), y_pred_x)
        fpr_roc, tpr_roc, thresholds_roc = roc_curve(y_true, y_pred_x)
        label_str = df_i + "-{0:.2f}".format(AUROC * 100) + "%"
        ax_roc.plot(fpr_roc, tpr_roc, color=colors[df_i], label=label_str)
    # Done with the page

    ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="AUROC - Training")
    ax_roc.set_xlabel("False positive rate (1-specificuty)")
    ax_roc.set_ylabel("True positive rate (sensitive)")
    ax_roc.legend(loc="lower right", fontsize=8, frameon=False)
    auroc_pdf.savefig(fig3)
    plt.close(fig3)


with PdfPages(folder + "/PRC_merged.pdf") as pr_pdf:
    ####################
    fig1, ax_pr = plt.subplots(figsize=(8, 8))
    ## train
    de_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/test_y.csv", index_col=0, header=0)
    de_prs_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/test_y_PRS.csv", index_col=0, header=0)
    de_prs_clinical_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/test_y_PRS_clinical.csv",
                                     index_col=0, header=0)
    #####
    dfs = {"DEGs": de_df, "DEGs_PRS": de_prs_df, "DEGs_PRS_Clinical": de_prs_clinical_df}
    colors = {"DEGs": '#f8a289', "DEGs_PRS": "#b365a6", "DEGs_PRS_Clinical": '#7cae01'}
    for df_i in dfs:
        y_pred_x = dfs[df_i].loc[:, "predY"].to_numpy()
        y_pred = y_pred_x > 0.5
        y_true = dfs[df_i].loc[:, "trueY"].to_numpy()
        pos_rate = y_true.mean()
        ax_pr.plot([0, 1], [pos_rate, pos_rate], linestyle='--', lw=1.2, color=colors[df_i], label='Baseline',alpha=.8)
        precision_pr, recall_pr, thresholds_pr = precision_recall_curve(y_true, y_pred_x)
        average_precision = average_precision_score(y_true, y_pred_x)
        recall_pr = np.insert(recall_pr, 0, 1.0)
        precision_pr = np.insert(precision_pr, 0, 0.0)
        label_str = df_i + " P-R (AUC = %0.2f)-{0:.2f}".format(average_precision * 100) + "%"
        ax_pr.plot(recall_pr,precision_pr, color=colors[df_i], label=label_str, lw=1.5, alpha=.8)

    ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="P-R - Testing" )
    ax_pr.set_xlabel("True positive rate (Recall)")
    ax_pr.set_ylabel("Positive predictive value (Precision)")
    ax_pr.legend(loc="lower left", fontsize=8, frameon=False)

    pr_pdf.savefig(fig1)
    plt.close(fig1)
    # Done with the page

    ####################
    fig1, ax_pr = plt.subplots(figsize=(8, 8))
    ## train
    de_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/val_y.csv", index_col=0, header=0)
    de_prs_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/val_y_PRS.csv", index_col=0, header=0)
    de_prs_clinical_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/val_y_PRS_clinical.csv",
                                     index_col=0, header=0)
    #####
    dfs = {"DEGs": de_df, "DEGs_PRS": de_prs_df, "DEGs_PRS_Clinical": de_prs_clinical_df}
    colors = {"DEGs": '#f8a289', "DEGs_PRS": "#b365a6", "DEGs_PRS_Clinical": '#7cae01'}
    for df_i in dfs:
        y_pred_x = dfs[df_i].loc[:, "predY"].to_numpy()
        y_pred = y_pred_x > 0.5
        y_true = dfs[df_i].loc[:, "trueY"].to_numpy()
        pos_rate = y_true.mean()
        ax_pr.plot([0, 1], [pos_rate, pos_rate], linestyle='--', lw=1.2, color=colors[df_i], label='Baseline',
                   alpha=.8)
        precision_pr, recall_pr, thresholds_pr = precision_recall_curve(y_true, y_pred_x)
        average_precision = average_precision_score(y_true, y_pred_x)
        recall_pr = np.insert(recall_pr, 0, 1.0)
        precision_pr = np.insert(precision_pr, 0, 0.0)
        label_str = df_i + " P-R (AUC = %0.2f)-{0:.2f}".format(average_precision * 100) + "%"
        ax_pr.plot(recall_pr, precision_pr, color=colors[df_i], label=label_str, lw=1.5, alpha=.8)

    ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="P-R - Validation")
    ax_pr.set_xlabel("True positive rate (Recall)")
    ax_pr.set_ylabel("Positive predictive value (Precision)")
    ax_pr.legend(loc="lower left", fontsize=8, frameon=False)

    pr_pdf.savefig(fig1)
    plt.close(fig1)
    # Done with the page

    ####################
    fig1, ax_pr = plt.subplots(figsize=(8, 8))
    ## train
    de_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/train_y.csv", index_col=0, header=0)
    de_prs_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/train_y_PRS.csv", index_col=0, header=0)
    de_prs_clinical_df = pd.read_csv("../results/only_PPMI/LASSO_pool_1692/logz/LASSO_LN/train_y_PRS_clinical.csv",
                                     index_col=0, header=0)
    #####
    dfs = {"DEGs": de_df, "DEGs_PRS": de_prs_df, "DEGs_PRS_Clinical": de_prs_clinical_df}
    colors = {"DEGs": '#f8a289', "DEGs_PRS": "#b365a6", "DEGs_PRS_Clinical": '#7cae01'}
    for df_i in dfs:
        y_pred_x = dfs[df_i].loc[:, "predY"].to_numpy()
        y_pred = y_pred_x > 0.5
        y_true = dfs[df_i].loc[:, "trueY"].to_numpy()
        pos_rate = y_true.mean()
        ax_pr.plot([0, 1], [pos_rate, pos_rate], linestyle='--', lw=1.2, color=colors[df_i], label='Baseline',
                   alpha=.8)
        precision_pr, recall_pr, thresholds_pr = precision_recall_curve(y_true, y_pred_x)
        average_precision = average_precision_score(y_true, y_pred_x)
        recall_pr = np.insert(recall_pr, 0, 1.0)
        precision_pr = np.insert(precision_pr, 0, 0.0)
        label_str = df_i + " P-R (AUC = %0.2f)-{0:.2f}".format(average_precision * 100) + "%"
        ax_pr.plot(recall_pr, precision_pr, color=colors[df_i], label=label_str, lw=1.5, alpha=.8)

    ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="P-R - Training")
    ax_pr.set_xlabel("True positive rate (Recall)")
    ax_pr.set_ylabel("Positive predictive value (Precision)")
    ax_pr.legend(loc="lower left", fontsize=8, frameon=False)

    pr_pdf.savefig(fig1)
    plt.close(fig1)
    # Done with the page


anno_df = pd.read_csv("../../data/RNAseq/genes_annotation_rh.tsv", sep="\t",index_col=0, header=0)
### barplot:
importance_df = pd.read_csv("../results/only_PPMI/DT_SVM_pool_1696/mm/LR_coef.txt", sep="\t",index_col=0, header=0)
# importance_df["CumSum"] = importance_df[["importance"]].cumsum()
# importance_df = importance_df.loc[importance_df["CumSum"]<0.91]

importance_df["abs_coef"] = np.abs(importance_df["coef"])
importance_df.sort_values(by=['abs_coef'],inplace=True,ascending=False)

importance_df = importance_df.iloc[0:20,:]

importance_df["feature"] = [i[4:] if i.startswith("DEG") else i for i in importance_df["feature"]]

importance_df = pd.merge(importance_df,anno_df["gene_name"],how="left",left_on="feature",right_index=True)
importance_df["gene_label"] = importance_df["feature"] + " (" + importance_df["gene_name"] + ")"

importance_df.loc[importance_df["gene_label"].isna(),"gene_label"] = importance_df["feature"]
importance_df.loc[importance_df["gene_name"].isna(),"gene_name"] = importance_df["feature"]

fig,ax = plt.subplots(figsize=(8,12))
colors = ['r' if c >= 0 else 'b' for c in importance_df["coef"]]
sns.barplot(data=importance_df, x="coef", y="gene_name",palette=colors,ax=ax)
plt.ylabel("Features")
plt.tight_layout()
plt.savefig("../results/only_PPMI/DT_SVM_pool_1696/mm/Top20_feature_importance_bar_Symbol.pdf")
plt.close()




# x_label = ["Known genes","Non-coding genes","eRNAs", "circRNAs"]
# data = [0.6288,0.6199,0.6174,0.5624]
# plt.bar(x_label,data,color='#f8a289')
# plt.ylim(0.5,0.7)
# plt.savefig(folder + "/Bar_genes_AUC.pdf")

