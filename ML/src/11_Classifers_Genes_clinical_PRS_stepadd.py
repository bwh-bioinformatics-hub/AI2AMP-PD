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

from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression,SGDClassifier,LinearRegression
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier,GradientBoostingClassifier
from sklearn.svm import SVC

from sklearn.metrics import roc_curve, precision_recall_curve, confusion_matrix, balanced_accuracy_score
from sklearn.metrics import auc, accuracy_score, precision_score, recall_score, roc_auc_score, average_precision_score
from sklearn.metrics import mean_squared_error

from matplotlib.backends.backend_pdf import PdfPages

from funcs import *

# ["SVM","LR","SGD","ABC","XGB","RF","KNN","MLP"]
# model = sys.argv[1]
feature_selection= sys.argv[1]
scale_method = sys.argv[2]
features = sys.argv[3] ## DEG, eRNA, circRNA
print("*******************************************")
########################
## train
train_data_file = "data/PPMIvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
train_data_df = pd.read_csv(train_data_file,sep="\t",index_col=0,header=0)

train_DEG = train_data_df.loc[:,train_data_df.columns.str.startswith("DEG")]
train_eRNA = train_data_df.loc[:,train_data_df.columns.str.startswith("eRNA")]
train_circRNA = train_data_df.loc[:,train_data_df.columns.str.startswith("circRNA")]
train_PRS = train_data_df.loc[:,train_data_df.columns.str.startswith("PRS")]
train_Clinical = train_data_df.loc[:,train_data_df.columns.str.startswith("Clin")]

PRS_features = ["PRS"]
clinicals_features = train_Clinical.columns.tolist() ## ["upsit_total_score","sex_int","age_at_baseline"]

## test
test_data_file = "data/PDBFvst_DEG_eRNA_circRNA_Clinical_PRS.tsv"
test_data_df = pd.read_csv(test_data_file,sep="\t",index_col=0,header=0)

dataset = "only_PPMI"
###################
for model in ["SVM","SVM_rbf","RG","LR","SGD","ABC","XGB","RF","KNN","MLP"]:
    print("=====================")
    output_dir = "../results/"+dataset+"/"+scale_method+"/"+features+"/" + model + "_" + feature_selection
    print(output_dir)
    # Create folder if the directory doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir,exist_ok=True)

    ######################
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
       
    kept_features = []
    # kept_features = DEG_features + eRNA_features + circRNA_features + PRS_features + clinicals_features
    if "DEG" in features: kept_features += list(DEG_features) 
    if "eRNA" in features: kept_features += list(eRNA_features) 
    if "circRNA" in features: kept_features += list(circRNA_features) 
    if "PRS" in features: kept_features += list(PRS_features) 
    if "Clinical" in features: kept_features += list(clinicals_features) 

    ######################
    train_X = train_data_df.loc[:,kept_features]
    train_X = train_X.loc[~train_X.isna().any(axis=1),:]
    train_y = train_data_df.loc[train_X.index.tolist(),["Y"]]
    
    print("Total feature number:", len(kept_features))

    #####################
    train_i_X, val_i_X, train_i_Y, val_i_Y = train_test_split(train_X, train_y, test_size=0.2, shuffle=True, random_state=42)

    train_index_ls =train_i_X.index.tolist()
    val_index_ls = val_i_X.index.tolist()

    print("Selected feature number:", len(kept_features))
    train_i_X =train_i_X.loc[:,kept_features]
    val_i_X =val_i_X.loc[:,kept_features]

    # if "Clin_sex_int" in kept_features:
    #     train_i_X.drop(['Clin_sex_int'], axis=1,inplace=True)
    #     val_i_X.drop(['Clin_sex_int'], axis=1,inplace=True)

    #########################
    ## log transform
    if scale_method.startswith("log"):
        print("log")
        train_i_X = np.log10(train_i_X+1)
        val_i_X = np.log10(val_i_X+1)

    ##################
    ## MinMax
    if scale_method.endswith("mm"):
        print("mm")
        train_i_X = minmax_scale(train_i_X, feature_range=(0, 1))
        val_i_X = minmax_scale(val_i_X, feature_range=(0, 1))
        # scaler = MinMaxScaler()
        # train_i_X = scaler.fit_transform(train_i_X)
        # val_i_X = scaler.transform(val_i_X)

    ## Standardization
    if scale_method.endswith("z"):
        print("z")
        train_i_X = stats.zscore(train_i_X)
        val_i_X = stats.zscore(val_i_X)
        # scaler = StandardScaler()
        # train_i_X = scaler.fit_transform(train_i_X)
        # val_i_X = scaler.transform(val_i_X)
    
    # if "Clin_sex_int" in kept_features:
    #     train_i_X = np.append(train_i_X, train_Clinical.loc[train_index_ls,["Clin_sex_int"]].values, axis=1)
    #     val_i_X=np.append(val_i_X, train_Clinical.loc[val_index_ls,["Clin_sex_int"]].values, axis=1)

    print("Final train size:", train_i_X.shape)
    print("Final val size:", val_i_X.shape)

    ###################################################
    ##### SVM
    if model == "SVM":
        classifier = SVC(kernel="linear",probability=True, C=1)
    
    if model == "SVM_rbf":
        classifier = SVC(probability=True,C=1)
    
    if model == "RG":
        classifier = LinearRegression()

    #### Logistic Regression
    if model == "LR":
        classifier = LogisticRegression(random_state=12, max_iter=300)

    #### SGDClassifier
    if model == "SGD":
        classifier = SGDClassifier(loss="log",max_iter=1000, tol=1e-3,random_state=12)

    #### AdaBoostClassifier
    if model == "ABC":
        classifier = AdaBoostClassifier(n_estimators=5, learning_rate=1,random_state=12)

    #### GradientBoostingClassifier
    if model == "XGB":
        classifier = GradientBoostingClassifier(n_estimators=5,min_samples_split=2,tol=1e-3,random_state=12)

    #### Random Forest
    if model == "RF":
        classifier = RandomForestClassifier(n_estimators=3,random_state=12)

    #### KNN
    if model == "KNN":
        classifier = KNeighborsClassifier(n_neighbors=5,algorithm='kd_tree')

    #### MLP
    if model == "MLP":
        classifier = MLPClassifier(hidden_layer_sizes=(32,),activation="relu",random_state=12)


    classifier.fit(train_i_X, train_i_Y)

    ##################
    ## val
    if model=="RG":
        pred_res = classifier.predict(val_i_X)
        val_y_pred_x = pred_res
        val_y_pred = val_y_pred_x > 0.5
    else:
        pred_res = classifier.predict_proba(val_i_X)
        val_y_pred_x = pred_res[:, 1]
        val_y_pred = val_y_pred_x > 0.5

    val_pred_df = pd.DataFrame(data=val_y_pred_x, index=val_index_ls, columns=["Pred_Y"])
    val_pred_df["True_Y"] = train_data_df.loc[val_index_ls, ["Y"]]
    val_pred_df.to_csv(output_dir + "/val_prediction_result.csv")
    sns.histplot(data=val_pred_df, x="Pred_Y", hue="True_Y", bins=50)
    sns.despine()
    plt.savefig(output_dir + "/val_prediction_distribution.pdf")
    plt.close()

    tn, fp, fn, tp = confusion_matrix(val_i_Y, val_y_pred).ravel()
    val_sensitivity=tp / (tp + fn)
    val_specificity=tn / (tn + fp)

    val_accuracy =accuracy_score(val_i_Y, val_y_pred)
    val_precision=precision_score(val_i_Y, val_y_pred)
    val_bal_acc=balanced_accuracy_score(val_i_Y, val_y_pred)
    val_recall = recall_score(val_i_Y, val_y_pred)

    ##AUC: cv-val
    val_fpr_roc, val_tpr_roc, val_thresholds_roc = roc_curve(val_i_Y, val_y_pred_x)
    val_roc = roc_auc_score(val_i_Y, val_y_pred_x)

    ##PR: cv-val
    val_precision_pr, val_recall_pr, val_thresholds_pr = precision_recall_curve(val_i_Y, val_y_pred_x)
    val_average_precision = average_precision_score(val_i_Y, val_y_pred_x)

    ######################
    ## train
    if model=="RG":
        train_pred_res = classifier.predict(train_i_X)
        train_y_pred_x = train_pred_res
        train_y_pred = train_y_pred_x > 0.5
    else:
        train_pred_res = classifier.predict_proba(train_i_X)
        train_y_pred_x = train_pred_res[:, 1]
        train_y_pred = train_y_pred_x > 0.5

    train_pred_df = pd.DataFrame(data=train_y_pred_x, index=train_index_ls, columns=["Pred_Y"])
    train_pred_df["True_Y"] = train_data_df.loc[train_index_ls, ["Y"]]
    train_pred_df.to_csv(output_dir + "/train_prediction_result.csv")
    sns.histplot(data=train_pred_df, x="Pred_Y", hue="True_Y", bins=50)
    sns.despine()
    plt.savefig(output_dir + "/train_prediction_distribution.pdf")
    plt.close()

    tn, fp, fn, tp = confusion_matrix(train_i_Y, train_y_pred).ravel()
    train_sensitivity=tp / (tp + fn)
    train_specificity=tn / (tn + fp)

    train_accuracy=accuracy_score(train_i_Y, train_y_pred)
    train_precision=precision_score(train_i_Y, train_y_pred)
    train_bal_acc=balanced_accuracy_score(train_i_Y, train_y_pred)
    recall_train = recall_score(train_i_Y, train_y_pred)

    ##AUC: cv-train
    train_fpr_roc, train_tpr_roc, train_thresholds_roc = roc_curve(train_i_Y, train_y_pred_x)
    train_roc = roc_auc_score(train_i_Y, train_y_pred_x)

    ##PR: cv-train
    train_precision_pr, train_recall_pr, train_thresholds_pr = precision_recall_curve(train_i_Y, train_y_pred_x)
    train_average_precision = average_precision_score(train_i_Y, train_y_pred_x)

    ###################
    # Plot
    with PdfPages(output_dir + "/Train_Val_AUROC_PR_"+model+".pdf") as pdf:
        fig1, ax_roc = plt.subplots(figsize=(8,8))
        fig2, ax_pr = plt.subplots(figsize=(8,8))
        ## Plots
        ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=1.2, color='grey', label='Chance', alpha=.8)
        ##1-1  CV-Train,AUC
        ax_roc.plot(train_fpr_roc, train_tpr_roc, color='#f8a289',label=r'Mean CV-Train ROC (AUC = %0.2f)' % (train_roc),lw=1.5, alpha=1)
        ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],  title="AUROC - "+model)

        ## 1-2  CV-val, AUC
        ax_roc.plot(val_fpr_roc, val_tpr_roc, color='#b365a6',label=r'Mean CV-Test ROC (AUC = %0.2f)' % (val_roc), lw=1.5, alpha=1)
        ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="AUROC - "+model)

        ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="AUROC - " + model)
        ax_roc.set_xlabel("False positive rate (1-specificuty)")
        ax_roc.set_ylabel("True positive rate (sensitive)")
        ax_roc.legend(loc="lower right", fontsize=8, frameon=False)

        ##################
        pos_rate = train_y.mean()
        ax_pr.plot([0, 1], [pos_rate, pos_rate], linestyle='--', lw=1.2, color='#f8a289', label='Baseline', alpha=.8)

        ## 2-1  CV-train,PR
        train_recall_pr = np.append(train_recall_pr, 1.0)
        train_precision_pr = np.append(train_precision_pr, 0.0)
        ax_pr.plot(train_recall_pr, train_precision_pr, color='#f8a289', label=r'Mean CV-Train P-R (AUC = %0.2f)' % (train_average_precision),lw=1.5, alpha=.8)
        ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],title="P-R - "+model)

        ## 2-2  CV-val,PR
        val_recall_pr = np.append(np.flip(val_recall_pr), 1.0)
        val_precision_pr = np.append(np.flip(val_precision_pr), 0.0)
        ax_pr.plot(val_recall_pr, val_precision_pr, color='#b365a6', label=r'Mean CV-Test P-R (AUC = %0.2f)' % (val_average_precision),lw=1.5, alpha=.8)
        ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="P-R - "+model)

        ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="P-R - " + model)
        ax_pr.set_xlabel("True positive rate (Recall)")
        ax_pr.set_ylabel("Positive predictive value (Precision)")
        ax_pr.legend(loc="lower left", fontsize=8, frameon=False)

        ######
        ax_roc.spines['right'].set_visible(False)
        ax_roc.spines['top'].set_visible(False)

        ax_pr.spines['right'].set_visible(False)
        ax_pr.spines['top'].set_visible(False)

        pdf.savefig(fig1)
        pdf.savefig(fig2)

        plt.close(fig1)
        plt.close(fig2)

    print(model)
    print("Acc:",train_accuracy,val_accuracy)
    print("Sen:",train_sensitivity,val_sensitivity)
    print("Spe:",train_specificity,val_specificity)
    print("Prec:",train_precision,val_precision)
    print("Bal_Acc:",train_bal_acc,val_bal_acc)
    print("ROC:",train_roc,val_roc)
    print("PR:",train_average_precision,val_average_precision)

    ###################################################################
    print("Run test...")
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

    #### Scale
    if scale_method.endswith("mm"):
        print("mm")
        test_X = minmax_scale(test_X, feature_range=(0, 1))

    if scale_method.endswith("z"):
        print("z")
        test_X = stats.zscore(test_X)

    # if "Clin_sex_int" in kept_features:
    #     test_X=np.append(test_X, test_data_df.loc[test_index_ls,["Clin_sex_int"]].values, axis=1)
    
    ## test
    if model=="RG":
        test_pred_res = classifier.predict(test_X)
        test_y_pred_x = test_pred_res
        test_y_pred = test_y_pred_x > 0.5

        coef = classifier.coef_[0,]
        coef_df = pd.DataFrame({"feature": kept_features,"coef":coef})
        coef_df.to_csv(output_dir+"/../rg_coef.txt",sep="\t")
    else:
        test_pred_res = classifier.predict_proba(test_X)
        test_y_pred_x = test_pred_res[:, 1]
        test_y_pred = test_y_pred_x > 0.5

    if model=="SVM":
        coef = classifier.coef_[0,]
        coef_df = pd.DataFrame({"feature": kept_features,"coef":coef})
        coef_df.to_csv(output_dir+"/../svm_coef.txt",sep="\t")

    test_pred_df = pd.DataFrame(data=test_y_pred_x,index=test_index_ls,columns=["Pred_Y"])
    test_pred_df["True_Y"]= test_data_df.loc[test_index_ls,["Y"]]

    test_pred_df.to_csv(output_dir + "/test_prediction_result.csv")
    sns.histplot(data=test_pred_df, x="Pred_Y", hue="True_Y",bins=50)
    sns.despine()
    plt.savefig(output_dir + "/test_prediction_distribution.pdf")
    plt.close()
    # test_mse = mean_squared_error(test_y, test_y_pred_x)

    tn, fp, fn, tp = confusion_matrix(test_y, test_y_pred).ravel()
    test_accuracy = accuracy_score(test_y, test_y_pred)
    test_bal_acc = balanced_accuracy_score(test_y, test_y_pred)
    test_precision = precision_score(test_y, test_y_pred)
    test_recall = recall_score(test_y, test_y_pred)
    test_sensitivity = tp / (tp+fn)
    test_specificity = tn / (tn+fp)

    ## AUC
    test_fpr_roc, test_tpr_roc, thresholds_roc = roc_curve(test_y, test_y_pred_x)
    test_roc = roc_auc_score(test_y, test_y_pred_x)

    ## PR
    test_precision_pr, test_recall_pr, test_thresholds_pr = precision_recall_curve(test_y, test_y_pred_x)
    test_average_precision = average_precision_score(test_y, test_y_pred_x)

    with PdfPages(output_dir + "/Test_AUROC_PR_"+model+".pdf") as pdf:
        fig1, ax_roc = plt.subplots(figsize=(8, 8))
        fig2, ax_pr = plt.subplots(figsize=(8, 8))
        ## Plots
        ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=1.2, color='grey', label='Chance', alpha=.8)

        ## 1-3
        ax_roc.plot(test_fpr_roc, test_tpr_roc, color='#7cae01', label=r'Testing ROC (AUC = %0.2f)' % (test_roc),
                    lw=1.5, alpha=1)
        ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="AUROC - " + model)
        ax_roc.set_xlabel("False positive rate (1-specificuty)")
        ax_roc.set_ylabel("True positive rate (sensitive)")
        ax_roc.legend(loc="lower right", fontsize=8, frameon=False)

        ##################
        pos_rate = test_y.mean()
        ax_pr.plot([0, 1], [pos_rate, pos_rate], linestyle='--', lw=1.2, color='#7cae01', label='Testing baseline',
                alpha=.8)

        ## 2-3 PDBF, PR
        test_recall_pr = np.append(np.flip(test_recall_pr), 1.0)
        test_precision_pr = np.append(np.flip(test_precision_pr), 0.0)
        ax_pr.plot(test_recall_pr, test_precision_pr, color='#7cae01',label=r'Testing P-R (AUC = %0.2f)' % (test_average_precision), lw=1.5, alpha=.8)
        ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="P-R - " + model)
        ax_pr.set_xlabel("True positive rate (Recall)")
        ax_pr.set_ylabel("Positive predictive value (Precision)")
        ax_pr.legend(loc="lower left", fontsize=8, frameon=False)

        # Hide the right and top spines
        ax_roc.spines['right'].set_visible(False)
        ax_roc.spines['top'].set_visible(False)

        ax_pr.spines['right'].set_visible(False)
        ax_pr.spines['top'].set_visible(False)

        pdf.savefig(fig1)
        pdf.savefig(fig2)

        plt.close(fig1)
        plt.close(fig2)

        # fig1.show()
        # fig2.show()

    print(model)
    print("Acc:",test_accuracy)
    print("Sen:",test_sensitivity)
    print("Spe:",test_specificity)
    print("Prec:",test_precision)
    print("Bal_Acc:",test_bal_acc)
    print("ROC:", test_roc)
    print("PR:",test_average_precision)

    f_perform = open(output_dir+"/../model_performances_"+feature_selection+".txt","a")
    f_perform.write("\t".join([model,str(train_accuracy),str(train_sensitivity),str(train_specificity),str(train_precision),str(train_bal_acc),str(train_roc),str(train_average_precision),
    str(val_accuracy),str(val_sensitivity),str(val_specificity),str(val_precision),str(val_bal_acc),str(val_roc),str(val_average_precision),
    str(test_accuracy),str(test_sensitivity),str(test_specificity),str(test_precision),str(test_bal_acc),str(test_roc),str(test_average_precision)])+"\n")
    f_perform.close()



