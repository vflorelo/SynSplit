#!/usr/bin/python3
import argparse
import json
import matplotlib.pyplot as plt
import numpy  as np
import pandas as pd
import random
import seaborn as sns
from codons                  import *
from scipy.stats             import ks_2samp
from sklearn.ensemble        import RandomForestClassifier
from sklearn.metrics         import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing   import StandardScaler
parser = argparse.ArgumentParser(prog="ml-splitter.py",
                                 description="A simple program to separate populations of genes based on compositional bias",
                                 epilog="""LICENSE:
                                 Distributed under the GNU General Public License (GPLv3). See LICENSE.md
                                 https://github.com/vflorelo/mlsplitter""")
parser.add_argument('-A',
                    '--train_list_a',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A plain text file containing a list of known genes belonging to group A')
parser.add_argument('-B',
                    '--train_list_b',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A plain text file containing a list of known genes belonging to group B')
parser.add_argument('-t',
                    '--test_list',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A plain text file containing a list of unknown genes')
parser.add_argument('-i',
                    '--info_file',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A tsv file containing the metrics to evaluate to separate the gene populations')
parser.add_argument('-p',
                    '--p_value',
                    type=float,
                    default=1e-6,
                    help='Adjust the p-value cutoff for including or excluding metrics from the two sets')
parser.add_argument('-x',
                    '--prefix',
                    type=str,
                    help='Name of the svg file with plots of the evaluated metrics')
args              = parser.parse_args()
train_list_a      = args.train_list_a
train_list_b      = args.train_list_b
test_list         = args.test_list
info_file         = args.info_file
p_value           = args.p_value
prefix            = args.prefix
train_list_a_fh   = train_list_a
train_list_a_data = train_list_a_fh.read()
train_list_a      = train_list_a_data.split("\n")
train_list_b_fh   = train_list_b
train_list_b_data = train_list_b_fh.read()
train_list_b      = train_list_b_data.split("\n")
test_list_fh      = test_list
test_list_data    = test_list_fh.read()
test_list         = test_list_data.split("\n")
full_list         = []
full_list_a       = train_list_a
full_list_b       = train_list_b
full_list         = full_list + full_list_a + full_list_b + test_list


if (len(train_list_a) > len(train_list_b) ):
    train_list_a = random.sample(train_list_a, k=len(train_list_b))
else:
    train_list_b = random.sample(train_list_b, k=len(train_list_a))


train_list = []
train_list = train_list + train_list_a
train_list = train_list + train_list_b
df = pd.read_csv(info_file,sep="\t")
def assign_set(gene_id,gene_list_a,gene_list_b):
    if (gene_id in gene_list_a):
        gene_set = "A"
    elif(gene_id in gene_list_b):
        gene_set = "B"
    else:
        gene_set = np.nan
    return gene_set
training_df = df.copy()
training_df = training_df[training_df["id"].isin(train_list)]
training_df["set"] = training_df["id"].apply(lambda x: assign_set(x,train_list_a,train_list_b))
training_df.dropna(subset=["set"],inplace=True)
var_set  = []
ks_res   = {}
var_list = ["%G+C"] + codon_list
for var in var_list:
    dist_a = list(training_df[training_df["set"]=="A"][var])
    dist_b = list(training_df[training_df["set"]=="B"][var])
    ks_test = ks_2samp(dist_a,dist_b)
    ks_stat = ks_test[0]
    ks_pval = ks_test[1]
    ks_dict = {}
    ks_dict["count"]  = len(train_list)
    ks_dict["stat"]   = ks_stat
    ks_dict["pvalue"] = ks_pval
    ks_res[var]       = ks_dict
    if ks_pval < p_value:
        var_set.append(var)
num_vars = len(var_set)
ks_res["pval_cutoff"] = p_value
ks_res["var_count"] = num_vars
ks_out_file = prefix + ".ks.json"
with open(ks_out_file, "w") as fp:
    json.dump(ks_res, fp)
X = training_df[var_set]
y = training_df["set"]
X_train, X_test, y_train, y_test = train_test_split(X, y , test_size=0.25, random_state=0)
scaler  = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test  = scaler.fit_transform(X_test)
model   = RandomForestClassifier()
model.fit(X_train,y_train)
predictions = model.predict(X_test)
acc_score   = accuracy_score(predictions,  y_test)
pre_a_score = precision_score(predictions, y_test, average="binary", pos_label="A")
rec_a_score = recall_score(predictions,    y_test, average="binary", pos_label="A")
pre_b_score = precision_score(predictions, y_test, average="binary", pos_label="B")
rec_b_score = recall_score(predictions,    y_test, average="binary", pos_label="B")
model_dict = {}
model_dict["name"]      = "Random Forest"
model_dict["count"]     = len(train_list)
model_dict["accuracy"]  = acc_score
model_dict["precision_A"] = pre_a_score
model_dict["recall_A"]    = rec_a_score
model_dict["precision_B"] = pre_a_score
model_dict["recall_B"]    = rec_a_score
model_out_file = prefix + ".model.json"
with open(model_out_file, "w") as fp:
    json.dump(model_dict, fp)
testing_X    = training_df[var_set]
testing_X_ss = scaler.fit_transform(testing_X)
testing_y    = model.predict(testing_X_ss)
training_df["pred_set"] = testing_y
tp_a_counts = len(training_df[(training_df["set"]=="A") & (training_df["pred_set"]=="A")]["id"].values.flatten().tolist())
fp_a_counts = len(training_df[(training_df["set"]=="A") & (training_df["pred_set"]=="B")]["id"].values.flatten().tolist())
tp_b_counts = len(training_df[(training_df["set"]=="B") & (training_df["pred_set"]=="B")]["id"].values.flatten().tolist())
fp_b_counts = len(training_df[(training_df["set"]=="B") & (training_df["pred_set"]=="A")]["id"].values.flatten().tolist())
testing_dict = {}
testing_dict["name"]  = "Random Forest"
testing_dict["TP_A"]  = tp_a_counts
testing_dict["FP_A"]  = fp_a_counts
testing_dict["TP_B"]  = tp_b_counts
testing_dict["FP_B"]  = fp_b_counts
model_eval_out_file = prefix + ".model.eval.json"
with open(model_eval_out_file, "w") as fp:
    json.dump(testing_dict, fp)
class_df   = df.copy()
class_df   = class_df[class_df["id"].isin(full_list)]
class_df["set"] = class_df["id"].apply(lambda x: assign_set(x,full_list_a,full_list_b))
class_X    = class_df[var_set]
class_X_ss = scaler.fit_transform(class_X)
class_y    = model.predict(class_X_ss)
class_df["pred_set"] = class_y
class_df["pred_set"] = "pred" + class_df["pred_set"].astype(str)
class_df["set"] = class_df["set"].fillna(value=class_df["pred_set"])
svg_file = prefix + ".svg"
var_mod  = num_vars % 6
if (var_mod != 0):
    num_rows   = int(((num_vars - var_mod) / 6 ) + 1)
else:
    num_rows = int(num_vars / 6)
fig_height = int(num_rows * 5)
print(str(num_rows))
print(str(fig_height))
fig, axes = plt.subplots(num_rows,6,figsize=(30,fig_height))
counter_1 = 0
counter_2 = 0
for var in var_set:
    sns.violinplot(
        ax           = axes[counter_1,counter_2],
        data         = class_df,
        x            = "set",
        y            = var,
        inner        = "quartile",
        order        = ["A","predA","predB","B"],
        density_norm = "count",
        hue          = "set",
        cut          = 0)
    axes[counter_1,counter_2].set_title(var)
    axes[counter_1,counter_2].set_ylim(-0.125,1.125)
    axes[counter_1,counter_2].set_xlabel('')
    axes[counter_1,counter_2].set_ylabel('')
    counter_2 += 1
    if counter_2%6 == 0:
        counter_2 = 0
        counter_1 += 1
plt.rcParams['svg.fonttype']    = 'none'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family']     = 'Arial'
plt.tight_layout()
plt.savefig(svg_file)
plt.close()
tsv_out_file = prefix + ".out.tsv"
class_df.drop(columns=["pred_set"],inplace=True)
class_df.to_csv(tsv_out_file,sep="\t",index=False)