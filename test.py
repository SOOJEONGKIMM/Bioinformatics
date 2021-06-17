import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import permutation_test_score
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme()
from tqdm import tqdm
import seaborn as sns


#from mlxtend.evaluate import permutation_test
def clinical_var_preprocessing(arr):
    value = arr[1]
    if value == "N/A":
        value = np.nan
    return value

'''
var_table = pd.read_excel("GSE39582_series_matrix_variation.xlsx", index_col="!Sample_title")
var_table = var_table.T

var_num = 1
for i in range(var_table.shape[1]):
    if var_table.columns[i][:-3] == "!Sample_characteristics_ch1":
        split_idx = var_table.iloc[0, i].find(":")
        var_name = var_table.iloc[0, i][:split_idx]
        var_table.iloc[:, i] = var_table.iloc[:, i].str.split(": ").str[1]
        var_table = var_table.rename(columns={var_table.columns[i]: "Var_{:02d}_{}".format(var_num, var_name)})
        var_num += 1
var_table.to_excel("GSE39582_var_table.xlsx", header=True, index=False)
'''
var_table = pd.read_excel("GSE39582_var_table.xlsx")
var_table_valid = var_table[var_table["!Sample_submission_date"] != "May 12 2015"]
'''
print(var_table_valid.shape)
'''
matrix_table = pd.read_excel("GSE39582_series_matrix_table.xlsx", index_col="ID_REF")

matrix_table = matrix_table.iloc[:, :566]
'''
print(matrix_table.shape)
'''
gpl = pd.read_excel("GPL570-55999.xlsx", index_col="ID")
'''
gpl_drop = gpl.drop_duplicates()
gene_id = gpl_drop["Gene Symbol"].values
unique_gene = set()
for i in gene_id:
    if type(i) == str:
        id_list = list(i.split(" /// "))
        unique_gene.add(id_list[0])
    else:
        pass
print("unique",len(unique_gene))
'''
#########################################################################
#data preprocessing
male = np.where(var_table_valid["Var_03_Sex"].values == "Male", True, False)
male_matrix = matrix_table.iloc[:, male]
female = np.where(var_table_valid["Var_03_Sex"].values == "Female", True, False)
female_matrix = matrix_table.iloc[:, female]

print("# of male: {0} / # of female: {1}".format(male.sum(), female.sum()))

sex_p_value_set = {}
for i in range(female_matrix.shape[0]):
    male_data, female_data = male_matrix.iloc[i, :].values, female_matrix.iloc[i, :].values
    static, p_value = stats.ttest_ind(male_data, female_data)
    sex_p_value_set[female_matrix.index[i]] = float(p_value)
print(len(sex_p_value_set)) #54676

sex_p_value_set_sort = sorted(sex_p_value_set.items(), key=lambda item: item[1])
print("top 10 unique genes:")
for i in range(10):
    print(sex_p_value_set_sort[i])
print("# of unique genes:")
sex_p = np.array(list(sex_p_value_set.values()), dtype=float)
print(np.sum(np.where(sex_p < 0.001, 1, 0)))#222
sig_genes = pd.DataFrame(np.where(sex_p < 0.001,1,0))
print("sig genes:{}", sig_genes)

significant_probe = list()
for dt in range(200):
    if np.where(sex_p[dt] < 0.001):
        significant_probe.append(sex_p_set_sort[dt][:])
pd.options.display.float_format = '{:.5f}'.format

df =pd.DataFrame(significant_probe,columns=['gene_id','p_value'])
df=pd.DataFrame(df, columns=['gene_id','p_value'])
df = pd.pivot_table(df, values='p_value',columns='gene_id',index='gene_id')
df = df.fillna(0)

from sklearn.preprocessing import MinMaxScaler
m = MinMaxScaler(feature_range=(-1,1))
data_num_minmax = m.fit_transform(significant_probe)
data_num_minmax = pd.DataFrame(data_num_minmax ,
                               columns=significant_probe.columns.tolist())

#df =pd.DataFrame(significant_probe, columns=significant_probe.columns.tolist())
rounded_df = df.round(decimals=5)
#pdf = pd.DataFrame(sig_genes)
#sig_genes['foo'] = sig_genes['foo'].astype('bool')
print("sig genes bool:{}", sig_genes)
#df = matrix_table.loc[:significant_probe, :]
#df = analys.get_data(sig_genes)
df =df.dropna()
g = sns.clustermap(df)
#g = sns.clustermap(rounded_df.loc[:,1])
plt.savefig('hierarchical_clustered_corr_heatmap_with_Seaborn_clustermap_python_2nd_try.png',dpi=150)
#visuz.gene_exp.hmap(df=sig_genes, dim=(3,6), tickfont=(6,4))

'''
#sex_p_value_sig = {}
sex_p_value_sig = np.where(sex_p_value<0.001, 1, 0)
sex_p_value_sig_sort = {}
print(sex_p_value_sig)
sex_p_value_sig_sort = sorted(sex_p_value_sig.items(), key=lambda item: item[1])
print(sex_p_value_sig_sort)
sex_p_value_sig_two=[]
for i in range(200):
    sex_p_value_sig_two.append(sex_p_value_sig_sort[i])
print("200 significance")
print(sex_p_value_sig)
#df = analys.get_data(sex_p_value_set).data
'''


print("#3#########################################################################")

kras_mutation = np.unique(var_table_valid["Var_23_kras.mutation"].dropna().values)
print(kras_mutation)
# only two class: "M" / "WT"
#data preprocessing
kras_M = np.where(var_table_valid["Var_23_kras.mutation"].values == "M", True, False)
kras_M_matrix = matrix_table.iloc[:, kras_M]
kras_WT = np.where(var_table_valid["Var_23_kras.mutation"].values == "WT", True, False)
kras_WT_matrix = matrix_table.iloc[:, kras_WT]
print("# of M: {0} / # of WT: {1} / # of Nan: {2}".format(kras_M.sum(), kras_WT.sum(),
                                                          matrix_table.shape[1] - kras_M.sum() - kras_WT.sum()))

kras_p_set = {}
for i in range(female_matrix.shape[0]):
    male_data, female_data = male_matrix.iloc[i, :].values, female_matrix.iloc[i, :].values
    static, p_value = stats.ttest_ind(male_data, female_data)
    kras_p_set[female_matrix.index[i]] = float(p_value)
print(len(kras_p_set))

kras_p_sort = sorted(kras_p_set.items(), key=lambda item: item[1])
for i in range(10):
    print("iter", i)
    print(kras_p_sort[i])

kras_p_value = np.array(list(kras_p_set.values()), dtype=np.float)

print(np.sum(np.where(kras_p_value < 0.001, 1, 0)))
clf = SVC(kernel='linear', random_state=7)
cv = StratifiedKFold(2, shuffle=True, random_state=0)

score_iris, perm_scores_iris, pvalue_iris = permutation_test_score(
    clf, male_data[:256], female_data, scoring="accuracy", cv=cv, n_permutations=5)