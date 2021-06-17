import statsmodels.api as sm
import numpy as np
import pandas as pd
import scipy.stats as stats

var_table = pd.read_excel("GSE39582_var_table.xlsx")
matrix_table = pd.read_excel("GSE39582_series_matrix_table.xlsx", index_col="ID_REF")
#matrix_table = matrix_table.iloc[:, :566]
gpl = pd.read_excel("GPL570-55999.xlsx", index_col="ID")

tumor_stage_na = np.where(var_table["Var_05_tnm.stage"].values == 'N/A', True, False)
tumor_stage_0 = np.where(var_table["Var_05_tnm.stage"].values == 0, True, False)
tumor_stage_1 = np.where(var_table["Var_05_tnm.stage"].values == 1, True, False)
t01_or = np.logical_or(tumor_stage_0, tumor_stage_1)

t01_or_X = matrix_table.iloc[:, t01_or]
t01_or_y = var_table["Var_05_tnm.stage"][t01_or]

print("# of stage0: {0} / # of stage1: {1}".format(tumor_stage_0.sum(), tumor_stage_1.sum()))
print("Total # of stage01_data: {}".format(t01_or.sum()))


tumor = {}
for i in range(t01_or_X.shape[0]):
    X, y = t01_or_X.iloc[i, :].values, list(t01_or_y)
    X = sm.add_constant(X)
    X[:,1] = 0 #remove NaN to 0
    model = sm.OLS(y, X)
    result = model.fit()
    tumor[t01_or_X.index[i]] = float(result.t_test([1, 0]).pvalue)

print(len(tumor))

tumor_sort = sorted(tumor.items(), key=lambda item: item[1])

tumor_p_value = np.array(list(tumor.values()), dtype=np.float)
print(np.sum(np.where(tumor_p_value < 0.001, 1, 0)))

significant_probe = list()
for dt in tumor_sort:
    if dt[1] < 0.001:
        significant_probe.append(dt[0])

gpl_drop = gpl.loc[significant_probe, "Gene Symbol"].values
unique_gene = set()
for i in gpl_drop:
    if type(i) == str:
        id_list = list(i.split(" /// "))
        unique_gene.add(id_list[0])
    else:
        pass
len(unique_gene)

top10 = list()
i = 0
while len(top10) < 10:
    probe_name = tumor_sort[i][0]
    i += 1
    candidate = gpl.loc[probe_name, "Gene Symbol"]
    if candidate in top10 or candidate is np.nan:
        continue
    else:
        top10.append(candidate)

print(top10)

#(used code) – t-test
benign = np.logical_or(np.where(var_table["Var_05_tnm.stage"].values == 0, True, False),
                       np.where(var_table["Var_05_tnm.stage"].values == 1, True, False))
benign = np.logical_or(benign, np.where(var_table["Var_05_tnm.stage"].values == 2, True, False))
benign_matrix = matrix_table.iloc[:, benign]
malignant = np.logical_or(np.where(var_table["Var_05_tnm.stage"].values == 3, True, False),
                          np.where(var_table["Var_05_tnm.stage"].values == 4, True, False))
malignant_matrix = matrix_table.iloc[:, malignant]

print("# of benign: {0} / # of malignant: {1}".format(benign.sum(), malignant.sum()))

malig_p_value_dict = {}
for i in range(benign_matrix.shape[0]):
    M_data, WT_data = benign_matrix.iloc[i, :].values, malignant_matrix.iloc[i, :].values
    static, p_value = stats.ttest_ind(M_data, WT_data)
    malig_p_value_dict[benign_matrix.index[i]] = float(p_value)

malig_p_value_dict_sort = sorted(malig_p_value_dict.items(), key=lambda item: item[1])

malig_p_value = np.array(list(malig_p_value_dict.values()), dtype=np.float)
print(np.sum(np.where(malig_p_value < 0.001, 1, 0)))


significant_probe = list()
for dt in malig_p_value_dict_sort:
    if dt[1] < 0.001:
        significant_probe.append(dt[0])

gpl_drop = gpl.loc[significant_probe, "Gene Symbol"].values
unique_gene = set()
for i in gpl_drop:
    if type(i) == str:
        id_list = list(i.split(" /// "))
        unique_gene.add(id_list[0])
    else:
        pass
len(unique_gene)

top10 = list()
i = 0
while len(top10) < 10:
    probe_name = malig_p_value_dict_sort[i][0]
    i += 1
    candidate = gpl.loc[probe_name, "Gene Symbol"]
    if candidate in top10 or candidate is np.nan:
        continue
    else:
        top10.append(candidate)

print(top10)



print("#5#######################################################")
entrez = pd.read_csv("c5.all.v7.4.entrez.gmt.csv", header=None, index_col=0)
entrez_p_value = {}

malig_entrez = list()
i = 0
while malig_p_value_dict_sort[i][1] < 0.001:
    probe_name = malig_p_value_dict_sort[i][0]
    i += 1
    candidate = gpl.loc[probe_name, "ENTREZ_GENE_ID"]
    if candidate is np.nan:
        continue
    else:
        if type(candidate) == str:
            candidate = candidate.split(" /// ")[0]
        malig_entrez.append(candidate)

malig_size = len(np.unique(malig_entrez))


gene_id = gpl["ENTREZ_GENE_ID"].values
unique_gene = set()
for i in gene_id:
    if type(i) == str:
        id_list = list(i.split(" /// "))
        unique_gene.add(id_list[0])
    elif type(i) == int:
        unique_gene.add(i)
    else:
        continue
entrez_size = len(unique_gene)
unique_gene = np.array(list(unique_gene))

for i in range(entrez.shape[0]):
    if entrez.index[i][:2] == "GO":
        info = entrez.iloc[i, :].values
        info = info[~np.isnan(info)].astype(np.int32)
        a = len(np.unique(np.intersect1d(malig_entrez, info)))
        b = len(np.unique(np.intersect1d(unique_gene, info)))
        #print(a, b)
        _, pvalue = stats.fisher_exact([[a, b], [malig_size - a, entrez_size - b]])
        entrez_p_value[entrez.index[i]] = pvalue


entrez_p_value_dict_sort = sorted(entrez_p_value.items(), key=lambda item: item[1])

entrez_p_value = np.array(list(entrez_p_value.values()), dtype=np.float)
print(np.sum(np.where(entrez_p_value < 0.001, 1, 0)))

for i in range(10):
    print(entrez_p_value_dict_sort[i])
'''
go_table = pd.read_csv("c5.all.v7.4.entrez.gmt.csv", header=None, index_col=0).T

go_table_pvalue = {}


#malig = tumor 3 4 stage
malig = list()
i =0
while malig_p_value_dict_sort[i][1] <0.001:
    probe_name = malig_p_value_dict_sort[i][0]
    i += 1
    candidate = gpl.loc[probe_name, "ENTREZ_GENE_ID"] #mapping gpl and GO_id
    if candidate is np.nan:
        continue #nan exception handling
    else:
        if type(candidate) == str:
            candidate = candidate.split(" /// ")[0] #remove useless genes after ///
        malig.append(candidate)

malig_size = len(np.unique(malig)) #unique gene size_malig go table
print("malig_unique gene size:{}",malig_size) #1633

gene_id = gpl["ENTREZ_GENE_ID"].values
unique_gene_go = set()
for i in gene_id:
    if type(i) == str:
        id_list = list(i.split(" /// "))
        unique_gene_go.add(id_list[0])
    elif type(i) == int:
        unique_gene_go.add(i)
    else:
        continue
go_table_size = len(np.unique(unique_gene_go)) #unique gene size in gpl #21180
unique_gene_go = np.array(list(unique_gene_go)) #len: 21712

for i in range(go_table.shape[0]):
    if go_table.index[i][:2] == "GO":
        info = go_table.iloc[i, :].values
        info = info[~np.isnan(info)].astype(np.int32)
        a = len(np.unique(np.intersect1d(malig, info)))
        b = len(np.unique(np.intersect1d(unique_gene, info)))
        print(a, b)
        _, pvalue = stats.fisher_exact([[a, b], [malig_size - a, go_table_size - b]])
        go_table_pvalue[go_table.index[i]] = pvalue
'''
