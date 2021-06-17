
# library
import numpy as np
import scipy as sp
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats as stats

from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

# read data
SAMPLE_ANNOT = pd.read_csv('sample_annotation.csv')
PS_ANNOT = pd.read_csv('probeset_annotation.csv')
EXPR = pd.read_csv('expression.csv',index_col=0).transpose()

# train / test separation
train_idx = np.array(SAMPLE_ANNOT.Dataset=='discovery')#[true:false]   #443
test_idx = np.array(SAMPLE_ANNOT.Dataset=='validation')#[false:true]   #123 
EXPR_TRAIN = EXPR[train_idx]#앞부분
SAMPLE_TRAIN = SAMPLE_ANNOT[train_idx]#뒷부분
EXPR_TEST = EXPR[test_idx]#앞부분
SAMPLE_TEST = SAMPLE_ANNOT[test_idx]#뒷부분

########################################################
# Inference of Significant Genes for Sex
########################################################

# two groups of samples, male vs. female
idx1 = np.array(SAMPLE_TRAIN.Sex=='Male')#[true,false,true,false...]
idx2 = np.array(SAMPLE_TRAIN.Sex=='Female')#[false,true,false,true...]
dat1 = EXPR_TRAIN[idx1]#male expression
dat2 = EXPR_TRAIN[idx2]#female expression


# t-test
pv_list = np.ones(EXPR_TRAIN.shape[1])
for i in range(len(pv_list)):#54675
    pv_list[i] = sp.stats.ttest_ind(dat1.iloc[:,i],dat2.iloc[:,i]).pvalue
'''
0.9857430473645499 
0.10481318440159079
0.07484513797465758
0.7868608418864361
0.7420495004806249
'''

# select top 100 significant probesets
pv_rank = pd.Series(pv_list).rank()
sig_ps_idx = np.array(pv_rank < 101)#[False False False ... False False False] 
sig_ps_list = EXPR_TRAIN.columns[sig_ps_idx].values #top 100 ranked probesets  ['1554447_at' '1556102_x_at' '1556677_at'....]
sig_gene_list = PS_ANNOT.Gene_Symbol[sig_ps_idx].values #top 100 mapped gene_symbols   ['JPX' 'LOC389906' nan 'TXLNG'...]
sig_gene_unique =  set(sig_gene_list) #unique gene symbols  {nan, 'LOC102724689', 'STS', 'KDM6A'....}


print(pd.Series(pv_list).rank())



########################################################
# Predidction for Sex
########################################################


# using only significant probesets
xtrain = EXPR_TRAIN[sig_ps_list]   #ID_REF     1554447_at  1556102_x_at  1556677_at  1557954_at  1558045_a_at  \
                                   # GSM971957    4.783981      5.266989 
    
ytrain = SAMPLE_TRAIN.Sex.values   #['Male' 'Female' 'Male' 'Female'....]

#1) feature selection
select = SelectKBest(score_func = chi2, k=4)
fit = select.fit(xtrain, ytrain)

#summarize scores
np.set_printoptions(precision=3)
arr = fit.scores_
print(arr)
features = fit.transform(xtrain)
print(features[0:5,:])


#printing-column-variable-names-after-feature-selection
idx = (-arr).argsort()[:2]
print (idx)

cols = EXPR_TRAIN.columns[idx]#feature selection에 변수 2개 뽑았는데, gene_symbol 이름 두 개 떠서 제대로 된건지 모르곘음. 


#2) prediction algorithm
# prediction using all significant probesets
f = LogisticRegression()
f.fit(xtrain,ytrain)
yhat_train = f.predict(xtrain)   #['Male' 'Female' 'Male' 'Female'...]
yhat_train_prob = f.predict_proba(xtrain)   #[[2.98681214e-04 9.99701319e-01]
                                             #[9.83788713e-01 1.62112875e-02]
                                                #[3.19379679e-03 9.96806203e-01]
pd.crosstab(yhat_train,ytrain)
cross_val_score(f,xtrain,ytrain) #array([0.96621622, 0.97297297, 0.97278912])  #Evaluate a score by cross-validation


print(pd.crosstab(yhat_train,ytrain))


#2) dimension reduction
# scaling & PCA(principal component analyis: in order to reduce dimenion, find another orthogonal axis (z1 z2 than x1 x2))
sca = StandardScaler() #Standardize features by removing the mean and scaling to unit variance
sca.fit(xtrain)
#transform(): Perform standardization by centering and scaling
xtrain_s = sca.transform(xtrain)  #[[ 1.06736805 -1.46122668 -0.75381845 ... -1.16661788  0.1356544 -0.73049722]...]
pca = PCA(n_components=10)
pca.fit(xtrain_s)
xtrain_t = pca.transform(xtrain_s)  #[[-7.55210571 -0.67877184 -0.25033022 ... -3.68675668  1.60155495  0.78502745]...]

# prediction after scaling and PCA
f = LogisticRegression()
f.fit(xtrain_t,ytrain)
yhat_train = f.predict(xtrain_t)   #['Male' 'Female' 'Male' 'Female' ...]
yhat_train_prob = f.predict_proba(xtrain_t)   #[[2.52270072e-03 9.97477299e-01] [9.97944916e-01 2.05508424e-03]...
pd.crosstab(yhat_train,ytrain)
cross_val_score(f,xtrain_t,ytrain)    #array([0.96621622, 0.93918919, 0.95918367])


print(pd.crosstab(yhat_train,ytrain))   #scaling&PCA 전후가 crosstab와 cross_val_score 결과값이 다르다는 것을 확인 할 수 있음. 


# parameter tuning using cross-validation
#GridSearchCV(): Exhaustive search over specified parameter values for an estimator.
f = GridSearchCV(LogisticRegression(),{'C': [0.0001,0.001,0.01,0.1,1,10,100]})
f.fit(xtrain_t,ytrain)
f_final = f.best_estimator_  #LogisticRegression()

# final test on test set
xtest = EXPR_TEST[sig_ps_list]   #top 100 probe sets_expressions for test set
ytest = SAMPLE_TEST.Sex.values   #['Female' 'Male' 'Male...]

xtest_s = sca.transform(xtest)  #[[-0.47004215  0.23569942 -0.60299454 ...  2.11320599 -1.02927915   0.44420555][...]...]
xtest_t = pca.transform(xtest_s)   #[[ 4.82380952 -2.64045732 -0.56515518 ...  1.09409492  0.21693289  -0.85316592] [...

yhat_test = f_final.predict(xtest_t)   #['Female' 'Male' 'Male' 'Female'...
yhat_test_prob = f_final.predict_proba(xtest_t)   #[[9.96827312e-01 3.17268809e-03] [2.03618731e-03 9.97963813e-01]...
pd.crosstab(yhat_test,ytest)

#####final prediction of gender#####
f_final.score(xtest_t,ytest)  #0.975609756097561

# roc curve
fpr,tpr,th = roc_curve(ytest,yhat_test_prob[:,1],pos_label='Male')
auc(fpr,tpr)
plt.plot(fpr,tpr)
plt.title("AUC: %.2f"%auc(fpr,tpr))
plt.show()


pd.crosstab(yhat_test,ytest)


f_final.score(xtest_t,ytest)


#####################mapping with C1 gene sets of MSigDB##############################
ENTREZ = pd.read_csv('c1.all.v7.4.entrez.gmt.csv',header=None, index_col=0)

var_table = pd.read_excel("GSE39582_var_table.xlsx")
matrix_table = pd.read_excel("GSE39582_series_matrix_table.xlsx", index_col="ID_REF")



gpl = pd.read_excel("GPL570-55999.xlsx",index_col="ID")


male = np.where(var_table["Var_03_Sex"].values == "Male", True, False)
female = np.where(var_table["Var_03_Sex"].values == "Female", True, False)
male_matrix = matrix_table.iloc[:, male]
female_matrix = matrix_table.iloc[:, female]
print("# of male: {0} / # of female: {1}".format(male.sum(), female.sum()))


p_value_set = {}
for i in range(female_matrix.shape[0]):
    male_data, female_data = male_matrix.iloc[i, :].values, female_matrix.iloc[i, :].values
    static, p_value = stats.ttest_ind(male_data, female_data)
    for col in EXPR_TRAIN[:0]: 
        if female_matrix.index[i]== col: #getting only train sets from predictor
            p_value_set[female_matrix.index[i]] = float(p_value)
print(len(p_value_set))


p_sort = sorted(p_value_set.items(), key=lambda item: item[1])
p_value = np.array(list(p_value_set.values()), dtype=np.float)
print(p_value)
print(np.sum(np.where(p_value < 0.001, 1, 0)))


significant_probe = list()
for dt in p_sort:
    if dt[1] < 0.001:
        significant_probe.append(dt[0])
print(significant_probe)  ##train set 에서의 significant probe 추출 
print(len(significant_probe))


entrez = list()
i = 0
while p_sort[i][1] < 0.001: #among significant (234개)
    probe_name = p_sort[i][0]
    i += 1
    candidate = gpl.loc[probe_name, "ENTREZ_GENE_ID"]
    if candidate is np.nan:
        continue
    else:
        if type(candidate) == str:
            candidate = candidate.split(" /// ")[0]
        entrez.append(candidate)

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
print(unique_gene)
print(len(unique_gene))


ENTREZ = pd.read_csv('c1.all.v7.4.entrez.gmt.csv',header=None, index_col=0)
entrez_p_value = {}



for i in range(ENTREZ.shape[0]):
    info = ENTREZ.iloc[i, :].values
    info = info[~np.isnan(info)].astype(np.int32)
    a = len(np.unique(np.intersect1d(malig_entrez, info)))
    b = len(np.unique(np.intersect1d(unique_gene, info)))
    print(a, b)
    _, pvalue = stats.fisher_exact([[a, b], [malig_size - a, entrez_size - b]])
    entrez_p_value[ENTREZ.index[i]] = pvalue


entrez_p_value_dict_sort = sorted(entrez_p_value.items(), key=lambda item: item[1])

entrez_p_value = np.array(list(entrez_p_value.values()), dtype=np.float)
print(np.sum(np.where(entrez_p_value < 0.001, 1, 0)))



print(pvalue)




