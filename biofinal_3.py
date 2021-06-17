#!/usr/bin/env python
# coding: utf-8

# In[3]:


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

var_table = pd.read_excel("GSE39582_var_table.xlsx")
matrix_table = pd.read_excel("GSE39582_series_matrix_table.xlsx", index_col="ID_REF")
gpl = pd.read_excel("GPL570-55999.xlsx", index_col="ID")

# train / test separation
train_idx = np.array(SAMPLE_ANNOT.Dataset=='discovery')#[true:false]   #443
test_idx = np.array(SAMPLE_ANNOT.Dataset=='validation')#[false:true]   #123 
EXPR_TRAIN = EXPR[train_idx]#앞부분
SAMPLE_TRAIN = SAMPLE_ANNOT[train_idx]#뒷부분
EXPR_TEST = EXPR[test_idx]#앞부분
SAMPLE_TEST = SAMPLE_ANNOT[test_idx]#뒷부분


tumor_stage_na = np.array(SAMPLE_TRAIN.TNM_Stage=='N/A')#[true,false,true,false...]
benign = np.logical_or(np.array(SAMPLE_TRAIN.TNM_Stage==0), np.array(SAMPLE_TRAIN.TNM_Stage==1))
benign = np.logical_or(benign, np.array(SAMPLE_TRAIN.TNM_Stage==2))
benign_matrix = EXPR_TRAIN[benign]

malignant = np.logical_or(np.array(SAMPLE_TRAIN.TNM_Stage==3), np.array(SAMPLE_TRAIN.TNM_Stage==4))
malignant_matrix = EXPR_TRAIN[malignant]


# In[160]:


benign_test= np.logical_or(np.array(SAMPLE_TEST.TNM_Stage==0), np.array(SAMPLE_TEST.TNM_Stage==1))
benign_test = np.logical_or(benign_test, np.array(SAMPLE_TEST.TNM_Stage==2))
benign_matrix_test = EXPR_TEST[benign_test]

malignant_test = np.logical_or(np.array(SAMPLE_TEST.TNM_Stage==3), np.array(SAMPLE_TEST.TNM_Stage==4))
malignant_matrix_test = EXPR_TEST[malignant_test]


# In[ ]:





# In[192]:


benign_matrix_test


# In[168]:


print(malignant_test)


# In[22]:


# t-test
pv_list = np.ones(EXPR_TRAIN.shape[1])
for i in range(len(pv_list)):#54675
    pv_list[i] = sp.stats.ttest_ind(benign_matrix.iloc[:,i],malignant_matrix.iloc[:,i]).pvalue


# In[9]:


# select top 100 significant probesets
pv_rank = pd.Series(pv_list).rank()
sig_ps_idx = np.array(pv_rank < 101)#[False False False ... False False False] 
sig_ps_list = EXPR_TRAIN.columns[sig_ps_idx].values #top 100 ranked probesets  ['1554447_at' '1556102_x_at' '1556677_at'....]
sig_gene_list = PS_ANNOT.Gene_Symbol[sig_ps_idx].values #top 100 mapped gene_symbols   ['JPX' 'LOC389906' nan 'TXLNG'...]
sig_gene_unique =  set(sig_gene_list) #unique gene symbols  {nan, 'LOC102724689', 'STS', 'KDM6A'....}


# In[69]:


print(SAMPLE_TRAIN.TNM_Stage[benign].index)


# In[146]:


##making a dict with binary tumor stage 

binary_tumor = {}
for i in range(SAMPLE_TRAIN.TNM_Stage.shape[0]):
    for j in range(len(SAMPLE_TRAIN.TNM_Stage[benign].index)):
        if i == SAMPLE_TRAIN.TNM_Stage[benign].index[j]:
            binary_tumor[i]="benign"
    for j in range(len(SAMPLE_TRAIN.TNM_Stage[malignant].index)):
        if i == SAMPLE_TRAIN.TNM_Stage[malignant].index[j]:
            binary_tumor[i]="malignant"
      
print(binary_tumor)


# In[147]:


print(pd.Series(binary_tumor).values)


# In[148]:



########################################################
# Predidction for TMNtumor binary 2 stage
########################################################


# using only significant probesets
xtrain = EXPR_TRAIN[sig_ps_list]   #ID_REF     1554447_at  1556102_x_at  1556677_at  1557954_at  1558045_a_at  \
                                   # GSM971957    4.783981      5.266989 
    
ytrain = pd.Series(binary_tumor).values   #[4. 4. 2. 1. 4. 3. 2. 3. 2. 3. 2. 2. 1. 1. 4. 1. 2. 4. 2. 

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


# In[149]:


print(cols)


# In[150]:


#2) prediction algorithm
# prediction using all significant probesets
f = LogisticRegression()
f.fit(xtrain,ytrain)
yhat_train = f.predict(xtrain)   #[2. 3. 2. 2. 3. 3. 3. 2. 2. 3. 2. 4. 2. 2. 2. 2. 2. 2. 2. 2. 3. 2. 2. 3....]]
yhat_train_prob = f.predict_proba(xtrain)   #[[2.98681214e-04 9.99701319e-01]
                                             #[9.83788713e-01 1.62112875e-02]
                                                #[3.19379679e-03 9.96806203e-01]
pd.crosstab(yhat_train,ytrain)
cross_val_score(f,xtrain,ytrain) #array([0.506, 0.663, 0.36 , 0.466, 0.602]) #Evaluate a score by cross-validation


# In[151]:


print(ytrain )


# In[152]:


print(pd.crosstab(yhat_train,ytrain))


# In[153]:


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
cross_val_score(f,xtrain_t,ytrain)    #array([0.517, 0.73 , 0.292, 0.5  , 0.58 ])


# In[154]:


print(pd.crosstab(yhat_train,ytrain))


# In[184]:


print(SAMPLE_TEST.TNM_Stage[benign_test].index)


# In[187]:


binary_tumor_test = {}

for i in range(len(SAMPLE_TEST.TNM_Stage)):
    for j in range(len(SAMPLE_TEST.TNM_Stage[benign_test].index)):
        if i+443 == SAMPLE_TEST.TNM_Stage[benign_test].index[j]:
            binary_tumor_test[i]="benign"
    for j in range(len(SAMPLE_TEST.TNM_Stage[malignant_test].index)):
        if i+443 == SAMPLE_TEST.TNM_Stage[malignant_test].index[j]:
            binary_tumor_test[i]="malignant"
      
print(binary_tumor_test)


# In[198]:



# parameter tuning using cross-validation
#GridSearchCV(): Exhaustive search over specified parameter values for an estimator.
f = GridSearchCV(LogisticRegression(),{'C': [0.0001,0.001,0.01,0.1,1,10,100]})
f.fit(xtrain_t,ytrain)
f_final = f.best_estimator_  #LogisticRegression()

# final test on test set
xtest = EXPR_TEST[sig_ps_list]   #top 100 probe sets_expressions for test set
ytest = pd.Series(binary_tumor_test).values   #['Female' 'Male' 'Male...]

xtest_s = sca.transform(xtest)  #[[-0.47004215  0.23569942 -0.60299454 ...  2.11320599 -1.02927915   0.44420555][...]...]
xtest_t = pca.transform(xtest_s)   #[[ 4.82380952 -2.64045732 -0.56515518 ...  1.09409492  0.21693289  -0.85316592] [...

yhat_test = f_final.predict(xtest_t)   #['Female' 'Male' 'Male' 'Female'...
yhat_test_prob = f_final.predict_proba(xtest_t)   #[[9.968271312e-0 3.17268809e-03] [2.03618731e-03 9.97963813e-01]...
pd.crosstab(yhat_test,ytest)

#####final prediction of binary tmn#####
print(f_final.score(xtest_t,ytest)) #0.975609756097561

# roc curve
fpr,tpr,th = roc_curve(ytest,yhat_test_prob[:,1],pos_label='benign')
auc(fpr,tpr)
plt.plot(fpr,tpr)
plt.title("AUC: %.2f"%auc(fpr,tpr))
plt.show()


# In[199]:


pd.crosstab(yhat_test,ytest)


# In[ ]:




