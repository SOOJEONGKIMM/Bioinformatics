#!/usr/bin/env python
# coding: utf-8

# In[1]:


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




# In[35]:


tumor_stage = np.array(SAMPLE_TRAIN.TNM_Stage!='N/A')
#print(tumor_stage)
dat = EXPR_TRAIN[tumor_stage]


print(len(tumor_stage))#443
print(len(dat))#443
print(dat.iloc[:,1])
print(dat)


# In[23]:


###########discrete (use it for question 3)
'''
tumor_stage_na = np.array(SAMPLE_TRAIN.TNM_Stage=='N/A')#[true,false,true,false...]
tumor_stage_0 = np.array(SAMPLE_TRAIN.TNM_Stage==0)#[false,true,false,true...]
tumor_stage_1 = np.array(SAMPLE_TRAIN.TNM_Stage==1)
tumor_stage_2 = np.array(SAMPLE_TRAIN.TNM_Stage==2)
tumor_stage_3 = np.array(SAMPLE_TRAIN.TNM_Stage==3)
tumor_stage_4 = np.array(SAMPLE_TRAIN.TNM_Stage==4)
dat0 = EXPR_TRAIN[tumor_stage_0]
dat1 = EXPR_TRAIN[tumor_stage_1]
dat2 = EXPR_TRAIN[tumor_stage_2]
dat3 = EXPR_TRAIN[tumor_stage_3]
dat4 = EXPR_TRAIN[tumor_stage_4]

print(dat0)'''


# In[37]:


# using only significant probesets
xtrain = EXPR_TRAIN[tumor_stage]   #ID_REF     1554447_at  1556102_x_at  1556677_at  1557954_at  1558045_a_at  \
                                   # GSM971957    4.783981      5.266989 
    
ytrain = SAMPLE_TRAIN.TNM_Stage.values


# In[41]:


f = LogisticRegression()
f.fit(xtrain,ytrain)
yhat_train = f.predict(xtrain)   #[12344123210]
yhat_train_prob = f.predict_proba(xtrain)   #[[2.98681214e-04 9.99701319e-01]
                                             #[9.83788713e-01 1.62112875e-02]
                                                #[3.19379679e-03 9.96806203e-01]
pd.crosstab(yhat_train,ytrain)
cross_val_score(f,xtrain,ytrain) #array([0.472, 0.494, 0.348, 0.511, 0.511]) #Evaluate a score by cross-validation


# In[42]:


print(pd.crosstab(yhat_train,ytrain))


# In[49]:


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
yhat_train = f.predict(xtrain_t)   #['12310123 ...]
yhat_train_prob = f.predict_proba(xtrain_t)   #[[2.52270072e-03 9.97477299e-01] [9.97944916e-01 2.05508424e-03]...
pd.crosstab(yhat_train,ytrain)
cross_val_score(f,xtrain_t,ytrain)    #array([0.449, 0.551, 0.404, 0.466, 0.557])
print(f.intercept_)
print(f.coef_)


# In[ ]:


#Interpreting model coefficients

'''
[-1.961 -1.135  1.728  1.542 -0.174]
a "unit" increase in TNMstage0 is associated with a -1.961 "unit" decrease in expression.
a "unit" increase in TNMstage1 is associated with a -1.135 "unit" decrease in expression.
a "unit" increase in TNMstage2 is associated with a 1.728 "unit" increase in expression.
a "unit" increase in TNMstage3 is associated with a 1.542 "unit" increase in expression.
a "unit" increase in TNMstage4 is associated with a -0.174 "unit" decrease in expression.
'''


# In[48]:


# parameter tuning using cross-validation
#GridSearchCV(): Exhaustive search over specified parameter values for an estimator.
f = GridSearchCV(LogisticRegression(),{'C': [0.0001,0.001,0.01,0.1,1,10,100]})
f.fit(xtrain_t,ytrain)
f_final = f.best_estimator_  #LogisticRegression()

tumor_test = np.array(SAMPLE_TEST.TNM_Stage!='N/A')
dat_test = EXPR_TEST[tumor_test]

# final test on test set
xtest = EXPR_TEST[tumor_test]   #probe sets_expressions for test set
ytest = SAMPLE_TEST.TNM_Stage.values   #['123401123...]

xtest_s = sca.transform(xtest)  #[[-0.47004215  0.23569942 -0.60299454 ...  2.11320599 -1.02927915   0.44420555][...]...]
xtest_t = pca.transform(xtest_s)   #[[ 4.82380952 -2.64045732 -0.56515518 ...  1.09409492  0.21693289  -0.85316592] [...

yhat_test = f_final.predict(xtest_t)
yhat_test_prob = f_final.predict_proba(xtest_t)   #[[9.96827312e-01 3.17268809e-03] [2.03618731e-03 9.97963813e-01]...
pd.crosstab(yhat_test,ytest)

#####final prediction of gender#####
f_final.score(xtest_t,ytest)  #0.5365853658536586


# In[ ]:




