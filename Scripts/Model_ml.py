
import os
import pandas as pd
from sklearn.model_selection import  train_test_split

from sklearn.model_selection import StratifiedKFold
from sklearn.base import clone

import numpy as np
from numpy import interp
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import *
from sklearn.metrics import confusion_matrix


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


##chapter 2
def getBPM(HCC_path, Control_path,  X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_combine_BPM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)


        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)

    for i in HCC:
        j = i.replace("_combine_BPM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)

        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)


    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getsBPM(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_sBPM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_sBPM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)


    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getEDM(HCC_path, Control_path,  X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_combine_EDM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_combine_EDM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)

    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getOJM(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_JM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_JM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)

    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)

def getOLM(HCC_path, Control_path,  X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_OLM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_OLM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)
    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getCNV(HCC_path, Control_path,  X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace(".cnv.txt", '')

        if j in X_train:
            tmp_data = pd.read_table(Control_path + i, sep='\t', header=None, names=['gene', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_table(Control_path + i, sep='\t', header=None, names=['gene', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace(".cnv.txt", '')

        if j in X_train:
            tmp_data = pd.read_table(HCC_path + i, sep='\t', header=None, names=['gene', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)

        if j in X_test:
            tmp_data = pd.read_table(HCC_path + i, sep='\t', header=None, names=['gene', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)


    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getOBR(HCC_path, Control_path,  X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_overlap_bases_ratio.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['base', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['base', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_overlap_bases_ratio.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['base', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['base', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            # motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)

    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)

##chapter 4

Control = os.listdir('YOUR/PATH/')
HCC = os.listdir('YOUR/PATH/')
all_name = []
n_label = []


for i in Control:
    final_name = os.path.splitext(os.path.basename(i))[0]
    final_name = final_name.replace("_overlap_bases_ratio", "")
    n_label.append(0)
    all_name.append(final_name)


for i in HCC:
    final_name = os.path.splitext(os.path.basename(i))[0]
    final_name = final_name.replace("_overlap_bases_ratio", "")
    n_label.append(1)
    all_name.append(final_name)

X_train, X_test, y_train, y_test = train_test_split(all_name, n_label, stratify=n_label,test_size=0.3, random_state=2)


##chapter 5
#control
BPM_Control_path = 'YOUR/PATH/'
sBPM_Control_path = 'YOUR/PATH/'
EDM_Control_path = 'YOUR/PATH/'
OLM_Control_path = 'YOUR/PATH/'
OJM_Control_path = 'YOUR/PATH/'
CNV_onco_Control_path = 'YOUR/PATH/'
CNV_im_Control_path ='YOUR/PATH/'
OBR_Control_path ='YOUR/PATH/'


#HCC
BPM_HCC_path = 'YOUR/PATH/'
EDM_HCC_path = 'YOUR/PATH/'
sBPM_HCC_path = 'YOUR/PATH/'
OLM_HCC_path = 'YOUR/PATH/'
CNV_onco_HCC_path ='YOUR/PATH/'
CNV_im_HCC_path ='YOUR/PATH/'
OJM_HCC_path ='YOUR/PATH/'
OBR_HCC_path ='YOUR/PATH/'

BPM_X_train,BPM_y_train,BPM_X_test,BPM_y_test = getBPM(BPM_HCC_path,BPM_Control_path,X_train,X_test)
EDM_X_train, EDM_y_train, EDM_X_test, EDM_y_test = getEDM(EDM_HCC_path, EDM_Control_path, X_train, X_test)
sBPM_X_train, sBPM_y_train, sBPM_X_test, sBPM_y_test = getsBPM(sBPM_HCC_path, sBPM_Control_path, X_train, X_test)
OLM_X_train, OLM_y_train, OLM_X_test, OLM_y_test = getOLM(OLM_HCC_path, OLM_Control_path, X_train, X_test)
CNV_onco_X_train, CNV_onco_y_train, CNV_onco_X_test, CNV_onco_y_test = getCNV(CNV_onco_HCC_path, CNV_onco_Control_path, X_train, X_test)
CNV_im_X_train, CNV_im_y_train, CNV_im_X_test, CNV_im_y_test = getCNV(CNV_im_HCC_path, CNV_im_Control_path, X_train, X_test)
OJM_X_train, OJM_y_train, OJM_X_test, OJM_y_test = getOJM(OJM_HCC_path, OJM_Control_path, X_train, X_test)
OBR_X_train, OBR_y_train, OBR_X_test, OBR_y_test = getOBR(OBR_HCC_path, OBR_Control_path, X_train, X_test)


clf1 = RandomForestClassifier(random_state=22)

clf_list = {
    'RF': clf1
}

BPM = RandomForestClassifier(random_state=32)
EDM = RandomForestClassifier(random_state=2)
sBPM = RandomForestClassifier(random_state=32)
OLM = RandomForestClassifier(random_state=32)
CNV_onco = RandomForestClassifier(random_state=32)
CNV_im = RandomForestClassifier(random_state=32)
OBR = RandomForestClassifier(random_state=32)
OJM = RandomForestClassifier(random_state=32)
##chapter 7
n_folds = 5
skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=3)
df_eval = pd.DataFrame(columns=['Accuracy', 'Precision', 'Recall', 'F1_score', 'auc'])


def training_BaseModel(X_train, y_train, disease_test, skf, clf, blend_train, j):
    blend_disease_j = np.zeros((disease_test.shape[0], 5))
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    for i, (train_index, cv_index) in enumerate(skf.split(X_train, y_train)):
        print('Fold [%s]' % (i))
        tr_X = X_train[train_index]
        tr_y = y_train[train_index]
        cv_X = X_train[cv_index]
        cv_y = y_train[cv_index]
        clf.fit(tr_X, tr_y)
        blend_train[[cv_index], j] = clf.predict(cv_X)

        blend_disease_j[:, i] = clf.predict(disease_test)

    return blend_train, blend_disease_j


tprs = []
r = []
mean_fpr1 = np.linspace(0, 1, 100)


blend_train = np.zeros((BPM_X_train.shape[0], 8))
blend_test = np.zeros((BPM_X_test.shape[0], 8))
a = EDM_X_test.shape[0]
j = 0
for clf, label in zip([BPM, EDM, sBPM, OLM, CNV_onco,CNV_im, OBR,OJM], ['BPM', 'EDM', 'sBPM', 'OLM', 'CNV_onco','CNV_im', 'OBR','OJM']):
    print("Training BaseModel [%s]" % (label))

    if label == 'BPM':
        blend_train, blend_test_j = training_BaseModel(BPM_X_train, BPM_y_train, BPM_X_test, skf, clf, blend_train, j)

    if label == 'EDM':
        blend_train, blend_test_j = training_BaseModel(EDM_X_train, EDM_y_train, EDM_X_test, skf, clf, blend_train,
                                                       j)
    if label == 'sBPM':
        blend_train, blend_test_j = training_BaseModel(sBPM_X_train, sBPM_y_train, sBPM_X_test, skf, clf, blend_train,
                                                       j)
    if label == 'OLM':
        blend_train, blend_test_j = training_BaseModel(OLM_X_train, OLM_y_train, OLM_X_test, skf, clf, blend_train,
                                                       j)
    if label == 'CNV_onco':
        blend_train, blend_test_j = training_BaseModel(CNV_onco_X_train, CNV_onco_y_train, CNV_onco_X_test, skf, clf, blend_train,
                                                       j)
    if label == 'CNV_im':
        blend_train, blend_test_j = training_BaseModel(CNV_im_X_train, CNV_im_y_train, CNV_im_X_test, skf, clf, blend_train,
                                                       j)
    if label == 'OBR':
        blend_train, blend_test_j = training_BaseModel(OBR_X_train, OBR_y_train, OBR_X_test, skf, clf, blend_train,
                                                       j)
    if label == 'OJM':
        blend_train, blend_test_j = training_BaseModel(OJM_X_train, OJM_y_train, OJM_X_test, skf, clf, blend_train,
                                                       j)

    blend_test[:, j] = blend_test_j.mean(1)

    j = j + 1

datasets = {
    'BPM': (BPM_X_train, BPM_y_train, BPM_X_test, BPM_y_test),
    'EDM': (EDM_X_train, EDM_y_train, EDM_X_test, EDM_y_test),
    'CNV_onco': (CNV_onco_X_train, CNV_onco_y_train, CNV_onco_X_test, CNV_onco_y_test),
    'CNV_im': (CNV_im_X_train, CNV_im_y_train, CNV_im_X_test, CNV_im_y_test),
    'sBPM': (sBPM_X_train, sBPM_y_train, sBPM_X_test, sBPM_y_test),
    'OLM': (OLM_X_train, OLM_y_train, OLM_X_test, OLM_y_test),
    'OBR': (OBR_X_train, OBR_y_train, OBR_X_test, OBR_y_test),
    'OJM': (OJM_X_train, OJM_y_train, OJM_X_test, OJM_y_test),
    'blend': (blend_train, BPM_y_train, blend_test, BPM_y_test),
}


results = []
results_disease = []
results_disease_matrix = []
roc = []
feature_importances_dict = {}
for clf_name, clf in clf_list.items():
    for dataset_name, (X_train, y_train, X_test, y_test) in datasets.items():
        print(f", Dataset: {dataset_name}")

        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

        accuracies, recalls, specificities, f1s, aucs, clf_mean_fprs, clf_mean_tprs = [], [], [], [], [], [],[]

        for fold, (train_idx, val_idx) in enumerate(skf.split(X_train, y_train), 1):

            fold_clf = clone(clf)


            X_fold_train, y_fold_train = X_train[train_idx], y_train[train_idx]

            # 训练与预测
            fold_clf.fit(X_fold_train, y_fold_train)
            y_pred = fold_clf.predict(X_test)

            # 计算指标
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

            accuracies.append(accuracy_score(y_test, y_pred))
            recalls.append(recall_score(y_test, y_pred))
            specificities.append(tn / (tn + fp) if (tn + fp) != 0 else 0)
            f1s.append(f1_score(y_test, y_pred))

            clf_prob = fold_clf.predict_proba(X_test)[:, 1]
            clf_fpr, clf_tpr, clf_thresholds = roc_curve(y_test, clf_prob)
            clf_mean_fpr1 = np.linspace(0, 1, 100)
            clf_tprs = []
            clf_tprs.append(interp(clf_mean_fpr1, clf_fpr, clf_tpr))
            clf_tprs[-1][0] = 0.0
            clf_mean_tpr1 = np.mean(clf_tprs, axis=0)
            clf_mean_tpr1[-1] = 1.0
            clf_mean_auc = auc(clf_mean_fpr1, clf_mean_tpr1)
            clf_tprs_mean = np.mean(clf_tprs, axis=0)
            clf_roc_auc = auc(clf_mean_fpr1, clf_mean_tpr1)
            clf_mean_fprs.append(clf_mean_fpr1)
            clf_mean_tprs.append(clf_mean_tpr1)
            aucs.append(clf_roc_auc)

        accuracy = np.mean(accuracies)
        recall = np.mean(recalls)
        specificity = np.mean(specificities)
        f1 = np.mean(f1s)
        Auc = np.mean(aucs)
        fpr = np.mean(clf_mean_fprs)
        tpr = np.mean(clf_mean_tprs)

        roc.append([clf_name, fpr, tpr, Auc])
        results_disease.append(
            [clf_name, dataset_name, accuracy, specificity, recall,
             f1, Auc])

# 创建数据框并设置列名
df_disease = pd.DataFrame(results_disease,
                          columns=['Classifier', 'Dataset', 'Accuracy', 'Specificity', 'Recall', 'F1', 'AUC'])
file_path = 'YOUR/PATH/'
df_disease.to_excel(file_path,index=False)
print(df_disease)

