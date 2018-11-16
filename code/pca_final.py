# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 23:50:38 2018

@author: Ayush Kumar Pathak
"""
from matplotlib import pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
import csv

try:
    gene = open("gene")
    
    
    
    
        
    gendat = []
    for rows in csv.reader(gene):
        rowdat = []
        for columns in rows:
            rowdat.append(columns)
        gendat.append(rowdat)
    
    for i in range(len(gendat)):
        gendat[i] = gendat[i]
        for j in range(len(gendat[i])):
            gendat[i][j] = float(gendat[i][j])
        
    
    
    gene.close()
    
    
    
    X = gendat
    X_std = StandardScaler().fit_transform(X)
    





    mean_vec = np.mean(X_std, axis=0)
    cov_mat = (X_std - mean_vec).T.dot((X_std - mean_vec)) / (X_std.shape[0]-1)
    
    
    
    cov_mat = np.cov(X_std.T)
    
    eig_vals, eig_vecs = np.linalg.eig(cov_mat)
    
    
    
    eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]
    
    eig_pairs.sort(key=lambda x: x[0], reverse=True)
    
    
    
    matrix_w = np.hstack((eig_pairs[0][1].reshape(len(cov_mat),1),eig_pairs[1][1].reshape(len(cov_mat),1)))
    
    
    
    Y = X_std.dot(matrix_w)
    
    
    
    pc1=list()
    pc2=list()
    
    for i in range (len(Y)):
        pc1.append(Y[i][0])
        pc2.append(Y[i][1])
    
    
    
    with plt.style.context('seaborn-whitegrid'):
        plt.scatter(pc1,pc2)
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend(loc='lower center')
        plt.tight_layout()
        plt.show()

except:
    print("wrong format error occured")



