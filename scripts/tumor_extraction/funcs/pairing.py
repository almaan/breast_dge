#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 10:04:53 2019

@author: alma
"""


import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.spatial import distance_matrix

import pandas as pd

import os
import os.path as osp



def make_rotation_matrix(theta):
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c,-s), (s, c)))
    return R

def align_spots(X,Y):
    
    def mini_fun(x0,x,y):
        theta = x0[0]
        v = x0[1::].reshape((1,2))
        A = make_rotation_matrix(theta)
        return np.sum(np.min(distance_matrix(np.dot(y,A)+v,x),axis=1))

   
    theta = np.random.random(1)
    v = np.random.random((1,2))
    
    x0 = np.zeros((3,))
    x0[0] = theta
    x0[1::] = v.reshape(-1,)
    
    res = minimize(mini_fun,x0, args=(X,Y),options=dict(maxiter = 10e18))
    
    theta_new = res['x'][0]
    v_new = res['x'][1::].reshape((1,2))
    
    R = make_rotation_matrix(theta_new)
    Y_new = np.dot(Y,R) + v_new

    return Y_new

def pair_multiple_sections(*sections,colname = 'tumor_id',xlab = "xcoord",ylab="ycoord"):
    tsections = [s.copy() for s in sections]
    sections = list(sections)
    
    n_sections = len(sections)
    min_foci = int(np.argmin([np.unique(s['tumor_id']).shape[0] for s in sections]))
    nidx = np.arange(n_sections)
    nidx = np.delete(nidx,np.where(nidx == min_foci))
    
    for secid in nidx:
        tsections[secid][[xlab,ylab]] = align_spots(sections[min_foci][[xlab,ylab]].values,
                                                 tsections[secid][[xlab,ylab]].values) 
        
        plist = pair_foci(tsections[min_foci][[xlab,ylab]].values,
                          tsections[secid][[xlab,ylab]].values,
                          tsections[min_foci][colname],
                          tsections[secid][colname])
        
        
        new_labels = -1*np.ones(tsections[secid][colname].shape[0])
        
        for i in range(plist.shape[0]):
            new_labels[sections[secid][colname] == plist[i,1]] = plist[i,0]

        sections[secid][colname] = new_labels
    
    return sections
        

def pair_foci(ocrd1,ocrd2,olab1,olab2):
    
    switch = [False if np.unique(olab1[olab1 != -1]).shape[0] >=
              np.unique(olab2[olab2 != -1]).shape[0] else True][0]
    
    
    if switch:
        crd1,crd2 = ocrd2,ocrd1
        lab1,lab2 = olab2,olab1
    else:
        crd1,crd2 = ocrd1,ocrd2
        lab1,lab2 = olab1,olab2
    
    pair_arr = np.zeros((np.unique(lab1[lab1 != -1]).shape[0],2))
    score = np.inf
    
    for (k,x) in enumerate(np.unique(lab1[lab1 != -1])):
        for y in np.unique(lab2[lab2 != -1]):
            min_dists = np.min(distance_matrix(crd1[lab1 == x,:],crd2[lab2 == y,:]),
                                axis = np.argmax([sum(lab1==x),sum(lab2 == y)]))
            
            tmp = np.sum(min_dists)
            
            if tmp < score:
                pair = [x,y]
                score = tmp
                
        pair_arr[k,:] = pair
        score = np.inf
    
    if switch:
        pair_arr = pair_arr[:,[1,0]]
    
    return pair_arr
