#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:56:21 2019

@author: alma
"""

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import distance_matrix

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

