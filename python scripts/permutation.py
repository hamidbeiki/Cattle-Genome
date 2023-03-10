#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 11:18:40 2020

@author: beiki
"""

""" please see http://people.duke.edu/~ccc14/sta-663-2016/15B_ResamplingAndSimulation.html """
import numpy as np
import matplotlib.pyplot as plt

x = np.r_[np.random.exponential(size=200),
          np.random.normal(0, 1, size=100)]
y = np.r_[np.random.exponential(size=250),
          np.random.normal(0, 1, size=50)]


n1, n2 = map(len, (x, y))
reps = 10000

data = np.r_[x, y]
ps = np.array([np.random.permutation(n1+n2) for i in range(reps)])

xp = data[ps[:, :n1]]
yp = data[ps[:, n1:]]
samples = np.percentile(xp, 7, axis=1) - np.percentile(yp, 7, axis=1)

#perm_res=[]
#for i in range(len(xp)):
#    perm_res.append(stats.ttest_ind(np.array(xp[i]),np.array(yp[i]))[0][0])

plt.hist(samples, 25, histtype='step', color='red')
test_stat = np.percentile(x, 7) - np.percentile(y, 7)
plt.axvline(np.percentile(samples, 2.5), linestyle='--')
plt.axvline(np.percentile(samples, 97.5), linestyle='--')
print("p-value =", 2*np.sum(samples >= np.abs(test_stat))/reps)