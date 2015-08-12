#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def kernel(r):
    if r > 1.0e-8:
        return r * r * np.log(r)
    else:
        return 0.0

control_point = np.loadtxt('/media/hub/project-automorphing/arap_dfm/source_voxel.dat')
exp_point = np.loadtxt('/media/hub/project-automorphing/arap_dfm/source_emd.dat')

a = 0.0

n_num = len(control_point)
print 'v num', n_num

L = np.empty([n_num+4, n_num+4])
for i in xrange(n_num):
    for j in xrange(i+1, n_num):
        b = kernel(np.linalg.norm(control_point[i] - control_point[j]))
        L[i, j] = L[i, j] = b
        a += 2 * b

        
smooth = 0.5
a /= np.float(n_num * n_num)
for i in xrange(n_num):
    L[i ,i] = smooth * a * a

plt.imshow(L)
