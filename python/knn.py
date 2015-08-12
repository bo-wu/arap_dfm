#! /usr/bin/env python
import numpy as np
from sklearn.neighbors import KDTree

class DistNeighbor(object):
    def __init__(self):
        source_dist = np.loadtxt('../target_dist.dat')
        self.target_dist = np.loadtxt('../source_dist.dat')
        self.kdtree = KDTree(source_dist, leaf_size=100, metric='euclidean')

    def find_nearest(self):
        source_voxel = np.loadtxt('../target_voxel.dat')
        source_index = self.kdtree.query(self.target_dist, k=4, return_distance=False)
        #target_corresp = (source_voxel[source_index[:,0]] + source_voxel[source_index[:,1]] + source_voxel[source_index[:,2]] + source_voxel[source_index[:,3]]) / 4.
        target_corresp = source_voxel[source_index[:,0]]# + source_voxel[source_index[:,1]] + source_voxel[source_index[:,2]] + source_voxel[source_index[:,3]]) / 4.
        np.savetxt('source_corresp.dat', target_corresp, fmt='%.6f')


if __name__ == '__main__':
    dn = DistNeighbor()
    dn.find_nearest()

