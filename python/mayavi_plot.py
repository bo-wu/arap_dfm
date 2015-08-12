#! /usr/bin/env python
import sys
import numpy as np
import mayavi.mlab as mmlab


"""
#source_dist = np.loadtxt("../source_dist.dat")
source_dist = np.loadtxt("../target_dist.dat")
sequence = []
#name = '../source_intermedium'
name = '../target_intermedium'
for i in xrange(11):
    s_name = name + str(i) + '.dat'
    source_voxel = np.loadtxt(s_name)
    mmlab.points3d(source_voxel[:,0]+(i-5)*0.5, source_voxel[:,1], source_voxel[:,2], source_dist[:,0], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')

"""
#data_path = '../output_data/382_396_not_filter_tps/'
data_path = '../'
source_emd_index = np.loadtxt(data_path + "source_control_index.dat", dtype=np.int)
target_emd_index = np.loadtxt(data_path + "target_control_index.dat", dtype=np.int)

source_voxel = np.loadtxt  (data_path + "source_voxel.dat")
source_dist = np.loadtxt   (data_path + "source_dist.dat")
source_emd = np.loadtxt    (data_path + "source_emd.dat")
source_corresp = np.loadtxt(data_path + "source_corresp.dat")
target_voxel = np.loadtxt  (data_path + "target_voxel.dat")
target_dist = np.loadtxt   (data_path + "target_dist.dat")
target_emd = np.loadtxt    (data_path + "target_emd.dat")
target_corresp = np.loadtxt(data_path + "target_corresp.dat")

#mmlab.figure(figure='vis', bgcolor=(1.0,1.0,1.0), size=(800,600))
mmlab.figure(figure='vis', size=(1600,900))
mmlab.points3d(  source_voxel[:,0]-0.5,   source_voxel[:,1],   source_voxel[:,2],     source_dist[:,               int(sys.argv[1])], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.points3d(source_corresp[:,0]+1.5, source_corresp[:,1], source_corresp[:,2],     source_dist[:,               int(sys.argv[1])], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.points3d(    source_emd[:,0]+0.5,     source_emd[:,1],     source_emd[:,2],     source_dist[source_emd_index,int(sys.argv[1])], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.points3d(  target_voxel[:,0]-0.5,   target_voxel[:,1],   target_voxel[:,2]+1.0, target_dist[:,               int(sys.argv[1])], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.points3d(target_corresp[:,0]+1.5, target_corresp[:,1], target_corresp[:,2]+1.0, target_dist[:,               int(sys.argv[1])], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.points3d(    target_emd[:,0]+0.5,     target_emd[:,1],     target_emd[:,2]+1.0, target_dist[target_emd_index,int(sys.argv[1])], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.show()


"""
for i in xrange(7):
    mmlab.points3d(source_voxel[:,0]+4*(3-i)*0.5, source_voxel[:,1], source_voxel[:,2], source_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
    mmlab.points3d(target_voxel[:,0]+4*(3-i)*0.5+1.0, target_voxel[:,1], target_voxel[:,2], target_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
    mmlab.points3d(  source_emd[:,0]+4*(3-i)*0.5+0.5,   source_emd[:,1],   source_emd[:,2], source_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
    mmlab.points3d(  target_emd[:,0]+4*(3-i)*0.5+1.5,   target_emd[:,1],   target_emd[:,2], target_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')

for i in xrange(7,15):
    mmlab.points3d(source_voxel[:,0]+4*(10-i)*0.5, source_voxel[:,1],     source_voxel[:,2]+2.0, source_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
    mmlab.points3d(target_voxel[:,0]+4*(10-i)*0.5+1.0, target_voxel[:,1], target_voxel[:,2]+2.0, target_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
    mmlab.points3d(  source_emd[:,0]+4*(10-i)*0.5+0.5,   source_emd[:,1],   source_emd[:,2]+2.0, source_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
    mmlab.points3d(  target_emd[:,0]+4*(10-i)*0.5+1.5,   target_emd[:,1],   target_emd[:,2]+2.0, target_dist[:,i], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
mmlab.show()

"""

#
#mmlab.points3d(source_voxel[:,0]-0.5, source_voxel[:,1], source_voxel[:,2], source_dist[:,1], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(target_voxel[:,0]+2.5, target_voxel[:,1], target_voxel[:,2], target_dist[:,1], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(  source_emd[:,0]+0.5,   source_emd[:,1],   source_emd[:,2], source_dist[:,1], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(  target_emd[:,0]+3.5,   target_emd[:,1],   target_emd[:,2], target_dist[:,1], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.show()

#mmlab.points3d(source_voxel[:,0]-0.5, source_voxel[:,1], source_voxel[:,2], source_dist[:,5], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(target_voxel[:,0]+2.5, target_voxel[:,1], target_voxel[:,2], target_dist[:,5], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(  source_knn[:,0]+0.5,   source_knn[:,1],   source_knn[:,2], source_dist[:,5], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(  target_emd[:,0]+3.5,   target_emd[:,1],   target_emd[:,2], target_dist[:,5], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.show()

#mmlab.points3d(source_voxel[:,0]-0.5, source_voxel[:,1], source_voxel[:,2], source_dist[:,3], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(target_voxel[:,0]+2.5, target_voxel[:,1], target_voxel[:,2], target_dist[:,3], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(  source_emd[:,0]+0.5,   source_emd[:,1],   source_emd[:,2], source_dist[:,3], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.points3d(  target_emd[:,0]+3.5,   target_emd[:,1],   target_emd[:,2], target_dist[:,3], colormap='jet', scale_mode='none', scale_factor=0.015, opacity=1.0, mode='cube')
#mmlab.show()

