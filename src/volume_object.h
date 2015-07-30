/*
 * =====================================================================================
 *
 *       Filename:  volume_object.h  Version:  1.0  Created:  01/16/2015 04:57:56 PM
 *
 *    Description:  volume representation
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef VOLUME_OBJECT_H_
#define VOLUME_OBJECT_H_
#include <string>
#include <openvdb/openvdb.h>
#include "kdtree.h"
#include "def_types.h"

struct VolumeObject
{
	TriMesh mesh;
	std::string mesh_name;
	openvdb::FloatGrid::Ptr grid;
    //all active are inside tiles and voxel
    openvdb::BoolGrid::Ptr interior_grid;
	std::vector<openvdb::Vec3s> points;
	std::vector<openvdb::Vec3I> triangles;
	std::vector<Vector3r> mAnchors;
    SpMat mLaplaceMatrix;
    Real transform_scale_;
    int voxel_num_;
    // constraint voxel index
    VectorXi constraint_index_;
    MatrixX3r mVoxelPosition;
    MatrixXr distance_vector_field;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<MatrixX3r, 3, nanoflann::metric_L2_Simple> kd_tree_type;

	VolumeObject(Real transform_scale=0.01);
	VolumeObject(std::string mesh_name, Real transform_scale=0.01);
	~VolumeObject();
	void initial_volume();
	//compute vector field on anchor points
	void calc_vector_field();
	void read_mesh(bool resize=true);
	void write_grid(std::string name);
    void test_volume();
    void set_anchors(std::vector<Vector3r>& anchors);
    void construct_laplace_matrix();
};

#endif

