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

class VolumeObject
{
public:
	TriMesh mesh;
	std::string mesh_name;
	openvdb::FloatGrid::Ptr grid, dense_grid;
    //all active are inside tiles and voxel
    openvdb::BoolGrid::Ptr interior_grid, interior_dense_grid;

	std::vector<openvdb::Vec3s> points;
	std::vector<openvdb::Vec3I> triangles;
    // anchors for distance field
	std::vector<Vector3r> mAnchors;
    //tetrahedron index
    std::vector<Vector4i> mTetIndex;
    std::vector<std::pair<Matrix3r, Matrix3r> > mTetTransform;

    SpMat mLaplaceMatrix;
    Real transform_scale_;
    Real dense_transform_scale_;

    // constraint voxel index
    VectorXi constraint_index_;
    int voxel_num_; //sparse interior voxel num
    int  dense_voxel_num_;

    // relative sparse grid for TPS (then use TPS )
    MatrixX3r mVoxelPosition;
    //Dense Grid Position used for interpolate distance field
    MatrixX3r mDenseVoxelPosition;
    Vector3r mass_center; // object center
    int mass_center_voxel_index;
    // keep solver stability
    std::vector<std::pair<Vector3r, int> > tet_anchor;

    MatrixXr distance_vector_field;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<MatrixX3r, 3, nanoflann::metric_L2_Simple> kd_tree_type;

    VolumeObject(Real transform_scale=0.02);
	VolumeObject(std::string mesh_name, Real transform_scale=0.02, Real dense_transform_scale=0.01);
	~VolumeObject();

    void initial(std::string name, Real transform_scale=0.02, Real dense_transform_scale=0.01);
	void initial_volume();
    void initial_dense_volume();

	//compute vector field on anchor points
    void construct_laplace_matrix();
	void calc_vector_field(); // harmonic field for distance 
	void read_mesh(bool resize=true, bool b_out=true); 
	void write_grid(std::string name); // output for test
    void test_volume(); //for debug
    void set_anchors(const std::vector<Vector3r>& anchors);

    void polar_decompose(const Matrix3r &rest, const Matrix3r &deform, Matrix3r &R, Matrix3r &S);
    void calc_tetrahedron_transform(const MatrixX3r &final_corresp_points);
    // find intermedium correspondence points during morphing
    void find_intermedium_points(MatrixX3r &inter_corresp_points, const Real t=0.5);

    std::vector<MyTriplet> triplet_with_constraint;

    bool first_time_;
    SpMat tet_matrix_, tet_normal_matrix;
    Eigen::ConjugateGradient<SpMat> cg;
};

#endif

