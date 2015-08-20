/*
 * =====================================================================================
 *
 *       Filename:  morph.h  Version:  1.0  Created:  01/16/2015 04:56:49 PM
 *
 *    Description:  morph two volume data
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */

#ifndef MORPH_H_
#define MORPH_H_
#include <vector>
#include <openvdb/openvdb.h>
#include "def_types.h"
#include "skeleton.h"
#include "volume_object.h"
#include "thin_plate_spline.h"
#include "dense_correspondence.h"

struct Morph
{
    Morph(std::string source_mesh_name, std::string target_mesh_name, CorrespType corresp_pairs, Real source_voxel_size=0.02, Real source_dense_voxel_size=0.01, Real target_voxel_size=0.02, Real target_dense_voxel_size=0.01);

	~Morph();
    Real source_voxel_size_, source_dense_voxel_size_;
    Real target_voxel_size_, target_dense_voxel_size_;
    //openvdb::FloatGrid::Ptr result_grid;
    std::vector<openvdb::FloatGrid::Ptr> grid_vec_;
    std::string source_mesh_name_, target_mesh_name_;
    Skeleton source_skel_, target_skel_;
    SkeletonPair skel_pair_;
	VolumeObject source_volume_, target_volume_;
    CorrespType corresp_pairs_;

    //MatrixX3r intermedia_from_source, intermedia_from_target;
    void initial();
	void start_basic_morph(Real step_size=0.01);
    void interpolate_grids(openvdb::FloatGrid::Ptr &morph_grid, MatrixX3r &grid_vertex, Real t);
    
    void morphing_by_skeleton_correspondence(Real step_size=0.01);
    void morphing_by_part(Real step_size=0.01);


    void write_sequence(std::string grid_name="");

};

#endif

