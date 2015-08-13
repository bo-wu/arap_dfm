/*
 * =====================================================================================
 *
 *       Filename:  morph.cpp  Version:  1.0  Created:  01/16/2015 04:56:28 PM
 *
 *    Description:  morph two volume data
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <openvdb/tools/Prune.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/SignedFloodFill.h>
#include "morph.h"
#include "util.h"

Morph::Morph(std::string source_mesh_name, std::string target_mesh_name, CorrespType corresp_pairs, Real voxel_size) : 
    source_mesh_name_(source_mesh_name),
    target_mesh_name_(target_mesh_name),
    corresp_pairs_(corresp_pairs),
    voxel_size_(voxel_size)
{
    source_volume_ = VolumeObject(source_mesh_name, voxel_size);
    target_volume_ = VolumeObject(target_mesh_name, voxel_size);

    /*
    source_volume_.initial(source_mesh_name, voxel_size);
    target_volume_.initial(target_mesh_name, voxel_size);
    */

    for(int i=0; i < corresp_pairs_.size(); ++i)
    {
        source_volume_.mAnchors.push_back(corresp_pairs_[i].first);
        target_volume_.mAnchors.push_back(corresp_pairs_[i].second);
    }
    std::cout << "source voxel num "<< source_volume_.voxel_num_<<std::endl;
    std::cout << "target voxel num "<< target_volume_.voxel_num_<<std::endl;
}

Morph::~Morph()
{

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initial
 *  Description:  
 * =====================================================================================
 */
void Morph::initial()
{
    std::cout <<"initialize morphing process, it will take minutes ... "<<std::endl;

    std::clock_t start;
    Real elapse;
    //correspondence after TPS, 
    //S_T: source-to-target, T_S: target-to-source
    MatrixX3r corresp_S_T_, corresp_T_S_;

    start = std::clock();

    //could be parallel
    source_volume_.calc_vector_field();  
    target_volume_.calc_vector_field();

    EMD emd_flow;
    emd_flow.construct_correspondence(source_volume_, target_volume_);

    ThinPlateSpline source_target_tps, target_source_tps;
    
    source_target_tps.compute_tps(emd_flow.source_control_points_, emd_flow.corresp_source_target_);
    source_target_tps.interpolate(source_volume_.mDenseVoxelPosition, corresp_S_T_);
    source_volume_.calc_tetrahedron_transform(corresp_S_T_);
    
    matrix_to_point_cloud_file(emd_flow.corresp_source_target_, "source_emd");
    matrix_to_point_cloud_file(emd_flow.corresp_target_source_, "target_emd");
    matrix_to_point_cloud_file(source_volume_.mDenseVoxelPosition, "source_voxel");
    matrix_to_point_cloud_file(target_volume_.mDenseVoxelPosition, "target_voxel");

    std::ofstream output_source_dist("source_dist.dat");
    output_source_dist << source_volume_.distance_vector_field;
    output_source_dist.close();
    std::ofstream output_target_dist("target_dist.dat");
    output_target_dist << target_volume_.distance_vector_field;
    output_target_dist.close();

    /*  
    std::ofstream output_source_emd("source_emd.dat");
    output_source_emd << emd_flow.corresp_source_target_;
    output_source_emd.close();

    std::ofstream output_source_voxel("source_voxel.dat");
    output_source_voxel<< source_volume_.mDenseVoxelPosition;
    output_source_voxel.close();
    std::ofstream output_source_corresp("source_corresp.dat");
    output_source_corresp << corresp_S_T_;
    output_source_corresp.close();
    */

    target_source_tps.compute_tps(emd_flow.target_control_points_, emd_flow.corresp_target_source_);
    target_source_tps.interpolate(target_volume_.mDenseVoxelPosition, corresp_T_S_);
    target_volume_.calc_tetrahedron_transform(corresp_T_S_);

    matrix_to_point_cloud_file(corresp_S_T_, "source_corresp");
    matrix_to_point_cloud_file(corresp_T_S_, "target_corresp");
    
    /*  
    std::ofstream output_target_emd("target_emd.dat");
    output_target_emd << emd_flow.corresp_target_source_;
    output_target_emd.close();

    std::ofstream output_target_voxel("target_voxel.dat");
    output_target_voxel<< target_volume_.mDenseVoxelPosition;
    output_target_voxel.close();
    std::ofstream output_target_corresp("target_corresp.dat");
    output_target_voxel << corresp_T_S_;
    output_target_voxel.close();
    */

    elapse = (std::clock() - start) / (Real)(CLOCKS_PER_SEC);
    std::cout << "done!     morph initial uses "<< elapse <<"s\n";
}		/* -----  end of function initial  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  start_morph
 *  Description:  
 * =====================================================================================
 */
void Morph::start_morph (Real step_size)
{
    openvdb::Coord xyz;
    int dim =  0.5 * 1.2 / voxel_size_;
    MatrixX3r grid_vertex(8*dim*dim*dim, 3);

    openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(voxel_size_);

    openvdb::FloatGrid::Ptr temp_grid = openvdb::FloatGrid::create(2.0);
    temp_grid->setTransform(grid_transform);

    for(int i=-dim; i < dim; ++i)
        for(int j=-dim; j < dim; ++j)
            for(int k=-dim; k < dim; ++k)
            {
                xyz.reset(i, j, k);
                auto voxel_pos = temp_grid->indexToWorld(xyz);
                grid_vertex.row(4*(i+dim)*dim*dim + 2*(j+dim)*dim + k+dim) << voxel_pos[0], voxel_pos[1], voxel_pos[2];
            }

    std::cout<<"start morphing \n";
    int steps = 1 / step_size;
    std::string grid_name;
    Real elapse;
    std::clock_t start;
    for(int i=0; i <= steps; ++i)
    {
        start = std::clock();
        openvdb::FloatGrid::Ptr morph_grid = openvdb::FloatGrid::create(2.0);
        morph_grid->setTransform(grid_transform);
        morph_grid->setGridClass(openvdb::GRID_LEVEL_SET);

        std::stringstream ss;
        ss << std::setw(4) << std::setfill('0') << i;
        grid_name = ss.str();
        morph_grid->setName(grid_name.c_str());
        
        interpolate_grids(morph_grid, grid_vertex, i*step_size);

        elapse = (std::clock() - start) / (Real)(CLOCKS_PER_SEC);
        std::cout<< i <<"/"<<steps << " totally use "<<elapse<<"s"<<std::endl;
        grid_vec_.push_back(morph_grid);
    }
}		/* -----  end of function start_morph  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  interpolate_grids
 *  Description:  
 * =====================================================================================
 */
void Morph::interpolate_grids(openvdb::FloatGrid::Ptr &morph_grid, MatrixX3r &grid_vertex, Real t)
{
    ThinPlateSpline source_tps, target_tps;
    MatrixX3r corresp_source_grid_points;
    MatrixX3r corresp_target_grid_points;

    t = std::min(std::max(0.0, t), 1.0);

    //backwards(source) mapping
    MatrixX3r source_intermedium;
    std::cout<<"source ";
    source_volume_.find_intermedium_points(source_intermedium, t);

    std::string source_inter_name = "./output_result_data/source_intermedium" + std::to_string(int(10*t));
    matrix_to_point_cloud_file(source_intermedium, source_inter_name);

    /*  
    std::ofstream output_source_tps("source_intermedium.dat");
    output_source_tps << source_intermedium;
    output_source_tps.close();
    std::ofstream output_source_voxel("source_voxel.dat");
    output_source_voxel << source_volume_.mDenseVoxelPosition;
    output_source_voxel.close();
    */

    std::cout<<"backwards to source ";
    source_tps.compute_tps(source_intermedium, source_volume_.mDenseVoxelPosition);
    source_tps.interpolate(grid_vertex, corresp_source_grid_points);

    //backwards(target) mapping
    MatrixX3r target_intermedium;
    std::cout <<"target ";
    target_volume_.find_intermedium_points(target_intermedium, 1-t);
    std::cout<<"backwards to target ";

    std::string target_inter_name = "./output_result_data/target_intermedium" + std::to_string(int(10*t));
    matrix_to_point_cloud_file(target_intermedium, target_inter_name);

    /*  
    std::ofstream output_target_tps("target_intermedium.dat");
    output_target_tps << target_intermedium;
    output_target_tps.close();
    std::ofstream output_target_voxel("target_voxel.dat");
    output_target_voxel << target_volume_.mDenseVoxelPosition;
    output_target_voxel.close();
    */

    target_tps.compute_tps(target_intermedium, target_volume_.mDenseVoxelPosition);
    target_tps.interpolate(grid_vertex, corresp_target_grid_points);

    openvdb::FloatGrid::ConstAccessor source_accessor = source_volume_.grid->getConstAccessor();
    openvdb::FloatGrid::ConstAccessor target_accessor = target_volume_.grid->getConstAccessor();

    openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> source_sampler(source_accessor, source_volume_.grid->transform());
    openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> target_sampler(target_accessor, target_volume_.grid->transform());

    openvdb::FloatGrid::Accessor accessor = morph_grid->getAccessor();
    openvdb::Coord xyz;

    Real value;
    int index;
    Vector3r source_vert, target_vert;
    int dim = 0.5 * 1.2 / voxel_size_;
    for(int i=-dim; i < dim; ++i)
        for(int j=-dim; j < dim; ++j)
            for(int k=-dim; k < dim; ++k)
            {
                xyz.reset(i, j, k);
                index = 4*(i+dim)*dim*dim + 2*(j+dim)*dim + k+dim;

                source_vert = corresp_source_grid_points.row(index);
                target_vert = corresp_target_grid_points.row(index);
                value = (1-t) * source_sampler.wsSample(openvdb::Vec3d(source_vert(0), source_vert(1), source_vert(2)))
                    + t * target_sampler.wsSample(openvdb::Vec3R(target_vert(0), target_vert(1), target_vert(2)));
                accessor.setValue(xyz, value);
            }

    openvdb::tools::signedFloodFill(morph_grid->tree());
    openvdb::tools::pruneInactive(morph_grid->tree());
}		/* -----  end of function interpolate_grids  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  write_sequence
 *  Description:  
 * =====================================================================================
 */
void Morph::write_sequence(std::string grid_name)
{
    if(grid_name.empty())
    {
        std::size_t found = target_mesh_name_.find_last_of("/\\");
        grid_name = source_mesh_name_.substr(0, source_mesh_name_.size()-4) + "_" 
            + target_mesh_name_.substr(found+1, target_mesh_name_.size()-found-5) + ".vdb";
    }

    openvdb::io::File file(grid_name);
    file.setCompression(openvdb::io::COMPRESS_ZIP|openvdb::io::COMPRESS_ACTIVE_MASK);

    openvdb::GridPtrVec grids;
    for(int i=0; i < grid_vec_.size(); ++i)
    {
        grids.push_back(grid_vec_[i]);
    }

    file.write(grids);
    file.close();
}		/* -----  end of function write_sequence  ----- */

