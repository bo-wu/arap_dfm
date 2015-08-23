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
#include "mass_transport.h"
#include "util.h"

#define BASIC_DEBUG_

Morph::Morph(std::string source_mesh_name, std::string target_mesh_name, CorrespType corresp_pairs, Real source_voxel_size, Real source_dense_voxel_size, Real target_voxel_size, Real target_dense_voxel_size) : 
    source_mesh_name_(source_mesh_name),
    target_mesh_name_(target_mesh_name),
    corresp_pairs_(corresp_pairs),
    source_voxel_size_(source_voxel_size),
    source_dense_voxel_size_(source_dense_voxel_size),
    target_voxel_size_(target_voxel_size),
    target_dense_voxel_size_(target_dense_voxel_size)
{
    /*
    source_volume_ = VolumeObject(source_mesh_name, voxel_size, dense_voxel_size);
    target_volume_ = VolumeObject(target_mesh_name, voxel_size, dense_voxel_size);
    */

    //======== initial skeleton info ========
    source_skel_.init(source_mesh_name);
    target_skel_.init(target_mesh_name);
    std::size_t found = target_mesh_name.find_last_of("/\\");
    std::string skel_pair_file_name = source_mesh_name.substr(0, source_mesh_name.size()-4) + "_" + target_mesh_name.substr(found+1, target_mesh_name.size()-found-5) + ".brp";
    skel_pair_.read_match_info(source_skel_, target_skel_, skel_pair_file_name);

    source_volume_.initial(source_mesh_name, source_voxel_size, source_dense_voxel_size);
    target_volume_.initial(target_mesh_name, target_voxel_size, target_dense_voxel_size);

    // ========================================
    // add more anchor points from skeleton
    corresp_pairs_.insert(corresp_pairs_.end(), skel_pair_.corresp_skel_points.begin(), skel_pair_.corresp_skel_points.end());

    for(int i=0; i < corresp_pairs_.size(); ++i)
    {
        source_volume_.mAnchors.push_back(corresp_pairs_[i].first);
        target_volume_.mAnchors.push_back(corresp_pairs_[i].second);
    }
    std::cout << "source EMD voxel num "<< source_volume_.voxel_num_ <<" DFI voxel num "<< source_volume_.dense_voxel_num_<<std::endl;
    std::cout << "target EMD voxel num "<< target_volume_.voxel_num_ <<" DFI voxel num "<< target_volume_.dense_voxel_num_<<std::endl;

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
void Morph::initial(bool part_morph)
{
    std::cout <<"initialize morphing process, it will take minutes ... "<<std::endl;

    std::clock_t start;
    Real elapse;
    //correspondence after TPS, 
    //S_T: source-to-target, T_S: target-to-source
    MatrixX3r corresp_S_T_, corresp_T_S_;

    start = std::clock();

    //could be parallel
//#pragma omp parallel sections
//{
//    #pragma omp section
//    {
    std::cout<<"calculating source harmonic field ... "<<std::flush;
    source_volume_.calc_vector_field();  
    std::cout<<"done!\n";
//    }
//    #pragma omp section
//    {
    std::cout<<"calculating target harmonic field ... "<<std::flush;
    target_volume_.calc_vector_field();
    std::cout<<"done!\n";
//    }
//}

    EMD emd_flow;
    PartMassTransport pmt_flow;
    ThinPlateSpline source_target_tps, target_source_tps;

    if(part_morph)
    {   
        source_volume_.segment_volume_voxel(source_skel_);
        target_volume_.segment_volume_voxel(target_skel_);
        pmt_flow.construct_correspondence(source_volume_, target_volume_, skel_pair_);
        
        source_target_tps.compute_tps(pmt_flow.source_control_points_, pmt_flow.corresp_source_target_);
        source_target_tps.interpolate(source_volume_.mDenseVoxelPosition, corresp_S_T_);
        source_volume_.calc_tetrahedron_transform(corresp_S_T_);

        target_source_tps.compute_tps(pmt_flow.target_control_points_, pmt_flow.corresp_target_source_);
        target_source_tps.interpolate(target_volume_.mDenseVoxelPosition, corresp_T_S_);
        target_volume_.calc_tetrahedron_transform(corresp_T_S_);
    }
    else
    {
            emd_flow.construct_correspondence(source_volume_, target_volume_);
        
        #pragma omp parallel sections
        {
            #pragma omp section
            {
            source_target_tps.compute_tps(emd_flow.source_control_points_, emd_flow.corresp_source_target_);
            source_target_tps.interpolate(source_volume_.mDenseVoxelPosition, corresp_S_T_);
            source_volume_.calc_tetrahedron_transform(corresp_S_T_);
            }
            
            #pragma omp section
            {
            target_source_tps.compute_tps(emd_flow.target_control_points_, emd_flow.corresp_target_source_);
            target_source_tps.interpolate(target_volume_.mDenseVoxelPosition, corresp_T_S_);
            target_volume_.calc_tetrahedron_transform(corresp_T_S_);
            }
        
        }
    }

    elapse = (std::clock() - start) / (Real)(CLOCKS_PER_SEC);
    std::cout << "done!     morph initial uses "<< elapse <<"s\n";

#ifdef BASIC_DEBUG_
    std::string path = source_mesh_name_.substr(0, source_mesh_name_.find_last_of("\\/")+1);
    path = path + "output/";

    std::ofstream output_source_dist(path+"source_dist.dat");
    output_source_dist << source_volume_.distance_vector_field;
    output_source_dist.close();
    std::ofstream output_target_dist(path+"target_dist.dat");
    output_target_dist << target_volume_.distance_vector_field;
    output_target_dist.close();

    std::ofstream output_source_fixed_point(path+"source_fixed_point.obj");
    output_source_fixed_point << "v "<< source_volume_.mass_center.transpose()<<std::endl;
    for(int i=0; i < source_volume_.tet_anchor.size(); ++i)
    {
        output_source_fixed_point <<"v "<<(source_volume_.tet_anchor[i].first).transpose()<<std::endl;
    }
    output_source_fixed_point.close();
    std::cout<<"fixed points num "<<source_volume_.tet_anchor.size()<<std::endl;

    std::ofstream output_target_fixed_point(path+"target_fixed_point.obj");
    output_target_fixed_point << "v "<< target_volume_.mass_center.transpose()<<std::endl;
    for(int i=0; i < target_volume_.tet_anchor.size(); ++i)
    {
        output_target_fixed_point <<"v "<<(target_volume_.tet_anchor[i].first).transpose()<<std::endl;
    }
    output_target_fixed_point.close();

    if(part_morph)
    {
        matrix_to_point_cloud_file(pmt_flow.corresp_source_target_, path+"source_pmt");
        matrix_to_point_cloud_file(pmt_flow.corresp_target_source_, path+"target_pmt");
        matrix_to_point_cloud_file(corresp_S_T_, path+"source_pmt_corresp");
        matrix_to_point_cloud_file(corresp_T_S_, path+"target_pmt_corresp");
    }
    else
    {
        matrix_to_point_cloud_file(emd_flow.corresp_source_target_, path+"source_emd");
        matrix_to_point_cloud_file(emd_flow.corresp_target_source_, path+"target_emd");
        matrix_to_point_cloud_file(corresp_S_T_, path+"source_emd_corresp");
        matrix_to_point_cloud_file(corresp_T_S_, path+"target_emd_corresp");
    }
    matrix_to_point_cloud_file(source_volume_.mDenseVoxelPosition, path+"source_voxel");
    matrix_to_point_cloud_file(target_volume_.mDenseVoxelPosition, path+"target_voxel");
#endif
    
}		/* -----  end of function initial  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  start_basic_morph
 *  Description:  
 * =====================================================================================
 */
void Morph::start_basic_morph (Real step_size)
{
    auto min_size = std::min(source_dense_voxel_size_, target_dense_voxel_size_);
    int dim =  0.5 * 1.2 / min_size;
    MatrixX3r grid_vertex(8*dim*dim*dim, 3);

    openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(min_size);

    openvdb::FloatGrid::Ptr temp_grid = openvdb::FloatGrid::create(2.0);
    temp_grid->setTransform(grid_transform);

//    openvdb::Coord xyz;
#pragma omp parallel for collapse(3)
    for(int i=-dim; i < dim; ++i)
        for(int j=-dim; j < dim; ++j)
            for(int k=-dim; k < dim; ++k)
            {
                openvdb::Coord xyz;
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
}		/* -----  end of function start_basic_morph  ----- */


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

///#pragma omp parallel sections
///{
///    #pragma omp section
///    {

    //backwards(source) mapping
    MatrixX3r source_intermedium;
    std::cout<<"source ";
    source_volume_.find_intermedium_points(source_intermedium, t);

    MatrixX3r target_intermedium;
    std::cout <<"target ";
    target_volume_.find_intermedium_points(target_intermedium, 1-t);

    //====== align intermedium points 
   // RowVector3r average_offset = RowVector3r::Zero(3);

   // for(int i=0; i < source_volume_.mAnchors.size(); ++i)
   // {
   //     average_offset += target_intermedium.row(target_volume_.constraint_index_(i)) - source_intermedium.row(source_volume_.constraint_index_(i));
   // }
   // average_offset /= source_volume_.mAnchors.size();

   // source_intermedium.rowwise() += average_offset;
    //====== end alignment 
    //
 
    std::cout<<"backwards to source ";
    source_tps.compute_tps(source_intermedium, source_volume_.mDenseVoxelPosition);
    source_tps.interpolate(grid_vertex, corresp_source_grid_points);

//    }
//
//    #pragma omp section
//    {

    //backwards(target) mapping
    std::cout<<"backwards to target ";

#ifdef BASIC_DEBUG_
    std::string path = source_mesh_name_.substr(0, source_mesh_name_.find_last_of("\\/")+1);
    path = path + "output/";
    std::string source_inter_name = path + "source_intermedium" + std::to_string(int(10*t));
    matrix_to_point_cloud_file(source_intermedium, source_inter_name);
    std::string target_inter_name = path + "target_intermedium" + std::to_string(int(10*t));
    matrix_to_point_cloud_file(target_intermedium, target_inter_name);
#endif

    target_tps.compute_tps(target_intermedium, target_volume_.mDenseVoxelPosition);
    target_tps.interpolate(grid_vertex, corresp_target_grid_points);

//    }
//
//}
    openvdb::FloatGrid::ConstAccessor source_accessor = source_volume_.dense_grid->getConstAccessor();
    openvdb::FloatGrid::ConstAccessor target_accessor = target_volume_.dense_grid->getConstAccessor();

    openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> source_sampler(source_accessor, source_volume_.dense_grid->transform());
    openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> target_sampler(target_accessor, target_volume_.dense_grid->transform());

    openvdb::FloatGrid::Accessor accessor = morph_grid->getAccessor();

    /*
    */
    openvdb::Coord xyz;
    Real value;
    int index;
    Vector3r source_vert, target_vert;

    auto min_size = std::min(source_dense_voxel_size_, target_dense_voxel_size_);
    int dim = 0.5 * 1.2 / min_size;

    int i, j, k;
////#pragma omp parallel for collapse(3) private(i, j, k)
    for(i=-dim; i < dim; ++i)
        for(j=-dim; j < dim; ++j)
            for(k=-dim; k < dim; ++k)
            {
                //openvdb::Coord xyz;
                xyz.reset(i, j, k);
                int index = 4*(i+dim)*dim*dim + 2*(j+dim)*dim + k+dim;

                source_vert = corresp_source_grid_points.row(index);
                target_vert = corresp_target_grid_points.row(index);
                value = (1-t) * source_sampler.wsSample(openvdb::Vec3d(source_vert(0), source_vert(1), source_vert(2)))
                    + t * target_sampler.wsSample(openvdb::Vec3R(target_vert(0), target_vert(1), target_vert(2)));
////#pragma omp critical
                {
                if(value < 0.1 && value > -0.1)
                    accessor.setValue(xyz, value);
                else
                    accessor.setValueOff(xyz);
                }
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

