/*
 * =====================================================================================
 *
 *       Filename:  mass_transport.cpp      Created:  08/19/2015 11:21:56 PM
 *
 *    Description:  part transport
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <nanoflann.hpp>
#include "mass_transport.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  network_simplex
 *  Description:  
 * =====================================================================================
 */
void PartMassTransport::network_simplex(const VolumeObject &source, const VolumeObject &target, SkeletonPair &skel_pair)
{
    flow_matrix_ = SpMat(source.voxel_num_, target.voxel_num_);

    int num_flow = 2 * int(std::max(Real(source.voxel_num_) / Real(target.voxel_num_), Real(target.voxel_num_) / Real(source.voxel_num_))  * std::max(source.voxel_num_, target.voxel_num_) );
    total_flow_triplet_.reserve(num_flow);

    for(int i=0; i < skel_pair.matched_branch.size(); ++i)
    {
        int source_part_index = skel_pair.matched_branch[i].first;
        int target_part_index = skel_pair.matched_branch[i].second;
        part_network_simplex(source, target, source_part_index, target_part_index);
    }

    flow_matrix_.setFromTriplets(total_flow_triplet_.begin(), total_flow_triplet_.end());
}		/* -----  end of function network_simplex  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  part_network_simplex
 *  Description:  
 * =====================================================================================
 */
void PartMassTransport::part_network_simplex(const VolumeObject &source, const VolumeObject &target, int source_part_idx, int target_part_idx)
{
    lemon::ListDigraph g;
    lemon::ListDigraph::NodeMap<int> supply(g);
    lemon::ListDigraph::ArcMap<double> cost(g);
    lemon::ListDigraph::ArcMap<int> flow(g);

    // part volume index
    std::vector<int> source_volume_index = source.volume_part_index_[source_part_idx];
    std::vector<int> target_volume_index = target.volume_part_index_[target_part_idx];
    int num_source_part_voxel = source_volume_index.size();
    int num_target_part_voxel = target_volume_index.size();

    for(int i=0; i < num_source_part_voxel; ++i)
    {
        lemon::ListDigraph::Node m = g.addNode();
        supply[m] = num_target_part_voxel;
    }

    for(int i=0; i < num_target_part_voxel; ++i)
    {
        lemon::ListDigraph::Node n = g.addNode();
        supply[n] = -1 * num_source_part_voxel;

        for(lemon::ListDigraph::NodeIt s(g); s!=lemon::INVALID; ++s)
        {
            if(lemon::ListDigraph::id(s) < num_source_part_voxel)
            {
                lemon::ListDigraph::Arc arc = g.addArc(s, n);
                cost[arc] = ( source.distance_vector_field.row(source_volume_index[lemon::ListDigraph::id(s)]) - target.distance_vector_field.row(target_volume_index[i]) ).norm();
            }
        }
    }

    lemon::NetworkSimplex<lemon::ListDigraph, int, double> ns(g);
    ns.costMap(cost).supplyMap(supply).run();
    ns.flowMap(flow);

    Real thres = 0.001;
    int num_thres = int(thres * num_target_part_voxel);
    for(lemon::ListDigraph::ArcIt s(g); s!=lemon::INVALID; ++s )
    {
        if(flow[s] > num_thres)
        {
            total_flow_triplet_.push_back(MyTriplet(source_volume_index[g.id(g.source(s))], target_volume_index[g.id(g.target(s)) - num_source_part_voxel], flow[s]/(Real)(num_source_part_voxel*num_target_part_voxel)) );
        }
        else if (flow[s] < 0)
        {
            std::cout<<"exists negative flow, sth is wrong!\n";
        }
    }
}		/* -----  end of function part_network_simplex  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_correspondence
 *  Description:  
 * =====================================================================================
 */
void PartMassTransport::find_correspondence(const VolumeObject &s, const VolumeObject &t, const Real threshold)
{
    VectorXi index_each_row = VectorXi::Ones(s.voxel_num_); //col with max value 
    VectorXr max_each_row = VectorXr::Zero(s.voxel_num_);

    VectorXi index_each_col = VectorXi::Ones(t.voxel_num_);
    VectorXr max_each_col = VectorXr::Zero(t.voxel_num_);
    index_each_row *= -1;
    index_each_col *= -1;

    // average coordinates
    MatrixX3r average_source_target = MatrixX3r::Zero(s.voxel_num_, 3);
    MatrixX3r average_target_source = MatrixX3r::Zero(t.voxel_num_, 3);

    for(int i=0; i < flow_matrix_.outerSize(); ++i)
    {
        for(SpMat::InnerIterator it(flow_matrix_, i); it; ++it)
        {
            if(max_each_row(it.row()) < it.value())
            {
                max_each_row(it.row()) = it.value();
                index_each_row(it.row()) = it.col();
            }

            if(max_each_col(it.col()) < it.value())
            {
                max_each_col(it.col()) = it.value();
                index_each_col(it.col()) = it.row();
            }

            average_source_target.row(it.row()) += it.value() * t.mVoxelPosition.row(it.col());
            average_target_source.row(it.col()) += it.value() * s.mVoxelPosition.row(it.row());
        }
    }

    //filter out bad flow
    std::vector<int> source_control, corresp_s_t;
    std::vector<int> target_control, corresp_t_s;
    source_control.reserve(s.voxel_num_);
    corresp_s_t.reserve(s.voxel_num_);
    target_control.reserve(t.voxel_num_);
    corresp_t_s.reserve(t.voxel_num_);

    VectorXr dist_vector_diff(s.mAnchors.size());
    Real real_threshold = threshold * s.mAnchors.size();

    for(int i=0; i < s.voxel_num_; ++i)
    {
        dist_vector_diff = s.distance_vector_field.row(i) - t.distance_vector_field.row(index_each_row(i));
        if(dist_vector_diff.squaredNorm() < real_threshold)
        {
            source_control.push_back(i);
            corresp_s_t.push_back(index_each_row(i));
        }
    }

    for(int i=0; i < t.voxel_num_; ++i)
    {
        dist_vector_diff = t.distance_vector_field.row(i) - s.distance_vector_field.row(index_each_col(i));
        if(dist_vector_diff.squaredNorm() < real_threshold)
        {
            target_control.push_back(i);
            corresp_t_s.push_back(index_each_col(i));
        }
    }

    source_control_points_ = MatrixX3r(source_control.size(), 3);
    corresp_source_target_ = MatrixX3r(source_control.size(), 3);
    for(int i=0; i < source_control.size(); ++i)
    {
        source_control_points_.row(i) = s.mVoxelPosition.row(source_control[i]);
        // the max voxel
        // corresp_source_target_.row(i) = t.mVoxelPosition.row(corresp_s_t[i]);
        // average voxel
        corresp_source_target_.row(i) = average_source_target.row(source_control[i]);
    }

    target_control_points_ = MatrixX3r(target_control.size(), 3);
    corresp_target_source_ = MatrixX3r(target_control.size(), 3);
    for(int i=0; i < target_control.size(); ++i)
    {
        target_control_points_.row(i) = t.mVoxelPosition.row(target_control[i]);
        // the max voxel
        // corresp_target_source_.row(i) = s.mVoxelPosition.row(corresp_t_s[i]);
        // average voxel
        corresp_target_source_.row(i) = average_target_source.row(target_control[i]);
    }

    corresp_source_target_ *= s.voxel_num_;
    corresp_target_source_ *= t.voxel_num_;

    std::cout << "part source control point num " << source_control.size()<<std::endl;
    std::cout << "part target control point num " << target_control.size()<<std::endl;

#ifdef BASIC_DEBUG_
    std::ofstream output_source_control_index("part_source_control_index.dat");
    std::ofstream output_source_control_target_index("part_source_control_target_index.dat");
    for(int i=0; i < source_control.size(); ++i)
    {
        output_source_control_index << source_control[i]<<std::endl;
        output_source_control_target_index << corresp_s_t[i]<<std::endl;
    }
    output_source_control_index.close();
    output_source_control_target_index.close();
    std::ofstream output_target_control_index("part_target_control_index.dat");
    std::ofstream output_target_control_source_index("part_target_control_source_index.dat");
    for(int i=0; i < target_control.size(); ++i)
    {
        output_target_control_index << target_control[i] <<std::endl;
        output_source_control_target_index << corresp_t_s[i] <<std::endl;
    }
    output_target_control_index.close();
    output_target_control_source_index.close();
#endif


}		/* -----  end of function find_correspondence  ----- */
