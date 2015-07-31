/*
 * =====================================================================================
 *
 *       Filename:  find_dense_correspondence.cpp      Created:  07/29/2015 10:20:38 PM
 *
 *    Description:  find dense correspondence from anchors
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include "find_dense_correspondence.h"
#include "full_bipartitegraph.h"
#include "network_simplex_simple.h"
#include "def_types.h"


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  construct_correspondence
 *  Description:  
 * =====================================================================================
 */
void EMD::construct_correspondence(VolumeObject &s, VolumeObject &t)
{
    min_cost_flow(s, t);
    find_correspondence(s, t);
}		/* -----  end of function construct_correspondence  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_correspondence
 *  Description:  find the correspondence between source and target
 * =====================================================================================
 */
void EMD::find_correspondence(VolumeObject &s, VolumeObject &t)
{
    corresp_source_target_ = MatrixX3r::Zero(s.voxel_num_, 3);
    corresp_target_source_ = MatrixX3r::Zero(t.voxel_num_, 3);
    for(int i=0; i < flow_matrix_.outerSize(); ++i)
    {
        for(SpMat::InnerIterator it(flow_matrix_, i); it; ++i)
        {
            corresp_source_target_.row(it.row()) += it.value() * t.mVoxelPosition.row(it.col());
            corresp_target_source_.row(it.col()) += it.value() * s.mVoxelPosition.row(it.row());
        }
    }
    corresp_source_target_ *= s.voxel_num_;
    corresp_target_source_ *= t.voxel_num_;
}		/* -----  end of function find_correspondence ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  min_cost_flow
 *  Description:  
 * =====================================================================================
 */
void EMD::min_cost_flow(VolumeObject &s, VolumeObject &t)
{
    using namespace lemon;
    typedef FullBipartiteDigraph Digraph;
    typedef NetworkSimplexSimple<Digraph, Real, Real, unsigned short int> MyNetwork;
    DIGRAPH_TYPEDEFS(FullBipartiteDigraph);
    std::cout<<"s.voxel_num " << s.voxel_num_<< std::endl;
    std::cout<<"t.voxel_num " << t.voxel_num_<< std::endl;
//    std::cout<<"s.voxel_num * t.voxel_num = "<<s.voxel_num_ * t.voxel_num_<<std::endl;
    Digraph di(s.voxel_num_, t.voxel_num_);
    MyNetwork net(di, false);
    int arc_id = 0;
    Real d;
    VectorXr source_dist_vector;
    for(int i=0; i < s.voxel_num_; ++i)
    {
        source_dist_vector = s.distance_vector_field.row(i);
        for(int j=0; j < t.voxel_num_; ++j)
        {
            d = (source_dist_vector - t.distance_vector_field.row(j)).norm();
            net.setCost(di.arcFromId(arc_id), d);
            ++arc_id;
        }
    }

    Digraph::NodeMap<Real> node_supplies_demand(di);
    for(int i=0; i < s.voxel_num_; ++i)
        node_supplies_demand[di.nodeFromId(i)] = 1.0 / Real(s.voxel_num_);

    for(int i=0; i < t.voxel_num_; ++i)
        node_supplies_demand[di.nodeFromId(i+s.voxel_num_)] = -1.0 / Real(t.voxel_num_);

    net.supplyMap(node_supplies_demand);

    auto ret = net.run(MyNetwork::BLOCK_SEARCH);
    int num_flow = 2 * int(std::max(Real(s.voxel_num_) / Real(t.voxel_num_), Real(t.voxel_num_) / Real(s.voxel_num_))  * std::max(s.voxel_num_, t.voxel_num_) );
//    std::cout << "num_flow: " << num_flow<<std::endl;
//    std::cout <<"std::vector max size " << result_flow_.max_size()<<std::endl;

    flow_matrix_ = SpMat(s.voxel_num_, t.voxel_num_);
    std::vector<Triplet> flow_trip;
    flow_trip.reserve(num_flow);
    Real amount;
    for(int i=0; i < s.voxel_num_; ++i)
        for(int j=0; j < t.voxel_num_; ++j)
        {
            amount = net.flow(di.arcFromId(i*t.voxel_num_ + j));
            if(fabs(amount) > 1e-8)
            {
                flow_trip.push_back(Triplet(i, j, amount));
            }
        }
    flow_matrix_.setFromTriplets(flow_trip.begin(), flow_trip.end());
/*  
    std::ofstream output_flow_row("flow_row.dat");
    std::ofstream output_flow_col("flow_col.dat");
    for(int i=0; i < s.voxel_num_; ++i)
    {
        output_flow_row << flow_matrix_.row(i).sum()<<std::endl;
    }
    for(int j=0; j < t.voxel_num_; ++j)
    {
        output_flow_col << flow_matrix_.col(j).sum()<<std::endl;
    }
    output_flow_row.close();
    output_flow_col.close();
*/
}		/* -----  end of function min_cost_flow  ----- */
