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
void EMD::construct_correspondence(const VolumeObject &s, const VolumeObject &t)
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
void EMD::find_correspondence(const VolumeObject &s, const VolumeObject &t)
{
    corresp_source_target_ = MatrixX3r::Zero(s.voxel_num_, 3);
    corresp_target_source_ = MatrixX3r::Zero(t.voxel_num_, 3);
    VectorXi index_each_row = VectorXi::Ones(s.voxel_num_); //store max value col index
    VectorXr max_each_row = VectorXr::Zero(s.voxel_num_);
    VectorXi index_each_col = VectorXi::Ones(t.voxel_num_); //store max vlaue row index
    VectorXr max_each_col = VectorXr::Zero(t.voxel_num_);
    index_each_row *= -1;
    index_each_col *= -1;

    for(int i=0; i < flow_matrix_.outerSize(); ++i)
    {
        for(SpMat::InnerIterator it(flow_matrix_, i); it; ++it)
        {
            if (max_each_row(it.row()) < it.value())
            {
                max_each_row(it.row()) = it.value();
                index_each_row(it.row()) = it.col();
            }

            if (max_each_col(it.col()) < it.value())
            {
                max_each_col(it.col()) = it.value();
                index_each_col(it.col()) = it.row();
            }
            /*  
            corresp_source_target_.row(it.row()) += it.value() * t.mVoxelPosition.row(it.col());
            corresp_target_source_.row(it.col()) += it.value() * s.mVoxelPosition.row(it.row());
            */
        }
    }

    for(int i=0; i < s.voxel_num_; ++i)
    {
        corresp_source_target_.row(i) = t.mVoxelPosition.row(index_each_row(i));
    }

    for(int i=0; i < t.voxel_num_; ++i)
    {
        corresp_target_source_.row(i) = s.mVoxelPosition.row(index_each_col(i));
    }

    /*  
    corresp_source_target_ *= s.voxel_num_;
    corresp_target_source_ *= t.voxel_num_;
    */

    /* //for debug
    std::ofstream output_st("corresp_st.dat");
    output_st << corresp_source_target_;
    output_st.close();
    std::ofstream output_ts("corresp_ts.dat");
    output_ts << corresp_target_source_;
    output_st.close();
    */
}		/* -----  end of function find_correspondence ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  min_cost_flow
 *  Description:  
 * =====================================================================================
 */
void EMD::min_cost_flow(const VolumeObject &s, const VolumeObject &t)
{
    using namespace lemon;
    typedef FullBipartiteDigraph Digraph;
    typedef NetworkSimplexSimple<Digraph, Real, Real, unsigned short int> MyNetwork;
    DIGRAPH_TYPEDEFS(FullBipartiteDigraph);

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

    flow_matrix_ = SpMat(s.voxel_num_, t.voxel_num_);
    std::vector<MyTriplet> flow_trip;
    flow_trip.reserve(num_flow);
    Real amount;
    for(int i=0; i < s.voxel_num_; ++i)
        for(int j=0; j < t.voxel_num_; ++j)
        {
            amount = net.flow(di.arcFromId(i*t.voxel_num_ + j));
            if(fabs(amount) > 1e-8)
            {
                flow_trip.push_back(MyTriplet(i, j, amount));
            }
        }
    flow_matrix_.setFromTriplets(flow_trip.begin(), flow_trip.end());


/*  // for debug
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
