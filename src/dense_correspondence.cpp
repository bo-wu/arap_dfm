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
#include <nanoflann.hpp>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include "dense_correspondence.h"
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
    /*  
    direct_correspondence(s, t);
    */
    //min_cost_flow(s, t); // something wrong
    network_simplex(s, t);
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
            corresp_source_target_.row(it.row()) +=  it.value() * t.mVoxelPosition.row(it.col());
            corresp_target_source_.row(it.col()) +=  it.value() * s.mVoxelPosition.row(it.row());
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
    /*  
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
 *  Description:  something wrong !!!!
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



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  direct_correspondence
 *  Description:  directly use harmonic field to get correspondence
 * =====================================================================================
 */
void EMD::direct_correspondence(const VolumeObject &s, const VolumeObject &t)
{
    int knn_num = 2;
    corresp_source_target_ = MatrixX3r::Zero(s.voxel_num_, 3);
    corresp_target_source_ = MatrixX3r::Zero(t.voxel_num_, 3);

    const int anchor_num = s.mAnchors.size();

    typedef nanoflann::KDTreeEigenMatrixAdaptor<MatrixXr, -1, nanoflann::metric_L2_Simple> KDTreeType;
    KDTreeType source_harmonic_kdtree(anchor_num, s.distance_vector_field);
    source_harmonic_kdtree.index->buildIndex();

    long int out_index[knn_num];
    Real out_distances_sq[knn_num];

    std::ofstream output_knn("knn_target.dat");
    Real w = 0;
    for(int i=0; i < t.voxel_num_; ++i)
    {
        w = 0;
        output_knn <<"voxel "<<i/*<<" ["<<t.distance_vector_field.row(i)<<"]*/ <<" neighboring source ";
        source_harmonic_kdtree.query(t.distance_vector_field.row(i).data(), knn_num, out_index, out_distances_sq);
        for(int j=0; j < knn_num; ++j)
        {
            w += 1.0 / out_distances_sq[j];
            corresp_target_source_.row(i) += (1.0 / out_distances_sq[j]) * s.mVoxelPosition.row(out_index[j]);
            output_knn<<out_index[j] <<": "<<out_distances_sq[j]<<" "; //" ["<<s.distance_vector_field.row(out_index[j])<<"] ";
        }
        corresp_target_source_.row(i) /= w;
        output_knn<<std::endl;
    }
    output_knn.close();

    KDTreeType target_harmonic_kdtree(anchor_num, t.distance_vector_field);
    target_harmonic_kdtree.index->buildIndex();

    for(int i=0; i < s.voxel_num_; ++i)
    {
        w = 0;
        target_harmonic_kdtree.query(s.distance_vector_field.row(i).data(), knn_num, out_index, out_distances_sq);
        for(int j=0; j < knn_num; ++j)
        {
            w += 1.0 / out_distances_sq[j];
            corresp_source_target_.row(i) += (1.0 / out_distances_sq[j]) * t.mVoxelPosition.row(out_index[j]);
        }
        corresp_source_target_.row(i) /= w;
    }

}		/* -----  end of function direct_correspondence  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  network_simplex
 *  Description:  
 * =====================================================================================
 */
void EMD::network_simplex (const VolumeObject &source, const VolumeObject &target)
{
    lemon::ListDigraph g;
    lemon::ListDigraph::NodeMap<int> supply(g);
    lemon::ListDigraph::ArcMap<double> cost(g);
    lemon::ListDigraph::ArcMap<int> flow(g);

    for(int i=0; i < source.voxel_num_; ++i)
    {
        lemon::ListDigraph::Node m = g.addNode();
        supply[m] = target.voxel_num_;
    }

    for(int i=0; i < target.voxel_num_; ++i)
    {
        lemon::ListDigraph::Node n = g.addNode();
        supply[n] = -1 * source.voxel_num_;
        for(lemon::ListDigraph::NodeIt s(g); s!=lemon::INVALID; ++s)
        {
            if(lemon::ListDigraph::id(s) < source.voxel_num_)
            {
                lemon::ListDigraph::Arc arc = g.addArc(s, n);
                cost[arc] = ( source.distance_vector_field.row(lemon::ListDigraph::id(s)) - target.distance_vector_field.row(i) ).norm();
            }
        }
    }

    lemon::NetworkSimplex<lemon::ListDigraph, int, double> ns(g);
    ns.costMap(cost).supplyMap(supply).run();
    ns.flowMap(flow);

    int num_flow = 2 * int(std::max(Real(source.voxel_num_) / Real(target.voxel_num_), Real(target.voxel_num_) / Real(source.voxel_num_))  * std::max(source.voxel_num_, target.voxel_num_) );

    flow_matrix_ = SpMat(source.voxel_num_, target.voxel_num_);
    std::vector<MyTriplet> flow_triplet;
    flow_triplet.reserve(num_flow);

    for(lemon::ListDigraph::ArcIt s(g); s != lemon::INVALID; ++s)
    {
        if(flow[s] > 10)
        {
            flow_triplet.push_back(MyTriplet(g.id(g.source(s)), g.id(g.target(s)) - source.voxel_num_, flow[s]/(Real)(source.voxel_num_*target.voxel_num_) ));
        }
        else if(flow[s] < 0)
        {
            std::cout<<"exists negative flow\n";
        }
    }

    flow_matrix_.setFromTriplets(flow_triplet.begin(), flow_triplet.end());

}		/* -----  end of function network_simplex  ----- */

