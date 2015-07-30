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
#include "find_dense_correspondence.h"
#include "full_bipartitegraph.h"
#include "network_simplex_simple.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compute_EMD
 *  Description:  compute earth movers' distance
 * =====================================================================================
 */
Real EMD::compute_EMD()
{
     return 0.0;
}		/* -----  end of function compute_EMD ----- */


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
    result_flow_.reserve(num_flow);
    for(int i=0; i < s.voxel_num_; ++i)
        for(int j=0; j < t.voxel_num_; ++j)
        {
            TsFlow f;
            f.amount = net.flow(di.arcFromId(i*t.voxel_num_ + j));
            f.from = i;
            f.to = j;
            if(fabs(f.amount) > 1e-8)
                result_flow_.push_back(f);
        }

}		/* -----  end of function min_cost_flow  ----- */
