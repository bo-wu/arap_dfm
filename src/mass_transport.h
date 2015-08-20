/*
 * =====================================================================================
 *
 *       Filename:  mass_transport.h      Created:  08/19/2015 11:21:29 PM
 *
 *    Description:  
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef MASS_TRANSPORT_H_
#define MASS_TRANSPORT_H_
#include "def_types.h"
#include "volume_object.h"
#include "skeleton.h"

class PartMassTransport
{
public:
    PartMassTransport() {}
    void construct_correspondence(const VolumeObject &s, const VolumeObject &t);
    void network_simplex(const VolumeObject &s, const VolumeObject &t, SkeletonPair &sp);
    void part_network_simplex(const VolumeObject &s, const VolumeObject &t, int source_part_index, int target_part_index);

    void find_correspondence(const VolumeObject &s, const VolumeObject &t, const Real threshold=0.001);
    
    MatrixX3r source_control_points_;
    MatrixX3r corresp_source_target_;

    MatrixX3r target_control_points_;
    MatrixX3r corresp_target_source_;

    SpMat flow_matrix_;
    std::vector<MyTriplet> total_flow_triplet_;
};

#endif

