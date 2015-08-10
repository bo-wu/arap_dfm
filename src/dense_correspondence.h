/*
 * =====================================================================================
 *
 *       Filename:  dense_correspondence.h      Created:  07/29/2015 10:34:35 PM
 *
 *    Description:  find dense correspondence from sparse anchors
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef DENSE_CORRESPONDENCE_H_
#define DENSE_CORRESPONDENCE_H_
#include <vector>
#include "def_types.h"
#include "volume_object.h"

class EMD
{
public:
    EMD() {}

    void construct_correspondence(const VolumeObject &s, const VolumeObject &t);
    //correspondence from source voxel point to target (some point)
    MatrixX3r corresp_source_target_;
    //correspondence from target voxel point to source (some point)
    MatrixX3r corresp_target_source_;
//private:
    SpMat flow_matrix_;
    void find_correspondence(const VolumeObject &s, const VolumeObject &t);
    void min_cost_flow(const VolumeObject &s, const VolumeObject &t);
};
#endif

