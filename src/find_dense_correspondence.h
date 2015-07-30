/*
 * =====================================================================================
 *
 *       Filename:  find_dense_correspondence.h      Created:  07/29/2015 10:34:35 PM
 *
 *    Description:  find dense correspondence from sparse anchors
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef FIND_DENSE_CORRESPONDENCE_H_
#define FIND_DENSE_CORRESPONDENCE_H_
#include <vector>
#include "def_types.h"
#include "volume_object.h"
class EMD
{
public:
    EMD() {}
    Real compute_EMD();
    void min_cost_flow(VolumeObject &s, VolumeObject &t);
    SpMat flow_matrix_;

};
#endif

