/*
 * =====================================================================================
 *
 *       Filename:  skeleton.h      Created:  08/17/2015 03:32:05 PM
 *
 *    Description:  skeleton: read zhou yang's skeleton
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef AUTO_MORPH_SKELTON_
#define AUTO_MORPH_SKELTON_
#include <string>
#include <vector>
#include <utility>
#include "def_types.h"

class Skeleton
{

public:
    Skeleton(){}

    void init(std::string mesh_name);

    int num_branch;
    // not matched skeleton, merge to matched
    // <not_matched, matched>
    std::vector<MatrixX3r> skel_branch_points; // points on skeleton
    std::vector<std::pair<int, int> > merged_branch; //used for merge voxels
    // each face point and it belongs to which skeleton branch
    MatrixX3r mesh_face_points;
    VectorXi mesh_face_point_tag;
};


class SkeletonPair
{
public:
    SkeletonPair(){}
    void read_match_info(Skeleton &ss, Skeleton &ts, std::string name);
    std::vector<std::pair<int, int> > matched_branch;
    // correspondence skeleton sampled points
    std::vector<std::pair<Vector3r, Vector3r> > corresp_skel_points;
};

#endif

