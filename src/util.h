/*
 * =====================================================================================
 *
 *       Filename:  util.h      Created:  08/09/2015 08:57:14 PM
 *
 *    Description:  helper function
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef AUTOMORPH_UTIL_H_
#define AUTOMORPH_UTIL_H_
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "def_types.h"

void matrix_to_point_cloud_file(const MatrixXr &points, std::string name);

struct NeighborIndex
{
    NeighborIndex();
    std::vector<Vector3i> neighbor18_index;
    std::vector<Vector3i> get_neighbor18_index();
};

#endif

