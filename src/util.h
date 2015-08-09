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
#include "def_types.h"

void matrix_to_point_cloud(const MatrixXr &points, std::string name)
{
    name = name + ".obj";
    std::ofstream output_point(name.c_str());
    int p_num = points.rows();

    for(int i=0; i < p_num; ++i)
    {
        output_point << "v "<<points.row(i)<<std::endl;
    }
    output_point.close();
}

#endif

