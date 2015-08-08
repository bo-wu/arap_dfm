/*
 * =====================================================================================
 *
 *       Filename:  thin_plate_spline.h  Created:  07/21/2015 08:34:56 PM 
 *       
 *    Description:  
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef THIN_PLATE_SPLINE_H_
#define THIN_PLATE_SPLINE_H_

#include "def_types.h"

//template<typename Real>
class ThinPlateSpline
{
public:
    ThinPlateSpline(Real lambda=0.5);
    ~ThinPlateSpline() { }
    // control_points: input points
    // expected_positions: expect control points to be
    void compute_tps(const MatrixX3r &control_points, const MatrixX3r &expected_positions);
    // @param r: length of radius
    inline Real kernel(Real r2);
    void interpolate(const MatrixX3r &input, MatrixX3r &output);

//private:
    Real m_lambda;
    MatrixX3r mControlPoints;
    // size (v_num + 4) * 3
    // first v_num for control point and the last 4 for regular term
    // * 3 for x, y, z
    MatrixX3r mCoeff;
};

#endif
