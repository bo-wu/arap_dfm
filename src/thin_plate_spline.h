/*
 * =====================================================================================
 *
 *       Filename:  thin_plate_spline.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/21/2015 08:34:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "def_types.h"

template<typename Real>
class ThinPlateSpline
{
public:
    ThinPlateSpline(MatrixX3r & control_points, Real lambda=0.0);
    ~ThinPlateSpline();
    void compute_tps();
    // @param r: length of radius
    inline Real kernel(Real r2);
    void interplate(MatrixX3r&, MatrixX3r&);

private:
    Real m_lambda;
    MatrixX3r m_control_points;
};

