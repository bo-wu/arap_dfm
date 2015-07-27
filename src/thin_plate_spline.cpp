/*
 * =====================================================================================
 *
 *       Filename:  thin_plate_spline.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/22/2015 12:25:19 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <cmath>
#include "thin_plate_spline.h"

template<typename Real>
ThinPlateSpline<Real>::ThinPlateSpline(MatrixX3r& control_points, Real lambda)
        :
        m_lambda(lambda),
        m_control_points(control_points)
{

}

template<typename Real>
inline Real ThinPlateSpline<Real>::kernel(Real r)
{
    if(r > (Real)0)
    {
        return r * r * log(r);
    }
    else
        return (Real)0.0;
}

template<typename Real>
void ThinPlateSpline<Real>::compute_tps()
{
    auto v_num = m_control_points.rows();
    MatrixXr L = MatrixXr::Zero(v_num+4, v_num+4);
    Real r, alpha = 0.0;
    for(int i=0; i < v_num; ++i)
        for(int j=i+1; j < v_num; ++j)
        {
            r = (m_control_points.row(i) - m_control_points.row(j)).norm();
            L(i, j) = L(j, i) = kernel(r);
            alpha += r * 2;
        }
    alpha /= (Real)(v_num * v_num);
    auto regular = alpha * alpha * m_lambda;
    for(int i=0; i < v_num; ++i)
    {
        L(i, i) = regular;
    }
    L.block(0, v_num, v_num, 1) = VectorXr::Ones(v_num);
    L.block(v_num, 0, 1, v_num) = RowVectorXr::Ones(v_num);
    L.block(0, v_num+1, v_num, 3) = m_control_points;
    L.block(v_num+1, 0, 3, v_num) = m_control_points.transpose();

}

template<typename Real>
void ThinPlateSpline<Real>::interplate(MatrixX3r& current, MatrixX3r& next_step)
{

}

