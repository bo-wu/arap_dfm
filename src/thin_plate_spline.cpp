/*
 * =====================================================================================
 *
 *       Filename:  thin_plate_spline.cpp    Created:  07/22/2015 12:25:19 AM
 *
 *    Description:  
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include "thin_plate_spline.h"

//template<typename Real>
//ThinPlateSpline<Real>::ThinPlateSpline(Real lambda)
ThinPlateSpline::ThinPlateSpline(Real lambda)
        : m_lambda(lambda)
{
}

//template<typename Real>
//inline Real ThinPlateSpline<Real>::kernel(Real r)
inline Real ThinPlateSpline::kernel(Real r)
{
    if(r > (Real)1e-05)
    {
        return r * r * log(r);
    }
    else
        return (Real)0.0;
}

//template<typename Real>
//void ThinPlateSpline<Real>::compute_tps(MatrixX3r &control_points, MatrixX3r &expected_positions)
void ThinPlateSpline::compute_tps(const MatrixX3r &control_points, const MatrixX3r &expected_positions)
{
    assert( control_points.rows() == expected_positions.rows() );
    mControlPoints = control_points;

    auto v_num = control_points.rows();
    mCoeff = MatrixX3r::Zero(v_num+4, 3);
    MatrixXr L = MatrixXr::Zero(v_num+4, v_num+4);
    MatrixX3r B = MatrixX3r::Zero(v_num+4, 3);
    B.block(0,0, v_num, 3) = expected_positions;

//#pragma omp parallel 
{
    Real r, alpha = 0.0;
    for(int i=0; i < v_num; ++i)
    {
        for(int j=i+1; j < v_num; ++j)
        {
            r = (control_points.row(i) - control_points.row(j)).norm();
            L(i, j) = L(j, i) = kernel(r);
            alpha += r * 2;
        }
    }
    alpha /= (Real)(v_num * v_num);
    auto smooth_term = alpha * alpha * m_lambda;
    for(int i=0; i < v_num; ++i)
    {
        L(i, i) = smooth_term;
    }
}
    L.block(0, v_num, v_num, 1) = VectorXr::Ones(v_num);
    L.block(v_num, 0, 1, v_num) = RowVectorXr::Ones(v_num);
    L.block(0, v_num+1, v_num, 3) = control_points;
    L.block(v_num+1, 0, 3, v_num) = control_points.transpose();

    /*
    //should NOT use svd (it is NOT least square !!!)
    // solve L*x = B, B=[b1, b2, b3] , SVD slow, robust
    Eigen::JacobiSVD<MatrixXr> svd(L, Eigen::ComputeThinU | Eigen::ComputeThinV);
    for(int i=0; i < 3; ++i)
    {
        mCoeff.col(i) = svd.solve(B.col(i));
    }
    */

    //LDLT, faster 
    Eigen::LDLT<MatrixXr> ldlt;
    ldlt.compute(L);
    if(ldlt.info() != Eigen::Success)
    {
        std::cerr << "thin plate LDLT solver error\n";
        exit(-1);
    }
    for(int i=0; i < 3; ++i)
    {
        mCoeff.col(i) = ldlt.solve(B.col(i));
    }

    m_L = L;
    /*  
    std::ofstream output_solve_error ("solve_error.dat");
    output_solve_error << L * mCoeff - B;
    output_solve_error.close();
    */
}

//template<typename Real>
//void ThinPlateSpline<Real>::interplate(MatrixX3r &input, MatrixX3r &output)
void ThinPlateSpline::interplate(const MatrixX3r &input, MatrixX3r &output)
{
    int input_rows = input.rows();
    int control_num = mControlPoints.rows();
    MatrixXr L_matrix(input_rows, control_num+4);

#pragma omp parallel for 
    for(int i=0; i < input_rows; ++i)
    {
        for(int j=0; j < control_num; ++j)
        {
            L_matrix(i, j) = kernel( (mControlPoints.row(j) - input.row(i)).norm() );
        }
    }

    L_matrix.block(0,control_num, input_rows, 1) = VectorXr::Ones(input_rows);
    L_matrix.block(0,control_num+1, input_rows, 3) = input;
    output = L_matrix * mCoeff;

}

