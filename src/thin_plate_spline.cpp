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
#include <ctime>
#include <armadillo>
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
    if(r == (Real)0.0)
    {
        return (Real)0.0;
    }
    else
        return r * r * log(r);
}

//template<typename Real>
//void ThinPlateSpline<Real>::compute_tps(MatrixX3r &control_points, MatrixX3r &expected_positions)
void ThinPlateSpline::compute_tps(const MatrixX3r &control_points, const MatrixX3r &expected_positions)
{
    assert( control_points.rows() > 3 );
    assert( control_points.rows() == expected_positions.rows() );
    mControlPoints = control_points;

    auto v_num = control_points.rows();
    mCoeff = MatrixX3r::Zero(v_num+4, 3);
    MatrixXr L = MatrixXr::Zero(v_num+4, v_num+4);
    MatrixX3r B = MatrixX3r::Zero(v_num+4, 3);
    B.block(0,0, v_num, 3) = expected_positions;

#ifdef PARALLEL_OMP_
#pragma omp parallel 
#endif
    {

    Real r, alpha = 0.0;
#ifdef PARALLEL_OMP_
#pragma omp for
#endif
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
#ifdef PARALLEL_OMP_
#pragma omp for
#endif
    for(int i=0; i < v_num; ++i)
    {
        L(i, i) = smooth_term;
    }

    }

    L.block(0, v_num, v_num, 1) = VectorXr::Ones(v_num);
    L.block(v_num, 0, 1, v_num) = RowVectorXr::Ones(v_num);
    L.block(0, v_num+1, v_num, 3) = control_points;
    L.block(v_num+1, 0, 3, v_num) = control_points.transpose();

    arma::mat arma_L = arma::mat(v_num+4, v_num+4);
    arma::mat arma_B = arma::mat(v_num+4, 3);
#ifdef PARALLEL_OMP_
#pragma omp parallel for
#endif
    for(int i=0; i < v_num+4; ++i)
    {
        for(int j=0; j < v_num+4; ++j)
        {
            arma_L(i, j) = L(i, j);
        }
        for(int k=0; k < 3; ++k)
        {
            arma_B(i, k) = B(i, k);
        }
    }

    /*  
    //LDLT, faster 
    //Eigen::LDLT<MatrixXr> ldlt; //wrong
    //Eigen::FullPivLU<MatrixXr> ldlt;
    Eigen::ColPivHouseholderQR<MatrixXr> ldlt;
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
    */

    // fastest solver
    std::clock_t start;
    start = std::clock();
    arma::mat X = arma::solve(arma_L, arma_B);
    Real elapse = (std::clock() - start) / (Real)(CLOCKS_PER_SEC);
    std::cout <<"solving equation elapse " << elapse <<std::endl;

#ifdef PARALLEL_OMP_
#pragma omp parallel for
#endif
    for(int i=0; i < v_num+4; ++i)
        for(int j=0; j < 3; ++j)
        {
            mCoeff(i, j) = X(i, j);
        }

    /*  
    std::ofstream output_coeff("coefficient.dat");
    output_coeff << mCoeff;
    output_coeff.close();
    std::ofstream output_solve_error ("matrix_L.dat");
    output_solve_error << L;
    output_solve_error.close();
    std::ofstream output_B("B.dat");
    output_B << B;
    output_B.close();
    */
}

//template<typename Real>
//void ThinPlateSpline<Real>::interplate(MatrixX3r &input, MatrixX3r &output)
void ThinPlateSpline::interplate(const MatrixX3r &input, MatrixX3r &output)
{
    int input_rows = input.rows();
    int control_num = mControlPoints.rows();
    MatrixXr L_matrix(input_rows, control_num+4);
#ifdef PARALLEL_OMP_
#pragma omp parallel for 
#endif
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

    /*
    std::ofstream output_coeff("coefficient.dat");
    output_coeff << mCoeff;
    output_coeff.close();
    std::ofstream output_result("l_matrix_result.dat");
    output_result << output;
    output_result.close();
    std::ofstream output_ml_result("mL_matrix_result.dat");
    output_ml_result << m_L.block(0,0, mControlPoints.rows(), mControlPoints.rows()+4) * mCoeff;
    output_ml_result.close();
    */

    /*  
    std::ofstream output_L_matrix("input_Lmatrix.dat");
    output_L_matrix << L_matrix;
    output_L_matrix.close();
    std::ofstream output_mL("m_L.dat");
    output_mL << m_L.block(0,0, mControlPoints.rows(), mControlPoints.rows()+4);
    output_mL.close();
    std::ofstream output_diffL("diff_L.dat");
    output_diffL << L_matrix - m_L.block(0,0, mControlPoints.rows(), mControlPoints.rows()+4);
    output_diffL.close();
    */
}

