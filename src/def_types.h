/*
 * =====================================================================================
 *
 *       Filename:  types.h  Version:  1.0  Created:  01/16/2015 04:58:40 PM
 *
 *    Description:  types used in program
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */

#ifndef DEF_TYPES_H_
#define DEF_TYPES_H_
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

//use igl static lib
#define IGL_STATIC_LIBRARY

typedef double Real;

typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 3> MatrixX3r;
typedef Eigen::Matrix<Real, 3, Eigen::Dynamic> Matrix3Xr;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 4> MatrixX4r;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorXr;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<Real, 3, 1> Vector3r;
typedef Eigen::Matrix<Real, 1, Eigen::Dynamic> RowVectorXr;
typedef Eigen::SparseMatrix<Real> SpMat;
typedef Eigen::Triplet<Real> Triplet;


struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::VectorT<Real, 3> Point;	
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TriMesh;

#endif
