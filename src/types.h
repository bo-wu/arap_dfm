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

#ifndef TYPES_H_
#define TYPES_H_
#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef double Real;

typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorXr;


struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::VectorT<Real, 3> Point;	
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TriMesh;

#endif
