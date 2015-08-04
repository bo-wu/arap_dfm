/*
 * =====================================================================================
 *
 *       Filename:  volume_object.cpp  Version:  1.0  Created:  01/16/2015 04:58:12 PM
 *
 *    Description:  volume representation of data
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <openvdb/Grid.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetUtil.h>
#include "volume_object.h"

//#define MIN_QUAD_WITH_FIXED_CPP_DEBUG
#include <igl/min_quad_with_fixed.h>

VolumeObject::VolumeObject(Real transform_scale)
{
    transform_scale_ = transform_scale;
	if (mesh_name.empty())
	{
		std::cerr<<"mesh_name without initialization \n";
		exit(-1);
	}
	initial_volume();
}

VolumeObject::VolumeObject(std::string name, Real transform_scale)
{
    transform_scale_ = transform_scale;
	mesh_name = name;
	initial_volume();
}

VolumeObject::~VolumeObject()
{
}


///////////////////////////////////////////
//generate volume from mesh
void VolumeObject::initial_volume()
{
	read_mesh();
	grid = openvdb::FloatGrid::create(10.0);
	openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(transform_scale_);
	grid->setTransform(grid_transform);
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setName("mesh_grid");
	grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(grid->transform(), points, triangles, float(openvdb::LEVEL_SET_HALF_WIDTH));

    //std::cout<<"before subdivide " << grid->tree().activeTileCount()<<std::endl;
    interior_grid = openvdb::tools::sdfInteriorMask(*grid);
    //make inside grid dense
    interior_grid->tree().voxelizeActiveTiles();
    //std::cout<<"after subdivide " << grid->tree().activeTileCount()<<std::endl;
} //end of initial_volume


void VolumeObject::test_volume()
{
    auto inside_grid = openvdb::tools::sdfInteriorMask(*grid);
	std::cout<<"leaf num "<<grid->tree().leafCount()<<std::endl;
	std::cout<<"inside leaf num "<<inside_grid->tree().leafCount()<<std::endl;

	std::cout<<"       active leaf voxel "<<grid->tree().activeLeafVoxelCount()<<" inactive leaf voxel "<<grid->tree().inactiveLeafVoxelCount()<<"\n";
	std::cout<<"inside_grid active leaf voxel "<<inside_grid->tree().activeLeafVoxelCount()<<" inactive leaf voxel "<<inside_grid->tree().inactiveLeafVoxelCount()<<"\n";
    std::cout<<"            active tile num "<<grid->tree().activeTileCount()<<std::endl;
    std::cout<<"inside_grid active tile num "<<inside_grid->tree().activeTileCount()<<std::endl;
    std::ofstream output_depth("depth.txt");
    auto accessor = grid->getAccessor();
    for(auto iter=inside_grid->cbeginValueOn(); iter; ++iter)
    {
        output_depth<<iter.getVoxelCount()<<" ";
        if(iter.isTileValue())
            output_depth<<"tile value ";
        if(accessor.getValue(iter.getCoord()) > 0.0)
            output_depth<<"bigger than 0 ";
        output_depth<< accessor.getValue(iter.getCoord())<<" coord ";
        output_depth<< iter.getCoord()<<" world "<<grid->indexToWorld(iter.getCoord())<<std::endl;
    }
    output_depth.close();
    /*  
    inside_grid->tree().voxelizeActiveTiles();
    std::cout<<"\nafter voxelize active tiles\n\n";

	std::cout<<"leaf num "<<grid->tree().leafCount()<<std::endl;
	std::cout<<"inside leaf num "<<inside_grid->tree().leafCount()<<std::endl;

	std::cout<<"       active leaf voxel "<<grid->tree().activeLeafVoxelCount()<<" inactive leaf voxel "<<grid->tree().inactiveLeafVoxelCount()<<"\n";
	std::cout<<"inside_grid active leaf voxel "<<inside_grid->tree().activeLeafVoxelCount()<<" inactive leaf voxel "<<inside_grid->tree().inactiveLeafVoxelCount()<<"\n";
    std::cout<<"            active tile num "<<grid->tree().activeTileCount()<<std::endl;
    std::cout<<"inside_grid active tile num "<<inside_grid->tree().activeTileCount()<<std::endl;
		*/
}

void VolumeObject::set_anchors(const std::vector<Vector3r>& anchors)
{
    mAnchors = anchors;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  construct_laplace_matrix
 *  Description:  laplace matrix for distance vector field (sparse grid) 
 * =====================================================================================
 */
void VolumeObject::construct_laplace_matrix()
{
    //for sparse grid, should be activeVoxel + activeTile
    auto voxelNum = interior_grid->tree().activeLeafVoxelCount();
    voxel_num_ = voxelNum;
    auto anchorNum = mAnchors.size();
    mLaplaceMatrix = SpMat(voxelNum, voxelNum);
    mVoxelPosition = MatrixX3r::Zero(voxelNum, 3);
    distance_vector_field = MatrixXr::Zero(voxelNum, anchorNum);
    int k = 0;
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        auto voxel_pos = grid->indexToWorld(iter.getCoord());
        mVoxelPosition.row(k) << voxel_pos[0], voxel_pos[1], voxel_pos[2];
    }
    // temporal for testing, need to change 
    mDenseVoxelPosition = mVoxelPosition;

    //find neighbor voxel index
    kd_tree_type voxelKDTree(3, mVoxelPosition);
    voxelKDTree.index->buildIndex();
    long int outIndex;
    Real outDistance;
    int degree;
    openvdb::Coord v_coord;
    Vector3r v_world_pos;
    // construct tetrahedron
    std::vector< std::vector<int> > neighbor_index_3d; //3 directions
    std::vector<int> neighbor_index_1d; // 1 direction

    std::vector<MyTriplet> laplace_triplet_list;
    laplace_triplet_list.reserve(7*voxelNum);
    k = 0;
    // laplace matrix
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        neighbor_index_3d.clear();
        v_coord = iter.getCoord();
        degree = 0;
        for(int i=0; i < 3; ++i)
        {
            neighbor_index_1d.clear();
            for(int j=-1; j <= 1; j+=2)
            {
                auto temp_coord = v_coord;
                temp_coord[i] = v_coord[i] + 1*j;
                if (interior_grid->tree().isValueOn(temp_coord))
                {
                    ++degree;
                    auto voxel_pos = grid->indexToWorld(temp_coord);
                    v_world_pos<< voxel_pos[0], voxel_pos[1], voxel_pos[2];
                    voxelKDTree.query(v_world_pos.data(), 1, &outIndex, &outDistance);
                    if(outDistance > 1.0e-5)
                    {
                        std::cerr<<"Distance "<<outDistance<<" should be 0.0\n";
                    }
                    laplace_triplet_list.push_back(MyTriplet(k, outIndex, -1));
                    neighbor_index_1d.push_back(outIndex);
                }
            }
            if (neighbor_index_1d.size() > 0)
            {
                neighbor_index_3d.push_back(neighbor_index_1d);
            }
        }
        laplace_triplet_list.push_back(MyTriplet(k, k, degree));

        // construct tetrahedron
        if(neighbor_index_3d.size() == 3)
        {
            Vector4i tet_index;
            tet_index(0) = k;
            for(int l=0; l < neighbor_index_3d[0].size(); ++l)
                for(int m=0; m < neighbor_index_3d[1].size(); ++m)
                    for(int n=0; n < neighbor_index_3d[2].size(); ++n)
                    {
                        tet_index(1) = neighbor_index_3d[0].at(l);
                        tet_index(2) = neighbor_index_3d[1].at(m);
                        tet_index(3) = neighbor_index_3d[2].at(n);
                        mTetIndex.push_back(tet_index);
                    }
        }
    }
    mLaplaceMatrix.setFromTriplets(laplace_triplet_list.begin(), laplace_triplet_list.end());
    // constraint part
    if (anchorNum > 0)
        constraint_index_ = VectorXi::Zero(anchorNum, 1);
    k = 0;
    for (auto a : mAnchors)
    {
        voxelKDTree.query(a.data(), 1, &outIndex, &outDistance);
        constraint_index_(k) = outIndex;
        ++k;
    }
}		/* -----  end of function construct_laplace_matrix  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calc_vector_field
 *  Description:  
 * =====================================================================================
 */
void VolumeObject::calc_vector_field()
{
//    test_volume();
    construct_laplace_matrix();
    /*
    std::ofstream output_laplace("laplace.dat");
    output_laplace << mLaplaceMatrix;
    output_laplace.close();
    */
    igl::min_quad_with_fixed_data<Real> mqwf;
    int num_row = constraint_index_.rows();
    VectorXr B = VectorXr::Zero(voxel_num_, 1);
    //Empyty constraints (except for constraint_index_/value)
    SpMat Aeq;
    VectorXr Beq;
    igl::min_quad_with_fixed_precompute(mLaplaceMatrix, constraint_index_, Aeq, true, mqwf);
    VectorXr D;
    /*  
    std::ofstream output_constraint_index("constraint_index.dat");
    output_constraint_index << constraint_index_;
    output_constraint_index.close();
    */
    // solve equation with constraint
    for(int i=0; i < num_row; ++i)
    {
        VectorXr constraint_value = VectorXr::Zero(num_row, 1);
        constraint_value(i) = 1.0;
        igl::min_quad_with_fixed_solve(mqwf, B, constraint_value, Beq, D);
        distance_vector_field.col(i) = D;
    }
}		/* -----  end of function calc_vector_field  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  polar_decompose
 *  Description:  
 * =====================================================================================
 */
void VolumeObject::polar_decompose(const Matrix3r &rest, const Matrix3r &deform, Matrix3r &R, Matrix3r &S)
{
    Matrix3r defrom_gradient = deform * rest.inverse();
    Eigen::JacobiSVD<Matrix3r> svd(defrom_gradient, Eigen::ComputeFullU|Eigen::ComputeFullV);
    Matrix3r eigen_values = Matrix3r::Zero();
    Vector3r singular = svd.singularValues();
    for(int i=0; i < 3; ++i)
        eigen_values(i,i) = singular(i);

    R = svd.matrixU() * svd.matrixV().transpose();
    S = svd.matrixV() * eigen_values * svd.matrixV().transpose();
}		/* -----  end of function polar_decompose  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calc_tetrahedron_transform
 *  Description:  fill  mTetTransform
 *  @param final_corresp_points: correspondence after warping  
 * =====================================================================================
 */
void VolumeObject::calc_tetrahedron_transform(const MatrixX3r &final_corresp_points)
{
    assert( mDenseVoxelPosition.rows() == final_corresp_points.rows() );
    Matrix3r R, S, rest, deform;
    Vector3r rest_v0, deform_v0;
    Vector4i tet_vertex_index;
    for(int i=0; i < mTetIndex.size(); ++i)
    {
        tet_vertex_index = mTetIndex[i];
        rest_v0 = mDenseVoxelPosition.row( tet_vertex_index(0) );
        deform_v0 = final_corresp_points.row(tet_vertex_index(0) );
        for(int j=1; j < 4; ++j)
        {
            rest.row(j-1) = mDenseVoxelPosition.row(j) - rest_v0;
            deform.row(j-1) = final_corresp_points.row(j) - deform_v0;
        }
        polar_decompose(rest, deform, R, S);

        mTetTransform.push_back(std::make_pair(R, S));
    }
}		/* -----  end of function calc_tetrahedron_transform  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_intermedium_points
 *  Description:  find intermedium points at time t [0,1]
 *  @param inter_corresp_points: correspondence by interpolating warping
 * =====================================================================================
 */
void VolumeObject::find_intermedium_points(MatrixX3r &inter_corresp_points, const Real t)
{
    int anchor_num = 1;
    // with one anchor point
    SpMat L(inter_corresp_points.rows()+anchor_num, inter_corresp_points.rows());
    MatrixX3r B(inter_corresp_points.rows()+anchor_num, 3);

    std::vector<MyTriplet> tet_triplet_list;
    int tet_num = mTetIndex.size();
    tet_triplet_list.reserve(4*tet_num);
    for(int i=0; i < tet_num; ++i)
    {
        //construct B
        Quaternionr quat(mTetTransform[i].first);
        //construct L, B
        tet_triplet_list.push_back(MyTriplet(mTetIndex[i](0), mTetIndex[i](0), 3));
        for(int j=1; j < 4; ++j)
        {
            tet_triplet_list.push_back(MyTriplet(mTetIndex[i](0), mTetIndex[i](j), -1));

        }
    }

}		/* -----  end of function find_intermedium_points  ----- */


///////////////////////////////////////////
//read_mesh
void VolumeObject::read_mesh(bool resize)
{
	if(!OpenMesh::IO::read_mesh(mesh, mesh_name))	
	{
		std::cerr<<"read mesh "<<mesh_name<<" error!\n";
		exit(-1);
	}

	// normalize the mesh
	OpenMesh::VertexHandle v_h = mesh.vertex_handle(0);
	MyTraits::Point bb_center, bb_min, bb_max;
	bb_min = bb_max = OpenMesh::vector_cast<MyTraits::Point>(mesh.point(v_h));
	TriMesh::VertexIter v_it ;
	//get bounding box
	for(v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
	{
		bb_min.minimize(OpenMesh::vector_cast<MyTraits::Point>(mesh.point(*v_it)));
		bb_max.maximize(OpenMesh::vector_cast<MyTraits::Point>(mesh.point(*v_it)));
	}
	bb_center = (bb_max + bb_min) / 2.0;
	bb_max = bb_max - bb_min;
	Real scalar_max = std::max(bb_max[0], bb_max[1]);
	scalar_max = std::max(scalar_max, bb_max[2]);
	scalar_max = 1.0 / scalar_max;
	openvdb::Vec3s v_pos;
	for(v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
	{
        if(resize)
        {
            mesh.point(*v_it) = (mesh.point(*v_it) - bb_center) * scalar_max;
        }
		//fill points
		for(int i=0; i < 3; ++i)
		{
			v_pos(i) = mesh.point(*v_it)[i];
		}
		points.push_back(v_pos);
	}

	//fill triangles
	TriMesh::FaceIter f_it;
	TriMesh::FaceVertexIter fv_it;
	openvdb::Vec3I fv_index;
	int counter;
	for(f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
	{
		counter = 0;
		for(fv_it=mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it, ++counter)
		{
			if(counter >= 3)
			{
				std::cerr<<"more that 3 points, not triangle mesh \n";
				exit(-1);
			}
			fv_index(counter) = fv_it->idx();
		}
		triangles.push_back(fv_index);
	}

}

////////////////////////////////////////////
//write grid
void VolumeObject::write_grid(std::string name)
{
	if(name.substr(name.length()-4, 4).compare(".vdb") != 0)
	{
		name = name + ".vdb";
	}
	openvdb::io::File file(name);
    /*  
	if(file.hasBloscCompression())
	{
		std::cout<<"my openvdb has blosc compression support\n";
	}
	std::cout<<"default compression flags "<<file.DEFAULT_COMPRESSION_FLAGS<<std::endl;
    */
	openvdb::GridPtrVec grids;
	file.setCompression(openvdb::io::COMPRESS_ZIP|openvdb::io::COMPRESS_ACTIVE_MASK);
	grids.push_back(grid);
	file.write(grids);
	file.close();
}

