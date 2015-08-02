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
#include <Eigen/Dense>
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

void VolumeObject::set_anchors(std::vector<Vector3r>& anchors)
{
    mAnchors = anchors;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  construct_laplace_matrix
 *  Description:  
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

    //find neighbor voxel index
    kd_tree_type voxelKDTree(3, mVoxelPosition);
    voxelKDTree.index->buildIndex();
    long int outIndex;
    Real outDistance;
    int degree;
    openvdb::Coord v_coord;
    Vector3r v_world_pos;
    std::vector<MyTriplet> laplace_triplet_list;
    laplace_triplet_list.reserve(7*voxelNum);
    // laplace matrix
    k = 0;
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        v_coord = iter.getCoord();
        degree = 0;
        for(int i=0; i < 3; ++i)
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
                }
            }
        laplace_triplet_list.push_back(MyTriplet(k, k, degree));
    }
    mLaplaceMatrix.setFromTriplets(laplace_triplet_list.begin(), laplace_triplet_list.end());
    // constrain part
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
//    std::ofstream output_laplace("laplace.dat");
//    output_laplace << mLaplaceMatrix;
//    output_laplace.close();
    igl::min_quad_with_fixed_data<Real> mqwf;
    int num_row = constraint_index_.rows();
    VectorXr B = VectorXr::Zero(voxel_num_, 1);
    //Empyty constraints (except for constraint_index_/value)
    SpMat Aeq;
    VectorXr Beq;
    igl::min_quad_with_fixed_precompute(mLaplaceMatrix, constraint_index_, Aeq, true, mqwf);
    VectorXr D;
//    std::ofstream output_constraint_index("constraint_index.dat");
//    output_constraint_index << constraint_index_;
//    output_constraint_index.close();
//    
    // solve equation with constraint
    for(int i=0; i < num_row; ++i)
    {
        VectorXr constraint_value = VectorXr::Zero(num_row, 1);
        constraint_value(i) = 1.0;
        igl::min_quad_with_fixed_solve(mqwf, B, constraint_value, Beq, D);
        distance_vector_field.col(i) = D;
    }
}		/* -----  end of function calc_vector_field  ----- */


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

