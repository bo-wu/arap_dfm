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

VolumeObject::VolumeObject()
{
	if (mesh_name.empty())
	{
		std::cerr<<"mesh_name without initialization \n";
		exit(-1);
	}
	initial_volume();
}

VolumeObject::VolumeObject(std::string name)
{
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
	openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(0.1);
	grid->setTransform(grid_transform);
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setName("mesh_grid");
//	openvdb::tools::MeshToVolume<openvdb::FloatGrid> mesh2volume(grid_transform);
//	mesh2volume.convertToLeveSet(points, triangles);
//	grid = mesh2volume.distGridPtr();
	grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(grid->transform(), points, triangles, float(openvdb::LEVEL_SET_HALF_WIDTH));

    interior_grid = openvdb::tools::sdfInteriorMask(*grid);
    //make inside grid dense
    interior_grid->tree().voxelizeActiveTiles();
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

	int inactive=0, active=0, total=0;
	for(openvdb::FloatGrid::ValueOnCIter iter=grid->cbeginValueOn(); iter.test(); ++iter)	
	{
		active++;
	}
	for(openvdb::FloatGrid::ValueOffIter iter=grid->beginValueOff(); iter.test(); ++iter)
	{
		inactive++;
	}
	for(openvdb::FloatGrid::ValueAllIter iter=grid->beginValueAll(); iter.test(); ++iter)
	{
		total++;
	}
	std::cout<<"my active is "<< active <<"\ninactive is "<<inactive<<std::endl;
	std::cout<<"total num "<<total<<std::endl;
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
    auto anchorNum = mAnchors.size();
    mLaplaceMatrix = MatrixXr::Zero(voxelNum+anchorNum, voxelNum);
    mVoxelPosition = MatrixX3r::Zero(voxelNum, 3);
    int k = 0;
    std::ofstream output_coord("coord.txt");
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        output_coord << iter.getCoord();
        auto voxel_pos = grid->indexToWorld(iter.getCoord());
        output_coord <<" "<<voxel_pos<<std::endl;
        mVoxelPosition.row(k) << voxel_pos[0], voxel_pos[1], voxel_pos[2];
    }
    output_coord.close();
    //find neighbor voxel index
    kd_tree_type voxelKDTree(mVoxelPosition);
    voxelKDTree.index->buildIndex();
    const size_t numClosest = 1;
    size_t outIndex;
    Real outDistance;
    int degree;
    openvdb::Coord v_coord;
    // laplace matrix
    k = 0;
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        v_coord = iter.getCoord();
        degree = 0;
        std::vector<int> neighborIndex;
        for(int i=0; i < 3; ++i)
            for(int j=-1; j <= 1; j+=2)
        {
            auto temp_coord = v_coord;
            temp_coord[i] = v_coord[i] + 1*j;
            if (interior_grid->tree().isValueOn(temp_coord))
            {
                ++degree;
                auto voxel_pos = grid->indexToWorld(temp_coord);
                Vector3r vPos(voxel_pos[0], voxel_pos[1], voxel_pos[2]);
                voxelKDTree.query(vPos.data(), numClosest, &outIndex, &outDistance);
                if(outDistance > 1.0e-5)
                {
                    //std::cerr<<"Distance "<<outDistance<<" should be 0.0\n";
                }
                neighborIndex.push_back(outIndex);
            }
        }
        mLaplaceMatrix(k, k) = degree;
        for(auto idx : neighborIndex)
        {
            mLaplaceMatrix(k, idx) = -1;
        }
    }
    std::ofstream output_laplace("laplace.txt");
    output_laplace << mLaplaceMatrix;
    output_laplace.close();
}		/* -----  end of function construct_laplace_matrix  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calc_vector_field
 *  Description:  
 * =====================================================================================
 */
void VolumeObject::calc_vector_field()
{
    test_volume();
    construct_laplace_matrix();
    for(auto anchor : mAnchors)
    {
        
    }
}		/* -----  end of function calc_vector_field  ----- */



///////////////////////////////////////////
//read_mesh
void VolumeObject::read_mesh()
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
		mesh.point(*v_it) = (mesh.point(*v_it) - bb_center) * scalar_max;
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


