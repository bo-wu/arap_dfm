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
#include <algorithm>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <openvdb/tools/MeshToVolume.h>
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
	openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(0.005);
	grid->setTransform(grid_transform);
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setName("mesh_grid");
//	openvdb::tools::MeshToVolume<openvdb::FloatGrid> mesh2volume(grid_transform);
//	mesh2volume.convertToLeveSet(points, triangles);
//	grid = mesh2volume.distGridPtr();
	grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(grid->transform(), points, triangles, float(openvdb::LEVEL_SET_HALF_WIDTH));


} //end of initial_volume


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
		for(int i=0; i<3; ++i)
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
	if(file.hasBloscCompression())
	{
		std::cout<<"my openvdb has blosc compression support\n";
	}
	std::cout<<"default compression flags "<<file.DEFAULT_COMPRESSION_FLAGS<<std::endl;
	openvdb::GridPtrVec grids;
	file.setCompression(openvdb::io::COMPRESS_ZIP|openvdb::io::COMPRESS_ACTIVE_MASK);
	grids.push_back(grid);
	file.write(grids);
	file.close();
}
