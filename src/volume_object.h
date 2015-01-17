/*
 * =====================================================================================
 *
 *       Filename:  volume_object.h  Version:  1.0  Created:  01/16/2015 04:57:56 PM
 *
 *    Description:  volume representation
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#ifndef VOLUME_OBJECT_H_
#define VOLUME_OBJECT_H_
#include <string>
#include <openvdb/openvdb.h>
#include "types.h"

struct VolumeObject
{
	TriMesh mesh;
	std::string mesh_name;
	openvdb::FloatGrid::Ptr grid;
	std::vector<openvdb::Vec3s> points;
	std::vector<openvdb::Vec3I> triangles;

	VolumeObject();
	VolumeObject(std::string name);
	void initial_volume();
	void read_mesh();
	void write_grid(std::string name);
	~VolumeObject();

};


#endif

