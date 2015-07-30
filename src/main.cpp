/*
 * =====================================================================================
 *
 *       Filename:  main.cpp  Version:  1.0  Created:  01/16/2015 12:51:08 AM
 *
 *    Description:  main function
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "def_types.h"
#include "volume_object.h"
#include "find_dense_correspondence.h"

int main(int argc, char** argv)
{
    if (argc <= 2)
    {
        std::cerr << "Not enough parameters, should run like\n ./morph source_mesh_name target_mesh_name \n";
        return 0;
    }
    // read positions from skeleton correspondence
     
    std::cout << "initializing source volume ... \n";
    // transform 0.02
    VolumeObject source_volume(argv[1], 0.02);
    Vector3r source_anchor1(-0.16, -0.4, -0.42);
    Vector3r source_anchor2(0.16, 0.08, 0.1);

    //transfrom 0.01
//    VolumeObject source_volume(argv[1]);
//    Vector3r source_anchor1(-0.09, -0.35, -0.33);
//    Vector3r source_anchor2(0.09, 0.16, 0.16);

    source_volume.mAnchors.push_back(source_anchor1);
    source_volume.mAnchors.push_back(source_anchor2);
	source_volume.calc_vector_field();
	source_volume.write_grid(argv[1]);

    std::cout << "initializing target volume ... \n";
    //transfrom 0.02
    VolumeObject target_volume(argv[2], 0.02);
    Vector3r target_anchor1(-0.12, -0.48, -0.06);
    Vector3r target_anchor2(0.0, 0.48, 0.12);
    
    // transform 0.01
//    VolumeObject target_volume(argv[2]);
//    Vector3r target_anchor1(-0.13, -0.49, -0.04);
//    Vector3r target_anchor2(0.09, 0.33, 0.24);

    target_volume.mAnchors.push_back(target_anchor1);
    target_volume.mAnchors.push_back(target_anchor2);
    target_volume.calc_vector_field();
    target_volume.write_grid(argv[2]);

    std::cout << "find dense correspondence between source and target ... \n";
    EMD emd_flow;
    emd_flow.min_cost_flow(source_volume, target_volume);
    std::cout << "all done !\n";

    /*
    std::ofstream output_dist_vector("distance_field.dat");
    output_dist_vector << source_volume.distance_vector_field;
    output_dist_vector.close();
    std::ofstream output_laplace ("laplace.dat");
    output_laplace << source_volume.mLaplaceMatrix;
    output_laplace.close();
    */
    /*
    std::ofstream output_voxel_pos("source_voxel_pos.dat");
    output_voxel_pos << source_volume.mVoxelPosition;
    output_voxel_pos.close();
    std::ofstream output_target_voxel("target_voxel_pos.dat");
    output_target_voxel << target_volume.mVoxelPosition;
    output_target_voxel.close();
    */
	return 0;
}
