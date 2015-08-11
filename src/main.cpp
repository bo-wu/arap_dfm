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
#include "dense_correspondence.h"
#include "thin_plate_spline.h"
#include "morph.h"

int main(int argc, char** argv)
{
    if (argc <= 2)
    {
        std::cerr << "Not enough parameters, should run like\n ./morph source_mesh_name target_mesh_name \n";
        return 0;
    }
    // read positions from skeleton correspondence
    // head, left front, right front, left behind, right behind, tail
    std::vector<std::pair<Vector3r, Vector3r> > corresp_pairs;
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.015, 0.29, 0.5), Vector3r(-0.2, 0.301, 0.248)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.127, -0.41, 0.157), Vector3r(-0.098, -0.487, 0.2)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.131, -0.407, 0.158), Vector3r(0.096, -0.48, 0.439)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.139, -0.387, -0.377), Vector3r(-0.1168, -0.48, 0.001)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.138, -0.39, -0.376), Vector3r(0.1514, -0.468, -0.3248)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.00, 0.139, -0.24), Vector3r(0.0449, 0.1, -0.304)));
    // brest, stamock(front, behind)
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.0248, 0.06, 0.2), Vector3r(-0.003, 0.008, 0.191)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.008, 0.08, 0.11), Vector3r(0.196, 0.0, -0.105)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.0134, 0.14, -0.153), Vector3r(0.036, 0.09, -0.22)));
    // leg top, left front, right front, left behind, right behind
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.13, -0.068, 0.16), Vector3r(-0.0174, -0.168, 0.16)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.106, -0.067, 0.17), Vector3r(0.0857, -0.1, 0.24)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.12, 0.017, -0.185), Vector3r(-0.064, -0.14, -0.135)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.11, 0.038, -0.177), Vector3r(0.134, -0.17, -0.23)));
    //ear left, right
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.06, 0.3815, 0.327), Vector3r(-0.03, 0.487, 0.115)));
//    corresp_pairs.push_back(std::make_pair(Vector3r(0.08, 0.4, 0.3), Vector3r(0.018, 0.477, 0.3)));

    Morph morph(argv[1], argv[2], corresp_pairs, 0.025);
    morph.initial();
    morph.start_morph(0.1);
    morph.write_sequence();
    
/*
    Real voxel_size = 0.02;
    std::cout << "initializing source volume ... ";
    // transform 0.02
    VolumeObject source_volume(argv[1], voxel_size);
    Vector3r source_anchor1(-0.1, -0.35, -0.42);
    Vector3r source_anchor2(0.10, 0.08, 0.0);

    //transfrom 0.01
//    VolumeObject source_volume(argv[1]);
//    Vector3r source_anchor1(-0.09, -0.35, -0.33);
//    Vector3r source_anchor2(0.09, 0.16, 0.16);

    source_volume.mAnchors.push_back(source_anchor1);
    source_volume.mAnchors.push_back(source_anchor2);
	source_volume.calc_vector_field();
	source_volume.write_grid(argv[1]);

    std::cout << "initializing target volume ... ";
    //transfrom 0.02
    VolumeObject target_volume(argv[2], voxel_size);
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
    std::cout<< "done!\n";
    std::cout<<"source voxel num " << source_volume.mVoxelPosition.rows()<<std::endl;
    std::cout << "target voxel num " << target_volume.mVoxelPosition.rows()<<std::endl;

    std::cout << "find dense correspondence between source and target ... "<<std::flush;
    EMD emd_flow;
    emd_flow.construct_correspondence(source_volume, target_volume);
    //emd_flow.min_cost_flow(source_volume, target_volume);
    //emd_flow.find_correspondence(source_volume, target_volume);
    std::cout << "done!"<<std::endl;

    std::cout<<"computing source_target thin plate spline ... "<<std::flush;
    ThinPlateSpline source_target_tps; 
    source_target_tps.compute_tps(source_volume.mVoxelPosition, emd_flow.corresp_source_target_);
    MatrixX3r final_corresp_points;
    source_target_tps.interpolate(source_volume.mVoxelPosition, final_corresp_points);
    std::cout<<"done!"<<std::endl;
    source_volume.calc_tetrahedron_transform(final_corresp_points);
    MatrixX3r inter_corresp_points;
    source_volume.find_intermedium_points(inter_corresp_points, 1.0);
*/


    /*  
    std::ofstream output_inter_corresp("inter_corresp.dat");
    output_inter_corresp << inter_corresp_points;
    output_inter_corresp.close();

    std::ofstream output_init_corresp("corresp_target_pos.dat");
    output_init_corresp << emd_flow.corresp_source_target_;
    output_init_corresp.close();
    std::ofstream output_tps_corresp("tps_corresp.dat");
    output_tps_corresp << final_corresp_points;
    output_tps_corresp.close();
    auto position_error = new_position - emd_flow.corresp_source_target_;
    std::ofstream output_error("error.dat");
    output_error<< position_error;
    output_error.close();
    */

    /*
    std::ofstream output_dist_vector("distance_field.dat");
    output_dist_vector << source_volume.distance_vector_field;
    output_dist_vector.close();
    std::ofstream output_laplace ("laplace.dat");
    output_laplace << source_volume.mLaplaceMatrix;
    output_laplace.close();
    std::ofstream output_voxel_pos("source_voxel_pos.dat");
    output_voxel_pos << source_volume.mVoxelPosition;
    output_voxel_pos.close();
    */

    /*
    std::ofstream output_target_voxel("target_voxel_pos.dat");
    output_target_voxel << target_volume.mVoxelPosition;
    output_target_voxel.close();
    */
	return 0;
}
