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
    std::vector<std::pair<Vector3r, Vector3r> > corresp_pairs;
    /*   // Model 389-399
    // head, left front, right front, left behind, right behind, tail
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
    */
    // model 382-396
    // head middle, head left&right, brest left&right
    corresp_pairs.push_back(std::make_pair(Vector3r(0.219851, 0.213655, 0.495029), Vector3r(0.00474127, 0.192, 0.499866)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.0381196, 0.382666, 0.46442), Vector3r(-0.0397595, 0.417716, 0.42786)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.107545, 0.401916, 0.363451), Vector3r(0.0369645, 0.42119, 0.422515)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.179393, 0.0734587, 0.215654), Vector3r(-0.116293, 0.0264923, 0.205875)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.0592673, 0.0569322, 0.169), Vector3r(0.119458, 0.0189597, 0.187779)));
    //four leg: left front, right front, left behind, right behind, 
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.165382, -0.408648, 0.101051), Vector3r(-0.070656, -0.46222, 0.120456)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.0142926, -0.413487, 0.225126), Vector3r(0.0763945, -0.378731, 0.13358)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.185868, -0.405328, -0.336726), Vector3r(-0.0713341, -0.41537, -0.203423)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.0238738, -0.411759, -0.449048), Vector3r(0.0535133, -0.462253, -0.397268)));
    //leg middle
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.152678, -0.172809, 0.144243), Vector3r(-0.0673469, -0.218646, 0.160231)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.000314611, -0.17719, 0.188939), Vector3r(0.0817346, -0.187616, 0.311013)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.186793, -0.161489, -0.298735), Vector3r(-0.0625119, -0.192748, -0.292893)));
    corresp_pairs.push_back(std::make_pair(Vector3r(0.026597, -0.151722, -0.39207), Vector3r(0.0677306, -0.211979, -0.335093)));


    // back, stamock, ass
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.0524261, 0.168286, 0.0228699), Vector3r(0.0168113, 0.182663, 0.012828)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.0647473, -0.0800051, -0.121843), Vector3r(-0.0242501, -0.0622054, -0.161122)));
    corresp_pairs.push_back(std::make_pair(Vector3r(-0.0945625, 0.206458, -0.266852), Vector3r(0.0225166, 0.19435, -0.295222)));
    

    Morph morph(argv[1], argv[2], corresp_pairs, 0.02);
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
