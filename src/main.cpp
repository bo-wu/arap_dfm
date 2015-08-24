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
#include <string>
#include <sstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "def_types.h"
#include "volume_object.h"
#include "dense_correspondence.h"
#include "thin_plate_spline.h"
#include "morph.h"
#include "skeleton.h"

std::vector<std::pair<Vector3r, Vector3r> > read_anchor_file(std::string anchor_file_name);

int main(int argc, char** argv)
{
    if (argc <= 2)
    {
        std::cerr << "Not enough parameters, should run like\n ./morph source_mesh_name target_mesh_name\n";
        return 0;
    }

    std::vector<std::pair<Vector3r, Vector3r> > corresp_pairs;

    std::string source_name = argv[1];
    std::string target_name = argv[2];
    std::size_t found = target_name.find_last_of("/\\");
    std::string target_base = target_name.substr(found+1, target_name.size()-found-5);

    std::string extra_anchor_name = source_name.substr(0, source_name.length()-4) + "_" + target_base + ".acr";
    std::cout<<"anchor name "<<extra_anchor_name<<std::endl;
    corresp_pairs = read_anchor_file(extra_anchor_name);

    //alladin_dengshen
    Morph morph(argv[1], argv[2], corresp_pairs, 0.0175, 0.008, 0.0175, 0.008);
    // 201-208 
    //Morph morph(argv[1], argv[2], corresp_pairs, 0.02, 0.01, 0.02, 0.01);
    //Dino-397
    //Morph morph(argv[1], argv[2], corresp_pairs, 0.015, 0.015, 0.0195, 0.015);
    // gecko-seal
    //Morph morph(argv[1], argv[2], corresp_pairs, 0.0145, 0.0145, 0.02, 0.02);
    // gorilla-horse
    //Morph morph(argv[1], argv[2], corresp_pairs, 0.018, 0.018, 0.018, 0.018);
    //
    //
    //morph.source_volume_.write_grid(argv[1]);
    //morph.target_volume_.write_grid(argv[2]);
    
    morph.initial(true);
    morph.start_basic_morph(0.1);
    morph.write_sequence();
//    
 	return 0;
}



std::vector<std::pair<Vector3r, Vector3r> > read_anchor_file(std::string anchor_file_name)
{
    std::ifstream input_anchor(anchor_file_name);

    std::vector<std::pair<Vector3r, Vector3r> > corresp_pairs;

    if(!input_anchor.is_open())
    {
        std::cerr << "no extra anchors\n";
        return corresp_pairs;
    }
    
    std::string line;
    Vector3r s_anchor, t_anchor;

    while(std::getline(input_anchor, line))
    {
        if (line.empty() || line.at(0) == '#' || line.find_first_not_of(' ') == std::string::npos)
            continue;
        else 
        {
            std::istringstream iss(line);
            for(int i=0; i < 3; ++i)
            {
                if( !(iss >> s_anchor(i)) ) 
                {
                    std::cerr << "input anchor file format wrong\n";
                    input_anchor.close();
                    exit(-1);
                }
            }

            for(int i=0; i < 3; ++i)
            {
                if( !(iss >> t_anchor(i)) ) 
                {
                    std::cerr << "input anchor file format wrong\n";
                    input_anchor.close();
                    exit(-1);
                }
            }
        }
        corresp_pairs.push_back(std::make_pair(s_anchor, t_anchor) );
    }

    input_anchor.close();

    return corresp_pairs;
}
