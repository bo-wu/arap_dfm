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

int main(int argc, char** argv)
{
    if (argc <= 2)
    {
        std::cerr << "Not enough parameters, should run like\n ./morph source_mesh_name target_mesh_name\n";
        return 0;
    }

    /*  
    std::string source_name = argv[1];
    std::string target_name = argv[2];
    std::size_t found = target_name.find_last_of("/\\");
    //std::string anchor_file_name = source_name.substr(0, source_name.size()-4) + "_" + target_name.substr(found+1, target_name.size()-found-5) + ".anchors";
    std::string skel_pair_file_name = source_name.substr(0, source_name.size()-4) + "_" + target_name.substr(found+1, target_name.size()-found-5) + ".brp";

    */
    
/*  
    std::ifstream input_anchor(argv[3]);
    if(!input_anchor.is_open())
    {
        std::cerr << "cannot open anchor file, exist!\n";
        exit(-1);
    }
    
    std::string line;
    Vector3r s_anchor, t_anchor;
    // read positions from skeleton correspondence
    std::vector<std::pair<Vector3r, Vector3r> > corresp_pairs;

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

*/
    std::vector<std::pair<Vector3r, Vector3r> > corresp_pairs;
    Morph morph(argv[1], argv[2], corresp_pairs, 0.018, 0.01, 0.009, 0.01);
    //Morph morph(argv[1], argv[2], corresp_pairs, 0.02, 0.01, 0.02, 0.01);
    //morph.source_volume_.write_grid(argv[1]);
    //morph.target_volume_.write_grid(argv[2]);
    
///    morph.initial();
///    morph.start_basic_morph(0.25);
///    morph.write_sequence();
    
 	return 0;
}
