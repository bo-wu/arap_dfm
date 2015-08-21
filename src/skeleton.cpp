/*
 * =====================================================================================
 *
 *       Filename:  skeleton.cpp      Created:  08/17/2015 03:46:40 PM
 *
 *    Description:  
 *
 *         Author:  Wu Bo (Robert), wubo.gfkd@gmail.com
 *		Copyright:	Copyright (c) 2015, Wu Bo
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "skeleton.h"


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  init
 *  Description:  
 * =====================================================================================
 */
void Skeleton::init(std::string mesh_name)
{
    TriMesh mesh;
    if(!OpenMesh::IO::read_mesh(mesh, mesh_name))
    {
        std::cerr <<"read mesh "<<mesh_name<<" error!\n";
        exit(-1);
    }
    
    std::string base_name = mesh_name.substr(0, mesh_name.length()-4);

    //============= set tag for each point =============
    mesh_face_point_tag = VectorXi(mesh.n_faces());
    mesh_face_points = MatrixX3r(mesh.n_faces(), 3);

    std::string face_list_name = base_name + "_allGC_0.tris";
    std::ifstream input_branch_tris(face_list_name);
    if(!input_branch_tris.is_open())
    {
        std::cerr <<"cannot open branch triangles list file "<<face_list_name<<std::endl;
        input_branch_tris.close();
        exit(-1);
    }

    int num_tris;
    int f_index;
    TriMesh::FaceHandle f_h;
    TriMesh::FaceVertexIter fv_it;
    OpenMesh::Vec3d face_center;
    int f_tag = 0;
    int p = 0;
    while(input_branch_tris >> num_tris)
    {
        std::ofstream output_part_verts(base_name+"_part"+std::to_string(p)+".obj");
        ++p;
        for(int i=0; i < num_tris; ++i)
        {
            face_center = OpenMesh::Vec3d(0,0,0);
            input_branch_tris >> f_index;
            f_h = mesh.face_handle(f_index);
            for(fv_it = mesh.fv_iter(f_h); fv_it.is_valid(); ++fv_it)
            {
                face_center += mesh.point(*fv_it);
            }
            face_center /= 3.0;

            mesh_face_points.row(f_index) << face_center[0], face_center[1], face_center[2];
            mesh_face_point_tag(f_index) = f_tag;


            output_part_verts <<"v "<< face_center[0]<<" "<<face_center[1]<<" "<<face_center[2]<<std::endl;
        }
        f_tag++; // group tag
        output_part_verts.close();
    }

    input_branch_tris.close();
    //=================================================
    //
    //===========read points of each branch=============
    std::string skel_point_name = base_name + "_allGC_0.obj";
    // to save some effort store as points
    TriMesh skel_point;
    if(!OpenMesh::IO::read_mesh(skel_point, skel_point_name))
    {
        std::cerr<<"read skel points "<<skel_point_name<<" error!\n";
        exit(-1);
    }
    
    int v_num = skel_point.n_vertices();
    const int num_each_branch = 10;
    MatrixX3r brach_points(num_each_branch, 3);

    for(int i=0; i < v_num; ++i)
    {
        TriMesh::VertexHandle v_h = skel_point.vertex_handle(i);
        auto p = skel_point.point(v_h);
        brach_points.row(i%num_each_branch) << p[0], p[1], p[2];

        if(i%num_each_branch == num_each_branch-1)
        {
            skel_branch_points.push_back(brach_points);
        }
    }
    //======================================================
    // ====== get which two part should be together
    std::string merge_part_file = base_name + ".merge";
    std::ifstream input_merged_parts(merge_part_file);
    if(input_merged_parts.is_open())
    {
        int a, b;
        while(input_merged_parts >> a >> b)
        {
            merged_branch.push_back(std::make_pair(a, b));
        }
    }

    input_merged_parts.close();

    //======================================================
}		/* -----  end of function init  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  SkeletonPair
 *  Description:  
 * =====================================================================================
 */
void SkeletonPair::read_match_info(Skeleton &ss, Skeleton &ts, std::string name)
{
    std::ifstream input_match_pairs(name);
    if(!input_match_pairs.is_open())
    {
        std::cerr << "open skeleton branch matched file error\n";
        input_match_pairs.close();
        exit(-1);
    }

    matched_branch.clear();
    int a, b;
    while(input_match_pairs >> a >> b)
    {
        matched_branch.push_back(std::make_pair(a, b));
    }

    const int num_corresp_each_brach = 4;
    Vector3r s_branch_point, t_branch_point;
    for(int i=0; i < matched_branch.size(); ++i)
    {
        for(int j=0; j < 10; j+=3)
        {
            s_branch_point = ss.skel_branch_points[matched_branch[i].first].row(j);
            t_branch_point = ts.skel_branch_points[matched_branch[i].second].row(j);
            corresp_skel_points.push_back(std::make_pair(s_branch_point, t_branch_point));
        }
    }

}		/* -----  end of function SkeletonPair  ----- */

