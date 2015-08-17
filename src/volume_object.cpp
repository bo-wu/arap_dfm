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
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <openvdb/Grid.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetUtil.h>
#include "volume_object.h"

//#define MIN_QUAD_WITH_FIXED_CPP_DEBUG
#include <igl/min_quad_with_fixed.h>

VolumeObject::VolumeObject(Real transform_scale) : transform_scale_(transform_scale)
{
}

VolumeObject::VolumeObject(std::string name, Real transform_scale, Real dense_transform_scale)
{
	mesh_name = name;
    transform_scale_ = transform_scale;
    dense_transform_scale_ = dense_transform_scale;
    first_time_ = true;
	initial_volume();
    initial_dense_volume();
}

void VolumeObject::initial(std::string name, Real transform_scale, Real dense_transform_scale)
{
    mesh_name = name;
    transform_scale_ = transform_scale;
    dense_transform_scale_ = dense_transform_scale;
    first_time_ = true;
    initial_volume();
    initial_dense_volume();
}

VolumeObject::~VolumeObject()
{
}


///////////////////////////////////////////
//generate volume from mesh
void VolumeObject::initial_volume()
{
    // if resize, if output mesh
	read_mesh(false, false);
	grid = openvdb::FloatGrid::create(2.0);
	openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(transform_scale_);
	grid->setTransform(grid_transform);
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setName("mesh_grid");
	grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(grid->transform(), points, triangles, float(openvdb::LEVEL_SET_HALF_WIDTH));

    interior_grid = openvdb::tools::sdfInteriorMask(*grid);
    //make inside grid dense
    interior_grid->tree().voxelizeActiveTiles();
    voxel_num_ = interior_grid->tree().activeVoxelCount();

    mVoxelPosition = MatrixX3r::Zero(voxel_num_, 3);

    int k = 0;
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        auto voxel_pos = grid->indexToWorld(iter.getCoord());
        mVoxelPosition.row(k) << voxel_pos[0], voxel_pos[1], voxel_pos[2];
    }

} //end of initial_volume


///generate dense voxel grid from mesh
// should call initial volume first
void VolumeObject::initial_dense_volume()
{
    dense_grid = openvdb::FloatGrid::create(2.0);
    openvdb::math::Transform::Ptr grid_transform = openvdb::math::Transform::createLinearTransform(dense_transform_scale_);
    dense_grid->setTransform(grid_transform);
    dense_grid->setGridClass(openvdb::GRID_LEVEL_SET);
    dense_grid->setName("dense_grid");
    dense_grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(dense_grid->transform(), points, triangles, float(openvdb::LEVEL_SET_HALF_WIDTH));

    interior_dense_grid = openvdb::tools::sdfInteriorMask(*dense_grid);
    interior_dense_grid->tree().voxelizeActiveTiles();
    dense_voxel_num_ = interior_dense_grid->tree().activeVoxelCount();

    // get dense point coordinates
    mDenseVoxelPosition = MatrixX3r(dense_voxel_num_, 3);

    int k = 0;
    for(auto iter=interior_dense_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        auto voxel_pos = dense_grid->indexToWorld(iter.getCoord());
        mDenseVoxelPosition.row(k) << voxel_pos[0], voxel_pos[1], voxel_pos[2];
    }

    for(int i=0; i < 3; ++i)
    {
        mass_center(i) = mDenseVoxelPosition.col(i).mean();
    }

    // construct tetrahedron
    kd_tree_type dense_voxel_kdtree(3, mDenseVoxelPosition);
    dense_voxel_kdtree.index->buildIndex();
    long int outIndex;
    Real outDistance;
    dense_voxel_kdtree.query(mass_center.data(), 1, &outIndex, &outDistance);

    for(int i=0; i < 3; ++i)
    {
        Vector3r v_voxel = mass_center;
        for(int j=-1; j<=1; j+=2)
        {
            v_voxel(i) += j * 2 * dense_transform_scale_;
            dense_voxel_kdtree.query(v_voxel.data(), 1, &outIndex, &outDistance);
            tet_anchor.push_back(std::make_pair(mDenseVoxelPosition.row(outIndex), outIndex));
        }
    }

    // keep fixed
    mass_center_voxel_index = outIndex;
    mass_center = mDenseVoxelPosition.row(outIndex);

    // construct tetrahedron
    std::vector< std::vector<int> > neighbor_index_3d; //3 directions
    std::vector<int> neighbor_index_1d; // 1 direction
    
    openvdb::Coord v_coord;
    Vector3r v_world_pos;

    k = 0;
    for(auto iter=interior_dense_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        v_coord = iter.getCoord();
        neighbor_index_3d.clear();
        for(int i=0; i < 3; ++i)
        {
            neighbor_index_1d.clear();
            for(int j=-1; j<=1; j+=2)
            {
                auto temp_coord = v_coord;
                temp_coord[i] = v_coord[i] + j;
                if(interior_dense_grid->tree().isValueOn(temp_coord))
                {
                    auto voxel_pos = dense_grid->indexToWorld(temp_coord);
                    v_world_pos << voxel_pos[0], voxel_pos[1], voxel_pos[2];
                    dense_voxel_kdtree.query(v_world_pos.data(), 1, &outIndex, &outDistance);
                    if(outDistance > 1.0e-10)
                    {
                        std::cerr <<"Dense Distance " << outDistance <<" should be 0.0\n";
                    }
                    neighbor_index_1d.push_back(outIndex);
                }
            }
            if(neighbor_index_1d.size() > 0)
            {
                neighbor_index_3d.push_back(neighbor_index_1d);
            }
        }

        // construct tetrahedron
        if(neighbor_index_3d.size() == 3)
        {
            Vector4i tet_index;
            tet_index(0) = k;
            for(int l=0; l < neighbor_index_3d[0].size(); ++l)
                for(int m=0; m < neighbor_index_3d[1].size(); ++m)
                    for(int n=0; n < neighbor_index_3d[2].size(); ++n)
                    {
                        tet_index(1) = neighbor_index_3d[0].at(l);
                        tet_index(2) = neighbor_index_3d[1].at(m);
                        tet_index(3) = neighbor_index_3d[2].at(n);
                        mTetIndex.push_back(tet_index);
                    }
        }

    }

}


void VolumeObject::set_anchors(const std::vector<Vector3r>& anchors)
{
    mAnchors = anchors;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  construct_laplace_matrix
 *  Description:  laplace matrix for distance vector field (sparse grid) 
 * =====================================================================================
 */
void VolumeObject::construct_laplace_matrix()
{
    //for sparse grid, should be activeVoxel + activeTile
    auto voxelNum = interior_grid->tree().activeLeafVoxelCount();
    auto anchorNum = mAnchors.size();
    mLaplaceMatrix = SpMat(voxelNum, voxelNum);
    distance_vector_field = MatrixXr::Zero(voxelNum, anchorNum);
    int k = 0;

    //find neighbor voxel index
    kd_tree_type voxelKDTree(3, mVoxelPosition);
    voxelKDTree.index->buildIndex();
    long int outIndex;
    Real outDistance;

    // 6 neighbor
    std::vector<MyTriplet> laplace_triplet_list;
    laplace_triplet_list.reserve(7*voxelNum);

    int degree;
    Vector3r v_world_pos;
    openvdb::Coord v_coord;
    k = 0;

    // laplace matrix
    for(auto iter=interior_grid->cbeginValueOn(); iter; ++iter, ++k)
    {
        v_coord = iter.getCoord();
        degree = 0;
        for(int i=0; i < 3; ++i)
        {
            for(int j=-1; j <= 1; j+=2)
            {
                auto temp_coord = v_coord;
                temp_coord[i] = v_coord[i] + j;
                if (interior_grid->tree().isValueOn(temp_coord))
                {
                    ++degree;
                    auto voxel_pos = grid->indexToWorld(temp_coord);
                    v_world_pos<< voxel_pos[0], voxel_pos[1], voxel_pos[2];
                    voxelKDTree.query(v_world_pos.data(), 1, &outIndex, &outDistance);
                    if(outDistance > 1.0e-8)
                    {
                        std::cerr<<"Distance "<<outDistance<<" should be 0.0\n";
                    }
                    laplace_triplet_list.push_back(MyTriplet(k, outIndex, -1));
                }
            }
        }
        
        laplace_triplet_list.push_back(MyTriplet(k, k, degree));

    }

    mLaplaceMatrix.setFromTriplets(laplace_triplet_list.begin(), laplace_triplet_list.end());
    
    //anchor points for harmonic field
    if (anchorNum > 0)
        constraint_index_ = VectorXi::Zero(anchorNum, 1);
    k = 0;
    for (auto a : mAnchors)
    {
        voxelKDTree.query(a.data(), 1, &outIndex, &outDistance);
        constraint_index_(k) = outIndex;
        //std::cout << "anchor "<<k<< " index is "<<outIndex<<std::endl;
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

    /*   // solving use igl
    */
    igl::min_quad_with_fixed_data<Real> mqwf;
    int num_row = constraint_index_.rows();
    VectorXr B = VectorXr::Zero(voxel_num_, 1);
    //Empyty constraints (except for constraint_index_/value)
    SpMat Aeq;
    VectorXr Beq;

    igl::min_quad_with_fixed_precompute(mLaplaceMatrix, constraint_index_, Aeq, true, mqwf);

    // solve equation with constraint
    VectorXr D;
    for(int i=0; i < num_row; ++i)
    {
        VectorXr constraint_value = VectorXr::Zero(num_row, 1);
        constraint_value(i) = 1.0;
        igl::min_quad_with_fixed_solve(mqwf, B, constraint_value, Beq, D);
        distance_vector_field.col(i) = D;
    }

    /*  
    std::ofstream output_constraint_index("constraint_index.dat");
    output_constraint_index << constraint_index_;
    output_constraint_index.close();
    */
}		/* -----  end of function calc_vector_field  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  polar_decompose
 *  Description:  
 * =====================================================================================
 */
void VolumeObject::polar_decompose(const Matrix3r &rest, const Matrix3r &deform, Matrix3r &R, Matrix3r &S)
{
    Matrix3r defrom_gradient = deform * rest.inverse();
    Eigen::JacobiSVD<Matrix3r> svd(defrom_gradient, Eigen::ComputeFullU|Eigen::ComputeFullV);
    Matrix3r eigen_values = Matrix3r::Zero();
    Vector3r singular = svd.singularValues();
    for(int i=0; i < 3; ++i)
        eigen_values(i,i) = singular(i);

    R = svd.matrixU() * svd.matrixV().transpose();
    S = svd.matrixV() * eigen_values * svd.matrixV().transpose();

    if(R.determinant() < 0)
    {
        Matrix3r I = Matrix3r::Identity();
        I(2,2) = -1.0;
        R = R * I;
        S = I * S;
    }
}		/* -----  end of function polar_decompose  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calc_tetrahedron_transform
 *  Description:  fill  mTetTransform
 *  @param final_corresp_points: correspondence after warping  
 * =====================================================================================
 */
void VolumeObject::calc_tetrahedron_transform(const MatrixX3r &final_corresp_points)
{
    if(mDenseVoxelPosition.rows() != final_corresp_points.rows())
    {
        std::cerr << "input points num wrong for tetrahedron transformation\n";
        exit(-1);
    }

    Matrix3r R, S, rest, deform;
    for(int i=0; i < mTetIndex.size(); ++i)
    {
        for(int j=1; j < 4; ++j)
        {
            rest.col(j-1) = mDenseVoxelPosition.row( mTetIndex[i](j) ) - mDenseVoxelPosition.row(mTetIndex[i](0));
            deform.col(j-1) = final_corresp_points.row( mTetIndex[i](j) ) - final_corresp_points.row(mTetIndex[i](0));
        }

        polar_decompose(rest, deform, R, S);

        mTetTransform.push_back(std::make_pair(R, S));
    }
}		/* -----  end of function calc_tetrahedron_transform  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_intermedium_points
 *  Description:  find intermedium points at time t [0,1]
 *  @param inter_corresp_points: correspondence by interpolating warping
 * =====================================================================================
 */
void VolumeObject::find_intermedium_points(MatrixX3r &inter_corresp_points, const Real t)
{
    int fixed_num = 1;
    fixed_num += tet_anchor.size();
    int dense_voxel_num = mDenseVoxelPosition.rows();
    int tet_num = mTetIndex.size();

    inter_corresp_points = MatrixX3r(dense_voxel_num, 3);

    Real weight = 10.0;
    Real center_neighbour_weight = 1.0;

    // with one anchor point
    // construct A
    if(first_time_)
    {
        tet_matrix_ = SpMat(3*tet_num+fixed_num, dense_voxel_num);

        std::vector<MyTriplet> tet_triplet_list;
        tet_triplet_list.reserve(6*tet_num+fixed_num);
        for(int i=0; i < tet_num; ++i)
        {
            for(int j=0; j < 3; ++j)
            {
                tet_triplet_list.push_back(MyTriplet(3*i+j, mTetIndex[i](0), 1));
                tet_triplet_list.push_back(MyTriplet(3*i+j, mTetIndex[i](j+1), -1));
            }
        }

        // mass_center is not real the anchor point's position(a little difference)
        tet_triplet_list.push_back(MyTriplet(3*tet_num, mass_center_voxel_index, weight));
        for(int i=0; i < tet_anchor.size(); ++i)
        {
            tet_triplet_list.push_back(MyTriplet(3*tet_num+1+i, tet_anchor[i].second, center_neighbour_weight));
        }

        tet_matrix_.setFromTriplets(tet_triplet_list.begin(), tet_triplet_list.end());

        tet_normal_matrix = tet_matrix_.transpose() * tet_matrix_;
        cg.compute(tet_normal_matrix);
        first_time_ = false;
    }

    MatrixX3r B = MatrixX3r::Zero(3*tet_num+fixed_num, 3);

    Quaternionr quat_I, quat_res;
    quat_I.setIdentity();

    Matrix3r M;
    Matrix3r mat_I = Matrix3r::Identity();
    Matrix43r tet_vert;

//    std::ofstream output_transform("transform.dat");
    // construct B
    for(int i=0; i < tet_num; i++)
    {
        //================  construct B  ====================
        Quaternionr quat(mTetTransform[i].first); // test if quat equal first
        quat_res = quat_I.slerp(t, quat);
        M = quat_res.toRotationMatrix() * ( (1-t)*mat_I + t*mTetTransform[i].second );

        for(int j=0; j < 4; ++j)
        {
            tet_vert.row(j) = M * mDenseVoxelPosition.row(mTetIndex[i](j)).transpose();
        }
        //===============  construct B  ===============
        for(int k=0; k < 3; ++k)
        {
            B.row(3*i+k) = tet_vert.row(0) - tet_vert.row(k+1);
        }

//        output_transform <<"R\n"<< mTetTransform[i].first <<"\nS\n"<<mTetTransform[i].second;
//        output_transform <<"\nR*S\n"<<mTetTransform[i].first * mTetTransform[i].second;
//        output_transform <<"\nM\n" << M <<"\n\n";

    }
//    output_transform.close();

    B.row(3*tet_num) = weight * mass_center;

    for(int i=0; i < tet_anchor.size(); ++i)
    {
        B.row(3*tet_num+1+i) = center_neighbour_weight * tet_anchor[i].first;
    }

    B = tet_matrix_.transpose() * B;

    std::clock_t start;
    start = std::clock();
    for(int i=0; i < 3; ++i)
    {
        inter_corresp_points.col(i) = cg.solve(B.col(i));
//        std::cout << "col("<<i<<") #iterations: " << cg.iterations()<<std::endl;
//        std::cout << "estimated error "<<cg.error() <<std::endl;
        if(cg.info() != Eigen::Success)
        {
            std::cout << "ConjugateGradient solver not converage\n";
            exit(-1);
        }
    }
    Real elapse = (std::clock() - start) / (Real)(CLOCKS_PER_SEC);
    std::cout<<"intermedium correspondence points elapse "<<elapse<<"s\n";

    /*  
    std::cout <<"mass center voxel index " << mass_center_voxel_index<<std::endl;
    std::cout<<"tet_triplet_list.size "<<tet_triplet_list.size()<<" reserve size " << 6*tet_num+fixed_num<<std::endl;
    */
    /*
    std::ofstream output_l_row("l_row.dat");
    std::ofstream output_l_col("l_col.dat");
    std::ofstream output_l_value("l_value.dat");
            output_l_row << 3*i+k <<" ";
            output_l_col << mTetIndex[i](0) <<" ";
            output_l_value << 1 <<" ";
            output_l_row << 3*i+k <<" ";
            output_l_col << mTetIndex[i](k+1) <<" ";
            output_l_value << -1 <<" ";
    //output_l_row << 3*tet_num;
    //output_l_col << mass_center_voxel_index;
    //output_l_value << weight;

    output_l_row.close();
    output_l_col.close();
    output_l_value.close();
    */

    /*
    std::ofstream output_L("L_matrix.dat");
    output_L << L;
    output_L.close();
    std::ofstream output_B("B_position.dat");
    output_B << B;
    output_B.close();

    std::ofstream output_L_col("L_col_sum.dat");
    std::ofstream output_L_row("L_row_sum.dat");

    VectorXr row_values = VectorXr::Zero(L.rows());
    VectorXr col_values = VectorXr::Zero(L.cols());

    for(int k=0; k < L.outerSize(); ++k)
    {
        for(SpMat::InnerIterator it(L, k); it; ++it)
        {
            row_values(it.row()) += it.value();
            col_values(it.col()) += it.value();
        }

    }
    output_L_col << col_values;
    output_L_row << row_values;
    output_L_col.close();
    output_L_row.close();
      */

    /*  
    std::ofstream output_res("inter_points.dat");
    output_res << inter_corresp_points;
    output_res.close();
    */
}		/* -----  end of function find_intermedium_points  ----- */


///////////////////////////////////////////
//read_mesh
void VolumeObject::read_mesh(bool resize, bool b_out)
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

    if(b_out)
    {
        if(!OpenMesh::IO::write_mesh(mesh, mesh_name.c_str()))
        {
            std::cerr << "write mesh on disk error\n";
            exit(-1);
        }
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
		name = name.substr(0, name.length()-4) + ".vdb";
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

