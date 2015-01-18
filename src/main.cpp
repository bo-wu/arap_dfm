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
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "volume_object.h"
int main(int argc, char** argv)
{
	VolumeObject vo(argv[1]);
	vo.calc_vector_field();
//	vo.write_grid(argv[1]);
	return 0;
}
