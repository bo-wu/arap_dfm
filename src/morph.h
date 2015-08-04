/*
 * =====================================================================================
 *
 *       Filename:  morph.h  Version:  1.0  Created:  01/16/2015 04:56:49 PM
 *
 *    Description:  morph two volume data
 *
 *         Author:  Bo Wu (Robert), wubo.gfkd@gmail.com
 *	    Copyright:  Copyright (c) 2015, Bo Wu
 *   Organization:  National University of Defense Technology
 *
 * =====================================================================================
 */

#ifndef MORPH_H_
#define MORPH_H_
#include "def_types.h"
#include "volume_object.h"

struct Morph
{
	Morph();
	~Morph();
	VolumeObject source, target;
    
    //MatrixX3r intermedia_from_source, intermedia_from_target;
	void start_morph();
};

#endif

