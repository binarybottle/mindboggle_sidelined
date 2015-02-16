// ------------------------ RicGraphNode.cpp --------------------------------
/*
@file
Implementation file for the RicGraphNode class to read sulcal node info
from a BrainVisa graph (.arg) file.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <RicUtil.h>
#include "RicGraph.h"
#include "RicGraphNode.h"

/*!
Constructor for a RicGraphNode
*/
RicGraphNode::RicGraphNode()
{
	type = "uninitialized";
	Tmtktri_filename = "";

}

/*!
Constructor for a RicGraphNode that reads a node from
a file stream.

@param gfile - stream to read node from
*/
RicGraphNode::RicGraphNode(ifstream *gfile)
{
	type = "uninitialized";
	Tmtktri_filename = "";

}

/*!
Destructor for a RicGraphNode
*/
RicGraphNode::~RicGraphNode()
{

}

/*!
Routine to read a RicGraphNode from a file stream. The stream is parsed on
a line by line basis assigning values in the file to member variables.

@param gfile - stream to read node from

@returns - 1 on success, 0 if there is a problem
*/
int RicGraphNode::Read(ifstream *gfile)
{
	char tmp[512];
	string tok[30];
	int n;

	// check to see that this really is a node
	gfile->getline(tmp,511);
	n=tokenize(tmp,tok);
	if ( tok[1] != "NODE") return 0;

	// get the type and number from the header
	type = tok[2];
	number = atoi(tok[3].c_str());

	gfile->getline(tmp,511);

	// in a tedious manner parse the header and assign values to members
	while ((n=tokenize(tmp,tok)) != 0 )
	{
		if ( tok[0] == "bottom_point_number" )
		{
			bottom_point_number = atoi(tok[1].c_str());
		}
		else if( tok[0] == "name" )
		{
			name = tok[1];
		}
		else if( tok[0] == "point_number" )
		{
			point_number = atoi(tok[1].c_str());
		}
		else if( tok[0] == "size" )
		{
			size = (float)atof(tok[1].c_str());
		}
		else if( tok[0] == "skeleton_label" )
		{
			skeleton_label = atoi(tok[1].c_str());
		}
		else if( tok[0] == "ss_point_number" )
		{
			ss_point_number = atoi(tok[1].c_str());
		}
		else if( tok[0] == "Tal_boundingbox_max" )
		{
			Tal_boundingbox_max.x = (float)atof(tok[1].c_str());
			Tal_boundingbox_max.y = (float)atof(tok[2].c_str());
			Tal_boundingbox_max.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "Tal_boundingbox_min" )
		{
			Tal_boundingbox_min.x = (float)atof(tok[1].c_str());
			Tal_boundingbox_min.y = (float)atof(tok[2].c_str());
			Tal_boundingbox_min.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "Tmtktri_label" )
		{
			Tmtktri_label = atoi(tok[1].c_str());
		}
		else if( tok[0] == "bottom_label" )
		{
			bottom_label = atoi(tok[1].c_str());
		}
		else if( tok[0] == "boundingbox_max" )
		{
			boundingbox_max.x = atoi(tok[1].c_str());
			boundingbox_max.y = atoi(tok[2].c_str());
			boundingbox_max.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "boundingbox_min" )
		{
			boundingbox_min.x = atoi(tok[1].c_str());
			boundingbox_min.y = atoi(tok[2].c_str());
			boundingbox_min.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "gravity_center" )
		{
			gravity_center.x = (float)atof(tok[1].c_str());
			gravity_center.y = (float)atof(tok[2].c_str());
			gravity_center.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "index" )
		{
			index = atoi(tok[1].c_str());
		}
		else if( tok[0] == "label" )
		{
			label = tok[1];
		}
		else if( tok[0] == "Tmtktri_filename" )
		{
			Tmtktri_filename = tok[1];
		}
		else if( tok[0] == "maxdepth" )
		{
			maxdepth = (float)atof(tok[1].c_str());
		}
		else if( tok[0] == "mindepth" )
		{
			mindepth = (float)atof(tok[1].c_str());
		}
		else if( tok[0] == "normal" )
		{
			normal.x = (float)atof(tok[1].c_str());
			normal.y = (float)atof(tok[2].c_str());
			normal.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "other_label" )
		{
			other_label = atoi(tok[1].c_str());
		}
		else if( tok[0] == "other_point_number" )
		{
			other_point_number = atoi(tok[1].c_str());
		}
		else if( tok[0] == "refgravity_center" )
		{
			refgravity_center.x = (float)atof(tok[1].c_str());
			refgravity_center.y = (float)atof(tok[2].c_str());
			refgravity_center.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "refnormal" )
		{
			refnormal.x = (float)atof(tok[1].c_str());
			refnormal.y = (float)atof(tok[2].c_str());
			refnormal.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "refsize" )
		{
			refsize = (float)atof(tok[1].c_str());
		}
		else if( tok[0] == "refsurface_area" )
		{
			refsurface_area = (float)atof(tok[1].c_str());
		}
		else if( tok[0] == "rootsbassin" )
		{
			rootsbassin = atoi(tok[1].c_str());
		}
		else if( tok[0] == "ss_label" )
		{
			ss_label = atoi(tok[1].c_str());
		}
		else if( tok[0] == "surface_area" )
		{
			surface_area = (float)atof(tok[1].c_str());
		}
		else if( tok[0] == "talcovar" )
		{
			talcovar[0] = (float)atof(tok[1].c_str());
			talcovar[1] = (float)atof(tok[2].c_str());
			talcovar[2] = (float)atof(tok[3].c_str());
			talcovar[3] = (float)atof(tok[4].c_str());
			talcovar[4] = (float)atof(tok[5].c_str());
			talcovar[5] = (float)atof(tok[6].c_str());
			talcovar[6] = (float)atof(tok[7].c_str());
			talcovar[7] = (float)atof(tok[8].c_str());
			talcovar[8] = (float)atof(tok[9].c_str());
		}

		// get another line
		gfile->getline(tmp,511);

	}
	return 1;
}
