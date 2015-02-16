// ------------------------------ RicGraph.cpp ---------------------------------
/*!
@file
Implementation file for the RicGraph class.
*/
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "RicGraphNode.h"
#include "RicGraph.h"

using namespace std;

/*!
Base constructor
*/
RicGraph::RicGraph()
{
	nnodes = 0;
	nodes = NULL;
	node_idx = 0;
}

/*!
Constructor to load graph from a file

@param fname - name of file to read
*/
RicGraph::RicGraph(string fname)
{
	nnodes = 0;
	nodes = NULL;
	node_idx = 0;
	this->Read(fname);
}

/*!
Destructor cleans up memory allocated for nodes
*/
RicGraph::~RicGraph()
{

	if ( nodes ) delete [] nodes;
}

/*!
Reads a graph from a file.

@param fname - name of file to read

@returns - 1 on success, 0 on failure
*/
int RicGraph::Read(string fname)
{
	ifstream gfile;
	gfile.open(fname.c_str());
	char tmp[512];
	string tok[30];

	if (!gfile.is_open())
	{
		cerr<<"Couldn't open file "<<fname<<" exiting"<<endl;
		return(0);
	}

	// check to see if really a graph file
	gfile.getline(tmp,511);
	string stmp=tmp;

	if (stmp.find("graph", 0)>500)
	{
		cerr<< "This file isn't a graph"<<endl;
		return 0;
	}


	// read the header

	// get a  line
	gfile.getline(tmp,511);

	// skip blanks
	int n;
	while( (n=tokenize(tmp,tok)) == 0 )
	{
		gfile.getline(tmp,511);
	}

	// we should now have the start of the graph header
	if ( n!=3 || tok[1].find("GRAPH")==string::npos)
	{
		cerr<< "Invalid graph header"<<endl;
		return 0;
	}

	this->gname = tok[2];

	// in a tedious manner parse the header and assign values to members
	gfile.getline(tmp,511);
	while ((n=tokenize(tmp,tok)) != 0 )
	{
		if ( tok[0] == "boundingbox_max" )
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
		else if( tok[0] == "voxel_size" )
		{
			voxel_size.x = (float)atof(tok[1].c_str());
			voxel_size.y = (float)atof(tok[2].c_str());
			voxel_size.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "CorticalFoldArg_VERSION" )
		{
			CorticalFoldArg_VERSION = tok[1];
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
		else if( tok[0] == "Talairach_rotation" )
		{
			Talairach_rotation[0] = (float)atof(tok[1].c_str());
			Talairach_rotation[1] = (float)atof(tok[2].c_str());
			Talairach_rotation[2] = (float)atof(tok[3].c_str());
			Talairach_rotation[3] = (float)atof(tok[4].c_str());
			Talairach_rotation[4] = (float)atof(tok[5].c_str());
			Talairach_rotation[5] = (float)atof(tok[6].c_str());
			Talairach_rotation[6] = (float)atof(tok[7].c_str());
			Talairach_rotation[7] = (float)atof(tok[8].c_str());
		}
		else if( tok[0] == "Talairach_scale" )
		{
			Talairach_scale.x = (float)atof(tok[1].c_str());
			Talairach_scale.y = (float)atof(tok[2].c_str());
			Talairach_scale.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "Talairach_translation" )
		{
			Talairach_translation.x = (float)atof(tok[1].c_str());
			Talairach_translation.y = (float)atof(tok[2].c_str());
			Talairach_translation.z = (float)atof(tok[3].c_str());
		}
		else if( tok[0] == "Tmtktri_label" )
		{
			Tmtktri_label.x = atoi(tok[1].c_str());
			Tmtktri_label.y = atoi(tok[2].c_str());
			Tmtktri_label.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "anterior_commissure" )
		{
			anterior_commissure.x = atoi(tok[1].c_str());
			anterior_commissure.y = atoi(tok[2].c_str());
			anterior_commissure.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "bottom_label" )
		{
			bottom_label.x = atoi(tok[1].c_str());
			bottom_label.y = atoi(tok[2].c_str());
			bottom_label.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "datagraph_VERSION" )
		{
			datagraph_VERSION = tok[1];
		}
		else if( tok[0] == "datagraph_compatibility_model_VERSION" )
		{
			datagraph_compatibility_model_VERSION = tok[1];
		}
		else if( tok[0] == "filename_base" )
		{
			filename_base = tok[1];
		}
		else if( tok[0] == "interhemi_point" )
		{
			interhemi_point.x = atoi(tok[1].c_str());
			interhemi_point.y = atoi(tok[2].c_str());
			interhemi_point.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "other_label" )
		{
			other_label.x = atoi(tok[1].c_str());
			other_label.y = atoi(tok[2].c_str());
			other_label.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "posterior_commissure" )
		{
			posterior_commissure.x = atoi(tok[1].c_str());
			posterior_commissure.y = atoi(tok[2].c_str());
			posterior_commissure.z = atoi(tok[3].c_str());
		}
		else if( tok[0] == "ss_label" )
		{
			ss_label.x = atoi(tok[1].c_str());
			ss_label.y = atoi(tok[2].c_str());
			ss_label.z = atoi(tok[3].c_str());
		}

		// get another line to parse
		gfile.getline(tmp,511);

	}

	// if here then time to read nodes
	nnodes = 0;
	if ( nodes ) delete [] nodes;
	nodes = new RicGraphNode[MAXNODE];

	// read nodes while they last
	while( nodes[nnodes].Read(&gfile) )
	{
		++nnodes;
		if ( nnodes >= MAXNODE )
		{
			cerr<< "Too many nodes"<<endl;
			return 0;
		}
	}

	return 1;
}

/*!
This method searches for a node with the label field containing
the passed string.

@param labelname - label name to search for

@returns - array index of node containing label if found, -1 if not found
*/
int RicGraph::FindNodeLabel(string labelname)
{
	node_idx = 0;

	for ( int i=0 ; i<nnodes ; ++i )
	{
		// see if labelname is contained in the label
		if ( nodes[i].label.find(labelname) != string::npos )
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for the next node with the label field containing
the passed string. The search is started after the last found node indicated
by member variable "node_index".

@param labelname - label name to search for

@returns - array index of node containing label if found, -1 if not found
*/
int RicGraph::FindNextNodeLabel(string labelname)
{
	for ( int i=node_idx+1 ; i<nnodes ; ++i )
	{
		if ( nodes[i].label.find(labelname) != string::npos )
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for a node with the label field exactly matching
the passed string.

@param labelname - label name to search for

@returns - array index of node containing label if found, -1 if not found
*/
int RicGraph::FindNodeLabelExact(string labelname)
{
	node_idx = 0;

	for ( int i=0 ; i<nnodes ; ++i )
	{
		if ( nodes[i].label == labelname)
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for the next node with the label field exactly matching
the passed string. The search is started after the last found node indicated
by member variable "node_index".

@param labelname - label name to search for

@returns - array index of node containing label if found, -1 if not found
*/
int RicGraph::FindNextNodeLabelExact(string labelname)
{
	for ( int i=node_idx+1 ; i<nnodes ; ++i )
	{
		if ( nodes[i].label == labelname)
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for a node with the name field containing
the passed string.

@param name - name to search for

@returns - array index of node containing name if found, -1 if not found
*/
int RicGraph::FindNodeName(string name)
{
	node_idx = 0;

	for ( int i=0 ; i<nnodes ; ++i )
	{
		// see if labelname is contained in the label
		if ( nodes[i].name.find(name) != string::npos )
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for the next node with the name field containing
the passed string. The search is started after the last found node indicated
by member variable "node_index".

@param name - name to search for

@returns - array index of node containing name if found, -1 if not found
*/
int RicGraph::FindNextNodeName(string name)
{
	for ( int i=node_idx+1 ; i<nnodes ; ++i )
	{
		if ( nodes[i].name.find(name) != string::npos )
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for a node with the name field exactly matching
the passed string.

@param name - name to search for

@returns - array index of node containing name if found, -1 if not found
*/
int RicGraph::FindNodeNameExact(string name)
{
	node_idx = 0;

	for ( int i=0 ; i<nnodes ; ++i )
	{
		if ( nodes[i].name == name)
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}

/*!
This method searches for the next node with the name field exactly matching
the passed string. The search is started after the last found node indicated
by member variable "node_index".

@param name - name to search for

@returns - array index of node containing name if found, -1 if not found
*/
int RicGraph::FindNextNodeNameExact(string name)
{
	for ( int i=node_idx+1 ; i<nnodes ; ++i )
	{
		if ( nodes[i].name == name)
		{
			node_idx = i;	// set node index to selected node
			return i;		// return node number
		}
	}

	// if it got here then it did not find the node so return error
	return -1;
}



// ------------------- function to break strings into tokens -----------------
/*!
This routine breaks a string into tokens (strings) separated by delimiters. A
string array is passed to the routine.

@param str - string to parse into tokens
@param tokens - array of strings to contain tokens
@param delimiters - string containing delimiters (default space)

@returns - number of tokens parsed or 0 if none found
*/
int tokenize(const string& str, string *tokens,	const string& delimiters)
{
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);
	int i=0;
	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens[i]=str.substr(lastPos, pos - lastPos);
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
		i++;
		if (i>=9)
			return i;
	}
	return i;
}

