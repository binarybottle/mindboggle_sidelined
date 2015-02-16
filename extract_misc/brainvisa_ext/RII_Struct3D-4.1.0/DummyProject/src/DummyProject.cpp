/// //////////////////////// DummyProject.cpp ///////////////////////////

/*! \file
Implementation file for DummyProject.
Copyright (C) 2008 by Bill Rogers - Research Imaging Center - UTHSCSA

This program is an empty project that can be built upon.

Command line switches

--gm	gray matter mesh file (required)

 */



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include <iostream>
#include <cstdlib>
#include "RicMesh.h"

using namespace std;

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string gm_name;				// GM mesh name
	
	///////// read in the command line - tclap ///////////////

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("Dummy Project", ' ', "1.0");
		
		// gray matter mesh
		TCLAP::ValueArg<string> gmname( "","gm", "Gray matter (outer) mesh file", true, "", "string");
		cmd.add(gmname);
		
		// parse the command line
		cmd.parse( argc, argv );
		
		// copy command line variable to program variables
		gm_name = gmname.getValue();
	
	}catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cerr << " Command Line Error" << e.error() << " for arg " << e.argId() << endl; 
	}

	////////////// Sanity checks on command line input ///////////

	// check for valid GM mesh file
	RicMesh gm_mesh((char*)gm_name.c_str());
	if ( gm_mesh.p_size != 0 )
	{
		cerr << "We just read in "<<gm_name<<endl;
		cerr << "Number of vertices read: "<<gm_mesh.v_size<<endl;
	}
	else
	{
		cerr << "Error reading file"<<gm_name<<endl;
		cerr << "Error reading file"<<gm_name<<endl;
		exit(1);
	}
	
  return EXIT_SUCCESS;
}
