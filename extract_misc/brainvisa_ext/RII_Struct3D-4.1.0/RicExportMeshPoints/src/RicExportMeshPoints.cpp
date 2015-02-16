// --------------------------- RicExportMeshPoints.cpp ------------------------
/*! \file
Implementation file for RicExportMeshPoints
Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA

This program reads a mesh then exports the mesh as a list of points. The 
output list can be a subset of the original mesh skipping a specified number
of mesh points for every output point. If no output file is specified, the name
will be the input file name with the extension .pnts. The output file is
tab seperated with a point on each line. The first line contains the number
of points output.


Command line switches

-i	input mesh

-o output points list

-n number of mesh points to skip per output point
 */


#include <iostream>
#include <cstdlib>
#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include <RicMesh.h>

using namespace std;

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string mesh_file_name;		// name of gyral texture file
	string points_file_name;	// name of output file
	int	nskip;					// number of points to skip
	
/// Read in the command line - tclap

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("RIC Export Mesh Points", ' ', "0.5");
		
		// input mesh file
		TCLAP::ValueArg<string> mfile( "i","i", "Mesh file", true, "", "string");
		cmd.add(mfile);
		
		// output points file
		TCLAP::ValueArg<string> pfile( "o","o", "Points file", false, "", "string");
		cmd.add(pfile);
		
		// number of points to skip
		TCLAP::ValueArg<int> ns( "n","n", "Number of points to skip", false, 2, "int");
		cmd.add(ns);
		
		// parse the command line
		cmd.parse( argc, argv );
		
		// copy command line variable to program variables
		mesh_file_name = mfile.getValue();
		points_file_name = pfile.getValue();
		nskip = ns.getValue();
	
	}catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cout << " Command Line Error " << e.error() << " for arg " << e.argId() << endl; 
	}

/// Check input
	if ( nskip <= 0 || nskip > 100 )
	{
		cout << "Error - nskip out of range " << nskip << endl;
		exit(1);
	}
	if ( points_file_name.length() == 0 ) // build a file name
	{
		int i=mesh_file_name.rfind(".");
		points_file_name = mesh_file_name.substr(0,i) + ".pnts";
	}
	
/// Read the mesh
	RicMesh mesh((char*)mesh_file_name.c_str());
	if ( mesh.v_size == 0 )
	{
		cout << "Error reading file"<<mesh_file_name<<endl;
		exit(1);
	}
	
/// Write the mesh to a text file
	ofstream fout(points_file_name.c_str(), ios::out);
	int nout = mesh.v_size/nskip;
	fout << nout << endl;
	for ( int i=0 ; i<mesh.v_size ; i+=nskip )
	{
		fout << mesh.vertices[i].pnt.x << "\t"
		<< mesh.vertices[i].pnt.y << "\t" 
		<< mesh.vertices[i].pnt.z << endl;
	}
	fout.close();

  return EXIT_SUCCESS;
}
