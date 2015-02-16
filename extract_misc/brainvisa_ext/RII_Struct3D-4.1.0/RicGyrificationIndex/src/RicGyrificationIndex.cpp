// ------------------------ RicGyrificationIndex.cpp -------------------------

/*! \file
Implementation file for RicGyrificationIndex
Copyright (C) 2008 by Bill Rogers - Research Imaging Center - UTHSCSA 

\mainpage
This program determines gyrification index of a cortical mesh by comparing the
area of the mesh to the area of its convex hull. An optional ventricle volume
file can be used to remove the ventricles from both the cortical mesh and the
convex hull. If a ventricle volume file is specified then it is useful to super
sample the meshes to make a closer match to the edges of the ventricles. There
is also an optional scale factor to scale the meshes. BIG Warning ... The 
ventricle volume voxel space must correspond to the mesh volume. There is also
an optional transformation matrix for PeterK to scale the meshes after they
have been culled by the ventricle volume.

Command line switches

--gm	gray matter mesh file (required)

--cm	convex hull mesh file (required)

--vv	ventricle volume

--ta	triangle area

--tm	transformation matrix

-o		output file base name

-s		subject name

-v		verbose output

--bv	output for BrainVisa
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include "RicMesh.h"
#include "RicTexture.h"
#include "RicVolumeSet.h"
#include "RicUtil.h"
#include "RicMatrix.h"

using namespace std;

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string gm_name;				// GM mesh name
	string cm_name;				// convex hull mesh name
	string vv_name;				// ventricle volume name
	string outname;				// base name for output files
	string matname;				// name for transformation matrix file
	string subjname;			// subject name
	float max_area=50;			// maximum triangle area
	bool verbose=false;			// if true output more stuff on stdout
	bool brainvisa=false;		// if true output stuff on stdout for brain visa
	
	///////// read in the command line - tclap ///////////////

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("Acme-Presto Gyrification Index Processor", ' ', "0.5");
		
		// gray matter mesh
		TCLAP::ValueArg<string> gmname( "","gm", "Gray matter (outer) mesh file", true, "", "string");
		cmd.add(gmname);
		
		// convex hull mesh
		TCLAP::ValueArg<string> cmname( "","cm", "Convex hull mesh file", true, "", "string");
		cmd.add(cmname);
		
		// ventricle volume file
		TCLAP::ValueArg<string> vvname( "","vv", "Ventricle volume file", false, "", "string");
		cmd.add(vvname);
		
		// triangle area
		TCLAP::ValueArg<float> ta( "","ta", "Maximum triangle area", false, 50.0, "float");
		cmd.add(ta);
		
		// transformation matrix name
		TCLAP::ValueArg<string> mname( "","tm", "Transformation matrix file name", false,\
			 "", "string");
		cmd.add(mname);
	
		// output file base name
		TCLAP::ValueArg<string> oname( "o","o", "Base output file name", false, "Gidx", "string");
		cmd.add(oname);
	
		// subject name name just for peterk
		TCLAP::ValueArg<string> sname( "s","s", "Subject name", false, "FBSubject", "string");
		cmd.add(sname);
		
		TCLAP::SwitchArg ver( "v","verbose", "Verbose output on stdout", false);
		cmd.add(ver);
		
		TCLAP::SwitchArg bvisa( "b","bv", "Brain Visa output on stdout", false);
		cmd.add(bvisa);
		
		// parse the command line
		cmd.parse( argc, argv );
		
		// copy command line variable to program variables
		gm_name = gmname.getValue();
		cm_name = cmname.getValue();
		vv_name = vvname.getValue();
		max_area = ta.getValue();
		outname = oname.getValue();
		matname = mname.getValue();
		subjname = sname.getValue();
		verbose = ver.getValue();
		brainvisa = bvisa.getValue();
	
	}catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cerr << " Command Line Error" << e.error() << " for arg " << e.argId() << endl; 
		if ( brainvisa ) cout << " Command Line Error" << e.error() << " for arg " << e.argId() << endl; 
	}

	////////////// Sanity checks on command line input ///////////

	/// Check for valid GM mesh file
	RicMesh gm_mesh((char*)gm_name.c_str());
	if ( gm_mesh.p_size != 0 )
	{
		if ( verbose ) cerr << "We just read in "<<gm_name<<endl;
		if ( verbose ) cerr << "Number of vertices read: "<<gm_mesh.v_size<<endl;
	}
	else
	{
		cerr << "Error reading file"<<gm_name<<endl;
		if ( brainvisa ) cout << "Error reading file"<<gm_name<<endl;
		exit(1);
	}
	
	/// Check for valid convex hull mesh file
	RicMesh cm_mesh((char*)cm_name.c_str());
	if ( cm_mesh.p_size != 0 )
	{
		if ( verbose ) cerr << "We just read in "<<cm_name<<endl;
		if ( verbose ) cerr << "Number of vertices read: "<<cm_mesh.v_size<<endl;
	}
	else
	{
		cerr << "Error reading file"<<cm_name<<endl;
		if ( brainvisa ) cout << "Error reading file"<<cm_name<<endl;
		exit(1);
	}

	/// Read the Transformation Matrix
	bool	transflag=false; // indicates presence of ventricle volume
	CMatrix tm;		// transformation matrix
	if ( matname.length() )
	{
		ifstream ifile;
		ifile.open (matname.c_str(), ifstream::in);
		if ( !ifile.is_open() )
		{
			cerr<<"Error reading transform matrix " << matname << endl;
			exit(1);
		}
		
		// read in the matrix to temp structure
		CMatrix tmat;
		for ( int i=0 ; i<16 ; ++i )
			ifile >> tmat.mf[i];	
		ifile.close();
		if ( verbose ) cerr << "We just read in "<<matname<<endl;
		
		// transpose matrix to get proper orientation
		tm = tmat.Transpose();
		transflag = true;
	}
	
	/// Read the ventricle volume
	RicVolumeSet *ivol;
	RicVolume	*v_vol;
	bool	vvflag=false; // indicates presence of ventricle volume
	if ( vv_name.length() )
	{
		ivol = new RicVolumeSet(vv_name);
		if ( ivol->filetype == RIC_NO_FILE )
		{
			cout << "RicGyrificationIndex - Unknown file type: " << vv_name << endl;
			exit(1);
		}
		
		ivol->VolSet[0].CalcMinMaxAvg();
		
		if ( verbose ) cerr << "We just read in "<<vv_name<<endl;
		
		// disclaimer for Analyze
		if ( ivol->filetype == RIC_ANALYZE_FILE )
		{
			cout << "Analyze file so all bets are off on orientation!" << endl;
		}
		v_vol = &ivol->VolSet[0];
		v_vol->CalcMinMaxAvg();
		vvflag = true;
	}

// //////////////////////////////////////////////////////////////////////////
/// If there is a ventricle volume file then we need to do a bunch of stuff:
	RicMesh *gm_mesh2=NULL, *gm_mesh3=NULL, *cm_mesh2=NULL, *cm_mesh3=NULL;
	RicMesh *gm_cull=NULL, *cm_cull=NULL;
	
	if ( vvflag ) // we have the ventricle volume so tune up the meshes
	{
	
/// Fix and super sample meshes

		
		gm_mesh2 = gm_mesh.triangle_fix(5);
		gm_mesh3 = gm_mesh2->adjust_triangle_size(10,max_area);
		
		cm_mesh2 = cm_mesh.triangle_fix(5);
		cm_mesh3 = cm_mesh2->adjust_triangle_size(10,max_area);
	
/// Labeling polygons in ventricles

		int nculled=0;	// number of polygons tossed out
		int nculled2 = 0;
		
		// check each polygon to see if in ventricles
		for ( int i=0 ; i< gm_mesh3->p_size ; ++i )
		{
			// check each vertex to see if it is in the ventricle mask (triangles only)
			for ( int j=0 ; j<3 ; ++j )
			{
				int idx = gm_mesh3->polygons[i].vidx[j];
				
				if ( v_vol->interp3D(gm_mesh3->vertices[idx].pnt.x,gm_mesh3->vertices[idx].pnt.y,
					gm_mesh3->vertices[idx].pnt.z ) )
				{
					gm_mesh3->polygons[i].labeled = 1;
					++nculled;
					break;
				}
			}
		}
		
		// check each polygon to see if in ventricles
		for ( int i=0 ; i< cm_mesh3->p_size ; ++i )
		{
			// check each vertex to see if it is in the ventricle mask (triangles only)
			for ( int j=0 ; j<3 ; ++j )
			{
				int idx = cm_mesh3->polygons[i].vidx[j];
				
				if ( v_vol->interp3D(cm_mesh3->vertices[idx].pnt.x,cm_mesh3->vertices[idx].pnt.y,
					cm_mesh3->vertices[idx].pnt.z ) )
				{
					cm_mesh3->polygons[i].labeled = 1;
					++nculled2;
					break;
				}
			}
		}
	
/// Tweak meshes to remove labeled polygons

		gm_cull = new RicMesh(gm_mesh3->v_size,gm_mesh3->n_size,gm_mesh3->p_size - nculled,gm_mesh3->p_dim);
		
		// fill up mesh with values from gm_mesh
		for ( int i=0 ; i < gm_cull->v_size ; ++i )
			gm_cull->vertices[i] = gm_mesh3->vertices[i];
			
		for ( int i=0 ; i < gm_cull->n_size ; ++i )
			gm_cull->normals[i] = gm_mesh3->normals[i];
			
		// only put the unlabled (non culled) polygons in the new mesh
		int j=0;
		for ( int i = 0 ; i < gm_mesh3->p_size ; ++i )
		{
			if ( gm_mesh3->polygons[i].labeled == 0 )
			{
				gm_cull->polygons[j] = gm_mesh3->polygons[i];
				++j;
			}
		}
	
		cm_cull = new RicMesh(cm_mesh3->v_size,cm_mesh3->n_size,cm_mesh3->p_size - nculled2,cm_mesh3->p_dim);
		
		// fill up mesh with values from cm_mesh
		j=0;
		for ( int i=0 ; i < cm_cull->v_size ; ++i )
			cm_cull->vertices[i] = cm_mesh3->vertices[i];
			
		for ( int i=0 ; i < cm_cull->n_size ; ++i )
			cm_cull->normals[i] = cm_mesh3->normals[i];
			
		// only put the unlabled (non culled) polygons in the new mesh
		for ( int i = 0 ; i < cm_mesh3->p_size ; ++i )
		{
			if ( cm_mesh3->polygons[i].labeled == 0 )
			{
				cm_cull->polygons[j] = cm_mesh3->polygons[i];
				++j;
			}
		}

/// Write tweaked meshes to output files (optional)

		if ( verbose || brainvisa )
		{
			string filename;
			filename = outname + subjname + "_gmcull.mesh";
			
			gm_cull->write_mesh_bin((char*)filename.c_str());
			
			filename = outname + subjname + "_cmcull.mesh";
			
			cm_cull->write_mesh_bin((char*)filename.c_str());
		}	
		
	} /// End of ventricle volume stuff

// ///////////////////////////////////////////////////////////////////////////
/// Optional transformation matrix stuff to correct for scaling

	if ( transflag )
	{
		// transform the original meshes
		for ( int i=0 ; i<gm_mesh.v_size ; ++i )
		{
			// first offset mesh coordinates to center of brain
			gm_mesh.vertices[i].pnt = tm*gm_mesh.vertices[i].pnt;
		}
		gm_mesh.calc_limits();
	
		for ( int i=0 ; i<cm_mesh.v_size ; ++i )
		{
			// first offset mesh coordinates to center of brain
			cm_mesh.vertices[i].pnt = tm*cm_mesh.vertices[i].pnt;
		}
		cm_mesh.calc_limits();
	
		if ( vvflag )
		{
			// transform the culled meshes
			for ( int i=0 ; i<gm_cull->v_size ; ++i )
			{
				// first offset mesh coordinates to center of brain
				gm_cull->vertices[i].pnt = tm*gm_cull->vertices[i].pnt;
			}
			gm_cull->calc_limits();
		
			for ( int i=0 ; i<cm_cull->v_size ; ++i )
			{
				// first offset mesh coordinates to center of brain
				cm_cull->vertices[i].pnt = tm*cm_cull->vertices[i].pnt;
			}
			cm_cull->calc_limits();
		}	
	}
	
// ///////////////////////////////////////////////////////////////////////////
/// Calculate Gyrification Index from mesh areas

	float areag,areac,gi;
	
	// calculate area and GI from original meshes
	areag = gm_mesh.area();
	areac = cm_mesh.area();
	gi = areag/areac;
	
	if ( vvflag ) // calculate and print gi from culled meshes
	{
		if ( verbose || brainvisa )
		{
			printf("Uncorrected for ventricle volume\n"); 
			printf("gray matter area %8.1f\n convex hull area %8.1f\n gyrification index %5.3f\n",areag,areac,gi);
		}
		areag = gm_cull->area();
		areac = cm_cull->area();
		gi = areag/areac;
		if ( verbose || brainvisa ) 
		{
			printf("Corrected for ventricle volume\n"); 
			printf("gray matter area %8.1f\n convex hull area %8.1f\n gyrification index %5.3f\n",areag,areac,gi);
		}
	}
	else // print GI from original meshes
	{
		if ( verbose || brainvisa ) 
			printf("gray matter area %8.1f\n convex hull area %8.1f\n gyrification index %5.3f\n",areag,areac,gi);
	}
	
	if ( !brainvisa )
		printf("%s\t%6.4f\n",subjname.c_str(),gi);
		
	if ( vvflag )
	{
		delete gm_mesh2;
		delete gm_mesh3;
		delete cm_mesh2;
		delete cm_mesh3;
		delete gm_cull;
		delete cm_cull;
	}
	
  return EXIT_SUCCESS;
}
