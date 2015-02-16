// -------------------- RicCorticalThicknessByNormal.cpp ----------------------
/*! \mainpage
Implementation file for RicCorticalThicknessByNormal.cpp
Copyright (C) 2007 by Bill Rogers - Research Imaging Institute - UTHSCSA

This program determines the thickness between two meshes using a
modified brute force method looking for the intersection of a
normal from a vertex from one mesh to a triangle in the other
mesh. There is a check to make sure that the normal from the
first mesh does not intersect itself before intersecting
the other mesh.

Limits on the thickness can be specified with min and max values. Thickness values
for a vertex outside the limits will be set to zero. In addition, input vertices can
be culled by curvature if a white matter curvature map is input as well.

This version will read in a transformation matrix to scale the meshes before calculation
of thickness.

Verbose output will output additional meshes to show input mesh surface normals
as well as a mesh showing connecting vectors between normal points.

Command line switches

--gm	gray matter mesh file (required)

--wm	white matter mesh file (required)

--tm	Transformation matrix file name

--gc	grey matter curvature file

--wc	white matter curvature file

--sm	surface to map on - white (default), grey, or both

-o		output file base name

--mind	minimum thickness allowed

--maxd	maximum thickness allowed

--minc	minimum curvature allowed

--maxc	maximum curvature allowed

--ft	flag for filling invalid nodes in texture

--fd	radial distance to average for filling texture nodes

-m		method of calculation (brute, sub, thread)

-v		verbose output

--bv	output for BrainVisa

--pk	special output for PeterK

 */

#include <string>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <time.h>
#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include "RicMesh.h"
#include "RicTexture.h"
#include "TexFill.h"
#include "GM_Normal.h"
#include "RicUtil.h"
#include "RicMatrix.h"
#include "RicVolumeSet.h"

// methods of calculation
#define BRUTE 1 	///< use brute force calculation
#define SUBDIVIDE 2	///< subdivide into smaller chunks
#define THREADS 3	///< subdivide with a thread for each chunk

// surface to map on
#define MAP_WHITE 1	///< map from gray to white
#define MAP_GRAY 2	///< map from white to gray
#define MAP_BOTH 3	///< map both directions

using namespace std;

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string gm_name;				// GM mesh name
	string wm_name;				// WM mesh name
	string gm_curve_name;		// GM curvature map name
	string wm_curve_name;		// WM curvature map name
	string surfmapstr;			// surface to map on
	int	   surfmap=MAP_WHITE;	// surface to map on, default is white matter
	string outname;				// base name for output files
	string matname;				// name for transformation matrix file
	float minc=0;				// minimum curvature value
	float maxc=0;				// maximum curvature value
	float mind=0.0f;			// minimum GM distance
	float maxd=10.0f;			// maximum GM distance
	float filld = 10.0f;		// radius to check for filling 	string methstr;
	bool filltex=false;			// if true then fill invald nodes
	int	method=THREADS;			// default method using threads
	string methstr;				// method
	bool transmat=false;		// if true then use transformation matrix
	bool verbose=false;			// if true output more stuff on stdout
	bool brainvisa=false;		// if true output for BrainVisa
	bool peterk=false;			// if true then add special processing for PeterK

	///////// read in the command line - tclap ///////////////

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("Acme-Presto Cortical Thickness Normal Distance Processor",\
			 ' ', "0.5");

		// gm mesh (outer mesh)
		TCLAP::ValueArg<string> gmname( "","gm", "Gray matter (outer) mesh file",\
			 true, "", "string");
		cmd.add(gmname);

		// wm mesh (outer mesh)
		TCLAP::ValueArg<string> wmname( "","wm", "White matter (inner) mesh file",\
			 true, "", "string");
		cmd.add(wmname);

		// gm curvature
		TCLAP::ValueArg<string> gmcurve( "","gc", "Grey matter (outer) curvature file",\
			 false, "", "string");
		cmd.add(gmcurve);

		// wm curvature
		TCLAP::ValueArg<string> wmcurve( "","wc", "White matter (inner) curvature file",\
			 false, "", "string");
		cmd.add(wmcurve);

		// surface to map on
		TCLAP::ValueArg<string> surmap( "","sm", "Surface to map on (grey,white,or both)",\
			 false, "", "string");
		cmd.add(surmap);

		// output file base name
		TCLAP::ValueArg<string> oname( "o","o", "Base output file name", false,\
			 "Normal_Thickness_Map", "string");
		cmd.add(oname);

		// transformation matrix name
		TCLAP::ValueArg<string> mname( "","tm", "Transformation matrix file name", false,\
			 "", "string");
		cmd.add(mname);

		TCLAP::ValueArg<float> mnc( "","minc", "Minimum curvature threshold",\
			 false, 0, "float");
		cmd.add(mnc);

		TCLAP::ValueArg<float> mxc( "","maxc", "Maximum curvature threshold",\
			 false, 0, "float");
		cmd.add(mxc);

		TCLAP::ValueArg<float> mnd( "","mind", "Minimum gray matter distance",\
			 false, 0, "float");
		cmd.add(mnd);

		TCLAP::ValueArg<float> mxd( "","maxd", "Maximum gray matter distance",\
			 false, 10, "float");
		cmd.add(mxd);

		TCLAP::SwitchArg ft( "","ft", "Fill invald nodes in thickness texture", false);
		cmd.add(ft);

		TCLAP::ValueArg<float> fd( "","fd", "Maximum distance to fill thickness texture",\
			false, 16, "float");
		cmd.add(fd);

		TCLAP::ValueArg<string> mth( "m","meth", \
			"Calculation method - brute, sub, or thread - default sub",
			false, "thread", "float");
		cmd.add(mth);

		TCLAP::SwitchArg ver( "v","verbose", "Verbose output on stdout", false);
		cmd.add(ver);

		TCLAP::SwitchArg bv( "","bv", "Brain Visa output on stdout", false);
		cmd.add(bv);

		TCLAP::SwitchArg pk( "","pk", "Special PeterK output on stdout", false);
		cmd.add(pk);

		// parse the command line
		cmd.parse( argc, argv );

		// copy command line variable to program variables
		gm_name = gmname.getValue();
		wm_name = wmname.getValue();
		gm_curve_name = gmcurve.getValue();
		wm_curve_name = wmcurve.getValue();
		surfmapstr = surmap.getValue();
		if ( surfmapstr.substr(0,4) == "gray" ) surfmap = MAP_GRAY;
		if ( surfmapstr.substr(0,4) == "both" ) surfmap = MAP_BOTH;
		outname = oname.getValue();
		matname = mname.getValue();
		minc = mnc.getValue();
		maxc = mxc.getValue();
		mind = mnd.getValue();
		maxd = mxd.getValue();
		filltex = ft.getValue();
		filld = fd.getValue();
		methstr = mth.getValue();
		method = THREADS;
		if ( methstr.substr(0,5) == "brute" ) method = BRUTE;
		if ( methstr.substr(0,3) == "sub" ) method = SUBDIVIDE;
		verbose = ver.getValue();
		brainvisa = bv.getValue();
		peterk = pk.getValue();

	}catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cout << " Command Line Error" << e.error() << " for arg " << e.argId() << endl;
	}

//////////////////////////////////////////////////////////////////////////////
/// Step 1 - Read Files

	////////////// Sanity checks on command line input ///////////
	RicTexture *wm_curve_map=NULL;	// WM curvature map
	RicTexture *gm_curve_map=NULL;	// GM curvature map

	// check for valid GM mesh file
	RicMesh gm_mesh((char*)gm_name.c_str());
	gm_mesh.calc_limits();
	if ( gm_mesh.p_size != 0 )
	{
		if ( verbose ) cout << "We just read in GM mesh "<<gm_name<<endl;
		if ( verbose ) cout << "Number of vertices read: "<<gm_mesh.v_size<<endl;
	}
	else
	{
		cout << "Error reading file"<<gm_name<<endl;
		exit(1);
	}

	// check for valid WM mesh file
	RicMesh wm_mesh((char*)wm_name.c_str());
	if ( wm_mesh.p_size != 0 )
	{
		if ( verbose ) cout << "We just read in WM mesh "<<wm_name<<endl;
		if ( verbose ) cout << "Number of vertices read: "<<wm_mesh.v_size<<endl;
	}
	else
	{
		cout << "Error reading file"<<wm_name<<endl;
		exit(1);
	}

	// check for GM curvature map
	if ( gm_curve_name[0] )
	{
		gm_curve_map = new RicTexture();
		gm_curve_map->read_texture((char*)gm_curve_name.c_str());
		gm_curve_map->CalcMinMaxAvg();
		if ( gm_curve_map->size != 0 )
		{
			if ( verbose )
			{
				cout << "We just read in GM curvature " <<gm_curve_name<<endl;
				cout << "Number of values read="<<gm_curve_map->size;
				cout << " Min="<<gm_curve_map->min<<" Max="<<gm_curve_map->max<<endl;
			}
		}
		else
		{
			cout << "Error reading file "<< gm_curve_name<<endl;
			exit(1);
		}
		// sanity check to see if the number of curvature values and the  number
		// of mesh vertices are the same
		if ( gm_curve_map->size != gm_mesh.v_size )
		{
			cout<<"GM curvature map and GM mesh difference sizes"<<endl;
			exit(1);
		}
	}

	// check for WM curvature map
	if ( wm_curve_name[0] )
	{
		wm_curve_map = new RicTexture();
		wm_curve_map->read_texture((char*)wm_curve_name.c_str());
		wm_curve_map->CalcMinMaxAvg();
		if ( wm_curve_map->size != 0 )
		{
			if ( verbose )
			{
				cout << "We just read in WM curvature " <<wm_curve_name<<endl;
				cout << "Number of values read="<<wm_curve_map->size;
				cout << " Min="<<wm_curve_map->min<<" Max="<<wm_curve_map->max<<endl;
			}
		}
		else
		{
			cout << "Error reading file "<< wm_curve_name<<endl;
			exit(1);
		}
		// sanity check to see if the number of curvature values and the  number
		// of mesh vertices are the same
		if ( wm_curve_map->size != wm_mesh.v_size )
		{
			cout<<"WM curvature map and WM mesh difference sizes"<<endl;
			exit(1);
		}
	}

	// make sure that the min and max curvatures are in the right order
	if ( minc != 0 || maxc != 0 )
	{
		if ( minc >= maxc )
		{
			cout<<"min_curve must be less than max_curve"<<endl;
			exit(1);
		}
	}

	// Read the Transformation Matrix  if there is one
	CMatrix tm;		// transformation matrix
	if ( matname.length() != 0 )
	{
		transmat = true; // we have a transformation matrix

		ifstream ifile;
		ifile.open (matname.c_str(), ifstream::in);
		if ( !ifile.is_open() )
		{
			cout<<"Error reading transform matrix " << matname << endl;
			exit(1);
		}


		// read in the matrix to temp structure
		CMatrix tmat;
		for ( int i=0 ; i<16 ; ++i )
			ifile >> tmat.mf[i];
		ifile.close();

		// transpose matrix to get proper orientation
		tm = tmat.Transpose();
	}

	//////////////////// Start the clock //////////////////////////

	// start a timer to see how long this all takes
	time_t starttime;
	time_t endtime;
	time(&starttime);

//////////////////////////////////////////////////////////////////////////////
///  Step 2 - Cull vertices with curvature map

	/// Vertices are culled by setting the vertex label value to 1. If it is
	/// set to one then that vertex will not be evaluated for thickness

	/// Cull GM mesh vertices with curvature map

	// zero the vertex labels
	for ( int i=0 ; i<gm_mesh.v_size ; ++i ) gm_mesh.vertices[i].label=0;

	// check for curvature in range
	int ngm_culled=0; 	// number of curvature values discarded
	int ngm_used=0;		// number of curvature values used
	float gm_avg_curvature=0;
	if ( gm_curve_map && (maxc!=0) )
	{
		gm_curve_map->CalcMinMaxAvg();

		// use vertex labels to cull vertices - if label set to 1 then do not use
		for ( int i=0 ; i<gm_mesh.v_size ; ++i )
		{
			if ( gm_curve_map->nodes[i] < minc || gm_curve_map->nodes[i] > maxc )
			{
				gm_mesh.vertices[i].label=1; // do not use
				++ngm_culled;
			}
			else // start added up the rectified curvature value
			{
				gm_avg_curvature += fabs(gm_curve_map->nodes[i]);
				++ngm_used;
			}
		}
	}

	// calculate average curvature
	if ( ngm_used ) // don't divide by zero
		gm_avg_curvature /= (float)ngm_used;
	else
		gm_avg_curvature = 0;

	/// Cull WM mesh vertices with curvature map

	// zero the vertex labels
	for ( int i=0 ; i<wm_mesh.v_size ; ++i ) wm_mesh.vertices[i].label=0;

	// check for curvature in range
	int nwm_culled=0; 	// number of curvature values discarded
	int nwm_used=0;		// number of curvature values used
	float wm_avg_curvature=0;
	if ( wm_curve_map && (maxc!=0) )
	{
		wm_curve_map->CalcMinMaxAvg();

		// use vertex labels to cull vertices - if label set to 1 then do not use
		for ( int i=0 ; i<wm_mesh.v_size ; ++i )
		{
			if ( wm_curve_map->nodes[i] < minc || wm_curve_map->nodes[i] > maxc )
			{
				wm_mesh.vertices[i].label=1; // do not use
				++nwm_culled;
			}
			else // start added up the rectified curvature value
			{
				wm_avg_curvature += fabs(wm_curve_map->nodes[i]);
				++nwm_used;
			}
		}
	}

	// calculate average curvature
	if ( nwm_used ) // don't divide by zero
		wm_avg_curvature /= (float)nwm_used;
	else
		wm_avg_curvature = 0;

//////////////////////////////////////////////////////////////////////////////
/// Step 3 - Calculate Thickness

	/// The thickness can be either calculated from the white to the gray
	/// (default) where the texture thickness texture is mapped on the white
	/// or from the gray to the white with the thickness texture mapped on
	/// the gray. Both directions is an option too.

	// find out how big the meshes are
	gm_mesh.calc_limits();
	wm_mesh.calc_limits();

	// it is assumed that the the GM mesh will be larger than the WM mesh
	// if they have been swapped then we need to swap the normal direction
	// to make this work
	int nflip = 1;
	if ( (gm_mesh.xmax-gm_mesh.xmin) < (wm_mesh.xmax-wm_mesh.xmin) ) nflip = -1;

	// textures
	RicTexture *thick_map_wm = NULL;
	RicTexture *thick_map_gm = NULL;

	// connection vector meshes
	int cvpoly_size_wm = wm_mesh.v_size;
	int cvpoly_size_gm = gm_mesh.v_size;
	RicMesh cv_mesh_wm(2*cvpoly_size_wm, 2*cvpoly_size_wm, cvpoly_size_wm, 2);
	RicMesh cv_mesh_gm(2*cvpoly_size_gm, 2*cvpoly_size_gm, cvpoly_size_gm, 2);

	if ( surfmap == MAP_WHITE || surfmap == MAP_BOTH) // map thickness on white
	{

		// create the texture that will hold the per vertex thickness values
		thick_map_wm = new RicTexture(wm_mesh.v_size);

		// make mesh to show connecting points between surfaces
		if ( verbose ) cout << "Making white matter texture ";

		// populate mesh with short vectors using gm mesh points
		for ( int i=0 ; i<cvpoly_size_wm ; ++i )
		{
			cv_mesh_wm.assign_node(2*i, wm_mesh.vertices[i]);
			cv_mesh_wm.assign_normal(2*i, wm_mesh.normals[i]);
			cv_mesh_wm.assign_node(2*i+1, wm_mesh.vertices[i]);
			cv_mesh_wm.assign_normal(2*i+1, wm_mesh.normals[i]);
			cv_mesh_wm.assign_polygon(i, 2*i, 2*i+1, 0);
			//cv_mesh_wm.vertices[2*i+1].pnt.x += 0.2f; // make tiny vector that looks like point
		}

		// pick method to use
		if ( method == BRUTE ) // check every wm
		{
			if ( verbose ) cout<<"using brute force"<<endl;
			FindNormalDistBrute(&wm_mesh, &gm_mesh, &cv_mesh_wm, thick_map_wm, nflip, mind, maxd);
		}
		else if ( method == SUBDIVIDE )
		{
			if ( verbose ) cout<<"using subdivision"<<endl;
			FindNormalDistSubdivide(&wm_mesh, &gm_mesh, &cv_mesh_wm, thick_map_wm, 2, maxd, nflip, mind, maxd);
		}
		else // the default
		{
			if ( verbose ) cout<<"using threads"<<endl;
			FindNormalDistThreads(&wm_mesh, &gm_mesh, &cv_mesh_wm, thick_map_wm, 2, maxd, nflip, mind, maxd);
		}

		// calculate the limits on the vector mesh
		cv_mesh_wm.calc_limits();
	}
	if ( surfmap == MAP_GRAY || surfmap == MAP_BOTH) // map on gray
	{

		// create the texture that will hold the per vertex thickness values
		thick_map_gm = new RicTexture(gm_mesh.v_size);

		// make mesh to show connecting points between surfaces
		if ( verbose ) cout << "Making gray matter texture ";

		// populate mesh with short vectors using gm mesh points
		for ( int i=0 ; i<cvpoly_size_gm ; ++i )
		{
			cv_mesh_gm.assign_node(2*i, gm_mesh.vertices[i]);
			cv_mesh_gm.assign_normal(2*i, gm_mesh.normals[i]);
			cv_mesh_gm.assign_node(2*i+1, gm_mesh.vertices[i]);
			cv_mesh_gm.assign_normal(2*i+1, gm_mesh.normals[i]);
			cv_mesh_gm.assign_polygon(i, 2*i, 2*i+1, 0);
			//cv_mesh_gm.vertices[2*i+1].pnt.x += 0.2f; // make tiny vector that looks like point
		}

		// pick method to use
		if ( method == BRUTE ) // check every wm
		{
			if ( verbose ) cout<<"using brute force"<<endl;
			FindNormalDistBrute(&gm_mesh, &wm_mesh, &cv_mesh_gm, thick_map_gm, -1*nflip, mind, maxd);
		}
		else if ( method == SUBDIVIDE )
		{
			if ( verbose ) cout<<"using subdivision"<<endl;
			FindNormalDistSubdivide(&gm_mesh, &wm_mesh, &cv_mesh_gm, thick_map_gm, 2, maxd, -1*nflip, mind, maxd);
		}
		else // the default
		{
			if ( verbose ) cout<<"using threads"<<endl;
			FindNormalDistThreads(&gm_mesh, &wm_mesh, &cv_mesh_gm, thick_map_gm, 2, maxd, -1*nflip, mind, maxd);
		}

		// calculate the limits on the vector mesh
		cv_mesh_gm.calc_limits();
	}

//////////////////////////////////////////////////////////////////////////////
/// Step 4 - Untransform mesh

	/// If there is a transformation matrix then use inverse transform to
	/// create untransformed meshes. Later we will also untransform the vectors
	/// connecting the normal points as well.

	// make untransformed meshes
	RicMesh wm_mesh2(wm_mesh.v_size,wm_mesh.n_size,wm_mesh.p_size,wm_mesh.p_dim);
	RicMesh gm_mesh2(gm_mesh.v_size,gm_mesh.n_size,gm_mesh.p_size,gm_mesh.p_dim);

	// make mesh to show connecting points between surfaces
	RicMesh cv_mesh_wm2(2*cvpoly_size_wm, 2*cvpoly_size_wm, cvpoly_size_wm, 2);
	RicMesh cv_mesh_gm2(2*cvpoly_size_gm, 2*cvpoly_size_gm, cvpoly_size_gm, 2);

	// textures for untransformed meshes
	RicTexture *thick_map_wm2=NULL;
	RicTexture *thick_map_gm2=NULL;

	if ( transmat )
	{
		// White
		// make a new mesh for untransformed data

		// copy vertices normals and polygons
		for ( int i=0 ; i<wm_mesh.v_size ; ++i )
			wm_mesh2.vertices[i] = wm_mesh.vertices[i];

		for ( int i=0 ; i<wm_mesh.n_size ; ++i )
			wm_mesh2.normals[i] = wm_mesh.normals[i];

		for ( int i=0 ; i<wm_mesh.p_size ; ++i )
			wm_mesh2.polygons[i] = wm_mesh.polygons[i];


		for ( int i=0 ; i<wm_mesh.v_size ; ++i )
		{
			// first offset mesh coordinates to center of brain
			wm_mesh2.vertices[i].pnt = tm*wm_mesh2.vertices[i].pnt;
		}
		wm_mesh2.calc_limits();

		// gray
		// make a new mesh for untransformed data

		// copy vertices normals and polygons
		for ( int i=0 ; i<gm_mesh.v_size ; ++i )
			gm_mesh2.vertices[i] = gm_mesh.vertices[i];

		for ( int i=0 ; i<gm_mesh.n_size ; ++i )
			gm_mesh2.normals[i] = gm_mesh.normals[i];

		for ( int i=0 ; i<gm_mesh.p_size ; ++i )
			gm_mesh2.polygons[i] = gm_mesh.polygons[i];


		for ( int i=0 ; i<gm_mesh.v_size ; ++i )
		{
			// first offset mesh coordinates to center of brain
			gm_mesh2.vertices[i].pnt = tm*gm_mesh2.vertices[i].pnt;
		}
		gm_mesh2.calc_limits();


		/// Now we untransform the normal vectors connecting the meshes.
		/// We are assuming that by transforming the vectors that connect
		/// the meshes they will match up to the untransformed meshes.

		if ( surfmap == MAP_WHITE || surfmap == MAP_BOTH) // map thickness on white
		{
			// create the texture that will hold the per vertex thickness values
			thick_map_wm2 = new RicTexture(wm_mesh2.v_size);

			// populate mesh with short vectors using wm mesh points
			for ( int i=0 ; i<cvpoly_size_wm ; ++i )
			{
				cv_mesh_wm2.assign_node(2*i, wm_mesh2.vertices[i]);
				cv_mesh_wm2.assign_normal(2*i, wm_mesh2.normals[i]);
				cv_mesh_wm2.assign_node(2*i+1, wm_mesh2.vertices[i]);
				cv_mesh_wm2.assign_normal(2*i+1, wm_mesh2.normals[i]);
				cv_mesh_wm2.assign_polygon(i, 2*i, 2*i+1, 0);
			}

			// transform the vertices
			for ( int i=0 ; i<cv_mesh_wm2.v_size ; ++i )
			{
				cv_mesh_wm2.vertices[i].pnt = tm*cv_mesh_wm.vertices[i].pnt;
			}
			cv_mesh_wm2.calc_limits();

			// Calculate distances between ends of transformed closest vectors
			// to get thickness
			for ( int i=0 ; i<wm_mesh2.v_size ; ++i )
			{
				if ( thick_map_wm->nodes[i] != ERRVAL )
				{
					Point p1,p2;
					p1 = cv_mesh_wm2.vertices[cv_mesh_wm2.polygons[i].vidx[0]].pnt;
					p2 = cv_mesh_wm2.vertices[cv_mesh_wm2.polygons[i].vidx[1]].pnt;
					thick_map_wm2->nodes[i] = dist(p1,p2);
				}
				else
				{
					thick_map_wm2->nodes[i] = ERRVAL;
				}
			}
		}
		if ( surfmap == MAP_GRAY || surfmap == MAP_BOTH) // map thickness on gray
		{
			// create the texture that will hold the per vertex thickness values
			thick_map_gm2 = new RicTexture(gm_mesh2.v_size);

			// populate mesh with short vectors using gm mesh points
			for ( int i=0 ; i<cvpoly_size_gm ; ++i )
			{
				cv_mesh_gm2.assign_node(2*i, gm_mesh2.vertices[i]);
				cv_mesh_gm2.assign_normal(2*i, gm_mesh2.normals[i]);
				cv_mesh_gm2.assign_node(2*i+1, gm_mesh2.vertices[i]);
				cv_mesh_gm2.assign_normal(2*i+1, gm_mesh2.normals[i]);
				cv_mesh_gm2.assign_polygon(i, 2*i, 2*i+1, 0);
			}

			// transform the vertices
			for ( int i=0 ; i<cv_mesh_gm2.v_size ; ++i )
			{
				cv_mesh_gm2.vertices[i].pnt = tm*cv_mesh_gm.vertices[i].pnt;
			}
			cv_mesh_gm2.calc_limits();

			// Calculate distances between ends of transformed closest vectors
			// to get thickness
			for ( int i=0 ; i<gm_mesh2.v_size ; ++i )
			{
				if ( thick_map_gm->nodes[i] != ERRVAL )
				{
					Point p1,p2;
					p1 = cv_mesh_gm2.vertices[cv_mesh_gm2.polygons[i].vidx[0]].pnt;
					p2 = cv_mesh_gm2.vertices[cv_mesh_gm2.polygons[i].vidx[1]].pnt;
					thick_map_gm2->nodes[i] = dist(p1,p2);
				}
				else
				{
					thick_map_gm2->nodes[i] = ERRVAL;
				}
			}
		}

	} // end of untransform

//////////////////////////////////////////////////////////////////////////////
/// Step 5 - Calculate texture stats

	// texture maps for stat image
	RicTexture *stat_map_wm=NULL,*stat_map_wm2=NULL;
	RicTexture *stat_map_gm=NULL,*stat_map_gm2=NULL;

	if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // stat map on white
	{
		stat_map_wm = new RicTexture(thick_map_wm->size);

		// copy values from thickness map flagging all nodes outside range
		// also zero corresponding vectors in connecting mesh
		for (int i = 0; i < thick_map_wm->size; ++i)
		{
			if (thick_map_wm->nodes[i] > mind && thick_map_wm->nodes[i] < maxd)
			{
				stat_map_wm->nodes[i] = thick_map_wm->nodes[i];
			}
			else
			{
				stat_map_wm->nodes[i] = ERRVAL;
			}
		}
		// calculate the stats
		stat_map_wm->CalcMinMaxAvg();

		// if there is a transformation matrix then untransform
		if (transmat)
		{
			stat_map_wm2 = new RicTexture(thick_map_wm2->size);

			for (int i = 0; i < thick_map_wm2->size; ++i)
			{
				if (thick_map_wm2->nodes[i] > mind && thick_map_wm2->nodes[i] < maxd)
				{
					stat_map_wm2->nodes[i] = thick_map_wm2->nodes[i];
				}
				else
				{
					stat_map_wm2->nodes[i] = ERRVAL;
				}
			}

			// calculate the stats
			stat_map_wm2->CalcMinMaxAvg();
		}
	}

	if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // stat map on gray
	{
		stat_map_gm = new RicTexture(thick_map_gm->size);

		// copy values from thickness map flagging all nodes outside range
		// also zero corresponding vectors in connecting mesh
		for (int i = 0; i < thick_map_gm->size; ++i)
		{
			if (thick_map_gm->nodes[i] > mind && thick_map_gm->nodes[i] < maxd)
			{
				stat_map_gm->nodes[i] = thick_map_gm->nodes[i];
			}
			else
			{
				stat_map_gm->nodes[i] = ERRVAL;
			}
		}

		// calculate the stats
		stat_map_gm->CalcMinMaxAvg();

		// if there is a transformation matrix then untransform
		if (transmat)
		{
			stat_map_gm2 = new RicTexture(thick_map_gm2->size);

			for (int i = 0; i < thick_map_gm2->size; ++i)
			{
				if (thick_map_gm2->nodes[i] > mind && thick_map_gm2->nodes[i] < maxd)
				{
					stat_map_gm2->nodes[i] = thick_map_gm2->nodes[i];
				}
				else
				{
					stat_map_gm2->nodes[i] = ERRVAL;
				}
			}

			// calculate the stats
			stat_map_gm2->CalcMinMaxAvg();
		}
	}

	// now do special stuff for the PeterK dude
	// combine the average and standard deviation for by directions
	double comb_avg=0, comb_std_dev=0;
	double comb_avg2=0, comb_std_dev2=0;
	if ( peterk && (surfmap==MAP_BOTH) )
	{
		double std_t=0, avg_t=0;
		int ccnt=0;	// count of number of values in average
		for (int i=0; i<stat_map_wm->size; ++i)
		{
			if (stat_map_wm->nodes[i] == ERRVAL) // skip invalid values
			{
				continue; // skip empty ones
			}
			avg_t += stat_map_wm->nodes[i];
			std_t += stat_map_wm->nodes[i]*stat_map_wm->nodes[i];

			++ccnt;
		}

		for (int i=0; i<stat_map_gm->size; ++i)
		{
			if (stat_map_gm->nodes[i] == ERRVAL) // skip invalid values
			{
				continue; // skip empty ones
			}
			avg_t += stat_map_gm->nodes[i];
			std_t += stat_map_gm->nodes[i]*stat_map_gm->nodes[i];

			++ccnt;
		}

		// average
		comb_avg = avg_t/(double)ccnt;

		// standard deviation
		comb_std_dev= sqrt(fabs(std_t-avg_t*avg_t/(double)ccnt)/((double)ccnt-1.0));

		// if there is an transform matrix then do the average again
		if ( transmat )
		{
			std_t = avg_t = 0;
			ccnt = 0;
			for (int i=0; i<stat_map_wm2->size; ++i)
			{
				if (stat_map_wm2->nodes[i] == ERRVAL) // skip invalid values
				{
					continue; // skip empty ones
				}
				avg_t += stat_map_wm2->nodes[i];
				std_t += stat_map_wm2->nodes[i]*stat_map_wm2->nodes[i];

				++ccnt;
			}

			for (int i=0; i<stat_map_gm2->size; ++i)
			{
				if (stat_map_gm2->nodes[i] == ERRVAL) // skip invalid values
				{
					continue; // skip empty ones
				}
				avg_t += stat_map_gm2->nodes[i];
				std_t += stat_map_gm2->nodes[i]*stat_map_gm2->nodes[i];

				++ccnt;
			}

			// average
			comb_avg2 = avg_t/(double)ccnt;

			// standard deviation
			comb_std_dev2= sqrt(fabs(std_t-avg_t*avg_t/(double)ccnt)/((double)ccnt-1.0));

		}
	}
//////////////////////////////////////////////////////////////////////////////
//// Step 5 - Calculate texture histogram

	// figure out histogram bin size
	int nbins = 20;
	float binsize = 10.0/(float)nbins; // set max to 10mm

	int *wm_bins = NULL, *wm_bins2 = NULL;
	int *gm_bins = NULL, *gm_bins2 = NULL;
	int num_wm_culled = 0, num_wm_culled2 = 0;
	int num_gm_culled = 0, num_gm_culled2 = 0;

	if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // stat map on white
	{
		wm_bins = new int[nbins];
		for (int i = 0; i < nbins; ++i)
			wm_bins[i] = 0; // zero bins

		// sort the thickness values into appropriate bins
		for (int i = 0; i < stat_map_wm->size; ++i)
		{
			if (stat_map_wm->nodes[i] != ERRVAL)
			{
				int j = (int) (stat_map_wm->nodes[i] / binsize);
				if (j > (nbins - 1)) // make sure it fits in bins
					j = nbins - 1;
				++wm_bins[j];
			}
			else // keep track of vertices not used
			{
				++num_wm_culled;
			}
		}

		// now do the untransformed values
		if (transmat)
		{
			wm_bins2 = new int[nbins];
			for (int i = 0; i < nbins; ++i)
				wm_bins2[i] = 0; // zero bins

			// sort the thickness values into appropriate bins for transformed data
			for (int i = 0; i < stat_map_wm2->size; ++i)
			{
				if (stat_map_wm2->nodes[i] != ERRVAL)
				{
					int j = (int) (stat_map_wm2->nodes[i] / binsize);
					if (j > (nbins - 1)) // make sure it fits in bins
						j = nbins - 1;
					++wm_bins2[j];
				}
				else // keep track of vertices not used
				{
					++num_wm_culled2;
				}
			}
		}
	}

	if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // stat map on gray
	{
		gm_bins = new int[nbins];
		for (int i = 0; i < nbins; ++i)
			gm_bins[i] = 0; // zero bins

		// sort the thickness values into appropriate bins
		for (int i = 0; i < stat_map_gm->size; ++i)
		{
			if (stat_map_gm->nodes[i] != ERRVAL)
			{
				int j = (int) (stat_map_gm->nodes[i] / binsize);
				if (j > (nbins - 1)) // make sure it fits in bins
					j = nbins - 1;
				++gm_bins[j];
			}
			else // keep track of vertices not used
			{
				++num_gm_culled;
			}
		}

		// now do the untransformed values
		if (transmat)
		{
			gm_bins2 = new int[nbins];
			for (int i = 0; i < nbins; ++i)
				gm_bins2[i] = 0; // zero bins

			// sort the thickness values into appropriate bins for transformed data
			for (int i = 0; i < stat_map_gm2->size; ++i)
			{
				if (stat_map_gm2->nodes[i] != ERRVAL)
				{
					int j = (int) (stat_map_gm2->nodes[i] / binsize);
					if (j > (nbins - 1)) // make sure it fits in bins
						j = nbins - 1;
					++gm_bins2[j];
				}
				else // keep track of vertices not used
				{
					++num_gm_culled2;
				}
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////
/// Step 6 - Fill Texture

	/// If Texfill is set then fill the holes in the texture map, first by
	/// interpolating. Two shots are given to interpolation, first fill the
	/// smallest holes then those up to filld distance. Any holes left over
	/// are set to zero.

	// Thickness map
	if ( filltex )
	{
		if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // white texture
		{
			// we need to remap the error values to something more appropriate.
			// first try to fill by looking at surrounding nodes
			// if that does not work then just map those values to zero
			int status=0;
			status = TexFillAvg(thick_map_wm, &wm_mesh, filld*0.5);
			status = TexFillAvg(thick_map_wm, &wm_mesh, filld);
			if ( status == 0 ) // map the bogus ones to zero
			{
				for ( int i=0 ; i<thick_map_wm->size ; ++i )
				{
					if ( thick_map_wm->nodes[i] == ERRVAL )
						thick_map_wm->nodes[i]=0;
				}
			}


			if ( transmat ) // then do the untransformed texture also
			{
				status = 0;
				if ( filltex )
				{
					status = TexFillAvg(thick_map_wm2, &wm_mesh2, filld*0.5);
					status = TexFillAvg(thick_map_wm2, &wm_mesh2, filld);
				}
				if ( status == 0 ) // map the bogus ones to zero
				{
					for ( int i=0 ; i<thick_map_wm2->size ; ++i )
					{
						if ( thick_map_wm2->nodes[i] == ERRVAL )
							thick_map_wm2->nodes[i]=0;
					}
				}
			}
		}

		if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // gray texture
		{
			// we need to remap the error values to something more appropriate.
			// first try to fill by looking at surrounding nodes
			// if that does not work then just map those values to zero
			int status=0;
			status = TexFillAvg(thick_map_gm, &gm_mesh, filld*0.5);
			status = TexFillAvg(thick_map_gm, &gm_mesh, filld);
			if ( status == 0 ) // map the bogus ones to zero
			{
				for ( int i=0 ; i<thick_map_gm->size ; ++i )
				{
					if ( thick_map_gm->nodes[i] == ERRVAL )
						thick_map_gm->nodes[i]=0;
				}
			}


			if ( transmat ) // then do the untransformed texture also
			{
				status = 0;
				if ( filltex )
				{
					status = TexFillAvg(thick_map_gm2, &gm_mesh2, filld*0.5);
					status = TexFillAvg(thick_map_gm2, &gm_mesh2, filld);
				}
				if ( status == 0 ) // map the bogus ones to zero
				{
					for ( int i=0 ; i<thick_map_gm2->size ; ++i )
					{
						if ( thick_map_gm2->nodes[i] == ERRVAL )
							thick_map_gm2->nodes[i]=0;
					}
				}
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////
/// Step 7 - Output Results

	/// Basic stats output to console - average, min, max, std dev

	if (verbose || brainvisa)
	{
		if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // stat map on white
		{
			if (transmat)
			{
				cout << endl
						<< "Cortical thickness calculations (white-gray) untransformed"
						<< endl;
				cout << "Min Cortical thickness, " << stat_map_wm->min << ", "
						<< stat_map_wm2->min << endl;
				cout << "Max Cortical thickness, " << stat_map_wm->max << ", "
						<< stat_map_wm2->max << endl;
				cout << "Cortical Std deviation, " << stat_map_wm->std_dev << ", "
						<< stat_map_wm2->std_dev << endl;
				cout << "Median Cortical thickness, " << stat_map_wm->med << ", "
						<< stat_map_wm2->med << endl;
				cout << "Average Cortical thickness, " << stat_map_wm->avg << ", "
						<< stat_map_wm2->avg << endl;
				cout << "Number of mesh vertices, " << stat_map_wm->size
						<< ", " << stat_map_wm2->size << endl;
				cout << "Number of WM vertices distance culled, " << num_wm_culled << ", "
						<< num_wm_culled2 << endl;
				cout << "Number of WM vertices curvature culled, "
						<< nwm_culled << endl;
			}
			else
			{
				cout << endl << "Cortical thickness calculations (white-gray)" << endl;
				cout << "Min Cortical thickness, " << stat_map_wm->min << endl;
				cout << "Max Cortical thickness, " << stat_map_wm->max << endl;
				cout << "Cortical Std deviation, " << stat_map_wm->std_dev << endl;
				cout << "Median Cortical thickness, " << stat_map_wm->med << endl;
				cout << "Average Cortical thickness, " << stat_map_wm->avg << endl;
				cout << "Number of mesh vertices, " << stat_map_wm->size
						<< endl;
				cout << "Number of WM vertices distance culled, " << num_wm_culled << endl;
				cout << "Number of WM vertices curvature culled, "
						<< nwm_culled << endl;
			}
		}

		if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // stat map on gray
		{
			if (transmat)
			{
				cout << endl
						<< "Cortical thickness calculations (gray-white) - untransformed"
						<< endl;
				cout << "Min Cortical thickness, " << stat_map_gm->min << ", "
						<< stat_map_gm2->min << endl;
				cout << "Max Cortical thickness, " << stat_map_gm->max << ", "
						<< stat_map_gm2->max << endl;
				cout << "Cortical Std deviation, " << stat_map_gm->std_dev << ", "
						<< stat_map_gm2->std_dev << endl;
				cout << "Median Cortical thickness, " << stat_map_gm->med << ", "
						<< stat_map_gm2->med << endl;
				cout << "Average Cortical thickness, " << stat_map_gm->avg << ", "
						<< stat_map_gm2->avg << endl;
				cout << "Number of mesh vertices, " << stat_map_gm->size
						<< ", " << stat_map_gm2->size << endl;
				cout << "Number of GM vertices distance culled, " << num_gm_culled << ", "
						<< num_gm_culled2 << endl;
				cout << "Number of GM vertices curvature culled, "
						<< ngm_culled << endl;
			}
			else
			{
				cout << endl << "Cortical thickness calculations (gray-white)" << endl;
				cout << "Min Cortical thickness, " << stat_map_gm->min << endl;
				cout << "Max Cortical thickness, " << stat_map_gm->max << endl;
				cout << "Cortical Std deviation, " << stat_map_gm->std_dev << endl;
				cout << "Median Cortical thickness, " << stat_map_gm->med << endl;
				cout << "Average Cortical thickness, " << stat_map_gm->avg << endl;
				cout << "Number of mesh vertices, " << stat_map_gm->size
						<< endl;
				cout << "Number of GM vertices distance culled, " << num_gm_culled << endl;
				cout << "Number of GM vertices curvature culled, "
						<< ngm_culled << endl;
			}
		}
	}
	else // not verbose or BrainVisa
	{
		if ( peterk && (surfmap==MAP_BOTH) )
		{
			// output is name, average, transformed average, std dev, transformed std dev
			if ( transmat )
			{
				cout << outname << "\t" << comb_avg << "\t" << comb_avg2 << "\t"
					<< comb_std_dev << "\t" << comb_std_dev2 << "\t" << endl;
			}
			else
			{
				cout << outname << "\t" << comb_avg << "\t"
					<< comb_std_dev << "\t" << endl;
			}
		}
		else // output values individually for white to gray and gray to white
		{
			if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // stat map on white
			{
				// output is name, average, transformed average, std dev, transformed std dev
				if ( transmat )
				{
					cout << outname << "_WM\t" << stat_map_wm->avg << "\t" << stat_map_wm2->avg << "\t"
						<< stat_map_wm->std_dev << "\t" << stat_map_wm2->std_dev << "\t" << endl;
				}
				else
				{
					cout << outname << "_WM\t" << stat_map_wm->avg << "\t"
						<< stat_map_wm->std_dev << "\t" << endl;
				}
			}

			if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // stat map on gray
			{
				// output is name, average, transformed average, std dev, transformed std dev
				if ( transmat )
				{
					cout << outname << "_GM\t" << stat_map_gm->avg << "\t" << stat_map_gm2->avg << "\t"
						<< stat_map_gm->std_dev << "\t" << stat_map_gm2->std_dev << "\t" << endl;
				}
				else
				{
					cout << outname << "_GM\t" << stat_map_gm->avg << "\t"
						<< stat_map_gm->std_dev << "\t" << endl;
				}
			}
		}
	}


	/// Stats output to files

	string filename;
	ofstream ofile;

	/// First write simple one line file to be read by Excel
	if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // stat map on white
	{
		// first write simple one line file to be read by Excel
		filename = outname + "_wm.tab";
		ofile.open (filename.c_str());

		// output is name, average, transformed average, std dev, transformed std dev
		if ( transmat )
		{
			ofile << outname << "\t" << stat_map_wm->avg << "\t" << stat_map_wm2->avg << "\t"
				<< stat_map_wm->std_dev << "\t" << stat_map_wm2->std_dev << "\t" << endl;
		}
		else
		{
			ofile << outname << "\t" << stat_map_wm->avg << "\t"
				<< stat_map_wm->std_dev << "\t" << endl;
		}
		ofile.close();
	}

	if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // stat map on gray
	{
		// first write simple one line file to be read by Excel
		filename = outname + "_gm.tab";
		ofstream ofile;
		ofile.open (filename.c_str());

		// output is name, average, transformed average, std dev, transformed std dev
		if ( transmat )
		{
			ofile << outname << "\t" << stat_map_gm->avg << "\t" << stat_map_gm2->avg << "\t"
				<< stat_map_gm->std_dev << "\t" << stat_map_gm2->std_dev << "\t" << endl;
		}
		else
		{
			ofile << outname << "\t" << stat_map_gm->avg << "\t"
				<< stat_map_gm->std_dev << "\t" << endl;
		}
		ofile.close();
	}

	/// Now write a more complete file with all the stuff (CSV)

	if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // stat map on white
	{
		filename = outname + "_comp_wm.csv";
		ofstream ofile;
		ofile.open (filename.c_str());

		if (transmat)
		{
			ofile << endl
					<< "Cortical thickness calculations (white-gray) - untransformed"
					<< endl;
			ofile << "Min Cortical thickness, " << stat_map_wm->min << ", "
					<< stat_map_wm2->min << endl;
			ofile << "Max Cortical thickness, " << stat_map_wm->max << ", "
					<< stat_map_wm2->max << endl;
			ofile << "Cortical Std deviation, " << stat_map_wm->std_dev << ", "
					<< stat_map_wm2->std_dev << endl;
			ofile << "Median Cortical thickness, " << stat_map_wm->med << ", "
					<< stat_map_wm2->med << endl;
			ofile << "Average Cortical thickness, " << stat_map_wm->avg << ", "
					<< stat_map_wm2->avg << endl;
			ofile << "Number of mesh vertices, " << stat_map_wm->size
					<< ", " << stat_map_wm2->size << endl;
			ofile << "Number of WM vertices distance culled, " << num_wm_culled << ", "
					<< num_wm_culled2 << endl;
			ofile << "Number of WM vertices curvature culled, "
					<< nwm_culled << endl;
			for ( int i=0 ; i<nbins ; ++i )
				ofile << (i*binsize) << "," << wm_bins2[i] << endl;
		}
		else
		{
			ofile << endl << "Cortical thickness calculations (white-gray)" << endl;
			ofile << "Min Cortical thickness, " << stat_map_wm->min << endl;
			ofile << "Max Cortical thickness, " << stat_map_wm->max << endl;
			ofile << "Cortical Std deviation, " << stat_map_wm->std_dev << endl;
			ofile << "Median Cortical thickness, " << stat_map_wm->med << endl;
			ofile << "Average Cortical thickness, " << stat_map_wm->avg << endl;
			ofile << "Number of mesh vertices, " << stat_map_wm->size
					<< endl;
			ofile << "Number of WM vertices distance culled, " << num_wm_culled << endl;
			ofile << "Number of WM vertices curvature culled, "
					<< nwm_culled << endl;
			for ( int i=0 ; i<nbins ; ++i )
				ofile << (i*binsize) << "," << wm_bins[i] << endl;
		}
		ofile.close();
	}

	if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // stat map on gray
	{
		filename = outname + "_comp_gm.csv";
		ofstream ofile;
		ofile.open (filename.c_str());

		if (transmat)
		{
			ofile << endl
					<< "Cortical thickness calculations (gray-white) - untransformed"
					<< endl;
			ofile << "Min Cortical thickness, " << stat_map_gm->min << ", "
					<< stat_map_gm2->min << endl;
			ofile << "Max Cortical thickness, " << stat_map_gm->max << ", "
					<< stat_map_gm2->max << endl;
			ofile << "Cortical Std deviation, " << stat_map_gm->std_dev << ", "
					<< stat_map_gm2->std_dev << endl;
			ofile << "Median Cortical thickness, " << stat_map_gm->med << ", "
					<< stat_map_gm2->med << endl;
			ofile << "Average Cortical thickness, " << stat_map_gm->avg << ", "
					<< stat_map_gm2->avg << endl;
			ofile << "Number of mesh vertices, " << stat_map_gm->size
					<< ", " << stat_map_gm2->size << endl;
			ofile << "Number of GM vertices distance culled, " << num_gm_culled << ", "
					<< num_gm_culled2 << endl;
			ofile << "Number of GM vertices curvature culled, "
					<< ngm_culled << endl;
			for ( int i=0 ; i<nbins ; ++i )
				ofile << (i*binsize) << "," << gm_bins2[i] << endl;
		}
		else
		{
			ofile << endl << "Cortical thickness calculations (gray-white)" << endl;
			ofile << "Min Cortical thickness, " << stat_map_gm->min << endl;
			ofile << "Max Cortical thickness, " << stat_map_gm->max << endl;
			ofile << "Cortical Std deviation, " << stat_map_gm->std_dev << endl;
			ofile << "Median Cortical thickness, " << stat_map_gm->med << endl;
			ofile << "Average Cortical thickness, " << stat_map_gm->avg << endl;
			ofile << "Number of mesh vertices, " << stat_map_gm->size
					<< endl;
			ofile << "Number of GM vertices distance culled, " << num_gm_culled << endl;
			ofile << "Number of GM vertices curvature culled, "
					<< ngm_culled << endl;
			for ( int i=0 ; i<nbins ; ++i )
				ofile << (i*binsize) << "," << gm_bins[i] << endl;
		}
		ofile.close();
	}

	/// Write textures to files

	if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // white texture map
	{
		if ( transmat )
		{
			// write the untransformed thickness map to a file
			filename = outname + "_wm_utrans_thick.tex";
			thick_map_wm2->write_texture((char*)filename.c_str());
		}

		// write the thickness map to a file
		filename = outname + "_wm_thick.tex";
		thick_map_wm->write_texture((char*)filename.c_str());
	}

	if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // gray texture map
	{
		if ( transmat )
		{
			// write the untransformed thickness map to a file
			filename = outname + "_gm_utrans_thick.tex";
			thick_map_gm2->write_texture((char*)filename.c_str());
		}

		// write the thickness map to a file
		filename = outname + "_gm_thick.tex";
		thick_map_gm->write_texture((char*)filename.c_str());
	}

	/// Write out diagnostic meshes

	if ( verbose )
	{
		if (surfmap == MAP_WHITE || surfmap == MAP_BOTH) // white diagnostic
		{
			// write out the vector mesh
			filename = outname + "_wm_vector.mesh";
			cv_mesh_wm.write_mesh_bin((char*)filename.c_str());

			if ( transmat )
			{
				// write out the vector mesh
				filename = outname + "_wm_vector2.mesh";
				cv_mesh_wm2.write_mesh_bin((char*)filename.c_str());

			}
		}

		if (surfmap == MAP_GRAY || surfmap == MAP_BOTH) // gray diagnostice
		{
			// write out the vector mesh
			filename = outname + "_gm_vector.mesh";
			cv_mesh_gm.write_mesh_bin((char*)filename.c_str());

			if ( transmat )
			{
				// write out the vector mesh
				filename = outname + "_gm_vector2.mesh";
				cv_mesh_gm2.write_mesh_bin((char*)filename.c_str());

			}
		}

		// for reference purposes, write out the untransformed meshes
		if ( transmat )
		{
			// write out the untransformed wm mesh
			filename = outname + "_utrans_wm.mesh";
			wm_mesh2.write_mesh_bin((char*)filename.c_str());

			// write out the untransformed gm mesh
			filename = outname + "_utrans_gm.mesh";
			gm_mesh2.write_mesh_bin((char*)filename.c_str());
		}
	}

	// clean up memory
	if ( thick_map_wm ) delete thick_map_wm;
	if ( thick_map_wm2 ) delete thick_map_wm2;
	if ( thick_map_gm ) delete thick_map_gm;
	if ( thick_map_gm2 ) delete thick_map_gm2;

	if ( stat_map_wm ) delete stat_map_wm;
	if ( stat_map_wm2 ) delete stat_map_wm2;
	if ( stat_map_gm) delete stat_map_gm;
	if ( stat_map_gm2 ) delete stat_map_gm2;

	if ( gm_bins ) delete gm_bins;
	if ( gm_bins2 ) delete gm_bins2;
	if ( wm_bins ) delete wm_bins;
	if ( wm_bins2 ) delete wm_bins2;

	// lets see how long this took us to do
	time(&endtime);
	double elapsed = difftime(endtime,starttime);
	if ( verbose ) cout << "Elapsed time: "<<(elapsed/60.0)<<endl;
	if ( verbose ) cout << "That's all folks"<<endl;

	return EXIT_SUCCESS;
}
