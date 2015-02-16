// --------------------------- RicGyralSpan.cpp ------------------------
/*! \file
Implementation file for Peter K's RicGyralSpan scheme.
Copyright (C) 2009 by Bill Rogers - Research Imaging Center - UTHSCSA

\mainpage
This program calculates the gyral span by gyrus. It uses a gyral skeletal
mesh for determination of the center of the gyrus. There are several
required inputs.

1) Gyral graph file

2) Name of gyrus to measure

3) Gyral skeletal mesh

The gyral name is used to find the name of the gyral mesh file listed in
the gyral graph file.

The gyral thickness is calculated from normals from the gyral skeleton mesh
to the gyral mesh.


Gyrus thickness values are output to stdout.

Command line switches

-g		gyral arg file

-n		gyral name

--sm	skeletal mesh

--tm	Transformation matrix file name

--sn	subject name

--mt	max gyral thickness allowed

--dp	max dot product for gyral normals

-o		output file base name

-v		verbose output

--bv	output for BrainVisa

 */


using namespace std;
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include <RicGraph.h>
#include <RicMesh.h>
#include <RicUtil.h>

bool FileExists(string strFilename);
bool IntersectMesh(Point p0, Point p1, int pnump, int pnumn, RicMesh *mesh);

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string gyral_name = "";		// name of gyrus to measure
	string gyral_arg_name = "";	// name of arg file for gyri
	string skel_name = "";		// name of  gyral skeleton mesh
	string matname = "";		// name for transformation matrix file
	string subj_name = "";		// subject name
	string outname = "";		// base name for output files
	float maxthick=0;			// max gyral thickness allowed
	float maxdot=0;				// max dot product for comparing normals
	bool verbose=false;			// if true output more stuff on stdout
	bool brainvisa=false;		// if true output stuff for BrainVisa
	string filename = "";

	///////// read in the command line - tclap ///////////////

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("Peter K's Gyral Thickness Measurement Program", ' ', "0.9");

		// gyrus name
		TCLAP::ValueArg<string> gname( "n","n", "Gyrus name", true, "", "string");
		cmd.add(gname);

		// gyral arg file
		TCLAP::ValueArg<string> garg( "g","g", "Gyral arg file", true, "", "string");
		cmd.add(garg);

		// transformation matrix name
		TCLAP::ValueArg<string> mname( "","tm", "Transformation matrix file name", false,\
			 "", "string");
		cmd.add(mname);

		// subject name
		TCLAP::ValueArg<string> sname( "s","sn", "Subject name", false, "", "string");
		cmd.add(sname);

		// gyral skeleton mesh name
		TCLAP::ValueArg<string> skname( "","sm", "Gyral skeleton mesh name", true, "", "string");
		cmd.add(skname);

		// max thickness
		TCLAP::ValueArg<float> mxthk( "","mt", "Max gyral thickness", false, 20.0, "float");
		cmd.add(mxthk);

		// max dot product
		TCLAP::ValueArg<float> mxdot( "","md", "Max dot product for comparing normals", false, 0.75, "float");
		cmd.add(mxdot);

		// output file base name
		TCLAP::ValueArg<string> oname( "o","o", "Base output file name", false, \
										"", "string");
		cmd.add(oname);


		TCLAP::SwitchArg bv( "","bv", "BrainVisa output", false);
		cmd.add(bv);

		TCLAP::SwitchArg ver( "v","verbose", "Verbose output on stdout", false);
		cmd.add(ver);

		// parse the command line
		cmd.parse( argc, argv );

		// copy command line variable to program variables
		gyral_name = gname.getValue();
		gyral_arg_name = garg.getValue();
		skel_name = skname.getValue();
		matname = mname.getValue();
		maxthick = mxthk.getValue();
		maxdot = mxdot.getValue();
		outname = oname.getValue();
		subj_name = sname.getValue();
		verbose = ver.getValue();
		brainvisa = bv.getValue();

	}
	catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cout << " Command Line Error " << e.error() << " for arg " << e.argId() << endl;
	}

///////////////////////////////////////////////////////////////////////////////
/// Step 1 - Read in files - check to see that specified gyrus is in the graph
///	file. The graph file will list a mesh file name for the gyrus. If that mesh
/// file is not in the current directory, it will check to see if it is a
/// subdirectory with the same name as the graph file.
///	Also do sanity checks to make sure we have valid data to work with.

	// read in the graph file
	RicGraph graph((char*)gyral_arg_name.c_str());
	if ( graph.nnodes < 1 )
	{
		cout << "Invalid graph file " << gyral_arg_name << endl;
		exit(1);
	}

	// check to see that the gyrus specified is in graph file
	int n;
	if ( (n=graph.FindNodeNameExact(gyral_name)) >= 0 )
	{
		if ( verbose )
			cout << "Found " << gyral_name << " at index " << n << endl;
	}
	else
	{
		cout << "Did not find " << gyral_name << " in " << gyral_arg_name << endl;
		exit(1);
	}

	// read in the mesh specified by the gyral name
	string gfile = graph.nodes[n].Tmtktri_filename + ".mesh";
	if ( verbose) cout << gfile << endl;

	// If it is not in the current directory then try a directory made from
	// the graph file name
	if ( FileExists(gfile) == false )
	{
		// add base part of graph file name as subdirectory name
		string gbase = gyral_arg_name.substr(0,gyral_arg_name.length()-4) + ".data/";
		gfile = gbase + gfile;
		if ( FileExists(gfile) == false )
		{
			cout << "Cannot find file " << gfile << endl;
			exit(1);
		}
	}

	RicMesh gyral_mesh((char*)gfile.c_str());
	if ( gyral_mesh.v_size == 0 )
	{
		cout << "Invalid gyral mesh file " << gfile << endl;
		exit(1);
	}
	if ( gyral_mesh.p_dim != 3 )
	{
		cout << "Gyral mesh file must only contain triangles " << gfile << endl;
		exit(1);
	}

	// read in gyral skeleton mesh
	RicMesh skel_mesh((char*)skel_name.c_str());
	if ( skel_mesh.v_size == 0 )
	{
		cout << "Invalid gyral skeleton mesh file " << skel_name << endl;
		exit(1);
	}
	if ( skel_mesh.p_dim != 3 )
	{
		cout << "Gyral skeleton mesh file must only contain triangles " << skel_name << endl;
		exit(1);
	}


	// Read the Transformation Matrix  if there is one
	CMatrix tm;		// transformation matrix
	bool transmat=false;		// if true then use transformation matrix
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


///////////////////////////////////////////////////////////////////////////////
/// Step 2 - Find appropriate vertices in gyral skeleton by only using those
/// that fit in the bounding box of the gyrus. Make two arrays, one for the
/// vertices and the other for the matching normals.

	// find bounding box of gyrus
	gyral_mesh.calc_limits();
	float xmin = gyral_mesh.xmin;
	float ymin = gyral_mesh.ymin;
	float xmax = gyral_mesh.xmax;
	float ymax = gyral_mesh.ymax;
	float zmin = gyral_mesh.zmin;
	float zmax = gyral_mesh.zmax;

	// find vertices in skeleton mesh that are in bounding box
	vertex *slist = new vertex[skel_mesh.v_size]; // list of possible vertices
	vertex *nlist = new vertex[skel_mesh.v_size]; // list of possible vertices
	int nsvert = 0;; // number of skeleton vertices in bounding box
	for ( int i=0 ; i<skel_mesh.v_size ; ++i)
	{
		Point p = skel_mesh.vertices[i].pnt;
		if ( p.x > xmin && p.x < xmax && p.y > ymin && p.y < ymax
				&& p.z > zmin && p.z < zmax )
		{
			nlist[nsvert] = skel_mesh.normals[i];
			slist[nsvert++] = skel_mesh.vertices[i];  // we found one
		}
	}

///////////////////////////////////////////////////////////////////////////////
/// Step 3 - Use normals from skeleton mesh polygons found in bounding box to
/// intersect opposing polygons in gyral mesh. A normal vector from the
///	skeleton mesh is extended in both directions to check for intersections
///	with gyral mesh polygons. A dot product is used to compare
///	the normal from a skeleton mesh points to the average normal of an
/// intersecting mesh polygon. If the normals are not aligned then the point is
/// not used. Lists are made of the the end points where the vector from the
///	skeleton mesh intersects the gyral mesh.

	float *tval = new float[nsvert];
	Point *pntsp = new Point[nsvert];
	Point *pntsn = new Point[nsvert];
	Point *norsp = new Point[nsvert];
	Point *norsn = new Point[nsvert];
	float maxthick2 = maxthick*maxthick;

	// populate the thickness map with ERRVAL as default value
	for ( int i=0 ; i<nsvert ; ++i ) tval[i] = ERRVAL;

	// work from skeleton to gyral mesh
	int nfound = 0;
	int cnt=0;
	int pnump=0,pnumn=0; // polygon numbers for intersected polygons in both directions
	for ( int i=0 ; i<nsvert ; ++i )
	{
		// skip if vertex labeled as not to use
		if ( slist[i].label == 1 )
			continue;

		Point p0,pp,pn;	// line normal to vertex
		Vector n0;	// normals from vertex and opposing polygon

		p0 = slist[i].pnt;
		n0 = nlist[i].pnt;

		// make a long line in the direction of the normal
		pp.x = p0.x + maxthick*n0.x; // add normal to vertex to get line
		pp.y = p0.y + maxthick*n0.y;
		pp.z = p0.z + maxthick*n0.z;

		// make a long line in the opposite direction of the normal
		pn.x = p0.x - maxthick*n0.x; // subtract normal to vertex to get line
		pn.y = p0.y - maxthick*n0.y;
		pn.z = p0.z - maxthick*n0.z;

		// test against all polygons in the gyral mesh
		// to find the closest inner polygon
		Point t0,t1,t2;
		Vector tn0,tn1,tn2,tna; // normals
		float d0,d1,d2;

		// find distance in positive direction
		float indist=ERRVAL;
		bool foundp = false;
		Point inpntp;
		Point innorp;
		for ( int j=0 ; j<gyral_mesh.p_size ; ++j )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			t0 = gyral_mesh.vertices[gyral_mesh.polygons[j].vidx[0]].pnt;
			t1 = gyral_mesh.vertices[gyral_mesh.polygons[j].vidx[1]].pnt;
			t2 = gyral_mesh.vertices[gyral_mesh.polygons[j].vidx[2]].pnt;

			tn0 = gyral_mesh.normals[gyral_mesh.polygons[j].vidx[0]].pnt;
			tn1 = gyral_mesh.normals[gyral_mesh.polygons[j].vidx[1]].pnt;
			tn2 = gyral_mesh.normals[gyral_mesh.polygons[j].vidx[2]].pnt;

			// get average normals
			tna = tn0 + tn1 + tn2;
			tna.Normalize();

			// make sure that our average normal is pointing the in the opposite
			// direction of our test point
			float dprod = tna.Dot(n0);
			if ( dprod < maxdot )
				continue;

			// skip if any of the triangle points are our current vertex
			if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
			if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
			if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>maxthick2 && d1>maxthick2 && d2>maxthick2 ) continue;

			// if no intersection of line with polygon then skip
			Point pint;
			if ( !line_intersect_triangle(t0,t1,t2,p0,pp,&pint) )
				continue;

			// if shorter than previous measurement then give it a try.
			float dint = dist(p0,pint);
			if ( dint < indist ) // we found a point
			{
				indist = dint;
				inpntp = pint;
				innorp = tna;
				pnump = j;
				foundp = true;
				++cnt;
			}

		} // end of loop for positive normals

		// find distance in positive direction
		if ( !foundp )
			continue;

		indist=ERRVAL;
		bool foundn = false;
		Point inpntn;
		Point innorn;
		for ( int j=0 ; j<gyral_mesh.p_size ; ++j )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			t0 = gyral_mesh.vertices[gyral_mesh.polygons[j].vidx[0]].pnt;
			t1 = gyral_mesh.vertices[gyral_mesh.polygons[j].vidx[1]].pnt;
			t2 = gyral_mesh.vertices[gyral_mesh.polygons[j].vidx[2]].pnt;

			tn0 = gyral_mesh.normals[gyral_mesh.polygons[j].vidx[0]].pnt;
			tn1 = gyral_mesh.normals[gyral_mesh.polygons[j].vidx[1]].pnt;
			tn2 = gyral_mesh.normals[gyral_mesh.polygons[j].vidx[2]].pnt;

			// get average normals
			tna = tn0 + tn1 + tn2;
			tna.Normalize();

			// make sure that our average normal is pointing the in the opposite
			// direction of our test point
			float dprod = tna.Dot(n0);
			if ( dprod > -maxdot )
				continue;

			// skip if any of the triangle points are our current vertex
			if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
			if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
			if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>maxthick2 && d1>maxthick2 && d2>maxthick2 ) continue;

			// if no intersection of line with polygon then skip
			Point pint;
			if ( !line_intersect_triangle(t0,t1,t2,p0,pn,&pint) )
				continue;

			// if shorter than previous measurement then give it a try.
			float dint = dist(p0,pint);
			if ( dint < indist ) // we found a point
			{
				indist = dint;
				inpntn = pint;
				innorn = tna;
				pnumn = j;
				foundn = true;
				++cnt;
			}

		} // end of loop for negative normals

		if ( foundp && foundn )
		{
			// assign distance for vertex
			float dd = dist(inpntp,inpntn);

			if (dd < maxthick)
			{
				// for a check, make sure that no polygons are in between
				if ( IntersectMesh(inpntp, inpntn, pnump, pnumn, &gyral_mesh) == false )
				{
					tval[nfound] = dd;
					pntsp[nfound] = inpntp;
					pntsn[nfound] = inpntn;
					norsp[nfound] = innorp;
					norsn[nfound] = innorn;
					++nfound;
				}
			}
		}

	}

///////////////////////////////////////////////////////////////////////////////
/// Step 4 - Make a mesh showing connecting lines between points spanning the
/// gyral mesh. This is at 2D mesh that displays as lines spanning the interior
/// of the gyrus.

	int svpoly_size = nfound;
	RicMesh sv_mesh(2*svpoly_size, 2*svpoly_size, svpoly_size, 2);

	// populate mesh with short vectors using the plus and minus points
	double avg_t = 0, std_t = 0;
	for ( int i=0 ; i<svpoly_size ; ++i )
	{
		sv_mesh.assign_node(2*i, pntsp[i]);
		sv_mesh.assign_normal(2*i, norsp[i]);
		sv_mesh.assign_node(2*i+1, pntsn[i]);
		sv_mesh.assign_normal(2*i+1, norsn[i]);
		sv_mesh.assign_polygon(i, 2*i, 2*i+1, 0);
		avg_t += tval[i];
		std_t += tval[i]*tval[i];
	}

	// if there is an output name then use it as the base file name
	if ( outname.length() > 0 )
		filename = outname + "_" +gyral_name + "_gyvec.mesh";
	else // just use the gyral name
		filename = gyral_name + "_gyvec.mesh";
	sv_mesh.Write((char*)filename.c_str());

///////////////////////////////////////////////////////////////////////////////
/// Step 5 - Calculate average and standard deviation of gyral thickness then
/// sent to stdout. Also send output to text file

	double avg, std_dev;

	avg = avg_t/svpoly_size;

	std_dev= sqrt(fabs(std_t-avg_t*avg_t/(float)svpoly_size)/((float)svpoly_size-1.0));

///////////////////////////////////////////////////////////////////////////////
/// Step 6 - Transform data and recalculate mean and std dev if transform matrix

	double trn_avg=0, trn_std_dev=0;
	avg_t = std_t = 0;
	if ( transmat)
	{
		for ( int i=0 ; i<svpoly_size ; ++i )
		{
			Point p0,p1;
			p0 = sv_mesh.vertices[sv_mesh.polygons[i].vidx[0]].pnt;
			p1 = sv_mesh.vertices[sv_mesh.polygons[i].vidx[1]].pnt;

			// transform points
			p0 = tm*p0;
			p1 = tm*p1;

			// distance between transformed points
			float d = dist(p0,p1);
			avg_t += d;
			std_t += d*d;
		}
		trn_avg = avg_t/svpoly_size;

		trn_std_dev= sqrt(fabs(std_t-avg_t*avg_t/(float)svpoly_size)/((float)svpoly_size-1.0));
	}

///////////////////////////////////////////////////////////////////////////////
/// Step 7 - Output the results to console and files

	if ( transmat ) // more info if there is a transformation matrix
		cout << gyral_name << " avg=" << avg << " std dev=" << std_dev
			<< " trn avg = " << trn_avg << " trn std dev = " << trn_std_dev
			<< " num pnts=" << nfound << endl;
	else
		cout << gyral_name << " avg=" << avg << " std dev=" << std_dev
			<< " num pnts=" << nfound << endl;

	// if there is an output name then use it as the base file name
	if ( outname.length() > 0 )
		filename = outname + "_" +gyral_name + ".txt";
	else // just use the gyral name
		filename = gyral_name + ".txt";
	FILE *ofile = fopen(filename.c_str(),"w");
	if ( transmat ) // more info if there is a transformation matrix
		fprintf(ofile,"%s\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%d\n",
				gyral_name.c_str(),avg,std_dev,trn_avg,trn_std_dev,nfound);
	else
		fprintf(ofile,"%s\t%6.3f\t%6.3f\t%d\n",
				gyral_name.c_str(),avg,std_dev,nfound);
	fclose(ofile);


	// lets see how long this took us to do
	time(&endtime);
	double elapsed = difftime(endtime,starttime);
	if ( verbose ) cout << "Elapsed time: "<<(elapsed/60.0)<<endl;
	if ( verbose ) cout << "Fin"<<endl;

	return EXIT_SUCCESS;
}

////////////////////////////////// FileExists //////////////////////////////
/*!
 * simple function to see if a file exists
 *
 * \param strFilename - filename to check for existance
 * \returns - true if file exists
 */
bool FileExists(string strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(), &stFileInfo);
	if (intStat == 0)
	{
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	}
	else
	{
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}

	return (blnReturn);
}

////////////////////////////////// IntersectMesh //////////////////////////////
/*!
 * This function checks to see if a line segment intersects a polygon
 * of the passed mesh. Mesh polygons corresponding to p0 and p1 are
 * ignored.
 *
 * \param p0 - first line point
 * \param p1 - second line point
 * \param pnump - mesh polygon to ignore for p0
 * \param pnumn - mesh polygon to ignore for p1
 * \param mesh - pointer to the mesh
 * \returns - true if p0-p1 line intersects a mesh polygon
 */
bool IntersectMesh(Point p0, Point p1, int pnump, int pnumn, RicMesh *mesh)
{
	// make reasonable bounding box
	float off=2;	// amount to enlarge bounding box
	float xmin,xmax,ymin,ymax,zmin,zmax;	// bounding box
	if ( p0.x < p1.x)
	{
		xmin = p0.x-off;
		xmax = p1.x+off;
	}
	else
	{
		xmin = p1.x-off;
		xmax = p0.x+off;
	}
	if ( p0.y < p1.y)
	{
		ymin = p0.y-off;
		ymax = p1.y+off;
	}
	else
	{
		ymin = p1.y-off;
		ymax = p0.y+off;
	}
	if ( p0.z < p1.z)
	{
		zmin = p0.z-off;
		zmax = p1.z+off;
	}
	else
	{
		zmin = p1.z-off;
		zmax = p0.z+off;
	}

	// check each mesh polygon to see if it intersects a line between the points
	int i;
	for ( i=0 ; i<mesh->p_size ; ++i )
	{
		if ( i==pnump || i==pnumn )
			continue;

		// see if any polygons are in a reasonable bounding box
		// made by the points
		bool found = false;
		for ( int j=0 ; j<3 ; ++j )
		{
			Point p = mesh->vertices[mesh->polygons[i].vidx[j]].pnt;

			// see if point in bounding box
			if ( p.x>xmin && p.x<xmax
					&& p.y>ymin && p.y<ymax
					&& p.z>zmin && p.z<zmax )
			{
				found = true;
				break;
			}
		}

		// if a point has been found in the bounding box then check intersection
		if ( found )
		{
			// point in box so check this polygon
			Point pp0,pp1,pp2; // polygon vertices
			pp0 = mesh->vertices[mesh->polygons[i].vidx[0]].pnt;
			pp1 = mesh->vertices[mesh->polygons[i].vidx[1]].pnt;
			pp2 = mesh->vertices[mesh->polygons[i].vidx[2]].pnt;

			// check intersection with this polygon
			Point pint;
			if ( line_intersect_triangle(pp0,pp1,pp2,p0,p1,&pint) )
				return true;
		}
	}

	// nothing found so return false
	return false;
}
