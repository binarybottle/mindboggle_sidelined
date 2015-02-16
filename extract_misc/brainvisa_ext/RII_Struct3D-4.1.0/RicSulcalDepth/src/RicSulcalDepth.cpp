//////////////////////////// RicSulcalDepth.cpp /////////////////////////////
// This program calculates the depth and thickness of a sulcus given the
// sulcal mesh and the lines defining the top and bottom of the sulcus.

/*! \file
Implementation file for RicSulcalDepth.cpp (aka RicSulcalDepth)
Copyright (C) 2007-2010 by Bill Rogers - Research Imaging Institute - UTHSCSA
*/

/*! \mainpage
This program determines depth of a sulcus based on the sulcal
mesh and the lines defining the top and bottom of the sulcus.

The top and bottom lines will be split into a set of segments. For each segment
a plane will be cut at right angles through the sulcus. A set of points will be
extracted from this mesh at the intersection of the plane. A line will be
traced through the center of this set of points, defining a sulcal depth line.
A fraction of the sulcal lines will be ignored at the ends.

The program output max and average values for sulcal depth.
In addition a mesh is output with the suffix _Depth.mesh that shows the
sulcal depth lines.


Command line switches:

--sm	sulcal mesh file (required)

--tm	sulcal top line mesh file (required)

--bm	sulcal bottom line mesh file (required)

--ns	number of segments to split sulcus into

--nw	number of width values for each segment

--fs	fraction of each end of the sulcal lines to skip

--fg	fill gaps by adding depth measurements where mesh is missing

-o		output file base name

-s		subject name

-v		verbose output

--bv	output for BrainVisa

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include "SulcalDepthUtil.h"
#include <RicMesh.h>
#include <RicCurve.h>
#include <RicUtil.h>
#include <RicPoint.h>
#include <RicMatrix.h>

#undef SDEPTHDEBUG

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string sm_name;				// sulcal mesh name
	string tm_name;				// sulcal top line mesh name
	string bm_name;				// sulcal bottom line mesh name
	string outname;				// base name for output files
	int	   nseg=10;				// number of segments to split sulcus into
	int    ndepth=6;			// number of depth measurements for each segment
	float percent_skip=0.05;	// fraction of end of sulcal lines to skip
	bool fillgaps=false;		// if true then measure depth in mesh gaps
	string subjname;			// subject name
	bool verbose=false;			// if true output more stuff on stdout
	bool brainvisa=false;		// if true output stuff on stdout for brain visa

	///////////////////////////////////////////////////////////////////////////
	/// Read in the command line using tclap

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("Atlas Sulcal Depth and Thickness Processor", ' ', \
			"0.5");

		// sulcal mesh
		TCLAP::ValueArg<string> smname( "","sm", "Sulcal mesh file", true, "",\
			"string");
		cmd.add(smname);

		// bottom line mesh
		TCLAP::ValueArg<string> bmname( "","bm", "Sulal bottom line mesh file",\
			true, "", "string");
		cmd.add(bmname);

		// top line mesh
		TCLAP::ValueArg<string> tmname( "","tm", "Sulal top line mesh file",\
			true, "", "string");
		cmd.add(tmname);

		// output file base name
		TCLAP::ValueArg<string> oname( "o","o", "Base output file name", false,\
			 "RicSulcalDepth", "string");
		cmd.add(oname);

		// number of segments
		TCLAP::ValueArg<int> ns( "","ns", \
			"nseg  is the number of segments to split the sulcus into (default 10)",\
			false, 10, "int");
		cmd.add(ns);

		// number of depth measures per segment
		TCLAP::ValueArg<int> nd( "","nd", \
			"ndepth  is the number of depth increments per segment (default 2)",\
			false, 2, "int");
		cmd.add(nd);

		// fraction of each end of the sulcal lines to skip
		TCLAP::ValueArg<float> fs( "","fs", \
			"fraction of each end of the sulcal lines to skip (default 0.05)",\
			false, 0.05, "float");
		cmd.add(fs);

		TCLAP::SwitchArg fg( "","fg", "Fill gaps with depth lines", false);
		cmd.add(fg);

		// subject name name just for peterk
		TCLAP::ValueArg<string> sname( "s",\
			"s", "Subject name", false, "FBSubject", "string");
		cmd.add(sname);

		TCLAP::SwitchArg ver( "v","verbose", "Verbose output on stdout", false);
		cmd.add(ver);

		TCLAP::SwitchArg bvisa( "b","bv", "Brain Visa output on stdout", false);
		cmd.add(bvisa);

		// parse the command line
		cmd.parse( argc, argv );

		// copy command line variable to program variables
		sm_name = smname.getValue();
		tm_name = tmname.getValue();
		bm_name = bmname.getValue();
		outname = oname.getValue();
		nseg = ns.getValue();
		ndepth = nd.getValue();
		percent_skip = fs.getValue();
		fillgaps = fg.getValue();
		subjname = sname.getValue();
		verbose = ver.getValue();
		brainvisa = bvisa.getValue();

	}catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cerr << " Command Line Error" << e.error() << " for arg " << e.argId()\
			<< endl;
		if ( brainvisa ) cout << " Command Line Error" << e.error() << \
			" for arg " << e.argId() << endl;
	}

	///////////////////////////////////////////////////////////////////////////
	/// Read in input mesh files

	// check for valid sulcal mesh file
	RicMesh sm_mesh((char*)sm_name.c_str());
	if ( sm_mesh.p_size != 0 )
	{
		if ( verbose ) cerr << "We just read in "<<sm_name<<endl;
		if ( verbose ) cerr << "Number of vertices read: " << \
			sm_mesh.v_size<<endl;
	}
	else
	{
		cerr << "Error reading file "<<sm_name<<endl;
		if ( brainvisa ) cout << "Error reading file "<<sm_name<<endl;
		exit(1);
	}

	// check for valid sulcal top line mesh file
	RicMesh tm_mesh((char*)tm_name.c_str());
	if ( tm_mesh.p_size != 0 )
	{
		if ( verbose ) cerr << "We just read in "<<tm_name<<endl;
		if ( verbose ) cerr << "Number of vertices read: "<<tm_mesh.v_size<<endl;
	}
	else
	{
		cerr << "Error reading file "<<tm_name<<endl;
		if ( brainvisa ) cout << "Error reading file "<<tm_name<<endl;
		exit(1);
	}

	// check for valid sulcal bottom line mesh file
	RicMesh bm_mesh((char*)bm_name.c_str());
	if ( bm_mesh.p_size != 0 )
	{
		if ( verbose ) cerr << "We just read in "<<bm_name<<endl;
		if ( verbose ) cerr << "Number of vertices read: "<<bm_mesh.v_size<<endl;
	}
	else
	{
		cerr << "Error reading file "<<bm_name<<endl;
		if ( brainvisa ) cout << "Error reading file "<<bm_name<<endl;
		exit(1);
	}

	///////////////////////////////////////////////////////////////////////////
	/// Step by Step Processing
	///


	///////////////////////////////////////////////////////////////////////////
	/// Step 1 - Sample the lines with equal spacing

	// pick a percentage of the points to skip and beginning and end of line
	int nskip;

	// top
	nskip = (int)(tm_mesh.v_size * percent_skip);
	RicCurve tcurve(tm_mesh.v_size-2*nskip, 6, 1);
	for (int i=0; i<tm_mesh.v_size-2*nskip;i++)
	{
		tcurve.cpnts[i].x = tm_mesh.vertices[nskip+i].pnt.x;
		tcurve.cpnts[i].y = tm_mesh.vertices[nskip+i].pnt.y;
		tcurve.cpnts[i].z = tm_mesh.vertices[nskip+i].pnt.z;
	}

	// initialize output curve
	tcurve.InitCurve();

	// interpolate a bunch of points
	tcurve.Interpolate(nseg+2);

	// bottom
	nskip = (int)(bm_mesh.v_size * percent_skip);
	RicCurve bcurve(bm_mesh.v_size-2*nskip, 6, 1);

	// include a check to see that the top and bottom start out at about
	// the same point - if not then reverse the order of the bottom
	if ( dist(tm_mesh.vertices[0].pnt,bm_mesh.vertices[0].pnt) <
		dist(tm_mesh.vertices[tm_mesh.v_size-1].pnt,bm_mesh.vertices[0].pnt) )
	{
		// order of points ok then just copy
		for (int i=0; i<bm_mesh.v_size-2*nskip;i++)
		{
			bcurve.cpnts[i].x = bm_mesh.vertices[nskip+i].pnt.x;
			bcurve.cpnts[i].y = bm_mesh.vertices[nskip+i].pnt.y;
			bcurve.cpnts[i].z = bm_mesh.vertices[nskip+i].pnt.z;
		}
	}
	else
	{
		// need to reverse the order of the points
		int nn = bm_mesh.v_size-nskip-1;
		for (int i=0; i<bm_mesh.v_size-2*nskip;i++)
		{
			bcurve.cpnts[i].x = bm_mesh.vertices[nn-i].pnt.x;
			bcurve.cpnts[i].y = bm_mesh.vertices[nn-i].pnt.y;
			bcurve.cpnts[i].z = bm_mesh.vertices[nn-i].pnt.z;
		}
	}

	// initialize output curve
	bcurve.InitCurve();

	// interpolate a bunch of points
	// two more than desires as we will discard first and last points
	bcurve.Interpolate(nseg+2);

	// Assign curve points to local array (float instead of double)
	Point *bpnts,*tpnts;
	bpnts = new Point[nseg+2];
	tpnts = new Point[nseg+2];
	for ( int i=0 ; i<nseg+2 ; ++i )
	{
		bpnts[i].x = bcurve.pnts[i].x;
		bpnts[i].y = bcurve.pnts[i].y;
		bpnts[i].z = bcurve.pnts[i].z;
		tpnts[i].x = tcurve.pnts[i].x;
		tpnts[i].y = tcurve.pnts[i].y;
		tpnts[i].z = tcurve.pnts[i].z;
	}

#ifdef SDEPTHDEBUG
	if ( verbose )
	{
		// Output a mesh of segment lines - skip first and last points
		RicMesh omesh(2*nseg,0,nseg,2);
		for ( int i=0,idx=0 ; i<nseg ; ++i,idx+=2 )
		{
			omesh.vertices[idx].pnt = bpnts[i+1];
			omesh.vertices[idx+1].pnt = tpnts[i+1];

			omesh.polygons[i].vidx[0] = idx;
			omesh.polygons[i].vidx[1] = idx+1;
			omesh.polygons[i].vidx[2] = 0;
		}

		string filename;
		filename = outname + "_slines.mesh";
		omesh.write_mesh_txt((char*)filename.c_str());


	}
#endif

	///////////////////////////////////////////////////////////////////////////
	///  Step 2 - Loop for each segment

	// make arrays to hold stuff
	float *depth;	// depth measurements for each segment
	Point **depthpnts;	// array of depth points for segments
	matrix(&depthpnts,nseg,ndepth+2);
	depth = new float[nseg];
	bool *gapflag;	// if true then fill gap in mesh with measurement line
	gapflag = new bool[nseg];
	Point *ipnts0,*ipnts; 	// intersection points
	ipnts0 = new Point[sm_mesh.v_size];	// don't think it can be bigger than this
	ipnts = new Point[sm_mesh.v_size];	// don't think it can be bigger than this
	int		nfnd=0;	// number of depth measurements found

	for ( int n=0 ; n<nseg ; ++n )
	{
		///////////////////////////////////////////////////////////////////////////
		/// Step 2a
		/// Create a plane for each pair of top and bottom points
		/// (skipping first and last)
		/// Average four normals, two each at top and bottom.
		/// Each normal is made from a point on the line, the comparable
		/// point on the other and the nearest neighbor

		Vector v,v1,v2,v3,v4;	// normal vectors
		Point pp1,pp2,pp3;		// points for plane
		gapflag[n] = false;		// don't fill by default

		// make four normals - try to make them point the same way
		// remember to ignore first and last - make i+1 the start
		v1 = normal_pnt(bpnts[n],bpnts[n+1],tpnts[n+1]);
		v2 = normal_pnt(bpnts[n+1],tpnts[n+1],tpnts[n]);
		v3 = normal_pnt(tpnts[n+2],tpnts[n+1],bpnts[n+1]);
		v4 = normal_pnt(tpnts[n+1],bpnts[n+1],bpnts[n+2]);

		// add the four together to make the average
		v = v1+v2+v3+v4;

		// use top and bottom points along with the normal added to the top
		// this specifies the plane through the mesh
		pp1 = bpnts[n+1];
		pp2 = tpnts[n+1];
		pp3 = tpnts[n+1]+v;

		// get equation of plane
		float a,b,c,d;
		equ_plane(pp1,pp2,pp3,&a,&b,&c,&d);

		///////////////////////////////////////////////////////////////////////////
		/// Step 2b
		/// find intersection points of plane with mesh
		int ni0;			// number of intersection points
		ni0 = PlaneIntersectMesh(&sm_mesh, &ipnts0[2], pp1, pp2, pp3);

		// Make first two points the top and bottom points from sulcal lines
		ipnts[0] = bpnts[n+1];
		ipnts[1] = tpnts[n+1];
		int ni=2; // start filling in after first two points

		// Do a little check to see if any points are a long way from our top
		// and bottom points (pp1 & pp2). This can happen where the plane intersects
		// another part of the sulcus away from our intended segment.
		Point avgpnt;
		avgpnt = pp1+pp2;
		avgpnt.x *= 0.5f;
		avgpnt.y *= 0.5f;
		avgpnt.z *= 0.5f;
		float tbdist = 0.7f * dist(pp1,pp2);// use dist between top and bot as gauge
		float td;	// test distance
		for ( int i=0; i<ni0 ; ++i )
		{
			td = dist(avgpnt,ipnts0[i]);
			if ( td < tbdist ) // point probably ok so include
			{
				ipnts[ni++] = ipnts0[i];
				if ( ipnts[ni-1].x == 0 )
					cout <<ni<<' '<<ipnts[ni].x <<' '<<ipnts[ni].y<<' '<<ipnts[ni].z << endl;
			}
		}


		// if there are not enough points then interpolate a couple of points
		// in middle if fillgaps in true
		if ( (ni < 4) )
		{
			if ( fillgaps ) // then create a depth line by interpolation
			{
				cerr << "RicSulcalDepth - ni < 4 - interpolating" << endl;
				float dx,dy,dz;
				dx = pp2.x-pp1.x;
				dy = pp2.y-pp1.y;
				dz = pp2.z-pp1.z;
				ipnts0[2].x = pp1.x + 0.333f*dx;
				ipnts0[2].y = pp1.y + 0.333f*dy;
				ipnts0[2].z = pp1.z + 0.333f*dz;
				ipnts0[3].x = pp1.x + 0.666f*dx;
				ipnts0[3].y = pp1.y + 0.666f*dy;
				ipnts0[3].z = pp1.z + 0.666f*dz;
				ni = 4;
			}
			else // skip this measurement due to gap in mesh
			{
				cerr << "RicSulcalDepth - skipping gap in mesh" << endl;
				continue;
			}
		}

#ifdef SDEPTHDEBUG
		// write these points as a mesh to see what we get
		RicMesh pmesh(ni,0,ni-1,2);
		for ( int i=0; i<ni ; ++i )
		{
			pmesh.vertices[i].pnt = ipnts[i];
		}

		for ( int i=0; i<ni-1 ; ++i )
		{
			pmesh.polygons[i].vidx[0] = i;
			pmesh.polygons[i].vidx[1] = i+1;
			pmesh.polygons[i].vidx[2] = 0;
		}

		pmesh.write_mesh_txt("Intersect.mesh");
#endif

		///////////////////////////////////////////////////////////////////////////
		/// Step 2c
		/// Translate and rotate so that we will be in xz plane on z axis
		/// which will greatly simplify the processing.
		float ang1, ang2, ang3;
		Point origin;
		origin = TranRotPntsXZ2(ipnts, ni, &ang1, &ang2, &ang3);

		///////////////////////////////////////////////////////////////////////////
		/// Step 2d
		/// Now make depth measurements
		/// The depth follows a line through the center of the sulcus from
		/// bottom z to top z.

		depth[nfnd] = SulcalDepthMeasures0(ipnts, ni, ndepth, depthpnts[nfnd],
				ipnts[0], ipnts[1]);

		///////////////////////////////////////////////////////////////////////////
		/// Step 2e
		/// Rotate and translate back to the original orientation
		CMatrix mat;

		mat.Rotate(-ang3,0,0,1);

		for ( int i=0 ; i<ndepth+2 ; ++i )
			depthpnts[nfnd][i] = mat * depthpnts[nfnd][i];


		mat.Identity();
		mat.Rotate(-ang2,0,1,0);

		for ( int i=0 ; i<ndepth+2 ; ++i )
			depthpnts[nfnd][i] = mat * depthpnts[nfnd][i];


		mat.Identity();
		mat.Rotate(-ang1,0,0,1);

		for ( int i=0 ; i<ndepth+2 ; ++i )
			depthpnts[nfnd][i] = mat * depthpnts[nfnd][i];


		// now translate back
		for ( int i=0 ; i<ndepth+2 ; ++i )
			depthpnts[nfnd][i] += origin;

		// increment counter for number of depth measurements
		++nfnd;
	}

	/// End of loop for each segment

	///////////////////////////////////////////////////////////////////////////
	/// Step 3
	/// Now make a nice mesh to show what we did.
	/// This mesh will show lines down the center of each slice through the
	/// sulcus.

	// Figure out how big a mesh we need
	// This will be a 2d mesh of line segments

	// vertices is nseg*(ndepth+2) for depth lines and nseg*(ndepth*2) for thickness
	int nverts = nfnd*(ndepth+2);

	// polygons (line segments) nseg*(ndepth+1) for depth lines
	int npoly = nfnd*(ndepth+1);

	RicMesh line_mesh(nverts,0,npoly,2);

	// fill the mesh a depth line for each segment
	int vidx=0,pidx=0;
	for ( int n=0 ; n<nfnd ; ++n )
	{

		// depth line vertices
		for ( int i=0 ; i<ndepth+2 ; ++i )
		{
			// add the depth line points
			line_mesh.vertices[vidx+i].pnt = depthpnts[n][i];
		}

		// depth line polygons (line segments)
		for ( int i=0 ; i<ndepth+1 ; ++i )
		{
			line_mesh.polygons[pidx].vidx[0] = vidx+i;
			line_mesh.polygons[pidx].vidx[1] = vidx+i+1;
			line_mesh.polygons[pidx].vidx[2] = 0;
			++pidx;
		}
		vidx += (ndepth+2);
	}

	// now write this sucker out as a mesh
	string filename;
	filename = outname + "_Depth.mesh";
	line_mesh.write_mesh_txt((char*)filename.c_str());

	// for Peter K write out the depth measurement for each segment
	filename = outname + "_SegDepth.txt";
	ofstream fout(filename.c_str());
	fout<< subjname;
	for ( int n=0 ; n<nfnd ; ++n )
		fout << '\t' << depth[n];
	fout << endl;
	fout.close();

#ifdef SDEPTHDEBUG
	// also write out first and last points of top curve
	if ( verbose )
	{
		filename = outname + "_EndPnts.txt";
		ofstream fout(filename.c_str());
		fout<< subjname << '\t';
		fout<< tpnts[0].x << '\t' << tpnts[0].y << '\t' << tpnts[0].z << '\t';
		fout<< tpnts[nseg+1].x << '\t' << tpnts[nseg+1].y << '\t' << tpnts[nseg+1].z << endl;
		fout.close();
	}
#endif

	///////////////////////////////////////////////////////////////////////////
	/// Step 4
	/// Finally, lets us figure out some stats like average depth and max depth.
	float maxdepth=0, avgdepth=0;

	for ( int n=0 ; n<nfnd ; ++n )
	{
		avgdepth += depth[n];
		if ( depth[n] > maxdepth ) maxdepth = depth[n];

	}
	avgdepth /= nfnd;


	if ( !brainvisa && !verbose ) cout << subjname << ", " <<maxdepth << ", "
		<< avgdepth << endl;

	if ( brainvisa || verbose )
	{
		cout << "Num segments requested = " << nseg << endl;
		cout << "Num segments measured = " << nfnd << endl;
		cout << "Max depth = " << maxdepth << endl;
		cout << "Avg depth = " << avgdepth << endl;
	}

	// do some memory housekeeping
	delete [] ipnts;
	delete [] tpnts;
	delete [] bpnts;
	delete [] depth;
	free_matrix(depthpnts);

	if ( verbose ) cout << "That is All!" << endl;

	return EXIT_SUCCESS;
}
