// --------------------------- ScalDepthUtil.cpp ----------------------------
/*! \file
This file contains utility routines for ricsulcaldepth.cpp (aka RicSulcalDepth)
Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA

*/

#include <RicMesh.h>
#include <RicCurve.h>
#include <RicUtil.h>
#include "RicPoint.h"
#include "RicMatrix.h"

#undef SDEPTHDEBUG

// finds the intersections of the plane with the mesh
/**
 * This routine finds the intersection of a plane and a mesh. A point is
 * calculated at each intersection of a mesh trangle and the plane.
 * @param mesh triangluar mesh
	* @param ipnts array of plane-mesh intersection points
 * @param p1 first point in plane
 * @param p2 second point in plane
 * @param p3 third point in plane
 * @return the number of intersection points found
 */
int PlaneIntersectMesh(RicMesh *mesh, Point *ipnts, Point p1, Point p2, Point p3)
{
	// Find the equation of the plane
	float a,b,c,d;
	equ_plane(p1,p2,p3,&a,&b,&c,&d);

	// Look for all the intersections between the mesh and each plane
	int ni=0;
	for ( int i=1 ; i<2 ; ++i )
	{
		for ( int j=0 ; j<mesh->p_size ; ++j )
		{
			// each polygon is really a triangle so we have three sides to try
			if ( line_intersect_plane(a,b,c,d,
					mesh->vertices[mesh->polygons[j].vidx[0]].pnt,
					mesh->vertices[mesh->polygons[j].vidx[1]].pnt , &ipnts[ni]) )
			{
				++ni;
			}
			if ( line_intersect_plane(a,b,c,d,
					mesh->vertices[mesh->polygons[j].vidx[1]].pnt,
					mesh->vertices[mesh->polygons[j].vidx[2]].pnt , &ipnts[ni]) )
			{
				++ni;
			}
			if ( line_intersect_plane(a,b,c,d,
					mesh->vertices[mesh->polygons[j].vidx[2]].pnt,
					mesh->vertices[mesh->polygons[j].vidx[0]].pnt , &ipnts[ni]) )
			{
				++ni;
			}
		}
		
	}
	return ni;
}

/**
 * Rotates a set of point in a plane to align with the xz plane. The origin
 * is set at the lower of the two farthest apart points in the points array.
 * @param pnts array of points in a plane
 * @param ni number of points in array
 * @param ang1 pointer to first rotation angle
 * @param ang2 pointer to second rotation angle
 * @param ang3 pointer to third rotation angle
 * @return the origin point
 */
Point TranRotPntsXZ(Point *pnts, int ni, float *ang1, float *ang2, float *ang3)
{
	// Find maximum distance between any two points
	// We will use this as the axis that we move around
	float pd,maxd=0;
	int ip1=0,ip2=0;
	for ( int i=0 ; i<ni ; ++i )
	{
		for ( int j=i ; j<ni ; ++j )
		{
			pd = dist(pnts[i],pnts[j]);
			if ( pd > maxd )
			{
				maxd = pd;
				ip1 = i;
				ip2 = j;
			}
		}
	}
	
	// Translate so that the lower point is 0,0,0
	Point dp1,dp2;
	if ( pnts[ip1].z < pnts[ip2].z )
	{
		dp1 = pnts[ip1];
		dp2 = pnts[ip2];
	}
	else
	{
		dp1 = pnts[ip2];
		dp2 = pnts[ip1];
	}
	
	// Set origin to lower point
	Point origin = dp1;
	
	// make a new offset array of points
	Point *opnts;
	opnts = new Point[ni];
	for ( int i=0 ; i<ni ; ++i )
		opnts[i] = pnts[i]-origin;
	
	// Find first axis to rotate about
	Point vpnt,vpnt2,vpnt3,vpnt4,apnt,dpnt;
	vpnt = dp2-dp1;
	
	// angle rotate around z
	*ang1 = atan2(vpnt.y,vpnt.x);
	
	// rotation matrix
	CMatrix mat;	
	mat.Rotate(*ang1,0,0,1);

	Point *rpnts1;
	rpnts1 = new Point[ni];
	for ( int i=0 ; i<ni ; ++i )
		rpnts1[i] = mat*opnts[i];
		
	// rotate around y
	vpnt2 = mat*vpnt;
	*ang2 = atan2(vpnt2.x,vpnt2.z);
	
	mat.Identity();
	mat.Rotate(*ang2,0,1,0);
	vpnt3 = mat*vpnt2;
	
	Point *rpnts2;
	rpnts2 = new Point[ni];
	for ( int i=0 ; i<ni ; ++i )
		rpnts2[i] = mat*rpnts1[i];
	
	// find angle between the plane of the points and the xz plane
	// this so that the points all end up in the xz plane
	
	// equation of plane for xz
	Point pp1,pp2,pp3;
	pp1.x=0; pp1.y=0;pp1.z=1;
	pp2.x=0; pp2.y=0;pp2.z=0;
	pp3.x=1; pp3.y=0;pp3.z=0;
	float a1,b1,c1,d1,a2,b2,c2,d2;
	
	equ_plane(pp1,pp2,pp3,&a1,&b1,&c1,&d1); // zx plane
	
	// pick three points from array for plane
	pp1 = rpnts2[0];
	pp2 = rpnts2[(int)(ni*0.25)];
	pp3 = rpnts2[(int)(ni*0.75)];
	
	equ_plane(pp1,pp2,pp3,&a2,&b2,&c2,&d2); // zx plane
		
	*ang3 = angle_between_planes(a2,b2,c2,a1,b1,c1);
	
	mat.Identity();
	mat.Rotate(*ang3,0,0,1);
	vpnt4 = mat*vpnt3;
	
	for ( int i=0 ; i<ni ; ++i )
		pnts[i] = mat*rpnts2[i];
	
#ifdef SDEPTHDEBUG
	// write these points as a mesh to see what we get
	RicMesh tmesh(ni,0,ni-1,2);
	for ( int i=0; i<ni ; ++i )
	{
		tmesh.vertices[i].pnt.x = rpnts2[i].x+dp1.x;
		tmesh.vertices[i].pnt.y = rpnts2[i].y+dp1.y;
		tmesh.vertices[i].pnt.z = rpnts2[i].z+dp1.z;
	}
		
	for ( int i=0; i<ni-1 ; ++i )
	{
		tmesh.polygons[i].vidx[0] = i;
		tmesh.polygons[i].vidx[1] = i+1;
		tmesh.polygons[i].vidx[2] = 0;
	}
	
	tmesh.write_mesh_txt("RotTrans.mesh");
#endif

	// clean up memory
	delete [] opnts;
	delete [] rpnts1;
	delete [] rpnts2;
	
	return origin;
}


/**
 * Rotates a set of point in a plane to align with the xz plane. The origin
 * is set at the first point in the points array with the z direction defined
 * by the second point in the array.
 * @param pnts array of points in a plane
 * @param ni number of points in array
 * @param ang1 pointer to first rotation angle
 * @param ang2 pointer to second rotation angle
 * @param ang3 pointer to third rotation angle
 * @return the origin point
 */
Point TranRotPntsXZ2(Point *pnts, int ni, float *ang1, float *ang2, float *ang3)
{

	// Translate so that the lower point is 0,0,0
	Point dp1,dp2;
	if ( pnts[0].z < pnts[1].z )
	{
		dp1 = pnts[0];
		dp2 = pnts[1];
	}
	else
	{
		dp1 = pnts[1];
		dp2 = pnts[0];
	}

	// Set origin to lower point
	Point origin = dp1;

	// make a new offset array of points
	Point *opnts;
	opnts = new Point[ni];
	for ( int i=0 ; i<ni ; ++i )
		opnts[i] = pnts[i]-origin;

	// Find first axis to rotate about
	Point vpnt,vpnt2,vpnt3,vpnt4,apnt,dpnt;
	vpnt = dp2-dp1;

	// angle rotate around z
	*ang1 = atan2(vpnt.y,vpnt.x);

	// rotation matrix
	CMatrix mat;
	mat.Rotate(*ang1,0,0,1);

	Point *rpnts1;
	rpnts1 = new Point[ni];
	for ( int i=0 ; i<ni ; ++i )
		rpnts1[i] = mat*opnts[i];

	// rotate around y
	vpnt2 = mat*vpnt;
	*ang2 = atan2(vpnt2.x,vpnt2.z);

	mat.Identity();
	mat.Rotate(*ang2,0,1,0);
	vpnt3 = mat*vpnt2;

	Point *rpnts2;
	rpnts2 = new Point[ni];
	for ( int i=0 ; i<ni ; ++i )
		rpnts2[i] = mat*rpnts1[i];

	// find angle between the plane of the points and the xz plane
	// this so that the points all end up in the xz plane

	// equation of plane for xz
	Point pp1,pp2,pp3;
	pp1.x=0; pp1.y=0;pp1.z=1;
	pp2.x=0; pp2.y=0;pp2.z=0;
	pp3.x=1; pp3.y=0;pp3.z=0;
	float a1,b1,c1,d1,a2,b2,c2,d2;

	equ_plane(pp1,pp2,pp3,&a1,&b1,&c1,&d1); // zx plane

	// pick three points from array for plane
	pp1 = rpnts2[0];
	pp2 = rpnts2[(int)(ni*0.25)];
	pp3 = rpnts2[(int)(ni*0.75)];

	equ_plane(pp1,pp2,pp3,&a2,&b2,&c2,&d2); // zx plane

	*ang3 = angle_between_planes(a2,b2,c2,a1,b1,c1);

	mat.Identity();
	mat.Rotate(*ang3,0,0,1);
	vpnt4 = mat*vpnt3;

	for ( int i=0 ; i<ni ; ++i )
		pnts[i] = mat*rpnts2[i];

#ifdef SDEPTHDEBUG
	// write these points as a mesh to see what we get
	RicMesh tmesh(ni,0,ni-1,2);
	for ( int i=0; i<ni ; ++i )
	{
		tmesh.vertices[i].pnt.x = rpnts2[i].x+dp1.x;
		tmesh.vertices[i].pnt.y = rpnts2[i].y+dp1.y;
		tmesh.vertices[i].pnt.z = rpnts2[i].z+dp1.z;
	}

	for ( int i=0; i<ni-1 ; ++i )
	{
		tmesh.polygons[i].vidx[0] = i;
		tmesh.polygons[i].vidx[1] = i+1;
		tmesh.polygons[i].vidx[2] = 0;
	}

	tmesh.write_mesh_txt("RotTrans.mesh");
#endif

	// clean up memory
	delete [] opnts;
	delete [] rpnts1;
	delete [] rpnts2;

	return origin;
}

/**
 * This routine measures the depth of an area enclosed by a set of points.
 * The set is already aligned in the xz plane along the z axis.  The depth
 * measurement is a line that goes from bottom z to top z following the
 * center of the set of points. For each z level, the average of the min
 * and max point x values is chosen as the center.
 *
 * @param pnts array of points aligned in the xz plane
 * @param n number of points
 * @param ndepth - number of intermediate measurements along center line for depth
 * @param dpnts depth line points
 * @return returns depth
 */
float SulcalDepthMeasures(Point *pnts, int n, int ndepth, Point *dpnts)
{
	// find min and max z points in the input points array
	float minz=10000,maxz=-10000;
	for ( int i=0 ; i<n ; ++i )
	{
		if ( pnts[i].z < minz )
		{
			minz = pnts[i].z;
			dpnts[0] = pnts[i];
		}
		if ( pnts[i].z > maxz )
		{
			maxz = pnts[i].z;
			dpnts[ndepth+1] = pnts[i];
		}
	}

	// increment to divide z space into for intermediate points
	float zinc = 0.5*(maxz-minz)/(ndepth+1);

	// find the intermediate points
	for ( int k=0 ; k<ndepth ; ++k )
	{
		float minx=10000,maxx=-10000;
		float z = minz+(k+1)*2*zinc;
		float min = z-zinc;
		float max = z+zinc;

		for ( int i=0 ; i<n ; ++i )
		{
			if ( pnts[i].z >= min && pnts[i].z <= max )
			{
				if ( pnts[i].x < minx )
				{
					minx = pnts[i].x;
				}
				if ( pnts[i].x > maxx )
				{
					maxx = pnts[i].x;
				}
			}
		}

		// check to make sure we actually found points
		if ( minx == 10000 || maxx == -10000 )
		{
			cerr << "SulcalDepthMeasures - cannot find z level\n";
			return 0;
		}

		// assign intermediate points with the x value the average
		dpnts[k+1].x = (minx+maxx)*0.5;
		dpnts[k+1].y = 0;
		dpnts[k+1].z = z;
	}

	float depth=0;
	for ( int i=0 ; i<ndepth+1 ; ++i )
		depth += dist(dpnts[i],dpnts[i+1]);

	return depth;
}

/**
 * This routine measures the depth of an area enclosed by a set of points.
 * The set is already aligned in the xz plane along the z axis.  The depth
 * measurement is a line that goes from bottom z to top z following the
 * center of the set of points. For each z level, the average of the min
 * and max point x values is chosen as the center.
 *
 * @param pnts array of points aligned in the xz plane
 * @param n number of points
 * @param ndepth - number of intermediate measurements along center line for depth
 * @param dpnts depth line points
 * @return returns depth
 */
float SulcalDepthMeasures0(Point *pnts, int n, int ndepth, Point *dpnts,
		Point bpnt, Point tpnt)
{
	float depth=0;

	// find min and max z points in the input points array
	float minz=10000,maxz=-10000;
/*	for ( int i=0 ; i<n ; ++i )
	{
		if ( pnts[i].z < minz )
		{
			minz = pnts[i].z;
			dpnts[0] = pnts[i];
		}
		if ( pnts[i].z > maxz )
		{
			maxz = pnts[i].z;
			dpnts[ndepth+1] = pnts[i];
		}
	}
*/
	// find min and max z points in the input points array
	minz=pnts[0].z;
	maxz=pnts[1].z;
	dpnts[0] = pnts[0];
	dpnts[ndepth+1] = pnts[1];
	if ( minz > maxz )
	{
		minz = pnts[1].z;
		maxz = pnts[0].z;
		dpnts[0] = pnts[1];
		dpnts[ndepth+1] = pnts[0];
	}

	// increment to divide z space into for intermediate points
	float zinc = 0.5*(maxz-minz)/(ndepth+1);

	// find the intermediate points
	for ( int k=1 ; k<ndepth+1 ; ++k )
	{
		float minx=10000,maxx=-10000;
		float z = minz+k*2*zinc;
		float min = z-zinc;
		float max = z+zinc;

		for ( int i=0 ; i<n ; ++i )
		{
			if ( pnts[i].z >= min && pnts[i].z <= max )
			{
				if ( pnts[i].x < minx )
				{
					minx = pnts[i].x;
				}
				if ( pnts[i].x > maxx )
				{
					maxx = pnts[i].x;
				}
			}
		}

		// check to make sure we actually found points
		if ( minx == 10000 || maxx == -10000 )
		{
			cerr << "SulcalDepthMeasures0 - cannot find z level\n";
			minx = maxx =  0; // put on midline
		}

		// assign intermediate points with the x value the average
		dpnts[k].x = (minx+maxx)*0.5;
		dpnts[k].y = 0;
		dpnts[k].z = z;

		// check for big outliers and put them on midline
		if ( dpnts[k].x > maxz*0.5 )
			dpnts[k].x = 0;
	}

	for ( int i=0 ; i<ndepth+1 ; ++i )
		depth += dist(dpnts[i],dpnts[i+1]);

	return depth;
}
