/* ------------------------ RicUtil.cpp -------------------------------- */
/*! \mainpage
The RicUtil library contains lots of useful functions for dealing with points,
arrays, and geometry. There are routines for multidimensional matrix allocation
and coordinate conversion.

These utilities have been collected for years, mostly during my tenure in
Rehabilitation Medicine. This as been extended for use with RIC projects.
Stuff has been added in random fashion over a long time. One of these days,
the obsolete functions will get deleted and the rest reorganized.

There are Numerical Recipes in C cubic spline interpolation routines. They
are made thoroughly obsolete by the RicCurve (SISL) library. The NR in
C routines are documented in the NR book.

The documentation for newer functions is in doxygen format. In all the older
functions, the parameter and return values are specified but not separated out
into doxygen params.
Bill Rogers
*/

/*! \file
This file contains the math utilities. There are matrix allocation
routines, geometry routines, and coordinate conversion routines.

*/
// update 1-4-08

#ifdef WINDOWS
#include "stdafx.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include "RicUtil.h"
#include "RicPoint.h"
#include <math.h>
#include <float.h>

// there must be an error message routine somewhere that looks like this
void sutil_error(char *msg);
int fcompare (const void * _a, const void * _b);

#undef MIN
#undef MAX
/// min of two numbers
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
/// max of two numbers
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

//#define TRIG_TABLE
#undef TRIG_TABLE // use for table lookups - math shortcut
#ifdef TRIG_TABLE

#define TRIG_INC 0.1 ///< 0.1 degree increments
#define TRIG_LEN 3601	///< 3600 entries for 0.1 increments - one extra for zero wrap
static int	trig_init = 0;	///< 1 if trig lookup table is initialized
static double cos_table[TRIG_LEN];	///< trig cosine function lookup table
static double sin_table[TRIG_LEN];	///< trig sine function lookup table

// ------------------------ init_trig_table.cpp ---------------------------
/*!

 * int init_trig_table(void)
 *
 * returns 1 on success or 0 on failure
 *
 * This routine sets up a trig table to be used with cart_to_cyl and
 * cyl_to_cart. It provides a way for very fast trig conversions
 * using a table lookup. Not really necessary for modern CPU but might
 * be usefull in an embedded application.
 */
int init_trig_table(void)
{
	int i;
	for ( i=0 ; i<(TRIG_LEN-1) ; ++i )
	{
		cos_table[i] = cos(i*TRIG_INC*DTOR);
		sin_table[i] = sin(i*TRIG_INC*DTOR);
	}

	// wrap around zero
	cos_table[TRIG_LEN-1] = cos_table[0];
	sin_table[TRIG_LEN-1] = sin_table[0];

	trig_init = 1;

	return 1;
}

// ------------------------ cart_to_cyl.cpp ---------------------------
/*!
 * Cylind cart_to_cyl(Point cart)
 *
 * cart - x y z cartisian coordinates
 *
 * This routine converts an xyz point to cylindrical coordinates.
 */
Cylind cart_to_cyl(Point cart)
{

	Cylind  cyl;	/* converted point */

	cyl.r = (float)sqrt(cart.x*cart.x+cart.y*cart.y);
		if ( cart.x == 0.0f ) cart.x = .0001f;
	cyl.theta = (float)atan2(cart.y,cart.x);
	cyl.z = cart.z;
	return(cyl);
}

// ------------------------ cyl_to_cart_deg.cpp ---------------------------
/*!
 * Point cyl_to_cart_deg(Cylind cyl)
 *
 * cyl - cylindrical coordinates of point
 *
 * This routine converts a cylindrical point to xyz coordinates.
 */
Point cyl_to_cart(Cylind cyl)
{
	if ( !trig_init ) init_trig_table();

	Point   cart;

	cart.x = cyl.r*(float)cos(cyl.theta);
	cart.y = cyl.r*(float)sin(cyl.theta);
	cart.z = cyl.z;
	return(cart);

}


// ------------------------ cart_to_cyl_deg.cpp ---------------------------
/*!
 * Cylind cart_to_cyl_deg(Point cart)
 *
 * cart - x y z cartisian coordinates
 *
 * This routine converts an xyz point to cylindrical coordinates with
 * theta in degrees.
 */
Cylind cart_to_cyl_deg(Point cart)
{

	Cylind  cyl;            /* converted point */

	if ( cart.x == ERRVAL ) // check for invalid point
	{
		cyl.r = cyl.theta = ERRVAL;
		cyl.z = cart.z;
	}
	else
	{
		cyl.r = (float)sqrt(cart.x*cart.x+cart.y*cart.y);
		if ( cart.x == 0.0f ) cart.x = .0001f;
		cyl.theta = RTOD * (float)atan2(cart.y,cart.x);
		if ( cyl.theta < 0.0f ) cyl.theta += 360.0f;
		if ( cyl.theta >= 360.0f ) cyl.theta =- 360.0f;
		cyl.z = cart.z;
	}

	return(cyl);
}

// ------------------------ cyl_to_cart_deg.cpp ---------------------------
/*!
 * Point cyl_to_cart_deg(Cylind cyl)
 *
 * cyl - cylindrical coordinates of point
 *
 * This routine converts a cylindrical point to xyz coordinates.  The
 * theta values is in degrees.
 */
Point cyl_to_cart_deg(Cylind cyl)
{

	if ( !trig_init ) init_trig_table();

	Point   cart;

	// check for error value
	if ( cyl.r == ERRVAL )
	{
		cart.x = cart.y = ERRVAL;
		cart.z = cyl.z;
	}
	else
	{
		if ( cyl.theta < 0 ) cyl.theta += 360;
		if ( cyl.theta > 360 ) cyl.theta -= 360;
		int idx = (int)(cyl.theta / TRIG_INC);
		cart.x = (float)(cyl.r*cos_table[idx]);
		cart.y = (float)(cyl.r*sin_table[idx]);
		cart.z = cyl.z;
	}
	return(cart);

}

#else

// ------------------------ cart_to_cyl.cpp ---------------------------
/*!
 * Cylind cart_to_cyl(Point cart)
 *
 * cart - x y z cartisian coordinates
 *
 * This routine converts an xyz point to cylindrical coordinates.
 */
Cylind cart_to_cyl(Point cart)
{

	Cylind  cyl;	/* converted point */

	cyl.r = (float)sqrt(cart.x*cart.x+cart.y*cart.y);
		if ( cart.x == 0.0f ) cart.x = .0001f;
	cyl.theta = (float)atan2(cart.y,cart.x);
	cyl.z = cart.z;
	return(cyl);
}

// ------------------------ cyl_to_cart.cpp ---------------------------
/*!
 * Point cyl_to_cart(Cylind cyl)
 *
 * cyl - cylindrical coordinates of point
 *
 * This routine converts a cylindrical point to xyz coordinates.
 */
Point cyl_to_cart(Cylind cyl)
{

	Point   cart;

	cart.x = cyl.r*(float)cos(cyl.theta);
	cart.y = cyl.r*(float)sin(cyl.theta);
	cart.z = cyl.z;
	return(cart);

}


// ------------------------ cart_to_cyl_deg.cpp ---------------------------
/*!
 * Cylind cart_to_cyl_deg(Point cart)
 *
 * cart - x y z cartisian coordinates
 *
 * This routine converts an xyz point to cylindrical coordinates with
 * theta in degrees.
 */
Cylind cart_to_cyl_deg(Point cart)
{

	Cylind  cyl;            /* converted point */

	if ( cart.x == ERRVAL ) // check for invalid point
	{
		cyl.r = cyl.theta = ERRVAL;
		cyl.z = cart.z;
	}
	else
	{
		cyl.r = (float)sqrt(cart.x*cart.x+cart.y*cart.y);
		if ( cart.x == 0.0f ) cart.x = .0001f;
		cyl.theta = RTOD * (float)atan2(cart.y,cart.x);
		if ( cyl.theta < 0.0f ) cyl.theta += 360.0f;
		if ( cyl.theta >= 360.0f ) cyl.theta =- 360.0f;
		cyl.z = cart.z;
	}

	return(cyl);
}

// ------------------------ cyl_to_cart_deg.cpp ---------------------------
/*!
 * Point cyl_to_cart_deg(Cylind cyl)
 *
 * cyl - cylindrical coordinates of point
 *
 * This routine converts a cylindrical point to xyz coordinates.  The
 * theta values is in degrees.
 */
Point cyl_to_cart_deg(Cylind cyl)
{

	Point   cart;

	// check for error value
	if ( cyl.r == ERRVAL )
	{
		cart.x = cart.y = ERRVAL;
		cart.z = cyl.z;
	}
	else
	{
		cart.x = cyl.r*(float)cos(cyl.theta*DTOR);
		cart.y = cyl.r*(float)sin(cyl.theta*DTOR);
		cart.z = cyl.z;
	}
	return(cart);

}
#endif

//////////////////////////////////////////////////////////////////////////////
////////////////////////// function cart_to_sphere ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// DSphere cart_to_sphere(DPoint cart)
//
// cart - x, y, z cartesian point
// Returns - a spherical point
//
// This routine converts a cartesian point to a point in spherical
// coordinates.  Angles are in radians.
*/
DSphere  cart_to_sphere(DPoint cart)
{
	DSphere  sphere;

	sphere.r = sqrt( cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
	if ( cart.x == 0.0f ) cart.x = 0.0000001;	/* divide by zero fudge */
	sphere.theta = atan2(cart.y,cart.x);
	sphere.phi = acos(cart.z/sphere.r);
	return(sphere);
}

/*!
// Sphere cart_to_sphere(Point cart)
//
// cart - x, y, z cartesian point
// Returns - a spherical point
//
// This routine converts a cartesian point to a point in spherical
// coordinates.  Angles are in radians.
*/
Sphere  cart_to_sphere(Point cart)
{
	Sphere  sphere;

	sphere.r = (float)sqrt( cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
	if ( cart.x == 0.0f ) cart.x = 0.0001f;	/* divide by zero fudge */
	sphere.theta = (float)atan2(cart.y,cart.x);
	sphere.phi = (float)acos(cart.z/sphere.r);
	return(sphere);
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////// function sphere_to_cart ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Point sphere_to_cart(Point cart)
//
// cart - a point in spherical points
// Returns - a cartesian point
//
// This routine converts a spherical point to a point in cartesian
// coordinates.  Angles are in radians.
*/
DPoint sphere_to_cart(DSphere sphere)
{
	DPoint   cart;

	cart.x = sphere.r * sin( sphere.phi ) * cos( sphere.theta );
	cart.y = sphere.r * sin( sphere.phi ) * sin( sphere.theta );
	cart.z = sphere.r * cos( sphere.phi );
	return(cart);
}

/*!
// Point sphere_to_cart(Point cart)
//
// cart - a point in spherical points
// Returns - a cartesian point
//
// This routine converts a spherical point to a point in cartesian
// coordinates.  Angles are in radians.
*/
Point sphere_to_cart(Sphere sphere)
{
	Point   cart;

	cart.x = sphere.r * (float)sin( sphere.phi ) * (float)cos( sphere.theta );
	cart.y = sphere.r * (float)sin( sphere.phi ) * (float)sin( sphere.theta );
	cart.z = sphere.r * (float)cos( sphere.phi );
	return(cart);
}


// ------------------------ dist.cpp ---------------------------
/*!
 * float dist(Point p1,Point p2)
 *
 * p1,p2 - points
 * return distance between points
 *
 * this routine finds the difference between two points.
 */
float dist(Point p1,Point p2)
{
	return (float)sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) +
		(p1.z-p2.z)*(p1.z-p2.z) );
}

// ------------------------ ddist.cpp ---------------------------
/*!
 * float ddist(DPoint p1,DPoint p2)
 *
 * p1,p2 - points
 * return distance between points
 *
 * this routine finds the difference between two points.
 */
double ddist(DPoint p1,DPoint p2)
{
	return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) +
		(p1.z-p2.z)*(p1.z-p2.z) );
}

// ------------------------ distsqu.cpp ---------------------------
/*!
 * float distsqu(Point p1,Point p2)
 *
 * p1,p2 - points
 * return distance between points
 *
 * this routine finds the square of the distance between two points.
 */
float distsqu(Point p1,Point p2)
{

	return (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) +
		(p1.z-p2.z)*(p1.z-p2.z);

}

// ------------------------ cyldist.cpp ---------------------------
/*!
 * float cyldist(Cylind p1,Cylind p2)
 *
 * p1,p2 - points in cylindrical coordinates
 *
 * this routine finds the difference between two Cylind points.
 */
float cyldist(Cylind p1,Cylind p2)
{
	float x1,x2,y1,y2;
	float d;

	x1 = p1.r * (float)cos( p1.theta*DTOR );
	x2 = p2.r * (float)cos( p2.theta*DTOR );
	y1 = p1.r * (float)sin( p1.theta*DTOR );
	y2 = p2.r * (float)sin( p2.theta*DTOR );

	d = (float)sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (p1.z-p2.z)*(p1.z-p2.z) );

	return d;
}

// ------------------------ cyldist2.cpp ---------------------------
/*!
 * float cyldist2(Cylind p1,Cylind p2, float *dval)
 *
 * p1,p2 - points in cylindrical coordinates
 * dval - pointer to distance between points
 *
 * this routine finds the difference between two Cylind points.
 */
float cyldist2(Cylind p1,Cylind p2,float *dval)
{
	float x1,x2,y1,y2;
	float d;

	x1 = p1.r * (float)cos( p1.theta*DTOR );
	x2 = p2.r * (float)cos( p2.theta*DTOR );
	y1 = p1.r * (float)sin( p1.theta*DTOR );
	y2 = p2.r * (float)sin( p2.theta*DTOR );

	d = (float)sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (p1.z-p2.z)*(p1.z-p2.z) );
	*dval = d;
	return d;
}

// ------------------------ cyldistsqu.cpp ---------------------------
/*!
 * float cyldistsqu(Cylind p1,Cylind p2)
 *
 * p1,p2 - points in cylindrical coordinates
 *
 * this routine returns the square of the distance between two
 * cylindrical points.
 */
float cyldistsqu(Cylind p1,Cylind p2)
{
	float x1,x2,y1,y2;
	float d;

	x1 = p1.r * (float)cos( p1.theta*DTOR );
	x2 = p2.r * (float)cos( p2.theta*DTOR );
	y1 = p1.r * (float)sin( p1.theta*DTOR );
	y2 = p2.r * (float)sin( p2.theta*DTOR );

	d = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (p1.z-p2.z)*(p1.z-p2.z);

	return d;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function theta_sort /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void theta_sort(Cylind a[],int n)
//
// a - array of cylindrical points to sort\n
// n - number of cylindrical points in array\n
// Returns - nothing
//
// This routine sorts an array of cylindrical points in ascending
// theta order.
 */
void theta_sort(Cylind a[],int n)
{
	int i,j;
	Cylind temp;

	for ( i=0 ; i<n-1 ; ++i )
	{
		for ( j=i+1 ; j<n ; ++j )
		{
			if ( a[i].theta > a[j].theta )
			{
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function z_sort /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void z_sort(Cylind a[],int n)
//
// a - array of cylindrical points to sort\n
// n - number of cylindrical points in array\n
// Returns - nothing\n
//
// This routine sorts an array of cylindrical points in ascending
// z order.
 */
void z_sort(Cylind a[],int n)
{
	int i,j;
	Cylind temp;

	for ( i=0 ; i<n-1 ; ++i )
	{
		for ( j=i+1 ; j<n ; ++j )
		{
			if ( a[i].z > a[j].z )
			{
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
//////////////////////////// function z_sort /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void z_sort(Point a[],int n)
//
// a - array of  points to sort\n
// n - number of  points in array\n
// Returns - nothing\n
//
// This routine sorts an array of  points in ascending
// z order.
 */
void z_sort(Point a[],int n)
{
	int i,j;
	Point temp;

	for ( i=0 ; i<n-1 ; ++i )
	{
		for ( j=i+1 ; j<n ; ++j )
		{
			if ( a[i].z > a[j].z )
			{
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function y_sort /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void y_sort(Point a[],int n)
//
// a - array of points to sort\n
// n - number of points in array\n
// Returns - nothing\n
//
// This routine sorts an array of  points in ascending
// y order.
 */
void y_sort(Point a[],int n)
{
	int i,j;
	Point temp;

	for ( i=0 ; i<n-1 ; ++i )
	{
		for ( j=i+1 ; j<n ; ++j )
		{
			if ( a[i].y > a[j].y )
			{
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function int_sort /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void int_sort(int a[],int n)
//
// a - array of ints to sort\n
// n - number of ints in array\n
// Returns - nothing\n
//
// This routine sorts an array of ints in ascending
// theta order.
 */
void int_sort(int a[],int n)
{
	int i,j,temp;

	for ( i=0 ; i<n-1 ; ++i )
	{
		for ( j=i+1 ; j<n ; ++j )
		{
			if ( a[i] > a[j] )
			{
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function float_sort /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void float_sort(float a[],int n)
//
// a - array of floats to sort\n
// n - number of floats in array\n
// Returns - nothing\n
//
// This routine sorts an array of floats in ascending
// theta order using qsort.
 */
void float_sort(float a[],int n)
{
	qsort (a, n, sizeof(float), fcompare);
}

///< comparison routine for float_sort
int fcompare (const void * _a, const void * _b)
{
      // you've got to explicitly cast to the correct type
	const float* a = (const float*) _a;
	const float* b = (const float*) _b;

	if( *a > *b)
		return 1;		// first item is bigger than the second one -> return 1
	else if( *a == *b)
		return  0;		// equality -> return 0
	else
		return -1;		// second item is bigger than the first one -> return -1

}


//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function cross_product ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Vector cross_product(Point pnt[])
//
// pnt - array of three points that define two vectors\n
// Returns - the vector cross product\n
//
// Three points that define two vectors are passed to the routine.  The
// cross product of the two vectors is returned as a vector.
 */
Vector	cross_product(Point pnt[]) {
	Vector	cross;
	float	x,y,z;
	float	ax,ay,az,bx,by,bz;
	float	m,mm;


	ax = pnt[1].x-pnt[0].x;
	ay = pnt[1].y-pnt[0].y;
	az = pnt[1].z-pnt[0].z;

	bx = pnt[1].x-pnt[2].x;
	by = pnt[1].y-pnt[2].y;
	bz = pnt[1].z-pnt[2].z;

	x = ay*bz-az*by;
	y = az*bx-ax*bz;
	z = ax*by-ay*bx;

	m = (float)sqrt(x*x+y*y+z*z);

	if ( m == 0.0f ) {
		cross.x = 0.0f;
		cross.y = 0.0f;
		cross.z = -1.0f;
	} else {
		mm = 1.0f/m;
		cross.x = x*mm;
		cross.y = y*mm;
		cross.z = z*mm;
	}
	return cross;
}

/*!
// DPoint cross_product(DPoint pnt[])
//
// pnt - array of three points that define two vectors\n
// Returns - the vector cross product\n
//
// Three points that define two vectors are passed to the routine.  The
// cross product of the two vectors is returned as a DPoint.
 */
DPoint	cross_product(DPoint pnt[]) {
	DPoint	cross;
	double	x,y,z;
	double	ax,ay,az,bx,by,bz;
	double	m,mm;


	ax = pnt[1].x-pnt[0].x;
	ay = pnt[1].y-pnt[0].y;
	az = pnt[1].z-pnt[0].z;

	bx = pnt[1].x-pnt[2].x;
	by = pnt[1].y-pnt[2].y;
	bz = pnt[1].z-pnt[2].z;

	x = ay*bz-az*by;
	y = az*bx-ax*bz;
	z = ax*by-ay*bx;

	m = sqrt(x*x+y*y+z*z);

	if ( m == 0.0 ) {
		cross.x = 0.0;
		cross.y = 0.0;
		cross.z = -1.0;
	} else {
		mm = 1.0/m;
		cross.x = x*mm;
		cross.y = y*mm;
		cross.z = z*mm;
	}
	return cross;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function dot_product ////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// float dot_product(Vector *v1, Vector *v2)
//
// v1, v2 - pointers to vectors to take dot product of\n
// Returns - dot product\n
//
// This routine returns the product of two vectors (V1,V2) in 3D space.
 */
float dot_product(Vector *v1, Vector *v2) {
	float dot;

	dot = v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
	return dot;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function normal_pnt /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Vector normal_pnt(Point p1, Point p2, Point p3)
//
// p1,p2,p3 - points - p2 is connecting point
// Returns - the normal vector
//
// Returns a normal vector created from the 3 connected points passed
// to the routine.
 */
Vector normal_pnt (Point p1, Point p2,Point p3)
{

	Point	ppnts[3];
	Vector	n1;
	Vector	norpnt;
	float	normag;

	ppnts[0] = p1;
	ppnts[1] = p2;
	ppnts[2] = p3;
	n1 = cross_product (ppnts);
	normag = 1.0f/(float)sqrt (n1.x*n1.x+n1.y*n1.y+n1.z*n1.z);
	norpnt.x = n1.x*normag;
	norpnt.y = n1.y*normag;
	norpnt.z = n1.z*normag;
	return norpnt;
}

/*!
// DPoint normal_pnt(DPoint p1, DPoint p2, DPoint p3)
//
// p1,p2,p3 - points - p2 is connecting point
// Returns - the normal vector
//
// Returns a DPoint normal vector created from the 3 connected points passed
// to the routine.
 */
DPoint normal_pnt (DPoint p1, DPoint p2,DPoint p3)
{

	DPoint	ppnts[3];
	DPoint	n1;
	DPoint	norpnt;
	double	normag;

	ppnts[0] = p1;
	ppnts[1] = p2;
	ppnts[2] = p3;
	n1 = cross_product (ppnts);
	normag = 1.0f/sqrt (n1.x*n1.x+n1.y*n1.y+n1.z*n1.z);
	norpnt.x = n1.x*normag;
	norpnt.y = n1.y*normag;
	norpnt.z = n1.z*normag;
	return norpnt;
}

// ------------------------ sutil_error.cpp ---------------------------
/*!
 * void sutil_error(char *msg)
 *
 * msg - error message
 *
 * This error function prints the error to stderr and exits
 */
void sutil_error(char *msg)
{
	fprintf(stderr,"sutil error - %s\n",msg);
	exit(1);
}

/* ------------------------- file_exist.c ---------------------------- */
/*!
 * int file_exist(char *filename)
 *
 * filename - name of file to check\n
 * returns 1 if file exists else 0\n
 *
 * checks to see if a file exists.  If it does then 1 is returned else
 * 0 is returned.
 */
int file_exist(char *filename)
{
	FILE	*infile;

	if ( (infile=fopen(filename,"r")) == NULL ) return(0);
	fclose(infile);
	return(1);
}

/* ------------------------------- initrotate -------------------------------------- */
/*!
Rotate::InitRotate

Sets up rotation matrix for Rotate

Computation of the rotation matrix

            | r11  r12  r13  0 |
       R =  | r21  r22  r23  0 |
            | r31  r32  r33  0 |
            | r41  r42  r43  1 |

   to be used as [x1  y1  z1  1] = [x  y  z  1] R,
   see function 'rotate'.
   Point (x1, y1, z1) is the image of (x, y, z).
   The rotation takes place about the axis
     (px, py, pz)+lambda(vx, vy, vz)
   and through the angle alpha.
 */
void Rotate::InitRotate(float px, float py, float pz, float vx, float vy,
		float vz, float angle)
{

	float rho, theta, cal, sal, cph, sph, cth, sth,cph2, sph2,
		cth2, sth2, pi, cal1;

	cal = (float)cos(angle); sal = (float)sin(angle);
	cal1 = 1.0f-cal;
	rho = (float)sqrt(vx*vx+vy*vy+vz*vz);
	pi = 4.0f * (float)atan(1.0);
	if (rho == 0.0) {theta=0.0; cph=1.0; sph=0.0;} else
	{
		if (vx == 0.0)
			theta = (vy >= 0.0 ? 0.5f*pi : 1.5f*pi);
		else
		{
			theta = (float)atan(vy/vx);
			if (vx < 0) theta += pi;
		}
		cph = vz/rho; sph = (float)sqrt(1.0f - cph*cph);
		/* cph = cos(phi), sph = sin(phi)  */
	}
	cth = (float)cos(theta); sth = (float)sin(theta);
	cph2 = cph*cph; sph2 = 1.0f - cph2;
	cth2 = cth*cth; sth2 = 1.0f - cth2;
	r11 = (cal*cph2+sph2)*cth2+cal*sth2;
	r12 = sal*cph+cal1*sph2*cth*sth;
	r13 = sph*(cph*cth*cal1-sal*sth);
	r21 = sph2*cth*sth*cal1-sal*cph;
	r22 = sth2*(cal*cph2+sph2)+cal*cth2;
	r23 = sph*(cph*sth*cal1+sal*cth);
	r31 = sph*(cph*cth*cal1+sal*sth);
	r32 = sph*(cph*sth*cal1-sal*cth);
	r33 = cal*sph2+cph2;
	r41 = px-px*r11-py*r21-pz*r31;
	r42 = py-px*r12-py*r22-pz*r32;
	r43 = pz-px*r13-py*r23-pz*r33;
}

/*!
Rotate::RotateIt
Rotates passed point.
 */
void Rotate::RotateIt(float *x, float *y, float *z)
{
	*x = (*x)*r11+(*y)*r21+(*z)*r31+r41;
	*y = (*x)*r12+(*y)*r22+(*z)*r32+r42;
	*z = (*x)*r13+(*y)*r23+(*z)*r33+r43;
}

/*!
Rotate::RotateIt
Rotates passed point.

returns rotated point
 */
Point Rotate::RotateIt(Point pnt)
{
	Point p;

	p.x = pnt.x*r11+pnt.y*r21+pnt.z*r31+r41;
	p.y = pnt.x*r12+pnt.y*r22+pnt.z*r32+r42;
	p.z = pnt.x*r13+pnt.y*r23+pnt.z*r33+r43;

	return p;
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////////// function catrom ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Point* catrom_curve (Point *pin,int nctrl, int step_size, float tension,
//      int *noutpnts)
// pin - array of input control points\n
// nctrl - number of control points\n
// step_size - sets num points to interpolate - see table below\n
// tension - curve tension 0-1, 0 is tight to 1 which is loose\n
// noutpnts - pointer to number of output points\n
// returns - pointer to array of interpolated points or NULL on failure\n
//
// This function is derived from a TurboGeometry routine.
//
// Catmull-Rom curves interpolate through control points. To connect
// curve to the end points, specify the end points twice.
//
// table of step_size values whick indicate number of points between each
// control point.
//  0 - 100
//  1 - 50
//  2 - 25
//  3 - 20
//  4 - 10
//  5 - 5
//  6 - 2
*/
Point* catrom_curve (Point *pin,int nctrl, int step_size, float tension,
int *noutpnts) {
	int 	i,m3;
	float 	t,t3,t2,tincr;
	float 	b01,b02,b10,b11,b12;
	float 	f1,f2,f3,f4;
	Point 	p1,p2,p3,p4;
	Point	*pout;		// pointer to output points
	int		nout;		// number of output points

	if (nctrl < 1) return NULL;

	// need at least 4 points if not 4 then just copy input points

	if (nctrl < 4) {
		pout = new Point[nctrl];
		for (i=0 ; i<nctrl ; ++i) pout[i] = pin[i];
		*noutpnts = nctrl;
		return pout;
	}

	// set the number of points between control points

	switch (step_size) {
		case 0:  tincr = 0.01f; break; // 100 points
		case 1:  tincr = 0.02f; break; //  50 points
		case 2:  tincr = 0.04f; break; //  25 points
		case 3:  tincr = 0.05f; break; //  20 points
		case 4:  tincr = 0.10f; break; //  10 points
		case 5:  tincr = 0.20f; break; //   5 points
		case 6:  tincr = 0.50f; break; //   2 points
		default: tincr = 0.05f;        //  20 points
	}

	// allocate memory for output points

	nout = (int)(( (nctrl-3)* (1.0f/tincr)+2.0f)+0.5f );
	pout = new Point[nout];
	if (pout == NULL) return 0;

	b01 = 2.0f - tension;
	b02 = tension - 2.0f;
	b10 = 2.0f * tension;
	b11 = tension - 3.0f;
	b12 = 3.0f - b10;
	m3 = (int) (nctrl - 3);
	nout = 0;
	for (i = 0; i < m3; i++) {
		t = 0.0f;
		p1 = pin[i];
		p2 = pin[i+1];
		p3 = pin[i+2];
		p4 = pin[i+3];
		do {
			t2 = t * t;
			t3 = t2 * t;
			f1 = (t3 * (-tension)) + (t2 * b10) + (t * (-tension));
			f2 = (t3 * b01) + (t2 * b11) + 1.0f;
			f3 = (t3 * b02) + (t2 * b12) + (t * tension);
			f4 = (t3 * tension) + (t2 * (-tension));
			pout[nout].x = (f1*p1.x) + (f2*p2.x) + (f3*p3.x) + (f4*p4.x);
			pout[nout].y = (f1*p1.y) + (f2*p2.y) + (f3*p3.y) + (f4*p4.y);
			pout[nout].z = (f1*p1.z) + (f2*p2.z) + (f3*p3.z) + (f4*p4.z);
			nout++;
			t += tincr;
		} while (t<0.999999);
	}

	// last point

	pout[nout] = p3;
	nout++;
	*noutpnts = nout;
	return pout;
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////// function line_thru_plane //////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Point line_thru_plane(float a,float b,float c,float d,Point p0,Point p1)
//
// a,b,c,d - coefficients of plane equation ax + by + cz = d\n
// p0,p1 - end points of vector intersecting plane\n
// Returns - the point of intersection between line and plane\n
//
// Finds the point of intersetion between and line and a plane.  If the
// line segment described by p0 and p1 does not go through the plane
// then an extension of the line segment is used.
 */
Point line_thru_plane (float a, float b, float c, float d, Point p0, Point p1) {
	float 	t,u;		// intermediate solution
	float   dx,dy,dz;	// xyz differences between end points
	Point	pout;		// point of intersection

	// see if p0=p1 - if so return p0

	if (p0.x==p1.x && p0.y==p1.y && p0.z==p1.z)
		return p0;

	dx = p1.x-p0.x;
	dy = p1.y-p0.y;
	dz = p1.z-p0.z;

	u = a*dx + b*dy + c*dz;
	if ( !u ) u = .000001f;
	t = ( d -a*p0.x -b*p0.y -c*p0.z ) / u;

	pout.x = dx*t + p0.x;
	pout.y = dy*t + p0.y;
	pout.z = dz*t + p0.z;

	return pout;
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////// function line_thru_triangle ///////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int line_thru_triangle(Point t0,Point t1,Point t2,Point p0,Point p1,
//		Point *pout)
// t0,t1,t2 - triangle points\n
// p0,p1 - line points\n
// pout - pointer to output point\n
// returns 1 one if line intersects triangle\n
//
// Finds the point of intersetion between and line and a triangle.
 */
int line_thru_triangle(Point t0,Point t1,Point t2,Point p0,Point p1,Point *pout)
{
	Point	pint;		// point of intersection
	float	a,b,c,d;	// coefficients of plane equation

	// get equation of plane of triangle
	equ_plane(t0,t1,t2,&a,&b,&c,&d);

	// see if there is an intersection of line with plane
	if ( line_intersect_plane(a,b,c,d,p0,p1,&pint) );
	{
		// see if intersection point lies within triangle
		if ( inside_triangle(t0,t1,t2,pint) )
		{
			*pout = pint;
			return 1;
		}
	}

	// assign pint to pout anyway just so there is a point that might be close
	*pout = pint;
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
////////////////////// function line_intersect_triangle //////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int line_intersect_triangle(Point t0,Point t1,Point t2,Point p0,Point p1,
//		Point *pout)
// t0,t1,t2 - triangle vertex points\n
// p0,p1 - line end points\n
// pout - pointer to output point\n
// returns 1 one if line intersects triangle\n
//
// Taken from Ben Trumbore's code from the Journal of Graphics Tools (JGT)
// This version Finds the point of intersection between and line segment and
// a triangle rather than the intersection of a ray with the triangle.
 */
int line_intersect_triangle(Point t0,Point t1,Point t2,Point p0,Point p1,Point *pout)
{

	// vector providing direction
	Vector d = p1-p0;

	//  vectors for two edges sharing vertex 0
	Vector e1,e2;
	e1 = t1-t0;
	e2 = t2-t0;

	Vector h = d.CrossProduct(e2);
	float a = e1.Dot(h);

	if (a > -0.00001 && a < 0.00001)
		return(0);

	float f = 1/a;

	//vector(s,p,v0);
	Vector s = p0-t0;

	float u = f * (s.Dot(h));

	if (u < 0.0 || u > 1.0)
		return(0);

	//crossProduct(q,s,e1);
	Vector q = s.CrossProduct(e1);

	float v = f * d.Dot(q);

	if (v < 0.0 || u + v > 1.0)
		return(false);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	//Point t;
	float t = f * e2.Dot(q);
	*pout = p0 + d*t;

	if (t < 0.00001) // ray intersection not on line
		return(false);

	// the triangle is on the ray, now check to see if is
	// on the line segment
	float dseg = distsqu(p0,p1);
	float dparts = distsqu(p0,*pout) + distsqu(p1,*pout);
	float w = dparts-dseg;
	if ( w < 0.00001 )
		return true;
	else // triangle intersection not on line segment
		 return (false);

}

//////////////////////////////////////////////////////////////////////////////
/////////////////////// function line_intersect_plane ////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int line_intersect_plane(float a,float b,float c,float d,Point p0,
//     Point p1,Point *pout)
//
// a,b,c,d - coefficients of plane equation ax + by + cz = d\n
// p0,p1 - end points of vector intersecting plane\n
// pout - the intersection of the line segment and the plane.\n
// Returns - 1 if line intersects plane, 0 if not\n
//
// Finds the point of intersetion between and line and a plane.  If the
// line between p1 and p2 intersects the plane then 1 is returned.  If
// the line segment does not intersect the plane then 0 is returned.
// Pout is the point of intersection.
 */
int line_intersect_plane(float a,float b,float c,float d,Point p0,
	Point p1,Point *pout) {
	float 	t,u;				/* intermediate solution */
	float   dx,dy,dz;			/* xyz differences between end points */
	float	dist,dist0,dist1;	/* distances squared between points */
	Point	p;

	// see if p0=p1 - if so return p0

	if ( (fabs(p0.x-p1.x)+fabs(p0.y-p1.y)+fabs(p0.z-p1.z)) < 0.01 ) return(0);

	dx = p1.x-p0.x;
	dy = p1.y-p0.y;
	dz = p1.z-p0.z;

	u = a*dx + b*dy + c*dz;
	if ( !u ) u = .000001f;
	t = ( d -a*p0.x -b*p0.y -c*p0.z ) / u;

	p.x = dx*t + p0.x;
	p.y = dy*t + p0.y;
	p.z = dz*t + p0.z;

	dist = dist_squ (p0, p1);
	dist0 = dist_squ (p0,p);
	dist1 = dist_squ (p1,p);

	// if line segment intersects plane return 1

	if ( dist0 <= dist && dist1 <= dist ) {
		*pout = p;
		return 1;
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
///////////////////////////// function equ_plane /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void equ_plane(Point p1,Point p2,Point p3,float *a,float *b,float *c,
//      float *d)
//
// p1,p2,p3 - three points through which the plane must pass\n
// *a,*b,*c,*d - pointers to coefficients of the plane equation\n
// Returns - nothing\n
//
// Calculates the equation of a plane from three points.  The coefficients
// of the equation ax + by + cz = d are passed as pointers to floats.
 */
void equ_plane(Point p1,Point p2,Point p3,float *a,float *b,float *c,float *d)
{
	*a = p1.y*p2.z +p1.z*p3.y +p2.y*p3.z -p3.y*p2.z -p3.z*p1.y -p2.y*p1.z;

	*b = -(p1.x*p2.z +p1.z*p3.x +p2.x*p3.z -p3.x*p2.z -p3.z*p1.x -p2.x*p1.z);

	*c = p1.x*p2.y +p1.y*p3.x +p2.x*p3.y -p3.x*p2.y -p3.y*p1.x -p2.x*p1.y;

	*d = p1.x*p2.y*p3.z +p1.y*p2.z*p3.x +p1.z*p2.x*p3.y -p3.x*p2.y*p1.z
	-p3.y*p2.z*p1.x -p3.z*p2.x*p1.y;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function equ_plane2 /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// void equ_plane2(Point p1,Vector n1,float *a,float *b,float *c,float *d)
//
// p1 - point on the plane\n
// n1 - normal vector to the plane\n
// *a,*b,*c,*d - pointers to coefficients of the plane equation\n
// Returns - nothing\n
//
// Calculates the equation of a plane from a point on the plane and a
// normal vector to the plane.  The coefficients of the equation
// ax + by + cz = d are passed as pointers to floats.
 */
void equ_plane2(Point p1,Vector n1,float *a,float *b,float *c,float *d)
{
	*a = n1.x;

	*b = n1.y;

	*c = n1.z;

	*d = n1.x*p1.x + n1.y*p1.y + n1.z*p1.z;
}

//////////////////////////////////////////////////////////////////////////////
///////////////////////////// function dist_squ //////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// float dist_squ(Point p1,Point p2)
//
// p1, p2 - The x,y,z data points to find the square of the distance between\n
// Returns - the square of the distance between the points\n
//
// This routine calculates the square of the distance between two
// x, y, z points
 */
float Dist_squ (Point *p1, Point *p2) {
	float	dx,dy,dz;

	dx = p1->x-p2->x;
	dy = p1->y-p2->y;
	dz = p1->z-p2->z;

	return( dx*dx + dy*dy + dz*dz );
}

// some variables for 2D in_polygon
#define MAXAOIPNTS		750		///< maximum number of points allowed in aoi

float	M2[MAXAOIPNTS]; ///< slope for in_polygon precalculation
float	B2[MAXAOIPNTS]; ///< intercept for in_polygon precalculation

//////////////////////////////////////////////////////////////////////////////
////////////////////// function initialize_in_polygon ////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int initialize_in_polygon (Point *v, int npnts)
//
// This function is used to initialize the in_polygon function.
// It precalculates some numbers used in the loop.
*/
void initialize_in_polygon (Point *v, int npnts) {
	int		i;

	// for all the edges...

	for (i=1; i<npnts; i++) {

		// check for vertical line segment; perterb

		if (v[i].x == v[i-1].x) v[i].x += 0.001f;

		// if the edge v_i-1..v_i intersects a horizontal line thru p
		// to the left of p, L++
		// find eqn. of line

		M2[i] = (v[i].y - v[i-1].y) / (v[i].x - v[i-1].x);
		B2[i] = v[i].y - (M2[i] * v[i].x);
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// function in_polygon /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int in_polygon (Point *v, int npnts, Point p)
//
// This function returns 1 if a given point is in a given 2D polygon
 */
int in_polygon (Point *v, int npnts, Point p) {
	float m1 = 0.0f;
	float b1 = p.y;
	int L = 0, i;

	// for all the edges...

	for (i=1; i<npnts; i++) {

		// (these values are precomputed)
		// if the edge v_i-1..v_i intersects a horizontal line thru p
		// to the left of p, L++
		// find eqn. of line

		float m2 = M2[i];
		float b2 = B2[i];

		// x is the x-coordinate of the intersection
		float x, min_x, max_x;
		if (m1 != m2) {
			x = (b2 - b1) / (m1 - m2);
			max_x = MAX (v[i].x, v[i-1].x);
			min_x = MIN (v[i].x, v[i-1].x);
			if (x <= max_x) if (x >= min_x) if (x < p.x) L++;
		}
	}

	// iff L is odd, p is inside the polygon

	return L & 1;
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function in_polygon2 /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int in_polygon2 (Point *v, int n, Point p, int add_pi)
//
// This function wraps some of the ugliness of using in_polygon
// inside aoi_lasso_vertices.
 */
int in_polygon2 (Point *v, int n, Point p, int add_pi) {
	Point	q;
	float	theta=0;

	// convert to cylindrical


	// !!!!!!!!!! look here - non functional !!!!!!!!!

//	theta = point_theta (p);

	// if the theta value wraps around zero, translate up

	if (theta < 0.0f) theta += 2*PI;
	if (add_pi) {
		theta += PI;
		if (theta >= 2.0f*PI) theta -= 2.0f*PI;
	}
	q.x = theta;
	q.y = p.z;
	return in_polygon (v, n, q);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////// function angle_between_planes ///////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// float angle_between_planes(float a1,float b1,float c1,float a2, float b2,
//								float c2)
//
// a1,b1,c1 - coefficients of the equation for plane 1\n
// a2,b2,c2 - coefficients of the equation for plane 2\n
// Returns - the angle between the planes in radians\n
//
// Calculates the angle between two planes. The dot product of the normals
// divided by the product of their magnitudes gives the cosine of the angle.
 */
float angle_between_planes(float a1,float b1,float c1,float a2, float b2,
						   float c2)
{
	// calculate dot product
	float dp = a1*a2 + b1*b2 + c1*c2;

	// find the magnitudes
	float m1 = (float)sqrt(a1*a1 + b1*b1 + c1*c1);
	float m2 = (float)sqrt(a2*a2 + b2*b2 + c2*c2);

	// check for divide by zero
	if ( m1 == 0.0f || m2 == 0.0f ) return 0.0f;

	// calculate the cosine of the angle
	float cosa = dp/(m1*m2);

	float a = (float)acos(cosa);

	return a;
}

//////////////////////////////////////////////////////////////////////////////
///////////////////////// function dist_point_plane //////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// float dist_point_plane(float a,float b,float c,float d,
//                        Point p0,Point *pout)
//
// a,b,c,d - coefficients of plane equation ax + by + cz = d\n
// p0 - point to check distance from\n
// pout - the point on the plane where the normal from the point intersects\n
// Returns - the normal distance from the point to the plane\n
//
// Finds the distance from a point to a plane along a normal.  The
// textbook Calculus and Analytic Geometry by Thomas was used as a
// reference.
*/
float dist_point_plane (float a,float b,float c,float d,Point p0, Point *pout)
{

	float 	t,u;		// intermediate solution
	float	delta;		// distance from point to plane
	u = a*a + b*b + c*c;
	if ( !u ) u = .000001f;
	t = ( d - a*p0.x - b*p0.y - c*p0.z ) / u;
	pout->x = a*t + p0.x;
	pout->y = b*t + p0.y;
	pout->z = c*t + p0.z;
	delta = dist (p0, *pout);
	return delta;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////// function DistPointLine ///////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// double dist_point_line(Point Pnt,Point LineStart, Point LineEnd,
//						float *Distance, Point &Pout )
//
// Pnt - input point
// LineStart - origin of line
// LineEnd - line endpoint
// Distance - distance between Pnt and line
// Pout - point on line nearest input point
//
// Calculates the perpendicular distance between a point and a line segment.
// It also calculates the point of on the line nearest the input point.
// If the point is not located perpendicular to the line, 0 is returned
// else 1 is returned on success.
 */

int dist_point_line(Point Pnt, Point LineStart, Point LineEnd, float *Distance, Point *Pout)
{
	float LineMag;
	float U;
	LineMag = dist( LineEnd, LineStart );
	U = ( ( ( Pnt.x - LineStart.x ) * ( LineEnd.x - LineStart.x ) ) +
			( ( Pnt.y - LineStart.y ) * ( LineEnd.y - LineStart.y ) ) +
			( ( Pnt.z - LineStart.z) * ( LineEnd.z - LineStart.z ) ) ) /
			( LineMag * LineMag );

	if( U < 0.0f || U > 1.0f )
		return 0; // closest point does not fall within the line segment

	Pout->x = LineStart.x + U * ( LineEnd.x - LineStart.x );
	Pout->y = LineStart.y + U * ( LineEnd.y - LineStart.y );
	Pout->z = LineStart.z + U * ( LineEnd.z - LineStart.z );
	*Distance = dist( Pnt, *Pout );
	return 1;
}

/*!
Distance between a point and a line
 */
double distpntlinespace(Point pnt,double xo,double yo,double zo,
								double dx,double dy,double dz)
{
	double dx1,dy1,dz1,dx2,dy2,dz2;

	dx1= pnt.x - xo;
	dy1 = pnt.y - yo;
	dz1 = pnt.z - zo;
	dx2 = dx * dy1 - dy * dx1;
     	dy2 = dx * dz1 - dz * dx1;
	dz2 = dy * dz1 - dz * dy1;
     	dx1 = dy * dx2 + dz * dy2;
   	dy1 = dz * dz2 - dx * dx2;
     	dz1 = -dx * dy2 - dy * dz2;
     return(sqrt((dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1)));
}


//////////////////////////////////////////////////////////////////////////////
/////////////////////////// Spline routines /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// function cubic spline routines
//
// This spline function is from NR in C. It has been converted from
// origin 1 to origin 0.
 */
void splie2(float x1a[],float x2a[],float **ya,int m,int n,float **y2a)
{
	int j;

	for (j=0;j<m;j++)
		spline(x2a,ya[j],n,1.0e30f,1.0e30f,y2a[j]);
}

/*!
// function spline routines
//
// This spline function is from NR in C. It has been converted from
// origin 1 to origin 0.
 */
void splin2(float x1a[],float x2a[],float **ya,float **y2a,int m,int n,
			   float x1,float x2,float *y)
{
	int j;
	float ytmp[100],yytmp[100];

	for (j=0;j<m;j++)
		splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	spline(x1a,yytmp,m,1.0e30f,1.0e30f,ytmp);
	splint(x1a,yytmp,ytmp,m,x1,y);
}

/*!
// function spline routines
//
// This spline function is from NR in C. It has been converted from
// origin 1 to origin 0.
 */
void spline(float x[],float y[],int n,float yp1,float ypn,float y2[])
{
	int i,k;
	float p,qn,sig,un,u[100];

	if (yp1 > 0.99e30f)
		y2[0]=u[0]=0.0f;
	else {
		y2[0] = -0.5f;
		u[0]=(3.0f/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<=n-2;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0f;
		y2[i]=(sig-1.0f)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0f*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30f)
		qn=un=0.0f;
	else {
		qn=0.5f;
		un=(3.0f/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0f);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}

/*!
// function spline routines
//
// This spline function is from NR in C. It has been converted from
// origin 1 to origin 0.
*/
void spline_d(double x[],double y[],int n,double yp1,double ypn,double y2[])
{
	int i,k;
	double p,qn,sig,un,u[100];

	if (yp1 > 0.99e30f)
		y2[0]=u[0]=0.0f;
	else {
		y2[0] = -0.5f;
		u[0]=(3.0f/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<=n-2;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0f;
		y2[i]=(sig-1.0f)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0f*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30f)
		qn=un=0.0f;
	else {
		qn=0.5f;
		un=(3.0f/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0f);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}

/*!
// function spline routines
//
// This spline function is from NR in C. It has been converted from
// origin 1 to origin 0.
 */
void splint(float xa[],float ya[],float y2a[],int n,float x,float *y)
{
	int klo,khi,k;
	float h,b,a;

	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	//if (h == 0.0f) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0f;
}

/*!
// function spline routines
//
// This spline function is from NR in C. It has been converted from
// origin 1 to origin 0.
*/
void splint_d(double xa[],double ya[],double y2a[],int n,double x,double *y)
{
	int klo,khi,k;
	double h,b,a;

	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	//if (h == 0.0f) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0f;
}


//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function circle_3pnts /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int circle_3pnts (Point *p1, Point *p2, Point *p3, Point *cp)\n
// p1,p2,p3 - 3 points to determine circle from\n
// cp - pointer to center point\n
// returns radius on success, 0.0 on failure\n
//
// This function uses three points to define a circle. It is assumed that the
// circle is in the x-y plane. The z value is ignored. The radius value is
// returned on success else 0.
*/
float circle_from_3pnts(Point p1, Point p2, Point p3, Point *cp)
{
   float xdiff12, ydiff12;	// difference in x and y between pnts P1 and P2
   float xdiff13, ydiff13;	// difference in x and y between pnts P1 and P3
   float sdist12, sdist13;	// Distance squared between pnts P1 and P2/P3
   float det;				// Determinate
   float rad;				// radius of circle

   xdiff12 = p2.x - p1.x;
   ydiff12 = p2.y - p1.y;
   xdiff13 = p3.x - p1.x;
   ydiff13 = p3.y - p1.y;
   det = ydiff13 * xdiff12 - xdiff13 * ydiff12;
   if ( det == 0.0f )
   {
	   cp->x = cp->y = cp->z = ERRVAL;
	   return 0.0f;
   }
   else
   {
     det = 1.0f / (2.0f * det);
     sdist12 = (xdiff12*xdiff12)+(ydiff12*ydiff12);
     sdist13 = (xdiff13*xdiff13)+(ydiff13*ydiff13);
     cp->x = det * (sdist12 * ydiff13 - sdist13 * ydiff12);
     cp->y = det * (sdist13 * xdiff12 - sdist12 * xdiff13);
     rad = (float)sqrt((cp->x * cp->x) + (cp->y * cp->y));
     cp->x = cp->x + p1.x;
     cp->y = cp->y + p1.y;
	 cp->z = p1.z;
   }
   return rad;
}

//////////////////////////////////////////////////////////////////////////////
////////////////////// function interpolate_along_z //////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Point interpolate_along_z (Point *p, z)\n
// p - array of two points for high and low\n
// z - value to interpolate to\n
// returns interpolated point\n
//
// This function interpolates between the two points passed to it to determine
// the x and y values of a point at the passed z level. The interpolated point is
// returned. It does not check to see that the z is between the two points so
// is is possible to extrapolate
 */
Point interpolate_along_z (Point *p, float z)
{
	Point pout;

	// check for divide by zero
	pout.x = ERRVAL;
	if ( p[0].z == p[1].z ) return p[0];

	float mfac = 1.0f/(p[1].z-p[0].z);
	pout.x = p[0].x + (p[1].x-p[0].x)*(z-p[0].z)*mfac;
	pout.y = p[0].y + (p[1].y-p[0].y)*(z-p[0].z)*mfac;
	pout.z = z;

	return pout;
}

//////////////////////////////////////////////////////////////////////////////
////////////////////// function scale_array //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int scale_array(float newMin, float newMax, float *array, int l, int invert)\n
// newMin,newMax - range to scale to
// array - array of data
// l - length of array
// invert - invert - make negative
// returns 1 on success
//
// This function scales an array between min and max values. It will make a
// negative image if invert is 1
*/
int scale_array(float newMin, float newMax, float *array, int l, int invert)
{
	float maxVal = -FLT_MAX;
	float minVal = FLT_MAX;
	float v,scale;
	int i;

	i = 0;
	while( i < l )
	{
		v = array[i];
		if( v > maxVal )
			maxVal = v;
		else if( v < minVal )
			minVal = v;
		i++;
	}

	scale = ( newMax - newMin ) / ( maxVal - minVal );

	if(invert)
	{
		i = 0;
		while( i < l )
		{
			array[i] = (newMax - scale*( array[i] - minVal ));
			i++;
		}
	}
	else
	{
		i = 0;
		while( i < l )
		{
			array[i] = (newMin + scale*( array[i] - minVal ));
			i++;
		}
	}
	return 1;
}
//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function line_intersect //////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// Point line_intersect (Point p11, Point p12, Point p21, Point p22)
//
// p11,p12 - defines first line
// p21,p22 - defines second line
//
// This function returns point of intersection between two lines. This is
// a 2d routine that only uses x and y
*/
Point line_intersect (Point p11, Point p12, Point p21, Point p22)
{
	float d1x,d1y,d2x,d2y;
	float m1,m2;	// slopes
	float b1,b2;	// intercepts
	Point p;	// intersection point

	d1x = p11.x - p12.x;
	d1y = p11.y - p12.y;
	d2x = p21.x - p22.x;
	d2y = p21.y - p22.y;

	if ( d1x == 0 ) d1x = 0.00000001f;
	m1 = d1y/d1x;

	if ( d2x == 0 ) d2x = 0.00000001f;
	m2 = d2y/d2x;

	b1 = p11.y - m1* p11.x;
	b2 = p21.y - m2* p21.x;

	p.z = p11.z;	// set z to passed point

	float mm = m1-m2;
	if ( mm == 0 ) mm = 0.00000001f;
	p.x = (b2-b1)/mm;
	p.y = m1*p.x + b1;

	return p;
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function area_triangle ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// float area_triangle (float a, float b, float c)
//
// a, b, c - lengths of sides
//
// This function returns the area of a triangle whose legs are length
// a, b, and c.  If the legs are colinear or close to colinear, 0.0 is
// returned.
*/
float area_triangle (float a, float b, float c) {
	float s, t;

	s = (a + b + c) * 0.5f;
	t = s * (s - a) * (s - b) * (s - c);

	/* hmm; colinear? */

	if (t < 0.0f) return 0.0f;
	return (float)sqrt (t);
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function area_triangle ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// float area_triangle (Point p1, Point p2, Point p3)
//
// a, b, c - triangle vertices
//
// This function returns the area of a triangle whose vertices are
// a, b, and c.  If the legs are colinear or close to colinear, 0.0 is
// returned.
 */
float area_triangle (Point p1, Point p2, Point p3) {
	float s, t;
	float a,b,c;

	a = dist(p1,p2);
	b = dist(p2,p3);
	c = dist(p3,p1);

	s = (a + b + c) * 0.5f;
	t = s * (s - a) * (s - b) * (s - c);

	/* hmm; colinear? */

	if (t < 0.0f) return 0.0f;
	return (float)sqrt (t);
}

//////////////////////////////////////////////////////////////////////////////
///////////////////////// function inside_triangle ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int inside_triange (Point t1, Point t2, Point t3, Point p)
//
// t1, t2, t3 - vertices of a triangle\n
// p          - point inside or outside triangle\n
//
// This function returns 1 if the triangle defined by t1, t2, and t3
// contains the point p.  This is done by adding the areas of the three
// triangles obtained by joining the vertices with the point and comparing
// the sum the the area of the main triangle; if they are equal, the point
// is inside; otherwise, it is outside.
*/
int inside_triangle (Point t1, Point t2, Point t3, Point p) {
	float a1, a2, t1_t2, t1_t3, t2_t3, t1_p, t2_p, t3_p;

	t1_t2 = dist (t1, t2);
	t1_t3 = dist (t1, t3);
	t2_t3 = dist (t2, t3);
	t1_p = dist (t1, p);
	t2_p = dist (t2, p);
	t3_p = dist (t3, p);
	a1 = area_triangle (t1_t2, t1_t3, t2_t3); /* big triangle  */

	/* if there's no area, something fishy is going on (e.g. colinear points) */

	if (a1 == 0.0) return 0;
	a2 = area_triangle (t3_p, t2_p, t2_t3) +
	area_triangle (t1_p, t2_p, t1_t2) +
	area_triangle (t1_p, t3_p, t1_t3);

	/* if the areas are pretty close... */

	if (fabs (a1 - a2) < (0.0001*a1)) return 1;
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
///////////////////// function sphere_line_intersection //////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int sphere_line_intersection (Point center, float r, Point l1, Point l2,
//								Point *p1 ,Point *p2)
// center - center of sphere\n
// r - radius of sphere\n
// 11,12 - points on line\n
// p1, p2 - pointers to intersection points\n
// returns number of intersections intersection\n
//
// This function checks to see if a line intersects a sphere. The intersection
// points are calculated. The return value is the number if intersections where
// 1 intersection indicates the line is tangent.
*/
int sphere_line_intersection (Point center, float r, Point l1, Point l2,
							  Point *p1, Point *p2)
{
	float a, b, c, mu, intv ;

	a = (l2.x-l1.x)*(l2.x-l1.x) + (l2.y-l1.y)*(l2.y-l1.y)
		+ (l2.z-l1.z)*(l2.z-l1.z);

	b =  2* ( (l2.x - l1.x )*(l1.x - center.x)
		+ (l2.y - l1.y)*(l1.y - center.y)
		+ (l2.z - l1.z)*(l1.z - center.z) ) ;

	c = center.x*center.x + center.y*center.y + center.z*center.z
		+ l1.x*l1.x + l1.y*l1.y + l1.z*l1.z
		- 2*( center.x*l1.x + center.y*l1.y + center.z*l1.z)
		- r*r;

	intv =   b * b - 4 * a * c ;

	if ( intv < 0.0 )
	{
		// no intersection
		return 0;
	}
	if ( intv == 0.0 )
	{
		// one intersection - tangent
		mu = -b/(2*a) ;
		p1->x = l1.x + mu*(l2.x-l1.x);
		p1->y = l1.y + mu*(l2.y-l1.y);
		p1->z = l1.z + mu*(l2.z-l1.z);
		return 1;
	}
	if ( intv > 0.0 )
	{
		// two intersections

		// first intersection
		mu = (-b + (float)sqrt( b*b - 4*a*c )) / (2*a);
		p1->x = l1.x + mu*(l2.x-l1.x);
		p1->y = l1.y + mu*(l2.y-l1.y);
		p1->z = l1.z + mu*(l2.z-l1.z);
		// second intersection
		mu = (-b - (float)sqrt( b*b - 4*a*c )) / (2*a);
		p2->x = l1.x + mu*(l2.x-l1.x);
		p2->y = l1.y + mu*(l2.y-l1.y);
		p2->z = l1.z + mu*(l2.z-l1.z);

		return 2;
	}

	return 3;	// should never do this
}


//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function int_circle_line /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// function int_circle_line
// **index int_circle_line
//
// int int_circle_line (Point c, float r, Point p1, Point p2)\n
// c - center of circle\n
// r - radius of circle\n
// p1,p2 - end points of line\n
// returns 1 if there is an intersection\n
//
// This function checks to see if a 2d line intersects a 2d circle. This
// is derived from TGL subroutine code. If there is no intersection or
// the line is tangent then 0 is returned. If there are two intersections
// the 1 is returned.
*/
int int_circle_line (Point c, float r, Point p1, Point p2)
{

	float dx,dy;	// Change in and x
	float dist;		// Distance between P1 and P2*/
	float invdist;	// Inverse of Dist
	float xo,yo;	// origin of line
	float xdiff,ydiff;
	float parm1,parm2;

	// convert line to parametric equation
	xo = p1.x;
	yo = p1.y;
	dx = (p1.x - p2.x);
	dy = (p1.y - p2.y);
	dist = (dx * dx) + (dy * dy);
	if ( dist == 0.0f ) return 0; // if points the same then return 0

	// normalize
	invdist = (float)(1.0 / sqrt(dist));
	dx *= invdist;
	dy *= invdist;

	/*Get the difference between XO and CP.X to determine
	if difference is less than Radius. If not no intersection*/
	xdiff = c.x - xo;
	ydiff = c.y - yo;
	parm1 = dx * ydiff - dy * xdiff;
	parm2 = (r * r) - (parm1 * parm1);
	if ( parm2 <= 0 ) return(0); /*Line does not intersect circle*/

	parm1 = dx * xdiff + dy * ydiff;

	Point	intp1,intp2;	// intersection point
	float t1,t2;
	if ( parm2 == 0)
	{  /*One Intersection */
		intp1.x = xo + dx * parm1;
		intp1.y = yo + dy * parm1;
		if (inbetween(p1,p2,intp1)) return 1;
	}
	else
	{  /*Two Intersections*/
		parm2 = (float)sqrt(parm2);
		t1 = parm1 - parm2;
		t2 = parm1 + parm2;
		intp1.x = xo + dx * t1;
		intp1.y = yo + dy * t1;
		intp2.x = xo + dx * t2;
		intp2.y = yo + dy * t2;
		if (inbetween(p1,p2,intp1)) return 1;
		if (inbetween(p1,p2,intp2)) return 1;
 	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////////// function inbetween /////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
// int inbetween(Point lp1, Point lp2, Point p)\n
// lp1,lp2 - end points on line\n
// p - point to check\n
// returns 1 if p is between lp1 and lp2\n
//
// This function checks to see if a 2D point is between the end points of
// a line segment. If so then 1 is returned.
*/
int inbetween(Point lp1, Point lp2, Point p)

 /* This function determines if a point is between the end points of
    line segment.*/

{ /*inbetween*/

   Point tp1,tp2;	// test points

   tp1.x = ( lp1.x < lp2.x ) ? lp1.x : lp2.x;
   tp1.y = ( lp1.y < lp2.y ) ? lp1.y : lp2.y;
   tp2.x = ( lp1.x > lp2.x ) ? lp1.x : lp2.x;
   tp2.y = ( lp1.y > lp2.y ) ? lp1.y : lp2.y;

   if ( (p.x >= tp1.x) && (p.x <= tp2.x) &&
       (p.y >= tp1.y) && (p.y <= tp2.y))
     return(1);
   return(0);
} /*inbetween*/
