// ------------------------------- RicCurve.cpp --------------------------
/*! \file
Implementation file for the RicCurve Class
 */

/*! \mainpage
The RicCurve class encapsulates routines from the SISL curve and surface library.
What is included is the ability to initialize and interpolate parametric curves.
Also curves can be written out as BrainVISA/Anatomist mesh files.
*/

#include <fstream> 
#include <iostream>
#include <cstdlib>
#include "RicCurve.h"
using namespace std;

/*!
Constructor for empty curve
 */
RicCurve::RicCurve(void)
{
	order = 4;
	kind = 1;
	ncntl = npnts = 0;
	pnts = NULL;
	cpnts = NULL;
	knots = NULL;
	scurve = NULL;
}

/*!
Constructor allocating control points

@param nc - number of control points
 */
RicCurve::RicCurve(int nc)
{
	order = 4;
	kind = 1;
	ncntl = nc;
	cpnts = new DPoint[ncntl];
	npnts = 0;
	pnts = NULL;
	nknots = ncntl+order;
	knots = new double[nknots];
	scurve = NULL;
}

/*!
Constructor allocating control points, selecting kind of curve and picking order

@param nc - number of control points
@param ord - order of curve
@param knd - kind of curve (1-4) 1-bspline 2-nurbs 3-bezier 4-rat bezier
 */
RicCurve::RicCurve(int nc, int ord, int knd)
{
	order = ord;
	kind = knd;
	ncntl = nc;
	cpnts = new DPoint[ncntl];
	npnts = 0;
	pnts = NULL;
	nknots = ncntl+order;
	knots = new double[nknots];
	scurve = NULL;
}

/*!
Constructor using passed float control points to build curve structure

@param nc - number of control points
@param cp - control points array
 */
RicCurve::RicCurve(int nc, Point *cp)
{
	order = 4;
	kind = 1;
	ncntl = nc;
	cpnts = new DPoint[ncntl];
	int i;
	for ( i=0 ; i<nc ; ++i )
	{
		cpnts[i].x = cp[i].x;
		cpnts[i].y = cp[i].y;
		cpnts[i].z = cp[i].z;
	}
	npnts = 0;
	pnts = NULL;
	
	// we have enough info to create a curve
	InitCurve();	
}

/*!
Constructor using passed double control points to build curve structure

 * @param nc - number of control points
 * @param cp - control points array
 */
RicCurve::RicCurve(int nc, DPoint *cp)
{
	order = 4;
	kind = 1;
	ncntl = nc;
	cpnts = new DPoint[ncntl];
	int i;
	for ( i=0 ; i<nc ; ++i )
	{
		cpnts[i].x = cp[i].x;
		cpnts[i].y = cp[i].y;
		cpnts[i].z = cp[i].z;
	}
	npnts = 0;
	pnts = NULL;
	
	// we have enough info to create a curve
	InitCurve();	
}

/*!
Destructor - clean up memory
 */
RicCurve::~RicCurve()
{
	if ( cpnts ) delete[] cpnts;
	if ( pnts ) delete [] pnts;
	if ( knots ) delete [] knots;
	if ( scurve ) freeCurve(scurve);
}


/*!
Initialize Sisl curve structure. It assumes that the control points have
been initialized already. The knots array is filled in and the Sisl
newCurve function is called to create a new curve

@return 1 if successful, 0 on failure
 */
int RicCurve::InitCurve()
{
	int i;
	
	// bail if no control points
	if ( !ncntl ) return 0;
	
	// create knot vector
	nknots = ncntl+order;
	knots = new double[nknots];
	double kval = 0;
	for ( i=0 ; i<order ; ++i ) knots[i] = kval;
	for ( i=order ; i<(nknots-order) ; ++i ) knots[i] = ++kval;
	++kval;
	for ( i=(nknots-order) ; i<nknots ; ++i ) knots[i] = kval;
	
	// use this info to initialze the sisl curve structure
	sislcpnts = (double*)&cpnts[0].x;
	scurve= newCurve(ncntl,  	// number of control points
					 order,  	// order of spline curve (degree + 1)
					 knots,  	// pointer to knot vector (parametrization)
					 sislcpnts,	// pointer to coefficient vector (control points)
					 kind,   	// 1-bspline 2-nurbs 3-bezier 4-rat bezier
					 3,      	// dimension
					 0);     	// no copying of information, 'borrow' arrays
	if ( !scurve ) return 0;
	return 1;
}

/*!
Interpolate points from control points. The points array is passed to
the function and is filled with interpolated points. The number of interpolated
points is returned.

@param n - number of points to interpolate
@return number of interpolated points
 */
int RicCurve::Interpolate(int n)
{
	// allocate memory for points array
	npnts = n;
	pnts = new DPoint[npnts];
	

	int left=0;
	int j; 
	DPoint	*tpnts = new DPoint[npnts];
	for (j=0; j<npnts; j++)
	{
		double t =	scurve->et[scurve->ik-1]+
					(scurve->et[scurve->in]-scurve->et[scurve->ik-1])*
					j/(npnts-1.0);
		int stat;
		    
		s1227(scurve, 0, t, &left, (double*)&tpnts[j], &stat);
		if (stat!=0)
			cerr << "s1227 returned status" << stat << endl;
	}
	
	// It turns out that these interpolated points are not evenly spaced.
	// The quick and dirty fix is to interpolate at equal increments between
	// the sisl interpolated points.
	
	// sum the point to point length
	double l = 0;
	for ( j=0 ; j<npnts-1; ++j )
	{
		l += ddist(tpnts[j],tpnts[j+1]);
	}
	double inc = l/(npnts-1);
	
	// first and last points are first and last points anyway
	pnts[0] = tpnts[0];
	pnts[npnts-1] = tpnts[npnts-1];
	
	int i;
	DPoint p1,p2;
	double d1,d2,d12;
	p1 = tpnts[0];
	p2 = tpnts[1];
	d1 = 0;
	d2 = d12 = ddist(p1,p2); // start with first two points
	int idx1,idx2;
	idx1 = 0;
	idx2 = 1;
	double tl;	// target length
	for ( i=1 ; i<npnts-1 ; ++i )
	{
		tl = i*inc; // target to this length
		
		// get the points surrounding the target distance
		do
		{
			if ( tl > d1 && tl < d2 ) // right place so break
				break;
			
			// increment 
			++idx1;
			++idx2;
			p1 = tpnts[idx1];
			p2 = tpnts[idx2];
			
			d1 += d12;
			d12 = ddist(p1,p2);
			d2 += d12;
			
		} while ( tl > d1 );
		
		// now interpolate the point
		double frac = (tl-d1)/d12;
		pnts[i].x = p1.x + (p2.x-p1.x)*frac;
		pnts[i].y = p1.y + (p2.y-p1.y)*frac;
		pnts[i].z = p1.z + (p2.z-p1.z)*frac;
	}
		
	delete tpnts;
	
	return npnts;
}

/*!
Write out the curve in SISL GO format
@param filename - output file name
 */
void RicCurve::WriteCurveGo(string filename)
{
	// make sure there is a valid file name
	if ( filename.find(".g2") == string::npos
				 & filename.find(".G2") == string::npos )
	{
		filename += ".g2";
	}

	// open output stream
	ofstream go_stream(filename.c_str());
	
    // write standard header
	go_stream << CURVE_INSTANCE_TYPE << " " << MAJOR_VERSION << ' ' << MINOR_VERSION << endl;
	
    // write basic curve properties
	const int& dim = scurve->idim;
	const int rational = (scurve->ikind % 2 == 0) ? 1 : 0;
	go_stream << dim << ' ' << rational << endl;

    // write bspline basis information
	write_go_basis(go_stream, scurve->in, scurve->ik, scurve->et);
   
    // write control points
	int coef_size = scurve->in * (rational ? (dim + 1) : dim);
	const double* coef_pointer = rational ? scurve->rcoef : scurve->ecoef;
	for (int i = 0; i < coef_size; ++i) {
		go_stream << coef_pointer[i] << ' ';
	}
	go_stream << endl;

	go_stream.close();
	
}

/*!
Write out the curve vertices in Anatomist mesh format
@param filename - output file name
 */
void RicCurve::WriteCurveAsMesh(string filename)
{
	// create the file in the BV mesh format
	ofstream fout(filename.c_str());
	fout<<"ascii"<<endl;
	fout<<"VOID"<<endl;
	fout<<2<<endl;	// dim 2 is a line
	fout<<1<<endl<<0<<endl; // only one timestep
	
	// output curve vertex points
	fout<<npnts<<endl;
	for (int i=0;i<npnts;i++)
		fout<<"("<<pnts[i].x<<", "<<pnts[i].y<<", "<<pnts[i].z<<") ";
	
	// zero is size for normals
	fout<<endl<<0<<endl; 
	
	// zero is size for texture
	fout<<0<<endl;
	
	// polygon size is npnts-1
	int p_size = npnts-1;
	fout<<p_size<<endl;
	
	for (int i=0;i<p_size;i++)
		fout<<"("<< i <<" ,"<< i+1 <<") ";

	fout << endl;
	
	fout.close();
}

/*!
Write out the control points in Anatomist mesh format
@param filename - output file name
 */
void RicCurve::WriteControlAsMesh(string filename)
{
	// create the file in the BV mesh format
	ofstream fout(filename.c_str());
	fout<<"ascii"<<endl;
	fout<<"VOID"<<endl;
	fout<<2<<endl;	// dim 2 is a line
	fout<<1<<endl<<0<<endl; // only one timestep
	
	// output curve vertex points
	fout<<ncntl<<endl;
	for (int i=0;i<ncntl;i++)
		fout<<"("<<cpnts[i].x<<", "<<cpnts[i].y<<", "<<cpnts[i].z<<") ";
	
	// zero is size for normals
	fout<<endl<<0<<endl; 
	
	// zero is size for texture
	fout<<0<<endl;
	
	// polygon size is npnts-1
	int p_size = ncntl-1;
	fout<<p_size<<endl;
	
	for (int i=0;i<p_size;i++)
		fout<<"("<< i <<" ,"<< i+1 <<") ";

	fout << endl;
	
	fout.close();
}

/*!
Sort control points along length. It is assumed that the two points farthest
apart are really the line endpoints.
return 1 is returned on success or 0 on failure.
*/
int RicCurve::SortAlongLength()
{
	int i,j;
	int e0,e1;	// endpoint indices
	float max;	// max distance
	float min;	// min distance
	float d;
	
	// check to see that we really have control points
	if ( !ncntl ) return 0;
	
	// pick a point to start with and find the point furthest away
	int st = ncntl/2;
	max = 0;
	e0=0;
	for ( i=0 ; i<ncntl ; ++i )
	{
		d = ddist(cpnts[st],cpnts[i]);
		if ( d > max )
		{
			max = d;
			e0 = i;
		}
	}

	// repeat using e0 as one end	
	max = 0;
	e1=0;
	for ( i=0 ; i<ncntl ; ++i )
	{
		d = ddist(cpnts[e0],cpnts[i]);
		if ( d > max )
		{
			max = d;
			e1 = i;
		}
	}
	
	// now sort the points finding the closest one to the previous one.
	DPoint *sarray = new DPoint[ncntl];	// temp array for sorting indices
	int *sidx = new int[ncntl]; // keeps track of array indices already used
	
	for ( i=0 ; i<ncntl ; ++i ) sidx[i] = 0;
			
	// first and last we already know
	sarray[0] = cpnts[e0];
	sarray[ncntl-1] = cpnts[e1];
	
	// look for the rest starting with e0
	DPoint spnt = cpnts[e0];
	int midx=0;
	for ( i=1 ; i<ncntl-1 ; ++i )
	{
		min = 10000000;
		for ( j=1 ; j<ncntl-1 ; ++j )
		{
			// skip if alreay used
			if ( sidx[j] ) continue;
			
			d = ddist(spnt,cpnts[j]);
			if ( d < min )
			{
				min = d;
				midx = j;
			}
		}
		
		// assign next array place to closest point and use as new start
		spnt = sarray[i] = cpnts[midx];
		
		// mark the closest point as used
		sidx[midx] = 1;
	}
	
	// copy temporary array back
	for ( i=0 ; i<ncntl ; ++i )
		cpnts[i] = sarray[i];
	
	delete [] sarray;
	delete [] sidx;
	
	return 1;
}

/*!
Calculate the length of a parametric curve. The tolerance is fixed at this point.
@return the length is returned or 0 on error.
 */
double RicCurve::CurveLength(void)
{
	// return if no control points
	if ( !ncntl ) return 0;
	
	// call the sisl routine for curve length
	double len;
	int stat;
	s1240(this->scurve,0.000005,&len,&stat);
	
	// check for error
	if ( stat ) return 0;
	return len;
}
