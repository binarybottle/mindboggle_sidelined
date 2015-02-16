// ---------------------------- SulcalRayTracer.cpp --------------------------

#include "SulcalRayTracer.h"

// local function
float InterpolateAcrossTriangle(Point tri[3], Point p, float thick[3]);

/*!
Constructor
*/
SulcalRayTracer::SulcalRayTracer()
{

	max_threshold=0;
	min_threshold=0;
	polygons_found=0;
	avg_thickness=0;
	boundary_thick = 0;
	number_steps=0;
	std_dev=0;
	std_tot=0;
	avg_tot=0;
	polygons_near_sulcus = NULL;
	EndPoints = NULL;
	ThickVal = NULL;
}

/*!
Constructor with with passed boundary mesh and distance thresholds

@param hemi_mesh - pointer to boundary mesh
@param ma_threshold - max distance from sulcal mesh point
@param mi_threshold - min distance from sulcal mesh point
*/
SulcalRayTracer::SulcalRayTracer(RicMesh *hemi_mesh, float ma_threshold, float mi_threshold)
{
	boundary_mesh=hemi_mesh;
	max_threshold=ma_threshold;
	min_threshold=mi_threshold;
	polygons_near_sulcus = new int[hemi_mesh->p_size ];
	polygons_found=0;
	avg_thickness=0;
	boundary_thick = 0;
	number_steps=0;
	std_dev=0;
	std_tot=0;
	avg_tot=0;
	ThickVal = NULL;

	// allocate memory for end points of sulcal span lines
	matrix(&EndPoints,3,hemi_mesh->v_size);

}

/*!
Constructor with with passed boundary mesh, distance thresholds, and thickness texture

@param hemi_mesh - pointer to boundary mesh
@param ma_threshold - max distance from sulcal mesh point
@param mi_threshold - min distance from sulcal mesh point
@param tex - texture mapping boundary mesh points to sulcal thickness
*/
SulcalRayTracer::SulcalRayTracer(RicMesh *hemi_mesh, float ma_threshold, float mi_threshold,
		RicTexture *tex)
{
	boundary_mesh=hemi_mesh;
	thick_tex=tex;
	max_threshold=ma_threshold;
	min_threshold=mi_threshold;
	polygons_near_sulcus = new int[hemi_mesh->p_size ];
	polygons_found=0;
	avg_thickness=0;
	boundary_thick = 0;
	number_steps=0;
	std_dev=0;
	std_tot=0;
	avg_tot=0;

	// allocate memory for end points of sulcal span lines
	matrix(&EndPoints,3,hemi_mesh->v_size);

	// allocate memory for thickness values
	ThickVal = new float[hemi_mesh->v_size];
}

/*!
Destructor
*/
SulcalRayTracer::~SulcalRayTracer()
{
	if ( polygons_near_sulcus ) delete [] polygons_near_sulcus;
	if ( EndPoints ) free_matrix(EndPoints);
}

/*!
 Calculates the min, max, average, and std deviation of the sulcal
 width vectors.

 returns - average sulcal width
 */
float SulcalRayTracer::CalcMinMaxAvg()
{
	min = 10000000;
	max = -10000000;
	std_tot=0;
	avg_tot=0;
	double d;
	float avg_bthick=0; // boundary thickness

	boundary_thick = 0;

	for (long i=0; i<number_steps; ++i)
	{

		d = dist(EndPoints[0][i],EndPoints[2][i]);

		if (d > max)
			max = d;
		if (d < min)
			min = d;
		avg_tot += d;
		std_tot += d*d;

		if ( ThickVal ) // if there is a boundary thickness array
		{
			avg_bthick += ThickVal[i];
		}
	}

	// average
	avg_thickness = avg_tot/(double)number_steps;
	boundary_thick = avg_bthick/(float)number_steps;

	// standard deviation
	std_dev = sqrt(fabs(std_tot-avg_tot*avg_tot/(double)number_steps)/((double)number_steps-1.0));

	return avg_thickness;
}

/*!
This member function calculates whether a point is within a box bounding
the sulcus

@param p - point to test
@param sulcus - pointer to sulcus
@param expand - distance to expand the sulcus bounding box

returns 1 if point in bounding box
*/
int SulcalRayTracer::is_point_inside_box(vertex p, RicMesh *sulcus, float expand)
{
	if ( ( (p.pnt.x<= sulcus->xmax+expand) && (p.pnt.x>= sulcus->xmin -expand) )
			&& ( (p.pnt.y<= sulcus->ymax+expand) && (p.pnt.y >= sulcus->ymin
					-expand) )&& ( (p.pnt.z<= sulcus->zmax+expand) && (p.pnt.z
			>= sulcus->zmin-expand)))
		return 1;
	else
		return 0;

}

/*!
This member function calculates whether a triangle is within a box bounding
the sulcus

@param p - triangle to test
@param sulcus - pointer to sulcus
@param expand - distance to expand the sulcus bounding box

returns 1 if triangle in bounding box
*/
int SulcalRayTracer::is_polygon_inside_box(triangle p, RicMesh *sulcus, float expand)
{
	if (is_point_inside_box(boundary_mesh->vertices[p.vidx[0]], sulcus, expand)
			|| is_point_inside_box(boundary_mesh->vertices[p.vidx[1]], sulcus,
					expand) || is_point_inside_box(
			boundary_mesh->vertices[p.vidx[2]], sulcus, expand)

	)
		return 1;
	else
		return 0;
}

/*!
This function finds all the polygons from the boundary mesh that are in a
bounding box about the sulcus

@param sulcus - pointer to sulcal mesh
@param max_distance - amount to expand bounding box around sulcus

returns - the number of polygons found
*/
int SulcalRayTracer::match_sulcus_to_mesh(RicMesh *sulcus, float max_distance)
{
	polygons_found=0;
	for (long i=0; i<boundary_mesh->p_size; i++)
	{
		if (is_polygon_inside_box(boundary_mesh->polygons[i], sulcus,
				max_distance)==1)
			polygons_near_sulcus[polygons_found++]=i;
	}
	return --polygons_found;
}

/*!
Writes sulcal width endpoints to 2D mesh file. This mesh is a set of lines
where each line corresponds to a sulcal with measurement.
@param name - name of output file

returns 1
*/
int SulcalRayTracer::WriteEndpointMesh(string name)
{
	RicMesh *EndPointMesh;
	int msize = number_steps;
	EndPointMesh = new RicMesh(msize*2,0,msize,2);

	// copy the endpoints to 2D mesh
	for ( int k=0 ; k<msize ; ++k )
	{
		EndPointMesh->assign_node(2*k, EndPoints[0][k]);
		EndPointMesh->assign_node(2*k+1, EndPoints[2][k]);
		EndPointMesh->assign_polygon(k, 2*k, 2*k+1, 0);
	}

	// write output file
	EndPointMesh->write_mesh_txt((char*)name.c_str());

	delete EndPointMesh;

	return 1;
}

/*!
Ths function filters the sulcal distance measurements by checking the ratio
of the distances of from a sulcal mesh point to the boundary mesh in the
positive and negative normal directions. If the ratio is greater than the
passed ratio value then that measurement is discarded.

@param ratio - reference ratio

returns 1
*/
int SulcalRayTracer::FilterEndpointsByRatio(float ratio)
{
	// check the ratio the distance between each endpoint and the center point
	int n=0;
	for ( int k=0 ; k<number_steps ; ++k )
	{
		float d1,d2;
		d1 = dist(EndPoints[0][k],EndPoints[1][k]);
		d2 = dist(EndPoints[2][k],EndPoints[1][k]);

		float r = d1/d2;
		if ( r < 1 ) r = d2/d1; // use the inverse

		// keep endpoints that have less than allowed ratio
		if ( r < ratio )
		{
			EndPoints[0][n] = EndPoints[0][k];
			EndPoints[1][n] = EndPoints[1][k];
			EndPoints[2][n] = EndPoints[2][k];

			if ( ThickVal ) // if there is a boundary thickness map
				ThickVal[n] = ThickVal[k];

			++n;
		}
	}

	// reassign number of steps to valid endpoints
	number_steps = n;

	return 1;
}

/*!
Ths function filters the sulcal distance measurements checking width
to see if they are in a range specified by a multiple of the std deviation
of the distance values.

@param mult - value to multiply std deviation by to set range values

returns 1
*/
int SulcalRayTracer::FilterEndpointsByStddev(float mult)
{
	// make sure we have the latest std dev
	this->CalcMinMaxAvg();

	// check the ratio the distace between each endpoint and the center point
	int n=0;
	float minval = avg_thickness - mult*std_dev;
	float maxval = avg_thickness + mult*std_dev;
	for ( int k=0 ; k<number_steps ; ++k )
	{
		float d = dist(EndPoints[0][k],EndPoints[2][k]);
		// keep endpoints that have less than allowed ratio

		if ( d > minval && d < maxval )
		{
			EndPoints[0][n] = EndPoints[0][k];
			EndPoints[1][n] = EndPoints[1][k];
			EndPoints[2][n] = EndPoints[2][k];

			if ( ThickVal ) // if there is a boundary thickness map
				ThickVal[n] = ThickVal[k];

			++n;
		}
	}

	// reassign number of steps to valid endpoints
	number_steps = n;

	return 1;
}

/*!
This multiplys all the end point values by a transformation matrix

@param tm - transformation matrix

returns - 1
*/
int SulcalRayTracer::TransformEndpoints(CMatrix tm)
{
	for ( int i=0 ; i<3 ; ++i )
	{
		for ( int j=0 ; j<number_steps ; ++j )
		{
			EndPoints[i][j] = tm*EndPoints[i][j];
		}
	}
	return 1;
}

/*!
This is the function that actually calculates the sulcal width. The first
step is to find all the polygons from the boundary mesh that are near the
sulcal mesh. The second step is to draw vectors from each sulcal mesh point
in the plus and minus normal directions to see if they hit a boundary mesh
polygon. If there are hits in both directions then the distance between the
intersection points is a sulcal width measurement for that sulcal mesh
vertex.

@param sulcus - pointer to sulcal mesh

returns number of distance measurements on success or 0 on failure
 */
int SulcalRayTracer::CalcSulcalSpan(RicMesh *sulcus)
{
	vertex node, norm, top_point;

	sulcus->calc_limits(); //will be safe to call it again

	// Step 1 find neighboring boundary mesh polygons
	polygons_found=match_sulcus_to_mesh(sulcus, 20);
	if (polygons_found<10)
	{
		cerr<<"run_csf_sulcal_tracing - Not enough polygons found"<<endl;
		return 0;
	}

	// Step 2 - check vectors from each sulcal mesh point to see if they
	// intersect nearby boundary mesh polygons. Find the distances to the
	// nearest boundary polygons in the plus and minus normal directions.
	int not_found=0;
	for (int i=0; i<sulcus->v_size; ++i)
	{
		norm=sulcus->get_normal(i); // get normal
		node=sulcus->get_node(i);	// get point
		Point p0;	// current mesh point
		Point pplus,pminus;	// end points in plus and minus directions
		Point npnt;	// normal vector
		Point nmin,pmin;	// closest cortex points along sulcal span line

		p0 = node.pnt;

		// multiply the unit normal times the max threshod
		npnt = norm.pnt;
		npnt.x *= max_threshold;
		npnt.y *= max_threshold;
		npnt.z *= max_threshold;

		// add the normal vector to the current point for plus directions
		// and subtract it for the minus directions
		pplus = p0 + npnt;
		pminus = p0 - npnt;

		float minplus=1000000,minminus=1000000,d,dc;
		Point t0,t1,t2,negpnt,pospnt;
		float d0,d1,d2;

		// check against every nearby boundary polygon
		for ( int j=0 ; j<polygons_found ; ++j )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			t0 = boundary_mesh->vertices[boundary_mesh->polygons[polygons_near_sulcus[j] ].vidx[0]].pnt;
			t1 = boundary_mesh->vertices[boundary_mesh->polygons[polygons_near_sulcus[j] ].vidx[1]].pnt;
			t2 = boundary_mesh->vertices[boundary_mesh->polygons[polygons_near_sulcus[j] ].vidx[2]].pnt;

			// skip if any of the triangle points are our current vertex
			if ( (d0=dist(p0,t0)) < 0.1 ) continue;
			if ( (d1=dist(p0,t1)) < 0.1 ) continue;
			if ( (d2=dist(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>max_threshold && d1>max_threshold && d2>max_threshold ) continue;

			// plus distance

			// if no intersection of line with plane then skip
			if ( line_intersect_triangle(t0, t1, t2,p0,pplus,&pospnt) )
			{
				// compare the distance to the current min distance
				dc = dist(pospnt,p0);
				if ( dc!=0 )
				{
					if ( dc < minplus )
					{
						minplus = dc;
						pmin = pospnt;
					}
				}
			}

			// minus distance

			// if no intersection of line with plane then skip
			if ( line_intersect_triangle(t0, t1, t2,p0,pminus,&negpnt) )
			{
				// compare the distance to the current min wm distance
				dc = dist(negpnt,p0);
				if ( dc!=0 )
				{
					if ( dc < minminus )
					{
						minminus = dc;
						nmin = negpnt;
					}
				}
			}
		}

		// skip if we have not found both plus and minus points
		if ( minminus < max_threshold && minplus < max_threshold )
		{
			d = minplus+minminus;
			if ( d > min_threshold && d < max_threshold )
			{
				EndPoints[0][number_steps] = nmin;
				EndPoints[1][number_steps] = p0;
				EndPoints[2][number_steps] = pmin;
				++number_steps;
			}
		}
		else
		{
			++not_found;
			continue;
		}
	}

	return number_steps;
}

/*!
This is the function that actually calculates the sulcal width. The first
step is to find all the polygons from the boundary mesh that are near the
sulcal mesh. The second step is to draw vectors from each sulcal mesh point
in the plus and minus normal directions to see if they hit a boundary mesh
polygon. If there are hits in both directions then the distance between the
intersection points is a sulcal width measurement for that sulcal mesh
vertex.

This version also matches a cortical thickness texture to the sulcal span
measurements. The thickness value for a given sulcal width measurement is
the average of the thickness values on the boundary mesh where the vector
intersects it in both directions.

@param sulcus - pointer to sulcal mesh

returns number of distance measurements on success or 0 on failure
 */
int SulcalRayTracer::CalcSulcalSpanThick(RicMesh *sulcus)
{
	vertex node, norm, top_point;

	sulcus->calc_limits(); //will be safe to call it again

	// Step 1 find neighboring boundary mesh polygons
	polygons_found=match_sulcus_to_mesh(sulcus, 20);
	if (polygons_found<10)
	{
		cerr<<"run_csf_sulcal_tracing - Not enough polygons found"<<endl;
		return 0;
	}

	// Step 2 - check vectors from each sulcal mesh point to see if they
	// intersect nearby boundary mesh polygons. Find the distances to the
	// nearest boundary polygons in the plus and minus normal directions.
	int not_found=0;
	for (int i=0; i<sulcus->v_size; ++i)
	{
		norm=sulcus->get_normal(i); // get normal
		node=sulcus->get_node(i);	// get point
		Point p0;	// current mesh point
		Point pplus,pminus;	// end points in plus and minus directions
		Point npnt;	// normal vector
		Point nmin,pmin;	// closest cortex points along sulcal span line

		p0 = node.pnt;

		// multiply the unit normal times the max threshod
		npnt = norm.pnt;
		npnt.x *= max_threshold;
		npnt.y *= max_threshold;
		npnt.z *= max_threshold;

		// add the normal vector to the current point for plus directions
		// and subtract it for the minus directions
		pplus = p0 + npnt;
		pminus = p0 - npnt;

		float minplus=1000000,minminus=1000000,d,dc;
		Point t0,t1,t2,negpnt,pospnt;
		float d0,d1,d2;
		int negpoly=0,pospoly=0;	// indices to nearest distance polygons

		// check against every nearby boundary polygon
		for ( int j=0 ; j<polygons_found ; ++j )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			t0 = boundary_mesh->vertices[boundary_mesh->polygons[polygons_near_sulcus[j] ].vidx[0]].pnt;
			t1 = boundary_mesh->vertices[boundary_mesh->polygons[polygons_near_sulcus[j] ].vidx[1]].pnt;
			t2 = boundary_mesh->vertices[boundary_mesh->polygons[polygons_near_sulcus[j] ].vidx[2]].pnt;

			// skip if any of the triangle points are our current vertex
			if ( (d0=dist(p0,t0)) < 0.1 ) continue;
			if ( (d1=dist(p0,t1)) < 0.1 ) continue;
			if ( (d2=dist(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>max_threshold && d1>max_threshold && d2>max_threshold ) continue;

			// plus distance

			// if no intersection of line with plane then skip
			if ( line_intersect_triangle(t0, t1, t2,p0,pplus,&pospnt) )
			{
				// compare the distance to the current min distance
				dc = dist(pospnt,p0);
				if ( dc!=0 )
				{
					if ( dc < minplus )
					{
						minplus = dc;
						pmin = pospnt;
						pospoly = j;	// keep track of the polygon we used here
					}
				}
			}

			// minus distance

			// if no intersection of line with plane then skip
			if ( line_intersect_triangle(t0, t1, t2,p0,pminus,&negpnt) )
			{
				// compare the distance to the current min wm distance
				dc = dist(negpnt,p0);
				if ( dc!=0 )
				{
					if ( dc < minminus )
					{
						minminus = dc;
						nmin = negpnt;
						negpoly = j;	// keep track of the polygon we used here
					}
				}
			}
		}

		// skip if we have not found both plus and minus points
		if ( minminus < max_threshold && minplus < max_threshold )
		{
			d = minplus+minminus;
			if ( d > min_threshold && d < max_threshold )
			{
				// store end points of vectors intersecting boundary mesh
				EndPoints[0][number_steps] = nmin;
				EndPoints[1][number_steps] = p0;
				EndPoints[2][number_steps] = pmin;

				// interpolate thickness from intersected polygons and
				// the thickness map
				Point tri[3];	// vertex points from triangle
				float tmap[3];	// thickness for vertices
				float negthick,posthick;	// thickness interpolated

				// positive thickness
				tri[0] = boundary_mesh->vertices[boundary_mesh->polygons[pospoly].vidx[0]].pnt;
				tri[1] = boundary_mesh->vertices[boundary_mesh->polygons[pospoly].vidx[1]].pnt;
				tri[2] = boundary_mesh->vertices[boundary_mesh->polygons[pospoly].vidx[2]].pnt;
				tmap[0] = thick_tex->nodes[boundary_mesh->polygons[pospoly].vidx[0]];
				tmap[1] = thick_tex->nodes[boundary_mesh->polygons[pospoly].vidx[1]];
				tmap[2] = thick_tex->nodes[boundary_mesh->polygons[pospoly].vidx[2]];

				posthick = InterpolateAcrossTriangle(tri,p0,tmap);

				// positive thickness
				tri[0] = boundary_mesh->vertices[boundary_mesh->polygons[negpoly].vidx[0]].pnt;
				tri[1] = boundary_mesh->vertices[boundary_mesh->polygons[negpoly].vidx[1]].pnt;
				tri[2] = boundary_mesh->vertices[boundary_mesh->polygons[negpoly].vidx[2]].pnt;
				tmap[0] = thick_tex->nodes[boundary_mesh->polygons[negpoly].vidx[0]];
				tmap[1] = thick_tex->nodes[boundary_mesh->polygons[negpoly].vidx[1]];
				tmap[2] = thick_tex->nodes[boundary_mesh->polygons[negpoly].vidx[2]];

				negthick = InterpolateAcrossTriangle(tri,p0,tmap);

				// average thickness values
				ThickVal[number_steps] = (posthick+negthick)*0.5f;

				++number_steps;
			}
		}
		else
		{
			++not_found;
			continue;
		}
	}

	return number_steps;
}

// Interpolate
float InterpolateAcrossTriangle(Point tri[3], Point p, float thick[3])
{
	// find distance from each vertex to point
	float d0,d1,d2;
	d0 = dist(tri[0],p);
	d1 = dist(tri[1],p);
	d2 = dist(tri[2],p);

	// normalization factor
	float nfac;
	nfac = 1.0f/(d0+d1+d2);

	// interpolated value
	float val;
	val = nfac*d0*thick[0] + nfac*d1*thick[1] + nfac*d2*thick[2];

	return val;
}
