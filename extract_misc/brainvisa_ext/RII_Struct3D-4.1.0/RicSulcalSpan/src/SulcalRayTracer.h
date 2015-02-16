// ------------------------------- SulcalRayTracer.h -------------------------
/*! \file
Header file for SulcalRayTracer class
Copyright (C) 2008 by Bill Rogers - Research Imaging Center - UTHSCSA
*/
#ifndef SULCALRAYTRACER_H_
#define SULCALRAYTRACER_H_

#include <RicUtil.h>
#include <RicMesh.h>
#include <RicTexture.h>

/*!
This class used used to determine the sulcal width of a sulcus. That sulcus can
be composed of one or more segments. The method is for each sulcal mesh point,
trace a normal vector in both directions from the sulcal mesh surface to where
it intersects the boundary mesh (usually cortical mesh). The output can be
filtered and a vector mesh showing the boundary mesh intersections can be
output.
*/
class SulcalRayTracer
{
public:
	RicMesh *boundary_mesh;		///< boundary mesh (usually cortical mesh)
	float avg_thickness;		///< average sulcal width
	double avg_tot;				///< running total for average
	float std_dev;				///< standard deviation of width
	double std_tot;				///< running total for std dev
	float min;					///< min width value
	float max;					///< max width value
	int number_steps;			///< number of vectors used for average
	float max_threshold;		///< max distance from sulcal mesh to boundary allowed
	float min_threshold;		///< min distance from sulcal mesh to boundary allowed
	int *polygons_near_sulcus;	///< this is the array of the indices where we're going to keep polygon ids that are found near the sulcus
	int polygons_found;			///< number of boundary mesh polygons near sulcus
	Point	**EndPoints;		///< matrix of end points for sulcal width line segments
	float	*ThickVal;			///< array of thickness values
	RicTexture *thick_tex;		///< thickness texture corresponding to boundary texture
	float	boundary_thick;		///< average thickness of boundary mesh around sulcus

	// constructors
	SulcalRayTracer();
	SulcalRayTracer(RicMesh *hemi_mesh, float ma_threshold, float mi_threshold);
	SulcalRayTracer(RicMesh *hemi_mesh, float ma_threshold, float mi_threshold, RicTexture *tex);
	~SulcalRayTracer();

private:
	// a member function to calculate whether a point is within a box
	int is_point_inside_box(vertex p, RicMesh *sulcus, float expand);
	// a member function to calculate whether a triangle is within a bounding box
	int is_polygon_inside_box(triangle p, RicMesh *sulcus, float expand);
	int match_sulcus_to_mesh(RicMesh *sulcus, float max_distance);

public:

	int CalcSulcalSpan(RicMesh *sulcus);
	int CalcSulcalSpanThick(RicMesh *sulcus);
	int WriteEndpointMesh(string name);
	int FilterEndpointsByRatio(float ratio);
	int FilterEndpointsByStddev(float mult);
	int TransformEndpoints(CMatrix tm);
	float CalcMinMaxAvg(void);

};

#endif /*SULCALRAYTRACER_H_*/
