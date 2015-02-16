// ---------------------------  GM_Normal.cpp ---------------------------------
/**
\file
Routines to find the distance between two meshes by intersecting
normals from one mesh with triangles of the other.

Copyright (C) 2007 by Bill Rogers, Research Imaging Center, UTHSCSA
rogers@uthscsa.edu
*/

#include "RicMesh.h"
#include "RicTexture.h"
#include "GM_Normal.h"
#include <math.h>
#include <pthread.h>
#include "RicUtil.h"

/// Structure to pass info to the threads for threaded version
typedef struct { int *inner_vlist;	///< array of indices for inner mesh vertices
		int n_in_verts; 		///< number of inner mesh vertices
		int *inner_plist;		///< array of indicies for inner mesh polygons
		int n_in_polys;			///< number of inner mesh polygon indices
		int *outer_plist;		///< array of indicies for outer mesh polygons
		int n_out_polys;		///< number of outer mesh polygon indices
		RicMesh *inner_mesh;	///< inner mesh
		RicMesh *outer_mesh;	///< outer mesh
		float	*inner_thick;	///< thickness array - one node for each inner mesh vertex
		RicMesh *closest_vects;	///< mesh of vectors connecting inner and outer meshes
		float mind;				///< minimum allowed distance from vertex
		float maxd;				///< maximum allowed distance from vertex
		float sfac;				///< scale factor for normal line
} NormThreadStuff;

/// thread function declaration
void *FindNormalDistThread( void *tstruct );


///////////////////////// FindNormalDistBrute ////////////////////////
/**
This routine determines the thickness at each inner mesh vertex by
finding the intersection of a line normal of that vertex with a triangle
in the outer mesh

@param inner_mesh - white matter mesh
@param outer_mesh - gray matter mesh
@param closest_vects - output mesh with vectors connecting nearest points
@param thick - output texture to put inner-outer distances in - one value for each inner vertex
@param nflip - normal direction, if -1 flip normals
@param mind - min distance to search for normal intersection with other mesh
@param maxd - max distance to search for normal intersection with other mesh
@return - 1 on success
 */
int FindNormalDistBrute(RicMesh *inner_mesh, RicMesh *outer_mesh, RicMesh *closest_vects,
							RicTexture *thick, int nflip, float mind, float maxd)
{
	// check to see that the gm mesh and the texture are the same size
	if ( inner_mesh->v_size != thick->size )
	{
		cout << "FindNormalDistSubdivide - error inner mesh and texture not same size"<<endl;
		return 0;
	}

	int i,j;
	Point minpnt;
	int ninfound, noutfound;
	float indist,outdist;
	Point pint;

	// scale factor for normal line end point - includes normal direction
	float sfac = nflip*maxd;

	float maxdsqu = maxd*maxd;	// square of distance for comparisons

	// populate the thickness map with ERRVAL as default value
	for ( i=0 ; i<thick->size ; ++i ) thick->nodes[i] = ERRVAL;

	// work from inner to outer
	for ( i=0 ; i<inner_mesh->v_size ; ++i )
	{
		// skip if vertex labeled as not to use
		if ( inner_mesh->vertices[i].label == 1 )
			continue;

		Point p0,p1;	// line normal to vertex
		Vector n0,n1;	// normals from vertex and opposing polygon

		p0 = inner_mesh->vertices[i].pnt;
		n0 = inner_mesh->normals[i].pnt;

		// make a long line in the direction of the normal
		p1.x = p0.x + sfac*n0.x; // add normal to vertex to get line
		p1.y = p0.y + sfac*n0.y;
		p1.z = p0.z + sfac*n0.z;

		// test against all polygons in the inner mesh
		// to find the closest inner polygon
		ninfound=0;
		indist=ERRVAL;
		Point t0,t1,t2;
		float d0,d1,d2;
		for ( j=0 ; j<inner_mesh->p_size ; ++j )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			t0 = inner_mesh->vertices[inner_mesh->polygons[j].vidx[0]].pnt;
			t1 = inner_mesh->vertices[inner_mesh->polygons[j].vidx[1]].pnt;
			t2 = inner_mesh->vertices[inner_mesh->polygons[j].vidx[2]].pnt;

			// skip if any of the triangle points are our current vertex
			if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
			if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
			if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>maxdsqu && d1>maxdsqu && d2>maxdsqu ) continue;

			if ( !line_intersect_triangle(t0, t1, t2,p0,p1,&pint) )
				continue;

			// compare the distance to the current min wm distance
			++ninfound;
			float d = dist(pint,p0);
			if ( d==0 )
				continue;

			if ( d < indist )
				indist = d;

		}

		// test against all polygons in the outer mesh
		noutfound=0;
		outdist = ERRVAL;
		for ( j=0 ; j<outer_mesh->p_size ; ++j )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			Point t0,t1,t2;
			t0 = outer_mesh->vertices[outer_mesh->polygons[j].vidx[0]].pnt;
			t1 = outer_mesh->vertices[outer_mesh->polygons[j].vidx[1]].pnt;
			t2 = outer_mesh->vertices[outer_mesh->polygons[j].vidx[2]].pnt;

			// skip if any of the triangle points are our current vertex
			if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
			if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
			if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>maxd && d1>maxd && d2>maxd ) continue;

			if ( !line_intersect_triangle(t0, t1, t2,p0,p1,&pint) )
				continue;

			// compare the distance to the current min wm distance
			++noutfound;
			float d = dist(pint,p0);
			if ( (d < indist) && (d < maxd) ) // must be closer than nearest wm polygon
			{
				if ( d < outdist ) // see if closer than last suitable distance
				{
					outdist = d;
					minpnt = pint;
				}
			}

		}

		// check against minimum distance
		if ( outdist < mind )
			outdist = mind;

		// see if a vertex found but too far away
		if ( outdist == ERRVAL && noutfound > 0 )
			outdist = maxd;

		// assign distance for vertex
		thick->nodes[i] = outdist;

		// assign vector for points connecting surfaces
		if ( outdist != ERRVAL && outdist > mind && outdist < maxd )
		{
			vertex cvert(minpnt.x,minpnt.y,minpnt.z);
			closest_vects->assign_node(2*i, inner_mesh->vertices[i]);
			closest_vects->assign_normal(2*i,inner_mesh->normals[i]);
			closest_vects->assign_node(2*i+1, cvert);
			closest_vects->assign_normal(2*i+1, inner_mesh->normals[i]);
			closest_vects->assign_polygon(i, 2*i, 2*i+1, 0);
		}
	}
	return 1;
}


//////////////////////// FindNormalDistSubdivide //////////////////////////////////
/**
Find the distance between two meshes by taking a normal from the vertex
of one mesh and intersecting the second mesh. The normal is extended
to the maximum distance to search for an intersection.
Care is taken so that the normal distance chosen does not intersect any
surfaces between end points.
This version subdivides the meshes into smaller portions so that the search will
be faster.

@param inner_mesh - grey matter mesh
@param outer_mesh - white matter mesh
@param closest_vects - output mesh with vectors connecting nearest points
@param thick - output texture to put inner-outer distances in - one value for each inner vertex
@param nsub - number of subdivisions for each axis
@param over - amount of overlap to add to subdivisions
@param nflip - normal direction, if -1 flip normals
@param mind - min distance to search for normal intersection with other mesh
@param maxd - max distance to search for normal intersection with other mesh
@return - 1 on success
 */
int FindNormalDistSubdivide(RicMesh *inner_mesh, RicMesh *outer_mesh, RicMesh *closest_vects,
				RicTexture *thick, int nsub, float over, int nflip, float mind, float maxd)
{
	int i,j,k,l,m,n;

	// check to see that the inner mesh and the texture are the same size
	if ( inner_mesh->v_size != thick->size )
	{
		cout << "FindNormalDistSubdivide - error inner mesh and texture not same size"<<endl;
		return 0;
	}

	// array pointers for vertex and polygon lists for subdividing the meshes
	int *inner_vlist;
	int *inner_plist;
	int *outer_plist;

	// allocate memory for the arrays
	inner_vlist = new int[inner_mesh->v_size];
	inner_plist = new int[inner_mesh->p_size];
	outer_plist = new int[outer_mesh->p_size];

	// scale factor for normal line end point - includes normal direction
	float sfac = nflip*maxd;

	float maxdsqu = maxd*maxd;	// square of distance for comparisons

	// make sure that limits have been calculated for the gm mesh
	inner_mesh->calc_limits();

	// find the step size for each axis for an iteration
	float xinc,yinc,zinc;
	xinc = (inner_mesh->xmax-inner_mesh->xmin)/(float)nsub;
	yinc = (inner_mesh->ymax-inner_mesh->ymin)/(float)nsub;
	zinc = (inner_mesh->zmax-inner_mesh->zmin)/(float)nsub;

	// initialize the thickness values to ERRVAL
	for ( i=0 ; i<thick->size ; ++i ) thick->nodes[i]=ERRVAL;

	// the number of actual searches will be the cube of the subdivision number
	float xstart,ystart,zstart;	// starting values for vertex search
	float xend,yend,zend;		// ending value for vertex search
	float xstart2,ystart2,zstart2;	// starting values for vertex search
	float xend2,yend2,zend2;		// ending value for vertex search
	int n_in_verts = 0;	// number of inner mesh vertices found in a subdivision
	int n_in_polys = 0;	// number of inner mesh polygons found in a subdivision
	int n_out_polys = 0;// number of outer mesh polygons found in a subdivision
	for ( i=0 ; i<nsub ; ++i )
	{
		xstart = inner_mesh->xmin + i*xinc;
		xend = xstart+xinc;

		for ( j=0 ; j<nsub ; ++j )
		{
			ystart = inner_mesh->ymin + j*yinc;
			yend = ystart+yinc;

			for ( k=0 ; k<nsub ; ++k )
			{
				///////////////// preprocess for current subdivision ///////////
				// find all the vertices and polygons in these bounds

				n_in_verts = n_in_polys = n_out_polys = 0;

				zstart = inner_mesh->zmin + k*zinc;
				zend = zstart + zinc;

				// find the appropriate gm vertices in this box
				for ( l=0 ; l<inner_mesh->v_size ; ++l )
				{
					if ( inner_mesh->vertices[l].pnt.x >= xstart && inner_mesh->vertices[l].pnt.x <= xend
						&&	inner_mesh->vertices[l].pnt.y >= ystart && inner_mesh->vertices[l].pnt.y <= yend
						&&	inner_mesh->vertices[l].pnt.z >= zstart && inner_mesh->vertices[l].pnt.z <= zend )
					{
						inner_vlist[n_in_verts++] = l;
					}
				}

				// allow for overlap
				xstart2 = xstart-over;
				ystart2 = ystart-over;
				zstart2 = zstart-over;
				xend2 = xend+over;
				yend2 = yend+over;
				zend2 = zend+over;

				// look for inner triangles that fit in this box
				for ( l=0 ; l<inner_mesh->p_size ; ++l )
				{
					// check each triangle vertex
					for ( m=0 ; m<3 ; ++m )
					{
						int vidx = inner_mesh->polygons[l].vidx[m];
						if ( inner_mesh->vertices[vidx].pnt.x >= xstart2
							   && inner_mesh->vertices[vidx].pnt.x <= xend2
							   &&	inner_mesh->vertices[vidx].pnt.y >= ystart2
							   && inner_mesh->vertices[vidx].pnt.y <= yend2
							   &&	inner_mesh->vertices[vidx].pnt.z >= zstart2
							   && inner_mesh->vertices[vidx].pnt.z <= zend )
						{
							inner_plist[n_in_polys++] = l;
							break; // no need to check the others
						}
					}
				}

								// look for outer triangles that fit in this box
				for ( l=0 ; l<outer_mesh->p_size ; ++l )
				{
					// check each triangle vertex
					for ( m=0 ; m<3 ; ++m )
					{
						int vidx = outer_mesh->polygons[l].vidx[m];
						if ( outer_mesh->vertices[vidx].pnt.x >= xstart2
								&& outer_mesh->vertices[vidx].pnt.x <= xend2
								&&	outer_mesh->vertices[vidx].pnt.y >= ystart2
								&& outer_mesh->vertices[vidx].pnt.y <= yend2
								&&	outer_mesh->vertices[vidx].pnt.z >= zstart2
											   && outer_mesh->vertices[vidx].pnt.z <= zend2 )
						{
							outer_plist[n_out_polys++] = l;
							break; // no need to check the others
						}
					}
				}

				//////////////////// Crunch the Subdivision ///////////////////

				// repeat for each inner vertex
				for ( l=0 ; l<n_in_verts ; ++l )
				{
					// skip if vertex labeled as not to use
					if ( inner_mesh->vertices[inner_vlist[l]].label == 1 )
						continue;

					Point p0,p1;	// line normal to vertex
					Vector n0,n1;	// normals from vertex and opposing polygon

					p0 = inner_mesh->vertices[inner_vlist[l]].pnt;
					n0 = inner_mesh->normals[inner_vlist[l]].pnt;

					// make a long line in the direction of the normal
					p1.x = p0.x + sfac*n0.x; // add normal to vertex to get line
					p1.y = p0.y + sfac*n0.y;
					p1.z = p0.z + sfac*n0.z;

					Point minpnt;
					int ninfound, noutfound;
					float indist,outdist;
					Point pint;
					ninfound=0;
					indist=ERRVAL;

					// First find the closest inner mesh triangle to our inner vertex
					// This is used to check to see if the normal intersects the
					// inner surface before the outer one
					int v0,v1,v2;	// vertex indices for current poly
					Point t0,t1,t2; // triangle vertex points
					float d0,d1,d2; // distance between triangle vertices and current point

					for ( m=0 ; m<n_in_polys ; ++m )
					{
						// see if the line intersects the plane of the current polygon
						// skip if it does not
						v0 = inner_mesh->polygons[inner_plist[m]].vidx[0];
						v1 = inner_mesh->polygons[inner_plist[m]].vidx[1];
						v2 = inner_mesh->polygons[inner_plist[m]].vidx[2];
						t0 = inner_mesh->vertices[v0].pnt;
						t1 = inner_mesh->vertices[v1].pnt;
						t2 = inner_mesh->vertices[v2].pnt;

						// skip if any of the triangle points are our current vertex
						if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
						if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
						if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

						// skip if all vertices are out of search distance
						if ( d0>maxdsqu && d1>maxdsqu && d2>maxdsqu ) continue;

						// see if intersection point lines not within triangle then skip
						if ( !line_intersect_triangle(t0, t1, t2,p0,p1,&pint) )
							continue;

						// compare the distance to the current min inner distance
						++ninfound;
						float d = dist(pint,p0);
						if ( d==0 )
							continue;

						if ( d < indist )
							indist = d;

					} // end loop to find distance to nearest inner surface

					// test against all polygons in the outer mesh
					noutfound=0;
					outdist = ERRVAL;
					for ( m=0 ; m<n_out_polys ; ++m )
					{
						// see if the line intersects the plane of the current polygon
						// skip if it does not
						v0 = outer_mesh->polygons[outer_plist[m]].vidx[0];
						v1 = outer_mesh->polygons[outer_plist[m]].vidx[1];
						v2 = outer_mesh->polygons[outer_plist[m]].vidx[2];
						t0 = outer_mesh->vertices[v0].pnt;
						t1 = outer_mesh->vertices[v1].pnt;
						t2 = outer_mesh->vertices[v2].pnt;

						// skip if any of the triangle points are our current vertex
						if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
						if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
						if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

						// skip if all vertices are out of search distance
						if ( d0>maxd && d1>maxd && d2>maxd ) continue;

						// see if intersection point lines not within triangle then skip
						if ( !line_intersect_triangle(t0, t1, t2,p0,p1,&pint) )
							continue;

						// compare the distance to the current min inner distance
						++noutfound;
						float d = dist(pint,p0);
						if ( (d < indist) && d < maxd ) // must be closer than nearest inner polygon
						{
							if ( d < outdist ) // see if closer than last suitable distance
							{
								outdist = d;
								minpnt = pint;
							}
						}

					}

					// check against minimum distance
					if ( outdist < mind )
						outdist = mind;

					// see if a vertex found but too far away
					if ( outdist == ERRVAL && noutfound > 0 )
						outdist = maxd;

					// assign distance for vertex
					thick->nodes[inner_vlist[l]] = outdist;

					// assign vector for points connecting surfaces
					if ( outdist != ERRVAL && outdist > mind && outdist < maxd)
					{
						vertex cvert(minpnt.x,minpnt.y,minpnt.z);
						n = inner_vlist[l];
						closest_vects->assign_node(2*n, inner_mesh->vertices[n]);
						closest_vects->assign_normal(2*n,inner_mesh->normals[n]);
						closest_vects->assign_node(2*n+1, cvert);
						closest_vects->assign_normal(2*n+1, inner_mesh->normals[n]);
						closest_vects->assign_polygon(n, 2*n, 2*n+1, 0);
					}

				} // end of subdivision

			} // z search
		} // y search
	} // x search

	// clean up memory allocation
	delete [] inner_vlist;
	delete [] inner_plist;
	delete [] outer_plist;

	return 1;

}

//////////////////////// FindNormalDistThreads //////////////////////////////////

/**
Find the distance between two meshes by taking a normal from the vertex
of one mesh and intersecting the second mesh. The normal is extended
to the maximum distance to search for an intersection.
Care is taken so that the normal distance chosen does not intersect any
surfaces between end points.
This version subdivides the meshes into smaller portions so that the search will
be faster. In addition, this version assigns a separate thread to each subdivision.

@param inner_mesh - white matter mesh
@param outer_mesh - gray matter mesh
@param closest_vects - output mesh with vectors connecting nearest points
@param thick - output texture to put inner-outer distances in - one value for each inner vertex
@param nsub - number of subdivisions for each axis
@param over - amount of overlap to add to subdivisions
@param nflip - normal direction, if -1 flip normals
@param mind - min distance to search for normal intersection with other mesh
@param maxd - max distance to search for normal intersection with other mesh
@return - 1 on success
 */
int FindNormalDistThreads(RicMesh *inner_mesh, RicMesh *outer_mesh, RicMesh *closest_vects,
					RicTexture *thick, int nsub, float over, int nflip, float mind, float maxd)
{
	int i,j,k,l,m;

	// check to see that the inner mesh and the texture are the same size
	if ( inner_mesh->v_size != thick->size )
	{
		cout << "FindNormalDistThreads - error inner mesh and texture not same size"<<endl;
		return 0;
	}

	// array pointers for vertex and polygon lists for subdividing the meshes
	int **inv_mat;	// inner mesh vertices
	int **inp_mat;	// inner mesh polygons
	int **outp_mat;	// outer mesh polygons

	// number of threads
	int nthreads = nsub*nsub*nsub;

	// allocate memory for the arrays for vertices and polygons
	matrix(&inv_mat,nthreads,inner_mesh->v_size); // inner vertices
	matrix(&inp_mat,nthreads,inner_mesh->p_size); // inner polygons
	matrix(&outp_mat,nthreads,outer_mesh->p_size); // outer polygons

	// allocate thread pointers
	pthread_t *bfThreadPtr;
	bfThreadPtr = new pthread_t[nthreads];

	// allocate array of structures passing data to threads
	NormThreadStuff *pstuff;
	pstuff = new NormThreadStuff[nthreads];

	// scale factor for normal line end point - includes normal direction
	float sfac = nflip*maxd;

	// make sure that limits have been calculated for the gm mesh
	inner_mesh->calc_limits();

	// find the step size for each axis for an iteration
	float xinc,yinc,zinc;
	xinc = (inner_mesh->xmax-inner_mesh->xmin)/(float)nsub;
	yinc = (inner_mesh->ymax-inner_mesh->ymin)/(float)nsub;
	zinc = (inner_mesh->zmax-inner_mesh->zmin)/(float)nsub;

	// initialize the thickness values to ERRVAL
	for ( i=0 ; i<thick->size ; ++i ) thick->nodes[i]=ERRVAL;

	// the number of actual searches will be the cube of the subdivision number
	float xstart,ystart,zstart;	// starting values for vertex search
	float xend,yend,zend;		// ending value for vertex search
	float xstart2,ystart2,zstart2;	// starting values for vertex search
	float xend2,yend2,zend2;		// ending value for vertex search
	int n_in_verts = 0;	// number of inner mesh vertices found in a subdivision
	int n_in_polys = 0;	// number of inner mesh polygons found in a subdivision
	int n_out_polys = 0;// number of outer mesh polygons found in a subdivision
	int cur_thread;

	////// repeat for each subdivision //////
	for ( i=0 ; i<nsub ; ++i )
	{
		xstart = inner_mesh->xmin + i*xinc;
		xend = xstart+xinc;

		for ( j=0 ; j<nsub ; ++j )
		{
			ystart = inner_mesh->ymin + j*yinc;
			yend = ystart+yinc;

			for ( k=0 ; k<nsub ; ++k )
			{
				///////////////// preprocess for current subdivision ///////////
				// find all the vertices and polygons in these bounds
				cur_thread = (i*nsub*nsub)+j*nsub+k;

				n_in_verts = n_in_polys = n_out_polys = 0;

				zstart = inner_mesh->zmin + k*zinc;
				zend = zstart + zinc;

				// find the appropriate gm vertices in this box
				for ( l=0 ; l<inner_mesh->v_size ; ++l )
				{
					if ( inner_mesh->vertices[l].pnt.x >= xstart && inner_mesh->vertices[l].pnt.x <= xend
						&&	inner_mesh->vertices[l].pnt.y >= ystart && inner_mesh->vertices[l].pnt.y <= yend
						&&	inner_mesh->vertices[l].pnt.z >= zstart && inner_mesh->vertices[l].pnt.z <= zend )
					{
						inv_mat[cur_thread][n_in_verts++] = l;
					}
				}

				// allow for overlap
				xstart2 = xstart-over;
				ystart2 = ystart-over;
				zstart2 = zstart-over;
				xend2 = xend+over;
				yend2 = yend+over;
				zend2 = zend+over;

				// look for inner triangles that fit in this box
				for ( l=0 ; l<inner_mesh->p_size ; ++l )
				{
					// check each triangle vertex
					for ( m=0 ; m<3 ; ++m )
					{
						int vidx = inner_mesh->polygons[l].vidx[m];
						if ( inner_mesh->vertices[vidx].pnt.x >= xstart2
							&& inner_mesh->vertices[vidx].pnt.x <= xend2
							&&	inner_mesh->vertices[vidx].pnt.y >= ystart2
							&& inner_mesh->vertices[vidx].pnt.y <= yend2
							&&	inner_mesh->vertices[vidx].pnt.z >= zstart2
							&& inner_mesh->vertices[vidx].pnt.z <= zend )
						{
							inp_mat[cur_thread][n_in_polys++] = l;
							break; // no need to check the others
						}
					}
				}

								// look for outer triangles that fit in this box
				for ( l=0 ; l<outer_mesh->p_size ; ++l )
				{
					// check each triangle vertex
					for ( m=0 ; m<3 ; ++m )
					{
						int vidx = outer_mesh->polygons[l].vidx[m];
						if ( outer_mesh->vertices[vidx].pnt.x >= xstart2
											   && outer_mesh->vertices[vidx].pnt.x <= xend2
											   &&	outer_mesh->vertices[vidx].pnt.y >= ystart2
											   && outer_mesh->vertices[vidx].pnt.y <= yend2
											   &&	outer_mesh->vertices[vidx].pnt.z >= zstart2
											   && outer_mesh->vertices[vidx].pnt.z <= zend2 )
						{
							outp_mat[cur_thread][n_out_polys++] = l;
							break; // no need to check the others
						}
					}
				}

				//////////////////// Crunch the Subdivision ///////////////////
				// fill up structure with data to pass to thread
				pstuff[cur_thread].inner_vlist = inv_mat[cur_thread];
				pstuff[cur_thread].n_in_verts = n_in_verts;
				pstuff[cur_thread].inner_plist = inp_mat[cur_thread];
				pstuff[cur_thread].n_in_polys = n_in_polys;
				pstuff[cur_thread].outer_plist = outp_mat[cur_thread];
				pstuff[cur_thread].n_out_polys = n_out_polys;
				pstuff[cur_thread].outer_mesh = outer_mesh;
				pstuff[cur_thread].inner_mesh = inner_mesh;
				pstuff[cur_thread].closest_vects = closest_vects;
				pstuff[cur_thread].inner_thick = thick->nodes;
				pstuff[cur_thread].mind = mind;
				pstuff[cur_thread].maxd = maxd;
				pstuff[cur_thread].sfac = sfac;

				// now create a thread for this subdivision
				pthread_create(&bfThreadPtr[cur_thread],NULL, FindNormalDistThread,
										   (void *)&pstuff[cur_thread]);

			} // z search
		} // y search
	} // x search

	// now wait for all the threads to finish
	for (i=0 ; i<nthreads ; ++i)
	{
		pthread_join(bfThreadPtr[i],NULL);
		//fprintf(stderr,"Thread %d finished\n",i);
	}

	// clean up memory allocation
	free_matrix(inv_mat);
	free_matrix(inp_mat);
	free_matrix(outp_mat);
	delete [] bfThreadPtr;
	delete [] pstuff;

	return 1;

}

////////////////////////// FindNormalDistThread ///////////////////////////
/**
This is the thread spun off for each subdivision in the threaded version.
Find the distance between two meshes by taking a normal from the vertex
of one mesh and intersecting the second mesh. The normal is extended
to the maximum distance to search for an intersection.

@param tstruct - pointer to NormThreadStuff structure
@return - null pointer
 */
void *FindNormalDistThread( void *tstruct )
{
	NormThreadStuff *stuff;
	stuff = (NormThreadStuff*) tstruct;

	// square of distance for comparisons
	float maxdsqu = stuff->maxd*stuff->maxd;

	// now do the brute force search for this subdivision
	int l,m,n;

	// repeat for each inner vertex
	for ( l=0 ; l<stuff->n_in_verts ; ++l )
	{
		// skip if vertex labeled as not to use
		if ( stuff->inner_mesh->vertices[stuff->inner_vlist[l]].label == 1 )
			continue;

		Point p0,p1;	// line normal to vertex
		Vector n0;		// normal from vertex

		p0 = stuff->inner_mesh->vertices[stuff->inner_vlist[l]].pnt;
		n0 = stuff->inner_mesh->normals[stuff->inner_vlist[l]].pnt;

		// make a long line in the direction of the normal
		p1.x = p0.x + stuff->sfac*n0.x; // add normal to vertex to get line
		p1.y = p0.y + stuff->sfac*n0.y;
		p1.z = p0.z + stuff->sfac*n0.z;

		Point minpnt;
		int ninfound, noutfound;
		float indist,outdist;
		Point pint;
		ninfound=0;
		indist=ERRVAL;

		// First find the closest inner mesh triangle to our inner vertex
		// This is used to check to see if the normal intersects the
		// inner surface before the outer one
		int v0,v1,v2;	// vertex indices for current poly
		Point t0,t1,t2; // triangle vertex points
		float d0,d1,d2; // distance between triangle vertices and current point

		for ( m=0 ; m<stuff->n_in_polys ; ++m )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			v0 = stuff->inner_mesh->polygons[stuff->inner_plist[m]].vidx[0];
			v1 = stuff->inner_mesh->polygons[stuff->inner_plist[m]].vidx[1];
			v2 = stuff->inner_mesh->polygons[stuff->inner_plist[m]].vidx[2];
			t0 = stuff->inner_mesh->vertices[v0].pnt;
			t1 = stuff->inner_mesh->vertices[v1].pnt;
			t2 = stuff->inner_mesh->vertices[v2].pnt;

			// skip if any of the triangle points are our current vertex
			if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
			if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
			if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>maxdsqu && d1>maxdsqu && d2>maxdsqu ) continue;

			// see if intersection point lines not within triangle then skip
			if ( !line_intersect_triangle(t0, t1, t2,p0,p1,&pint) )
				continue;

			// compare the distance to the current min inner distance
			++ninfound;
			float d = dist(pint,p0);
			if ( d==0 )
				continue;

			if ( d < indist )
				indist = d;

		} // end loop to find distance to nearest inner surface

		// test against all polygons in the outer mesh
		noutfound=0;
		outdist = ERRVAL;
		for ( m=0 ; m<stuff->n_out_polys ; ++m )
		{
			// see if the line intersects the plane of the current polygon
			// skip if it does not
			v0 = stuff->outer_mesh->polygons[stuff->outer_plist[m]].vidx[0];
			v1 = stuff->outer_mesh->polygons[stuff->outer_plist[m]].vidx[1];
			v2 = stuff->outer_mesh->polygons[stuff->outer_plist[m]].vidx[2];
			t0 = stuff->outer_mesh->vertices[v0].pnt;
			t1 = stuff->outer_mesh->vertices[v1].pnt;
			t2 = stuff->outer_mesh->vertices[v2].pnt;

			// skip if any of the triangle points are our current vertex
			if ( (d0=distsqu(p0,t0)) < 0.1 ) continue;
			if ( (d1=distsqu(p0,t1)) < 0.1 ) continue;
			if ( (d2=distsqu(p0,t2)) < 0.1 ) continue;

			// skip if all vertices are out of search distance
			if ( d0>maxdsqu && d1>maxdsqu && d2>maxdsqu ) continue;

			// see if intersection point lines not within triangle then skip
			if ( !line_intersect_triangle(t0, t1, t2,p0,p1,&pint) )
				continue;

			// compare the distance to the current min inner distance
			++noutfound;
			float d = dist(pint,p0);
			if ( (d < indist) && d < stuff->maxd ) // must be closer than nearest inner polygon
			{
				if ( d < outdist ) // see if closer than last suitable distance
				{
					outdist = d;
					minpnt = pint;
				}
			}

		}

		// check against minimum distance
		if ( outdist < stuff->mind )
			outdist = stuff->mind;

		// see if a vertex found but too far away
		if ( outdist == ERRVAL && noutfound > 0 )
			outdist = stuff->maxd;

		// assign distance for vertex
		stuff->inner_thick[stuff->inner_vlist[l]] = outdist;

		// assign vector for points connecting surfaces
		if ( outdist != ERRVAL && outdist > stuff->mind && outdist < stuff->maxd)
		{
			vertex cvert(minpnt.x,minpnt.y,minpnt.z);
			n = stuff->inner_vlist[l];
			stuff->closest_vects->assign_node(2*n, stuff->inner_mesh->vertices[n]);
			stuff->closest_vects->assign_normal(2*n,stuff->inner_mesh->normals[n]);
			stuff->closest_vects->assign_node(2*n+1, cvert);
			stuff->closest_vects->assign_normal(2*n+1, stuff->inner_mesh->normals[n]);
			stuff->closest_vects->assign_polygon(n, 2*n, 2*n+1, 0);
		}

	} // end of subdivision

	return NULL;
}

