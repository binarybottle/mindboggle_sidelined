// --------------------------  TexFill.cpp -----------------------------------
/*!
\file
This file contains routines fill voids in meshes due to missing or culled vertices.

Copyright (C) 2007 by Bill Rogers, Research Imaging Center, UTHSCSA   *
rogers@uthscsa.edu
 */

#include "RicMesh.h"
#include "RicTexture.h"
#include "TexFill.h"
#include <math.h>
#include "RicUtil.h"

/*!
 * This function fills in empty nodes in the texture with values averaged
 * from surrounding nodes. All valid nodes within a specificed distance
 * of the current empty node are averaged to get the node value.
 *
 * @param tex - texture to fill
 * @param mesh - mesh corresponding to tex
 * @param maxd - max distance to look for averaging
 * @return 1 if all mesh values filled
 */
int TexFillAvg(RicTexture *tex, RicMesh *mesh, float maxd)
{

	// sanity check to see if texture and mesh are the same size
	if ( tex->size != mesh->v_size )
	{
		cerr << "TexFillAvg - texture and mesh not same size "<<endl;
		return 0;
	}

	// create new mesh for filled values
	RicTexture *NewTex = new RicTexture(tex->size);

	// use distance squared for comparisons for efficiency
	float maxd2 = maxd*maxd;

	// look at every texture value for errors
	int i,numnotfix=0;
	for ( i=0 ; i<tex->size ; ++i )
	{
		// just copy to new mesh if ok and move on to next node
		if ( tex->nodes[i] != ERRVAL )
		{
			NewTex->nodes[i] = tex->nodes[i];
			continue;
		}

		// now check all surrounding nodes for valid values
		int nfound=0;	// number of surrounding values found
		float avg=0;	// average of surrounding valid values

		for ( int j=0 ; j<tex->size ; ++j )
		{
			// skip invalid entries
			if ( tex->nodes[j] == ERRVAL ) continue;

			if ( dist_squ(mesh->vertices[i].pnt,mesh->vertices[j].pnt) < maxd2 )
			{
				avg += tex->nodes[j];
				++nfound;
			}
		}

		// there must be at least one for the average
		if ( nfound == 0 )
		{
			NewTex->nodes[i] = ERRVAL;
			++numnotfix;
			continue;
		}

		// woo hoo - assign the average to texture
		NewTex->nodes[i] = avg/nfound;
		nfound = 0;
	}

	// copy back to original texture
	for ( i=0 ; i<tex->size ; ++i )
		tex->nodes[i] = NewTex->nodes[i];
	delete NewTex;

	// if all the nodes were filled then return 1 indicating success
	if ( numnotfix == 0 ) return 1;
	else return 0;
}
