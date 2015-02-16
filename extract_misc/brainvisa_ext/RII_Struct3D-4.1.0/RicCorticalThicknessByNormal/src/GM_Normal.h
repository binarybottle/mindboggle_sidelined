// ------------------------- GM_Normal.h ---------------------------------
/**
Header file for routines to find the distance between two meshes by
intersecting normals from one mesh with triangles of the other.

Copyright (C) 2007 by Bill Rogers, Research Imaging Center, UTHSCSA
rogers@uthscsa.edu
*/

#ifndef _GM_NORMAL_H
#define _GM_NORMAL_H
using namespace std;

#include "RicMesh.h"
#include "RicTexture.h"

int FindNormalDistBrute(RicMesh *gm_mesh, RicMesh *wm_mesh, RicMesh *closest_vects,
				RicTexture *thick, int nflip, float mind, float maxd);
int FindNormalDistSubdivide(RicMesh *gm_mesh, RicMesh *wm_mesh, RicMesh *closest_vects,
				RicTexture *tex, int nsub, float over, int nflip, float mind, float maxd);
int FindNormalDistThreads(RicMesh *gm_mesh, RicMesh *wm_mesh, RicMesh *closest_vects,
				RicTexture *tex, int nsub, float over, int nflip, float mind, float maxd);

#endif // _GM_NORMAL_H
