// ------------------------- TexFill.h ---------------------------------
/**
\file
Header file for routines to fill in voids in textures due to culling or bad values

Copyright (C) 2008 by Bill Rogers, Research Imaging Center, UTHSCSA
rogers@uthscsa.edu
*/

#ifndef _TEX_FILL_H
#define _TEX_FILL_H
using namespace std;

#include "RicMesh.h"
#include "RicTexture.h"

int TexFillAvg(RicTexture *tex, RicMesh *mesh, float maxd);

#endif // _MESH_FILL_H
