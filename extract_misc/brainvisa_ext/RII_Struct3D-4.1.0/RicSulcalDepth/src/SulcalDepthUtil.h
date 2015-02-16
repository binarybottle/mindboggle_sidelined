// -------------------------- SulcalDepthUtil.h --------------------------------
/*! \file
This file declarations for routines for RicSulcalDepth.cpp
Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA
*/

#include <RicUtil.h>
#include <RicMesh.h>

int PlaneIntersectMesh(RicMesh *mesh, Point *ipnts, Point p1, Point p2, Point p3);
Point TranRotPntsXZ(Point *pnts, int ni, float *ang1, float *ang2, float *ang3);
Point TranRotPntsXZ2(Point *pnts, int ni, float *ang1, float *ang2, float *ang3);
float SulcalDepthMeasures(Point *pnts, int n, int ndepth, Point *dpnts);
float SulcalDepthMeasures0(Point *pnts, int n, int ndepth, Point *dpnts,
		Point bpnt, Point tpnt);
