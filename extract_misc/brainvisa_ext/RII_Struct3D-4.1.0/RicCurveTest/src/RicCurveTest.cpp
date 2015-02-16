/***************************************************************************
 *   Copyright (C) 2007 by Bill Rogers                                     *
 *   rogers@uthscsa.edu                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "RicCurve.h"
#include "RicUtil.h"

using namespace std;

int main(int argc, char *argv[])
{
  float coef[] = {0.0,  0.0,   0,
	  1,  0, 0.5,
	  1,  1,   1,
	  1,  1,   3,
	  0,  0,   2,
	  1,  0, 2.5,
	  0,  1, 1.5,
	  0,  1, 3.5,
	  0,  0,   4,
	  1,  0, 6.5};

	Point *ptr;
	
	ptr = (Point*)coef;
	
	  cout << "RicCurve Class Creator!" << endl;
	  
	  // create a test curve
	  RicCurve rcurve(10,ptr);
	  
	  // interpolate a bunch of points
	  rcurve.Interpolate(51);
	  
	  // write to a file
	  rcurve.WriteCurveGo("RicCurveGo.g2");
	  rcurve.WriteCurveAsMesh("RicCurveMesh.mesh");
	  
	  // sort the data
//	  rcurve.SortAlongLength();
	  rcurve.Interpolate(51);
	  
	  double len = rcurve.CurveLength();
	  cout << "Curve Length " << len << endl;
	  
	  float d1 = ddist(rcurve.pnts[0],rcurve.pnts[1]);
	  float d2 = ddist(rcurve.pnts[1],rcurve.pnts[2]);
	  float d3 = ddist(rcurve.pnts[2],rcurve.pnts[3]);
	  float e1 = ddist(rcurve.pnts[50],rcurve.pnts[49]);
	  float e2 = ddist(rcurve.pnts[49],rcurve.pnts[48]);
	  float e3 = ddist(rcurve.pnts[48],rcurve.pnts[47]);
	  
	  cout << "Beginning Distances " << d1 << d2 << d3 <<endl;
	  cout << "Ending Distances " << e1 << e2 << e3 <<endl;
	  
	  // write a second file
	  rcurve.WriteCurveAsMesh("RicCurveMesh2.mesh");
	  
	  cout << "All done, sir" << endl;
	  return EXIT_SUCCESS;
}
