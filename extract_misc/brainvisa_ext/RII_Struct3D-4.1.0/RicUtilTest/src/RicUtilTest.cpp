//============================================================================
//   Copyright (C) 2007 by Bill Rogers  - RII-UTHSCSA
//   rogers@uthscsa.edu
//============================================================================

#include <iostream>
#include "RicUtil.h"
using namespace std;

int main()
{
	Point p;
	p.x = 2;
	Cylind c;
	c = cart_to_cyl(p);
	IPoint ip1(2,3,4),ip2(2,2,2),ip3;
	ip3 = ip1-ip2;
	DPoint dp;
	dp.x = 2.24;
	
	Point p1(0,0,0),p2(1,1,0),p3(1,-1,0.1);
	Point l1(0.6,-0.6,-1),l2(0.5,0.5,1);
	Point pout;
	
	int stat=line_thru_triangle(p1,p2,p3,l1,l2,&pout);
	
	stat=line_intersect_triangle(p1,p2,p3,l1,l2,&pout);

	
	cout << "RicUtil Test" << endl;
	return 0;
}
