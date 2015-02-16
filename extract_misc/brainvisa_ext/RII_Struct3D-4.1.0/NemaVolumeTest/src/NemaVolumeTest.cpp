/***************************************************************************
 *   Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA *
 *   rogers@uthscsa.edu                                                    *
 *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "nemardr.h"

using namespace std;

int main(int argc, char *argv[])
{
	cout << "Reading file gray.des" << endl;
	
	NEMAFILE *nfile = new NEMAFILE("gray.des");
	if ( nfile->tot_scans() == 0 )
	{
		cout << "Oops, cannot read file" << endl;
		exit(1);
	}
	
	cout << "Read file ok" << endl;
	
	// Now write the file out
	nfile->set_desfname("NemaWrite.des");
	
	delete nfile;
	
	return EXIT_SUCCESS;
}
