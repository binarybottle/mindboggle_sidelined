// -------------------------- RicTextureSet.cpp ---------------------------------
/*! \file
Implementation file of the RicTextureSet class
This class encapsulates a set of RicTextures for a series of time steps.This
mainly consists of initialization and reading and writing of RicTextures.

Bill Rogers - November 2008
 */


#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "RicUtil.h"
#include "RicTextureSet.h"

using namespace std;

//////////////////////////// constructors //////////////////////////////
/*!
Constructor for empty texture
 */
RicTextureSet::RicTextureSet(void)
{
	size=0;
	tex = NULL;
	ntstep = 1;
}

/*!
@brief
Constructor to allocate memory for a texture

@param numnodes
number of nodes generally one per mesh vertex

@param numtstep
number of time steps which is the number of textures in texture set
 */
RicTextureSet::RicTextureSet(int numnodes, int numtstep)
{
	size = numnodes;
	ntstep = numtstep;

	// allocate memory for textures
	tex=new RicTexture[ntstep];
	for ( int i=0 ; i<ntstep ; ++i)
		tex[i].init(size);
}

/*!
@brief
Constructor to allocate memory for a texture

@param filename
Name of file containing one or more textures
 */
RicTextureSet::RicTextureSet(string filename)
{
	size=0;
	tex = NULL;
	ntstep = 0;
	Read(filename);
}

/*!
Destructor - free memory
 */
RicTextureSet::~RicTextureSet()
{
	// clean up memory
	if ( tex ) delete [] tex;
}

///////////////////////////// member functions ///////////////////////////////

/*!
Initializes the texture array to passed number of nodes. If memory is already
allocated for the texture it is first freed.

@param numnodes
number of nodes generally one per mesh vertex

@param numtstep
number of time steps which is the number of textures in texture set
*/
void RicTextureSet::init(int numnodes, int numtstep)
{
	// delete textures if already there
	if ( tex ) delete tex;

	size = numnodes;
	ntstep = numtstep;

	// allocate a new set of textures
	tex=new RicTexture[ntstep];
	for ( int i=0 ; i<ntstep ; ++i)
		tex[i].init(size);
}

/*!
Read a texture or set of textures from the passed file name. POINT2DF is not
supported. Data types of FLOAT, S15, and U32 are supported. Both ascii and
binary binarDCBA file types are supported.

@param filename
Name of file containing one or more textures

@returns
1 on success, 0 if fails to read the file
 */
int RicTextureSet::Read(string filename)
{
	int i;
	string dtype;		// data type - FLOAT S15 U32;
	string openMode;	// file mode - ascii binary
	string tval;
	char param[8];

	// open our file
	ifstream fin;
	fin.open(filename.c_str(),ifstream::in);
	if (!fin.is_open())
	{
		cerr<<"Couldn't open file "<<filename<<" exiting"<<endl;
		return(0);
	}
	if ( !fin.good() )
		cerr << "Error" << endl;

	// look at header to determine file type
	fin >> openMode;

	// ok now, let us read in an ascii file
	if ( openMode.compare(0,5,"ascii") == 0 )
	{
		// get unit type
		fin >> dtype;
		if ( dtype == "POINT2DF" )
		{
			cerr<<"We do not read texture files of type POINT2DF"<<endl;
			return(0);
		}

		// get number of textures in file
		fin >> this->ntstep;

		// reallocate memory for textures
		this->init(10, ntstep);

		// now read the textures
		for ( int n=0 ; n<this->ntstep ; ++n )
		{
			// read current mode/time step
			fin>>tex[n].t_step;

			// get the size of each texture - only first counts
			fin >> tex[n].size;

			// allocate memory for texture
			tex[n].init(tex[n].size);

			// read in all the values
			for ( i=0 ; i<tex[n].size ; ++i )
			{
				fin >> tval;
				tex[n].nodes[i] = atof(tval.c_str());
			}

			// calc stats
			tex[n].CalcMinMaxAvg();
		}
		fin.close();
	}
	else if ( openMode.compare(0,9,"binarDCBA") == 0 ) // read in a DCBA binary file
	{
		// close file and reopen as binary
		fin.close();

		ifstream bfin;

		bfin.open(filename.c_str(),ifstream::binary|ifstream::in);
		if (!bfin.is_open())
		{
			cerr<<"Couldn't open file "<<filename<<" exiting"<<endl;
			return(0);
		}

		if ( !bfin.good() )
			cerr << "Error" << endl;

		// get the data type
		bfin.seekg(13,ifstream::beg);
		bfin.get(param,6);

		if ( !bfin.good() )
			cerr << "Error" << endl;

		param[3]=0;
		dtype = param;

		// the offset to the next field will depend on length of data type string
		if ( dtype.compare(0,3,"FLO") == 0 )
			bfin.seekg(18);
		else
			bfin.seekg(16);

		// get the number of textures
		bfin.read((char*)&this->ntstep,sizeof(this->ntstep));

		// reallocate memory for textures
		this->init(10, ntstep);

		// now read the textures
		for ( int n=0 ; n<ntstep ; ++n )
		{
			// get the current texture number
			unsigned ts;
			bfin.read((char*)&ts,sizeof(ts));
			tex[n].t_step = ts;

			// read the size
			bfin.read((char*)&tex[n].size,4);

			// allocate memory
			tex[n].init(tex[n].size);

			// reading depends on data type
			if ( dtype.compare(0,3,"S16") == 0 ) // 16 bit integer
			{
				// now read in the node values
				short stmp;
				for ( i=0 ; i<tex[n].size ; ++i )
				{
					bfin.read((char*)&stmp,2);
					tex[n].nodes[i] = stmp;
				}
			}
			else if ( dtype.compare(0,3,"U32") == 0 ) // 32 bit unsigned integer
			{
				// now read in the node values
				unsigned int utmp;
				for ( i=0 ; i<tex[n].size ; ++i )
				{
					bfin.read((char*)&utmp,4);
					tex[n].nodes[i] = utmp;
				}
			}
			else // must be float
			{
				// now read in the node values
				for ( i=0 ; i<tex[n].size ; ++i )
				{
					bfin.read((char*)&tex[n].nodes[i],4);
				}
			}

			// calc stats
			tex[n].CalcMinMaxAvg();
		}
		bfin.close();
	}
	else
	{
		cerr<<"Cannot read texture file mode " << openMode <<endl;
		return(0);
	}

	// assign size from first texture
	this->size = tex[0].size;



	return 1; // ok if we made it here
}

/*!
Writes the texture to the passed file name as an ascii file.

@param filename
File name to write one or more textures to

@returns
1 on success, 0 if fails to write file
 */
int RicTextureSet::Write_Ascii(string filename)
{
	// open up the file
	ofstream fout(filename.c_str());

	// check to see if file opens ok
	if ( fout.fail() )
		return 0;

	// lets make the header
	fout<<"ascii"<<endl<<"FLOAT"<<endl<<ntstep<<endl;

	// repeat for every time step

	for ( int n=0 ; n<ntstep ; ++n)
	{
		fout << n << " " << tex[n].size << " " << endl;

		// write out the nodes
		for (int i=0;i<tex[n].size;  i ++)
			fout<<tex[n].nodes[i]<<" ";
		fout << endl;
	}

	// close the file
	fout.flush();
	fout.close();

	return 1; // ok if it got here
}

/*!
Writes the texture to the passed file name as a binary file. This version
uses the binarDCBA format.

@param filename
File name to write one or more textures to

@returns
1 on success, 0 if fails to write file
 */
int RicTextureSet::Write(string filename)
{
	// open up the file
	ofstream fout(filename.c_str());

	// check to see if file opens ok
	if ( fout.fail() )
	{
		cerr<<"RicTextureSet - file open error"<<endl;
		return 0;
	}

	// file format
	fout<<"binarDCBA";

	// texture type
	int A = 5; // number of characters in float
	fout.write((char*) &A,sizeof (A));
	fout<<"FLOAT";

	// number of time steps
	fout.write((char*)&ntstep, sizeof (ntstep));

	// repeat for every time step

	for ( int n=0 ; n<ntstep ; ++n)
	{
		// time step
		unsigned tstep = n;
		fout.write((char*)&tstep, sizeof (tstep));

		// number of nodes in time step
		fout.write((char*) &tex[n].size, sizeof(tex[n].size) );

		// write out the nodes
		for (int i=0;i<tex[n].size;  i ++)
			fout.write((char*)&tex[n].nodes[i], sizeof(float));
	}

	// close the file
	fout.flush();
	fout.close();

	return 1; // ok if it got here
}


