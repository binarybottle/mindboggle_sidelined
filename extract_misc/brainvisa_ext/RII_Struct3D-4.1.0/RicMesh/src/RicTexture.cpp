// -------------------------- RicTexture.cpp ---------------------------------
/*! \file 
 Implementation file of the RicTexture class
 This was derived from the texture class developed by Tom Arnow but was 
 renamed in order to avoid confusion with the rest of the world.
 What a texture really is, is a per node value aligned with the vertices
 of a mesh. The value can be any useful parameter such a curvature, 
 thickness, etc.
 Bill Rogers - November 2008
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include "RicUtil.h"
#include "RicTexture.h"
#include "RicTextureSet.h"

using namespace std;

//////////////////////////// constructors //////////////////////////////
/*!
 Constructor for empty texture
 */
RicTexture::RicTexture(void)
{
	size=0;
	nodes = NULL;
	t_step = 0;
	avg = med = std_dev = min = max = 0;
}

/*!
 Constructor to allocate memory for a texture
 
 @param numnodes
 number of elements in the texture - generally one per mesh vertex
 */
RicTexture::RicTexture(int numnodes)
{
	size=numnodes;
	nodes=new float[size];
	t_step = 0;
	avg = med = std_dev = min = max = 0;
}

/*!
 Destructor - free memory
 */
RicTexture::~RicTexture()
{
	if (nodes)
		delete nodes;
}

///////////////////////////// member functions ///////////////////////////////

/*!
 Initializes the texture array to passed number of nodes
 
 @param nnodes
 number of elements in the texture - generally one per mesh vertex
 */
void RicTexture::init(int nnodes)
{
	if (nodes)
		delete [] nodes;
	size = nnodes;
	nodes=new float[size];
	avg = med = std_dev = min = max = 0;
}

/*!
 Sets a value for a texture node
 
@param index
index into texture node array

@param value
value to assign texture node		
 */
void RicTexture::set_node(int index, float value)
{
	nodes[index] = value;
}

/*!
 Calculates the min, max, median, standard deviation and average 
 values of the texture array.
 
@returns
returns the average value
 */
float RicTexture::CalcMinMaxAvg(void)
{
	min = 10000000;
	max = -10000000;
	avg = 0;
	long cnt=0;
	long nzero=0;
	double std_t=0, avg_t=0;
	float *medarray = new float[size]; // array to calculate median

	for (long i=0; i<size; ++i)
	{
		if (nodes[i] == ERRVAL) // skip invalid values
		{
			++nzero;
			continue; // skip empty ones
		}
		if (nodes[i] > max)
			max = nodes[i];
		if (nodes[i] < min)
			min = nodes[i];
		avg_t += nodes[i];
		std_t += nodes[i]*nodes[i];

		// put in array to calculate median
		medarray[cnt] = nodes[i];

		++cnt;
	}

	// average
	avg = avg_t/(double)cnt;
	
	// standard deviation
	std_dev= sqrt(fabs(std_t-avg_t*avg_t/(double)cnt)/((double)cnt-1.0));
	
	// median
	float_sort(medarray, cnt);
	med = medarray[cnt/2];

	delete [] medarray;

	return avg;
}

/*!
This function reads a single texture from a file. It uses the read function
from RicTextureSet and copies the first time point to the node array.

@param f_name
name of file to read texture from

@returns
1 on success or 0 if failed to read texture
*/
int RicTexture::read_texture(char *f_name)
{
	// create a one texture set from file
	RicTextureSet Tset(f_name);
	
	// check to see if texture read ok
	if (Tset.ntstep == 0) 
		return 0;

	// initialize the mesh to the size read in
	init(Tset.tex[0].size);

	// copy to the first mesh
	int i;
	for (i=0; i<size; i++)
		nodes[i] = Tset.tex[0].nodes[i];

	return 1;
}

/*!
This function writes the texture to a ascii file.

@param fname
file name to write texture to

@returns
1 on success or 0 if write failed
*/
int RicTexture::write_texture(char *fname)
{
	// open up the file
	ofstream fout(fname);

	// check to see if opens ok
	if ( fout.fail() )
		return 0;

	// lets make the header                 
	fout<<"ascii"<<endl<<"FLOAT"<<endl<<"1"<<endl;
	fout << t_step << " " << size << " " << endl;

	// write out the nodes
	for (int i=0; i<size; i ++)
		fout<<nodes[i]<<" ";
	fout << endl;

	// close the file
	fout.flush();
	fout.close();

	return 1;
}

