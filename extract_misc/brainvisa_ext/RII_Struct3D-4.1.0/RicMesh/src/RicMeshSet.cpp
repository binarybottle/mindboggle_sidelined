// ------------------------- RicMeshSet.cpp ----------------------------------
/*! \file
Implementation file for the RicMeshSet class. The class can open and write
binary and ascii BrainVisa/Anatomist file formats. This version will
handle polygons of dimension 2 and 3 - lines and triangles.
 */

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include "RicMeshSet.h"
#include <algorithm> //required for std::swap

using namespace std;

/*!
 * This function swaps order of the bytes of a floating point number
 * @param fNumber - the float to byte swap
 * @return - the swapped float
 */
float ByteSwap (float fNumber)
{
	union
	{
	float f;
	unsigned long l;
	}u;

	u.f = fNumber;
	u.l = (((u.l&0x000000FF)<<24)+((u.l&0x0000FF00)<<8)+
			   ((u.l&0x00FF0000)>>8)+((u.l&0xFF000000)>>24));
   return u.f;
}

/*!
 * This function swaps order of the bytes of a 4 byte integer
 * @param nLongNumber - the integer to byte swap
 * @return - the swapped integer
 */
int ByteSwap (int nLongNumber)
{

   return (int)((((unsigned)nLongNumber&0x000000FF)<<24)+(((unsigned)nLongNumber&0x0000FF00)<<8)+
   (((unsigned)nLongNumber&0x00FF0000)>>8)+(((unsigned)nLongNumber&0xFF000000)>>24));
}

/*!
 * This function swaps order of the bytes of a 4 byte long integer
 * @param nLongNumber - the long integer to byte swap
 * @return - the swapped integer
 */
unsigned long ByteSwap (unsigned long nLongNumber)
{
   return (((nLongNumber&0x000000FF)<<24)+((nLongNumber&0x0000FF00)<<8)+
   ((nLongNumber&0x00FF0000)>>8)+((nLongNumber&0xFF000000)>>24));
}


/*!
 * This function swaps order of the bytes of a short integer
 * @param nLongNumber - the short integer to byte swap
 * @return - the swapped short integer
 */
unsigned short ByteSwap (unsigned short nValue)
{
   return (((nValue>> 8)) | (nValue << 8));

}

///////////////////////////////// constructors //////////////////////////////

/*!
Construct an empty mesh setting values to zero and the mesh pointer to null.
 */
RicMeshSet::RicMeshSet (void)
{
	v_size = n_size = p_size = p_dim = ntstep = 0;
	mesh = NULL;
}

/*!
 * Construct a mesh to size.
 *
 * @param vs - number of vertices
 * @param ns - number of normals
 * @param ps - number of polygons
 * @param pd - polygon dimension (must be 2 or 3)
 * @param nt - number of time steps (1 mesh per step) or the number of meshes
 */
RicMeshSet::RicMeshSet (int vs, int ns, int ps ,int pd, int nt)
{
	mesh = NULL;

	if ( this->Init(vs, ns, ps, pd, nt) == 0 )
	{
		// allocation error - make sure mesh is null and sizes set to zero
		if ( mesh != 0 )delete [] mesh;
		mesh = NULL;
		v_size = n_size = p_size = p_dim = ntstep = 0;
	}

}

/*!
 * Construct mesh from file.
 * @param fname - name of file to read mesh from
 */
RicMeshSet::RicMeshSet(string fname)
{
	v_size = n_size = p_size = p_dim = ntstep = 0;
	mesh = NULL;

	struct stat file_size;
	int size=0;
	if(stat(fname.c_str(), &file_size)==0)
	{
		size=file_size.st_size;
	}
	else
	{
		cerr<<"No such file, exiting"<<endl;
		return;
	}

	// time to open this file
	ifstream meshfile(fname.c_str(), ios::in | ios::binary);

	// read the header info for file type
	char type[6];
	meshfile.read(type,5);
	meshfile.close();

	//  check to see if ascii or binary
	type[5] = 0; // terminate string
	if ( strcmp(type,"ascii") == 0 )
		this->ReadText(fname);
	else
		this->ReadBinary(fname);

}

/*!
Destructor - clean up memory by deleting the mesh
 */
RicMeshSet::~RicMeshSet()
{
	if ( mesh ) delete [] mesh;
}

/*!
 * Initialize mesh to size. If memory has been allocated for the mesh then
 * the existing arrays will be deleted first.
 *
 * @param vs - number of vertices
 * @param ns - number of normals
 * @param ps - number of polygons
 * @param pd - polygon dimension (must be 2 or 3)
 * @param nt - number of time steps (1 mesh per step) or the number of meshes
 * @return - 1 on success, 0 on failure
 */
int RicMeshSet::Init(int vs, int ns, int ps, int pd, int nt)
{
	// assign to member variables
	v_size = vs;
	n_size = ns;
	p_size = ps;
	p_dim = pd;
	ntstep = nt;

	// delete the meshes if allocated
	if ( mesh ) delete [] mesh;

	// create an array of meshes
	mesh = new RicMesh[ntstep] ;
	if ( !mesh )
	{
		cerr<<"RicMeshSet - allocation error"<<endl;
		return 0;
	}

	// allocate memory for arrays
	for ( int i=0 ; i<ntstep ; ++i )
	{
		if ( mesh[i].Init(v_size, n_size, p_size, p_dim) == 0 )
		{
			cerr<<"RicMeshSet - mesh allocation error"<<endl;
			return 0;
		}
	}

	return 1;
}

/*!
This function writes a mesh to a  file

@param file_out - name of file to write
@return - returns 1 on success, 0 on failure
 */
int RicMeshSet::Write(string file_out)
{
	// default write to binary
	return this->WriteBinary(file_out);
}

/*!
This function writes a mesh to a text file

@param file_out - name of file to write
@return - returns 1 on success, 0 on failure
 */
int RicMeshSet::WriteText(string file_out)
{
	ofstream fout(file_out.c_str());
	if ( fout.fail() )
	{
		cerr<<"RicMeshSet - file open error"<<endl;
		return 0;
	}

	fout<<"ascii"<<endl;
	fout<<"VOID"<<endl;
	fout<<p_dim<<endl;
	fout<<ntstep<<endl;

	for ( int t=0 ; t<ntstep ; ++t )
	{
		// output time step
		fout<<t<<endl;

		// output polygons
		fout<<mesh[t].v_size<<endl<<flush;
		for (int i=0;i<mesh[t].v_size;i++)
			fout<<"("<<mesh[t].vertices[i].pnt.x<<" ,"<<mesh[t].vertices[i].pnt.y
			<<" ,"<<mesh[t].vertices[i].pnt.z<<") ";


		// output normals
		fout<<endl<<mesh[t].n_size<<endl;
		for (int i=0;i<mesh[t].n_size;i++)
			fout<<"("<<mesh[t].normals[i].pnt.x<<" ,"<<mesh[t].normals[i].pnt.y
			<<" ,"<<mesh[t].normals[i].pnt.z<<") ";
		fout<<endl<<0<<endl<<mesh[t].p_size<<endl;

		// output polygons - only 2 or 3 dim supported
		if ( mesh[t].p_dim == 2 ) // lines
		{
			for (int i=0;i<mesh[t].p_size;i++)
				fout<<"("<<mesh[t].polygons[i].vidx[0]<<" ,"
				<<mesh[t].polygons[i].vidx[1]<<") ";
		}
		else // p_dim is 3
		{
			for (int i=0;i<mesh[t].p_size;i++)
				fout<<"("<<mesh[t].polygons[i].vidx[0]<<" ,"
				<<mesh[t].polygons[i].vidx[1]<<" ,"
				<<mesh[t].polygons[i].vidx[2]<<") ";
		}

		fout << endl;
	}
	fout.close();
	return 1;
}

/*!
This function writes a mesh to a binary file

@param file_out - name of file to write
@return - returns 1 on success, 0 on failure
 */
int RicMeshSet::WriteBinary(string file_out)
{
	ofstream fout(file_out.c_str(), ios::out|ios::binary);
	if ( fout.fail() )
	{
		cerr<<"RicMeshSet - file open error"<<endl;
		return 0;
	}

	fout<<"binarDCBA";
	int A=4	; // texture type - not used
	fout.write((char*) &A,sizeof (A));
	fout<<"VOID";


	A=p_dim;  fout.write((char*) &A,sizeof (A)); // polygon dimension
	A=ntstep;  fout.write((char*) &A,sizeof (A)); // number of time steps

	// write a mesh for each time step
	for (int t=0; t<ntstep; ++t)
	{
		A=t;
		fout.write((char*) &A, sizeof (A)); // current time step

		A=mesh[t].v_size;
		fout.write((char*) &A, sizeof (A));

		for ( int i=0; i<mesh[t].v_size; i++)
		{
			Point tmp= mesh[t].vertices[i].pnt;
			fout.write( (char*) &tmp.x, sizeof(float));
			fout.write( (char*) &tmp.y, sizeof(float));
			fout.write( (char*) &tmp.z, sizeof(float));
		}

		// write polygons
		A=mesh[t].n_size;
		fout.write((char*) &A, sizeof (A));
		for ( int i=0; i<mesh[t].n_size; i++)
		{
			Point tmp=mesh[t].normals[i].pnt;
			fout.write((char*) &tmp.x, sizeof(float));
			fout.write((char*) &tmp.y, sizeof(float));
			fout.write((char*) &tmp.z, sizeof(float));
		}

		// empty texture vector
		A=0;
		fout.write((char*) &A, sizeof (A));

		// write out polygon indices
		A=mesh[t].p_size;
		fout.write((char*) &A, sizeof (A));

		for ( int i=0; i<mesh[t].p_size; i++)
			for (int j=0; j<mesh[t].p_dim; j++)
				fout.write((char*) &mesh[t].polygons[i].vidx[j], sizeof(int));
	}

	fout.close();
	return 1;
}

/*!
This function reads a mesh from a file.

@param file_in - name of file to read mesh from
@return - 1 is returned on success, 0 on failure
 */
int RicMeshSet::Read(string file_in)
{

	struct stat file_size;
	int size=0;
	if(stat(file_in.c_str(), &file_size)==0)
	{
		size=file_size.st_size;
	}
	else
	{
		cerr<<"RicMeshSet - No such file, exiting"<<file_in<<endl;
		return 0;
	}

	// time to open this file
	ifstream meshfile(file_in.c_str(), ios::in | ios::binary);

	// read the header info for file type
	char type[6];
	meshfile.read(type,5);
	meshfile.close();

	//  check to see if ascii or binary
	type[5] = 0; // terminate string
	if ( strcmp(type,"ascii") == 0 )
		this->ReadText(file_in);
	else
		this->ReadBinary(file_in);

	return 1;
}

/*!
This function reads a mesh from a text file

@param file_in - name of file to read mesh from
@return - 1 is returned on success, 0 on failure
 */
int RicMeshSet::ReadText(string file_in)
{
	int i,j;
	v_size=0;n_size=0;p_size=0;

	// time to open this file
	ifstream meshfile(file_in.c_str(),ifstream::in);
	if ( meshfile.fail() )
	{
		cerr<<"RicMeshSet - file open error"<<endl;
		return 0;
	}

	// get the header info
	string type, junk;
	meshfile >> type;
	if ( type != "ascii" ) return 0;	// must be ascii
	meshfile >> junk;	// should be "VOID"
	meshfile >> p_dim;  // polygon dimension
	if(p_dim!=2 && p_dim!=3) {cerr<<"Error - polygon dimension must be 2 or 3"<<endl;}
	meshfile >> ntstep;

	// create a new set of meshes
	if ( mesh ) delete [] mesh;
	mesh = new RicMesh[ntstep];

	// read a mesh for each time step
	for (int t=0; t<ntstep; ++t)
	{
		// assign polygon dim for individual mesh
		mesh[t].p_dim = mesh[t].p_dim;

		// get the time step for this mesh
		meshfile >> mesh[t].t_step;

		// read vertices
		meshfile >> mesh[t].v_size;

		// now get all of those vertices using ')' as delimiter for end of vertex
		char tvertex[100];
		float x, y, z;
		this->mesh[t].init_vertices(mesh[t].v_size);
		for (i=0; i<mesh[t].v_size; ++i)
		{

			meshfile.getline(tvertex, 99, '('); // read to start of vertex
			meshfile.getline(tvertex, 99, ')'); // now read to end of vertex

			// remove commas
			j = tvertex[99] = 0;
			do
			{
				if (tvertex[j] == ',')
					tvertex[j] = ' ';
				++j;
			} while (tvertex[j] !=0);

			// parse the vector
			sscanf(tvertex, "%f %f %f", &x, &y, &z);

			mesh[t].assign_node(i, x, y, z);
		}

		// read normals
		meshfile >> mesh[t].n_size;
		if ( (mesh[t].n_size != mesh[t].v_size) && (mesh[t].n_size != 0))
		{
			cerr<<"Mesh size and normals size are not equal exiting"<<endl;
			exit(1);
		}

		if (mesh[t].n_size > 0)
			mesh[t].init_normals(mesh[t].n_size);

		for (i=0; i<n_size; ++i)
		{

			meshfile.getline(tvertex, 99, '('); // read to start of vertex
			meshfile.getline(tvertex, 99, ')'); // now read to end of vertex

			// remove commas
			j = tvertex[99] = 0;
			do
			{
				if (tvertex[j] == ',')
					tvertex[j] = ' ';
				++j;
			} while (tvertex[j] !=0);

			// parse the vector
			sscanf(tvertex, "%f %f %f", &x, &y, &z);

			mesh[t].assign_normal(i, x, y, z);
		}

		// read polygon indices
		meshfile >> mesh[t].p_size; // dummy read for texture (see file format doc)
		meshfile >> mesh[t].p_size; // read it one more time ??

		if (mesh[t].p_size > 0)
			mesh[t].init_polygons(mesh[t].p_size);

		int id1, id2, id3;
		for (i=0; i<mesh[t].p_size; ++i)
		{

			meshfile.getline(tvertex, 99, '('); // read to start of vertex
			meshfile.getline(tvertex, 99, ')'); // now read to end of vertex

			// remove commas
			j = tvertex[99] = 0;
			do
			{
				if (tvertex[j] == ',')
					tvertex[j] = ' ';
				++j;
			} while (tvertex[j] !=0);

			// parse the vector
			if (mesh[t].p_dim == 2)
			{
				id3 = 0;
				sscanf(tvertex, "%d %d", &id1, &id2);
				mesh[t].assign_polygon(i, id1, id2, id3);
			}
			else // p_dim == 3
			{
				sscanf(tvertex, "%d %d %d", &id1, &id2, &id3);
				mesh[t].assign_polygon(i, id1, id2, id3);
			}
		}

		mesh[t].calc_limits();
	}
	meshfile.close();

	// use size of first mesh for mesh set
	v_size = mesh[0].v_size;
	n_size = mesh[0].n_size;
	p_size = mesh[0].p_size;

	return 1;
}

/*!
This function reads a mesh from a binary file

@param file_in - name of file to read mesh from
@return - 1 is returned on success, 0 on failure
 */
int RicMeshSet::ReadBinary(string file_in)
{
	v_size=0;
	n_size=0;
	p_size=0;

	// time to open this file
	ifstream meshfile(file_in.c_str(), ios::in | ios::binary);
	if ( meshfile.fail() )
	{
		cerr<<"RicMeshSet - file open error"<<endl;
		return 0;
	}

	// read byte order
	char OrderStr[10];
	int ByteOrder;
	meshfile.read((char*) OrderStr, sizeof (OrderStr));
	if ( strncmp(OrderStr,"binarDCBA",9) == 0 )
	{
		ByteOrder = 0;	// i86
	}
	else if ( strncmp(OrderStr,"binarABCD",9) == 0 )
	{
		ByteOrder = 1;	// mac
	}
	else
	{
		cerr<<"RicMeshSet - unknown byte order\n";
		return 0;
	}

	meshfile.seekg(17);
	meshfile.read((char*) &p_dim, sizeof (p_dim));
	if ( ByteOrder ) p_dim = ByteSwap(p_dim);
	if (p_dim!=2 && p_dim!=3)
	{
		cerr<<"Error - polygon dimension must be 2 or 3"<<endl;
		return 0;
	}
	meshfile.read((char*) &ntstep, 4);
	if ( ByteOrder ) ntstep = ByteSwap(ntstep);

	// create a new set of meshes
	if ( mesh ) delete [] mesh;
	mesh = new RicMesh[ntstep];

	for (int t=0; t<ntstep; ++t)
	{
		// assign polygon dim to individual mesh
		mesh[t].p_dim = p_dim;

		// get the time step value
		meshfile.read((char*)&mesh[t].t_step, 4);
		if ( ByteOrder ) mesh[t].t_step = ByteSwap(mesh[t].t_step);

		// read vertices
		meshfile.read((char*) &v_size, 4);
		if ( ByteOrder ) v_size = ByteSwap(v_size);
		mesh[t].init_vertices(v_size);
		float x, y, z;
		for (int i=0; i<this->v_size; i++)
		{
			meshfile.read((char*) &x, 4);
			if ( ByteOrder ) x = ByteSwap(x);
			meshfile.read((char*) &y, 4);
			if ( ByteOrder ) y = ByteSwap(y);
			meshfile.read((char*) &z, 4);
			if ( ByteOrder ) z = ByteSwap(z);
			mesh[t].assign_node(i, x, y, z);
		}

		// read normals
		meshfile.read((char*) &n_size, 4) ;
		if ( ByteOrder ) n_size = ByteSwap(n_size);
		if (n_size!= this->v_size)
		{
			cerr<<"Mesh size and normals size are not equal exiting"<<endl;
			return 0;
		}
		mesh[t].init_normals(n_size);

		for (int i=0; i<this->n_size; i++)
		{
			meshfile.read((char*) &x, 4);
			if ( ByteOrder ) x = ByteSwap(x);
			meshfile.read((char*) &y, 4);
			if ( ByteOrder ) y = ByteSwap(y);
			meshfile.read((char*) &z, 4);
			if ( ByteOrder ) z = ByteSwap(z);
			mesh[t].assign_normal(i, x, y, z);
		}

		// read polygon indices
		meshfile.read((char*) &p_size, 4) ;
		if ( ByteOrder ) p_size = ByteSwap(p_size);
		if (p_size==0)
		{
			meshfile.read((char*) &p_size, 4) ; // read it one more time ??
			if ( ByteOrder ) p_size = ByteSwap(p_size);
		}

		mesh[t].init_polygons(p_size);
		int id1, id2, id3;
		if (this->p_dim == 2)
		{
			for (int i=0; i<this->p_size; i++)
			{
				meshfile.read((char*) &id1, 4);
				if ( ByteOrder ) id1 = ByteSwap(id1);
				meshfile.read((char*) &id2, 4);
				if ( ByteOrder ) id2 = ByteSwap(id2);
				mesh[t].assign_polygon(i, id1, id2, 0);
			}
		}
		else // must be dimension 3
		{
			for (int i=0; i<this->p_size; i++)
			{
				meshfile.read((char*) &id1, 4);
				if ( ByteOrder ) id1 = ByteSwap(id1);
				meshfile.read((char*) &id2, 4);
				if ( ByteOrder ) id2 = ByteSwap(id2);
				meshfile.read((char*) &id3, 4);
				if ( ByteOrder ) id3 = ByteSwap(id3);
				mesh[t].assign_polygon(i, id1, id2, id3);
			}
		}
		mesh[t].calc_limits();
	}
	meshfile.close();

	// use size of first mesh for mesh set
	v_size = mesh[0].v_size;
	n_size = mesh[0].n_size;
	p_size = mesh[0].p_size;

	return 1;
}

/*!
This function returns the index of the mesh in a mesh set with
the passed time step value.

@param tstep - tstep value to look for in array of meshes

@returns - index of the mesh with the time step value or -1 if not found
*/
int RicMeshSet::FindTimeStep(int tstep)
{
	for ( int i=0 ; i<ntstep ; ++i )
	{
		if ( mesh[i].t_step == tstep )
			return i;	// we found it
	}

	// we did not find it so return -1
	return -1;
}
