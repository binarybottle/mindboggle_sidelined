// ----------------------- RicMesh.cpp ---------------------------------
/*! \mainpage
The RicMesh and RicTexture classes included in this library are used to
manipulate BrainVISA meshes and the textures that are mapped on them.

The RicMeshSet and RicTextureSet classes are used to read and write
files containing multiple meshes or textures.

The meshes can be made from polygons of dimension 2 or 3 - lines and
triangles.

What a texture really is, is a per node value aligned with the vertices
of a mesh. The value can be any useful parameter such a curvature,
thickness, etc.

*/

/*! \file
Implementation file for the RicMesh class. The class can open and write
binary and ascii BrainVisa/Anatomist file formats. This version will
handle polygons of dimension 2 and 3 - lines and triangles.

This is derived from the mesh class developed by Tom Arnow
Bill Rogers - Research Imaging Center - UTHSCSA - Nov 2008
rogers@uthscsa.edu
 */

#include "RicMesh.h"
#include "RicMeshSet.h"

using namespace std;

///////////////////////////////// constructors //////////////////////////////

// ---------------------------------------------------------------------------
/*!
Construct an empty mesh. All unused pointers are set to NULL.
 */
RicMesh::RicMesh (void)
{
	v_size = n_size = p_size = p_dim = t_step = 0;
	vertices = NULL;
	normals = NULL;
	polygons = NULL;
}

// ---------------------------------------------------------------------------
/*!
Contruct mesh to size.  All unused pointers are set to NULL.

@param vsize
number of vertices in mesh

@param nsize
number of normals - generally one per vertex

@param psize
number of polygons in mesh

@param dim
polygon dimension - either 2 or 3
 */
RicMesh::RicMesh (int vsize, int nsize, int psize ,int dim)
{
	this->p_dim = dim;
	this->v_size = vsize;
	this->n_size = nsize;
	this->p_size = psize;
	this->t_step = 0;

	vertices = normals = NULL;
	polygons = NULL;

	// allocate memory for arrays
	this->init_vertices(v_size);
	this->init_normals(n_size);
	this->init_polygons(p_size);

}

// ---------------------------------------------------------------------------
/*!
Construct mesh from a file.  All unused pointers are set to NULL.

@param fname
file name containing mesh
 */
RicMesh::RicMesh(char *fname)
{
	// zero or NULL everything
	v_size=0;n_size=0;p_size=0;
	vertices = NULL;
	normals = NULL;
	polygons = NULL;

	this->Read(fname);
}

// ---------------------------------------------------------------------------
/*!
Destructor - clean up memory
 */
RicMesh::~RicMesh()
{
	if (v_size) delete [] vertices;
	if (n_size)	delete [] normals;
	if (p_size) delete [] polygons;
}

// ---------------------------------------------------------------------------
/*!
Allocate memory for mesh

@param vs
number of vertices in mesh

@param ns
number of normals - generally one per vertex

@param ps
number of polygons in mesh

@param pd
polygon dimension - either 2 or 3

@returns
1 if successful, 0 if failure
*/
int RicMesh::Init(int vs, int ns, int ps, int pd)
{
	// assign to member variables
	v_size = vs;
	n_size = ns;
	p_size = ps;
	p_dim = pd;
	t_step = 0;

	if ( p_dim > 3 )
	{
		cerr << "RicMesh - Polygon size > 3" << endl;
		return 0;
	}

	// delete arrays if allocated
	if ( vertices ) delete [] vertices;
	if ( normals ) delete [] normals;
	if ( polygons ) delete [] polygons;

	// allocate memory for data
	vertices = new vertex[v_size];
	if ( vertices == NULL )
	{
		cerr << "RicMesh - vertices allocation error" << endl;
		return 0;
	}

	normals = new vertex[n_size];
	if ( normals == NULL )
	{
		cerr << "RicMesh - normals allocation error" << endl;
		return 0;
	}

	polygons = new triangle[p_size];
	if ( polygons == NULL )
	{
		cerr << "RicMesh - polygons allocation error" << endl;
		return 0;
	}

	// initialize arrays
	for (int i=0;i<v_size;i++)   vertices[i].assign (0,0,0);
	for (int i=0;i<n_size;i++)   normals[i].assign (0,0,0);
	for (int i=0;i<p_size;i++)   polygons[i].assign (0,0,0);

	return 1;

}

// ---------------------------------------------------------------------------
/*!
Allocate memory for vertices and zeros them

@param vs
number of vertices to allocate

@returns 1 on success, 0 if a problem
 */
int RicMesh::init_vertices(int vs)
{
	// delete vertices if existing
	if (vertices) delete [] vertices;
	vertices = NULL;

	v_size=vs;
	if ( v_size < 1 ) return 0; // needs to be at least one vertex
	vertices = new vertex[v_size];

	// if error in allocation return failure
	if ( vertices == NULL ) return 0;

	// zero the vertices
	for (int i=0;i<v_size;i++)   vertices[i].assign (0,0,0);

	return 1;
}

// ---------------------------------------------------------------------------
/*!
Allocate memory for normals and zeros them

@param ns
number of normals to allocate

@returns 1 on success, 0 if a problem
*/
int RicMesh::init_normals(int ns)
{
	// delete normals if existing
	if (normals) delete [] normals;
	normals = NULL;

	n_size=ns;
	if ( n_size < 1 ) return 0; // needs to be at least 1 normal
	normals=new vertex[n_size];

	// if error in allocation return failure
	if ( normals == NULL ) return 0;

	// zero the normals
	for (int i=0;i<v_size;i++)   normals[i].assign (0,0,0);

	return 1;
}

// ---------------------------------------------------------------------------
/*!
Allocate memory for polygons and zero them. It assumes a polygon size of
2 or 3.

@param ps
number of polygons to allocate

@returns 1 on success, 0 if a problem
 */
int RicMesh::init_polygons(int ps)
{
	// delete polygons if existing
	if (polygons) delete [] polygons;
	polygons = NULL;

	p_size=ps;
	if ( p_size < 1 ) return 0; // needs to be at least one polygon
	polygons=new triangle[p_size];

	// if allocation error return failure
	if ( polygons == NULL ) return 0;

	// zero polygon indices
	for (int i=0;i<p_size;i++)   polygons[i].assign (0,0,0);

	return 1;
}

// ---------------------------------------------------------------------------
/*!
Find the spatial limits of the mesh vertices.
 */
void RicMesh::calc_limits(void)
{
	xmin = ymin = zmin = 1000000;
	xmax = ymax = zmax = -1000000;
	for ( int i=0 ; i<v_size ; ++i )
	{
		if ( vertices[i].pnt.x < xmin ) xmin = vertices[i].pnt.x;
		if ( vertices[i].pnt.y < ymin ) ymin = vertices[i].pnt.y;
		if ( vertices[i].pnt.z < zmin ) zmin = vertices[i].pnt.z;
		if ( vertices[i].pnt.x > xmax ) xmax = vertices[i].pnt.x;
		if ( vertices[i].pnt.y > ymax ) ymax = vertices[i].pnt.y;
		if ( vertices[i].pnt.z > zmax ) zmax = vertices[i].pnt.z;
	}

}

// ---------------------------------------------------------------------------
/*!
 Check to see if a vertex exists in the vertex array

 @param v - vertex to check for existance
 @returns - index of vertex in array if found or -1 if not
 */
 int RicMesh::check_vertex(vertex v)
 {
 	for ( int i=0 ; i<v_size ; ++i )
 	{
 		if ( vertices[i].pnt.x < (v.pnt.x+FP_FUDGE) &&  vertices[i].pnt.x > (v.pnt.x-FP_FUDGE)
 				&& vertices[i].pnt.y < (v.pnt.y+FP_FUDGE) &&  vertices[i].pnt.y > (v.pnt.y-FP_FUDGE)
 				&& vertices[i].pnt.z < (v.pnt.z+FP_FUDGE) &&  vertices[i].pnt.z > (v.pnt.z-FP_FUDGE) )
 		{
 			return i;
 		}
 	}

 	// if here then not found
 	return -1;
 }

// ---------------------------------------------------------------------------
/*!
Assign a passed xyz value to a vertex node

@param index
index in node array

@param x
x value

@param y
y value

@param z
z value
 */
void RicMesh::assign_node(int index, float x,float y, float z)
{
	vertices[index].assign(x,y,z);
}

// ---------------------------------------------------------------------------
/*!
Assign a passed vertex to a vertex node

@param index
index in node array

@param v
vertex to assign to node
 */
void RicMesh::assign_node(int index, vertex v)
{
	vertices[index].set(v);
}

// ---------------------------------------------------------------------------
/*!
Assign a passed xyz value to a normal

@param index
index in normal array

@param x
x value

@param y
y value

@param z
z value
*/
void RicMesh::assign_normal(int index, float x,float y, float z)
{
	normals[index].assign(x,y,z);
}

// ---------------------------------------------------------------------------
/*!
Assign a passed vertex value to a normal

@param index
index in vertex array

@param v
vertex to assign to normal
*/
void RicMesh::assign_normal(int index, vertex v)
{
	normals[index].set(v);
}

// ---------------------------------------------------------------------------
/*!
Assign vertex node indices to a polygon. Assumes that there
are three polygons. If only a 2d polygon just use the first
two indices.

@param index
index in polygon array

@param v1
first vertex index

@param v2
second index

@param v3
third index
 */
void RicMesh::assign_polygon(int index, int v1,int v2, int v3)
{
	polygons[index].assign(v1,v2,v3);
}

// ---------------------------------------------------------------------------
/*!
Search for a vertex by comparing the distance between the passed
vertex and the vertices in the mesh specified by start and end index. If
a mesh vertex is close (less than 0.1) the index of the mesh vertex is
returned.

@param start_index
index in mesh vertex array to start searching from

@param end_index
last index in mesh vertex array to search

@param vtx
vertex to check against mesh

@returns
index in mesh array or 0 if not found
 */
int RicMesh::search_for_vertex(int start_index, int end_index, vertex vtx)
{
	for (int i=start_index; i<end_index;i++)
		if (vtx.e_distance(vertices[i])<0.1) return i;
	return 0;
}

// ---------------------------------------------------------------------------
/*!
This function writes a mesh to a file (default binary). It plays a trick as
all the actual writing is done in RicMeshSet. A RicMeshSet is declared and
the mesh is copied to the first mesh in the mesh set. The mesh set then does
the writing.

@param outfile
name of output file

@returns
1 on success, 0 on failure
 */
int RicMesh::Write(string outfile)
{
	// create a one mesh set
	RicMeshSet Mset(v_size,n_size,p_size,p_dim,1);
	if ( Mset.mesh == NULL )
	{
		cerr << "Allocation in RicMesh::Write" << endl;
		return 0;
	}

	// copy to the first mesh
	int i;
	for ( i=0; i<v_size; i++) Mset.mesh[0].vertices[i] = vertices[i];
	for ( i=0; i<n_size; i++) Mset.mesh[0].normals[i] = normals[i];
	for ( i=0; i<p_size; i++) Mset.mesh[0].polygons[i] = polygons[i];


	return Mset.Write(outfile);

}

// ---------------------------------------------------------------------------
/*!
This function reads a mesh from a file. It uses the RicMeshSet read routine
and copies only the first mesh read in to the mesh. Reading from a file only
occurs in RicMeshSet.

@param infile
file to read file from

@returns
1 on success, 0 on failure
 */
int RicMesh::Read(string infile)
{
	// create a one mesh set
	RicMeshSet Mset(infile);
	if ( Mset.mesh == NULL )
		return 0; // cannot read mesh

	// initialize the mesh to the size read in
	int stat = Init(Mset.mesh[0].v_size,Mset.mesh[0].n_size,
			Mset.mesh[0].p_size,Mset.mesh[0].p_dim);
	if ( stat == 0 ) return 0;

	// copy to the first mesh
	int i;
	for ( i=0; i<v_size; i++) vertices[i] = Mset.mesh[0].vertices[i];
	for ( i=0; i<n_size; i++) normals[i] = Mset.mesh[0].normals[i];
	for ( i=0; i<p_size; i++) polygons[i] = Mset.mesh[0].polygons[i];
	p_dim = Mset.p_dim;

	return 1;
}

// ---------------------------------------------------------------------------
/*!
This function writes a mesh to a text file.
XXX - obsolete version. Only around for historical purposes.
 */
void RicMesh::write_mesh_txt(char *file_out)
{
	ofstream fout(file_out);
	fout<<"ascii"<<endl;
	fout<<"VOID"<<endl;
	fout<<p_dim<<endl;
	fout<<1<<endl<<0<<endl; // only one timestep
// output polygons
	fout<<v_size<<endl<<flush;
	for (int i=0;i<v_size;i++)
		fout<<"("<<vertices[i].pnt.x<<" ,"<<vertices[i].pnt.y<<" ,"<<vertices[i].pnt.z<<") ";
	fout<<endl<<n_size<<endl;
	for (int i=0;i<n_size;i++)
		fout<<"("<<normals[i].pnt.x<<" ,"<<normals[i].pnt.y<<" ,"<<normals[i].pnt.z<<") ";
	fout<<endl<<0<<endl<<p_size<<endl;
	if ( p_dim == 2 ) // lines
	{
		for (int i=0;i<p_size;i++)
			fout<<"("<<polygons[i].vidx[0]<<" ,"<<polygons[i].vidx[1]<<") ";
	}
	else // p_dim is 3
	{
		for (int i=0;i<p_size;i++)
			fout<<"("<<polygons[i].vidx[0]<<" ,"<<polygons[i].vidx[1]<<" ,"<<polygons[i].vidx[2]<<") ";
	}

	fout << endl;
	fout.close();
}

// ---------------------------------------------------------------------------
/*!
This function writes a mesh to a binary file
XXX - obsolete version. Only around for historical purposes.
 */
void RicMesh::write_mesh_bin(char *file_out)
{
	ofstream fout(file_out);
	fout<<"binarDCBA";
	int A=4	;
	fout.write((char*) &A,sizeof (A));
	fout<<"VOID";
	A=p_dim;  fout.write((char*) &A,sizeof (A)); // polygon dimension
	A=1;  fout.write((char*) &A,sizeof (A)); // number of time steps
	A=0;  fout.write((char*) &A,sizeof (A)); // current time step

	A=v_size; fout.write((char*) &A,sizeof (A));

	for (int i=0;i<v_size;i++)
	{
		Point tmp= vertices[i].pnt;
		fout.write( (char*) & tmp.x, sizeof (float) );
		fout.write( (char*) & tmp.y, sizeof (float) );
		fout.write( (char*) & tmp.z, sizeof (float) );
	}

	A=n_size; fout.write((char*) &A,sizeof (A));
	for (int i=0;i<n_size;i++)
	{
		Point tmp=normals[i].pnt;
		fout.write((char*) &tmp.x, sizeof(float));
		fout.write((char*) &tmp.y, sizeof(float));
		fout.write((char*) &tmp.z, sizeof(float));
	}
	A=0; fout.write((char*) &A,sizeof (A));
	A=p_size; fout.write((char*) &A,sizeof (A));

	for (int i=0;i<p_size;i++)
		for (int j=0;j<p_dim;j++)
			fout.write((char*) &polygons[i].vidx[j], sizeof(int));

	fout.close();
}

// ---------------------------------------------------------------------------
/*!
This function reads a mesh from a text file
XXX - obsolete version. Only around for historical purposes.
1 is returned on success, 0 on failure
 */
int RicMesh::read_mesh_txt(char *file_in)
{
	int i,j;
	v_size=0;n_size=0;p_size=0;

	int current; // current time step

	// time to open this file
	ifstream meshfile(file_in,ifstream::in);

	// get the header info
	string type, junk;
	meshfile >> type;
	if ( type != "ascii" ) return 0;	// must be ascii
	meshfile >> junk;	// should be "VOID"
	meshfile >> p_dim;  // polygon dimension
	if(p_dim!=2 && p_dim!=3) {cerr<<"Error - polygon dimension must be 2 or 3"<<endl;}
//	meshfile >> ntstep;
	meshfile >> current;

	// read vertices
	meshfile >> v_size;

	// now get all of those vertices using ')' as delimiter for end of vertex
	char tvertex[100];
	float x,y,z;
 	this->init_vertices(v_size);
	for ( i=0 ; i<v_size ; ++i )
	{

		meshfile.getline(tvertex,99,'(');	// read to start of vertex
		meshfile.getline(tvertex,99,')');	// now read to end of vertex

		// remove commas
		j = tvertex[99] = 0;
		do
		{
			if ( tvertex[j] == ',' ) tvertex[j] = ' ';
			++j;
		}while ( tvertex[j] !=0 );

		// parse the vector
		sscanf(tvertex,"%f %f %f",&x,&y,&z);

		this->assign_node(i,x,y,z);
	}

	// read normals
	meshfile >> n_size;
	if ( (n_size !=  this->v_size) && (n_size != 0) )
	{
		cerr<<"Mesh size and normals size are not equal exiting"<<endl;
		return(0);
	}

	if ( n_size > 0 ) this->init_normals(n_size);

	for ( i=0 ; i<n_size ; ++i )
	{

		meshfile.getline(tvertex,99,'(');	// read to start of vertex
		meshfile.getline(tvertex,99,')');	// now read to end of vertex

		// remove commas
		j = tvertex[99] = 0;
		do
		{
			if ( tvertex[j] == ',' ) tvertex[j] = ' ';
			++j;
		}while ( tvertex[j] !=0 );

		// parse the vector
		sscanf(tvertex,"%f %f %f",&x,&y,&z);

		this->assign_normal(i,x,y,z);
	}

	// read polygon indices
	meshfile >> p_size; // dummy read for texture (see file format doc)
	meshfile >> p_size ; // read it one more time ??

	if ( p_size > 0 ) this->init_polygons(p_size);

	int id1,id2,id3;
	for ( i=0 ; i<p_size ; ++i )
	{

		meshfile.getline(tvertex,99,'(');	// read to start of vertex
		meshfile.getline(tvertex,99,')');	// now read to end of vertex

	// remove commas
		j = tvertex[99] = 0;
		do
		{
			if ( tvertex[j] == ',' ) tvertex[j] = ' ';
			++j;
		}while ( tvertex[j] !=0 );

	// parse the vector
		if ( this->p_dim == 2 )
		{
			id3 = 0;
			sscanf(tvertex,"%d %d",&id1,&id2);
			this->assign_polygon(i,id1,id2,id3);
		}
		else // p_dim == 3
		{
			sscanf(tvertex,"%d %d %d",&id1,&id2,&id3);
			this->assign_polygon(i,id1,id2,id3);
		}
	}

	meshfile.close();
	return 1;
}

// ---------------------------------------------------------------------------
/*!
This function reads a mesh from a binary file
XXX - obsolete version. Only around for historical purposes.
1 is returned on success, 0 on failure
 */
int RicMesh::read_mesh_bin(char *file_in)
{
	v_size=0;n_size=0;p_size=0;

	int current; // current time step

	// time to open this file
	ifstream meshfile(file_in, ios::in | ios::binary);

	meshfile.seekg(17);
	meshfile.read((char*) &p_dim,sizeof (p_dim));
	if(p_dim!=2 && p_dim!=3) {cerr<<"Error - polygon dimension must be 2 or 3"<<endl;}
//	meshfile.read((char*) &ntstep,4);
	meshfile.read((char*)& current,4);

	// read vertices
	meshfile.read((char*) &v_size,4);
	this->init_vertices(v_size);
	float x,y,z;
	for (int i=0;i<this->v_size;i++){
		meshfile.read((char*) &x,4);
		meshfile.read((char*) &y,4);
		meshfile.read((char*) &z,4);
		this->assign_node(i,x,y,z);
	}

	// read normals
	meshfile.read((char*) &n_size,4) ;
	if (n_size!=  this->v_size)
	{cerr<<"Mesh size and normals size are not equal exiting"<<endl; return(0);}
	this->init_normals(n_size);

	for (int i=0;i<this->n_size;i++){
		meshfile.read((char*) &x,4);
		meshfile.read((char*) &y,4);
		meshfile.read((char*) &z,4);
		this->assign_normal(i,x,y,z);
	}

	// read polygon indices
	meshfile.read((char*) &p_size,4) ;
	if (p_size==0)  meshfile.read((char*) &p_size,4) ; // read it one more time ??

	this->init_polygons(p_size);
	int id1,id2,id3;
	if ( this->p_dim == 2 )
	{
		for (int i=0;i<this->p_size;i++){
			meshfile.read((char*) &id1,4);
			meshfile.read((char*) &id2,4);
			this->assign_polygon(i,id1,id2,0);
		}
	}
	else // must be dimension 3
	{
		for (int i=0;i<this->p_size;i++){
			meshfile.read((char*) &id1,4);
			meshfile.read((char*) &id2,4);
			meshfile.read((char*) &id3,4);
			this->assign_polygon(i,id1,id2,id3);
		}
	}

	meshfile.close();
	return 1;
}

// ---------------------------------------------------------------------------
/*!
This function merges the mesh with the input mesh to make one larger
mesh. Both meshes must have the same polygon dimension.

@param InMesh - input mesh to be merged with current mesh.

@returns - Combined mesh on success or NULL on failure
 */
RicMesh* RicMesh::Merge(RicMesh *InMesh)
{
	// check to see if both meshes have the same dimension
	if ( p_dim != InMesh->p_dim )
		return NULL;

	// create new mesh size
	RicMesh *CombMesh;
	CombMesh = new RicMesh(v_size+InMesh->v_size, n_size+InMesh->n_size,
			p_size+InMesh->p_size, p_dim);

	// copy vertices to combined mesh
	for ( int i=0 ; i<v_size ; ++i )
		CombMesh->vertices[i] = vertices[i];
	for ( int i=v_size, j=0 ; i<CombMesh->v_size ; ++i,++j )
		CombMesh->vertices[i] = InMesh->vertices[j];

	// copy normals
	for ( int i=0 ; i<n_size ; ++i )
		CombMesh->normals[i] = normals[i];
	for ( int i=n_size, j=0 ; i<CombMesh->n_size ; ++i,++j )
		CombMesh->normals[i] = InMesh->normals[j];

	// copy polygons
	for ( int i=0 ; i<p_size ; ++i )
		CombMesh->polygons[i] = polygons[i];
	for ( int i=p_size, j=0 ; i<CombMesh->p_size ; ++i,++j )
		CombMesh->polygons[i] = InMesh->polygons[j];

	// now offset the polygons from the second mesh to line
	// up with the proper vertex
	for ( int i=p_size, j=0 ; i<CombMesh->p_size ; ++i,++j )
	{
		CombMesh->polygons[i].vidx[0] += v_size;
		CombMesh->polygons[i].vidx[1] += v_size;
		CombMesh->polygons[i].vidx[2] += v_size;
	}

	return CombMesh;
}

// ---------------------------------------------------------------------------
/*!
This function finds the area of a triangular mesh by adding the
areas of all the triangles.

@returns
area of mesh
 */
float RicMesh::area(void)
{
	// if mesh not triangles then return 0
	if ( this->p_dim != 3 ) return 0;
	int nz=0;

	// sum the area of all the triangles
	float area=0;
	for ( int i=0 ; i<p_size ; ++i )
	{
		Point P0,P1,P2;
		P0 = vertices[polygons[i].vidx[0]].pnt;
		P1 = vertices[polygons[i].vidx[1]].pnt;
		P2 = vertices[polygons[i].vidx[2]].pnt;

		float a = area_triangle(P0,P1,P2);
		if ( a == 0 ) ++nz;
		area += a;
	}
	return area;
}

// ---------------------------------------------------------------------------
/*!
This function finds the volume of a closed 3D triangular mesh
by a 3D planimeter method. It will malfunction horribly if the mesh
is not closed. This came from Jack Lancaster.

@returns
volume of mesh
 */
float RicMesh::volume(void)
{
	// This method projects the area of each mesh triangle on a plane
	int i;
	float xvol=0,yvol=0,zvol=0; // volume on each plane
	float vol;		// average volume
	Vector n;		// normal
	float a;		// triangle area
	Point cp[3];	// points used to determine cross product
	Point c;		// triangle center

	// x direction
	for ( i=0 ; i<p_size ; ++i )
	{
		cp[0] = vertices[polygons[i].vidx[0]].pnt;
		cp[1] = vertices[polygons[i].vidx[1]].pnt;
		cp[2] = vertices[polygons[i].vidx[2]].pnt;

		// area of triangle
		a = area_triangle(cp[0],cp[1],cp[2]);

		// normal of triangle
		n = cross_product(cp);

		// center of triangle
		c.x = (cp[0].x+cp[1].x+cp[2].x)*0.333333;
		c.y = (cp[0].y+cp[1].y+cp[2].y)*0.333333;
		c.z = (cp[0].z+cp[1].z+cp[2].z)*0.333333;

		// multiply the area of the triangle by the distance to the axis
		// by the contribution of the unit normal in that direction
		xvol += a * c.x * n.x;
		yvol += a * c.y * n.y;
		zvol += a * c.z * n.z;
	}

	vol = fabs((xvol+yvol+zvol)*0.333333);

	return vol;
}

// ---------------------------------------------------------------------------
/*!
This function super samples a mesh to produce a new mesh. Each super sampled
mesh triangle is turned into 4 smaller triangles. If a triangle is below the
minumum area, it is just copied to the new mesh.

@param min_area - mininum triangle size to resample
@return new super sampled mesh
 */
RicMesh* RicMesh::super_sample(float min_area)
{
	// needs to be a 3d mesh
	if ( p_dim != 3 ) return NULL;

	RicMesh *smesh;
	smesh = new RicMesh(v_size*4, n_size*4, p_size*4, p_dim);
	int maxverts = smesh->v_size-4; // max allowed vertices in new mesh

	// just copy over the original vertices and normals to begin with
	for ( int i=0 ; i<v_size ; ++i )
		smesh->vertices[i] = vertices[i];
	for ( int i=0 ; i<n_size ; ++i )
		smesh->normals[i] = normals[i];

	// make four new polygons for each original polygon
	vertex p0,p1,p2,p01,p12,p02;
	vertex n0,n1,n2,n01,n12,n02;
	int	i0,i1,i2;
	int pc = 0;	// new polygon count
	int vc = v_size;	// new vertex
	float area;
	int i;
	for ( i=0 ; i<p_size ; ++i )
	{
		i0 = polygons[i].vidx[0];
		i1 = polygons[i].vidx[1];
		i2 = polygons[i].vidx[2];

		// original triangle vertices
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		// find area of polygon
		area = area_triangle(p0.pnt,p1.pnt,p2.pnt);
		if ( area < min_area ) // just copy this one and continue on
		{
			smesh->polygons[pc++] = polygons[i];
			continue;
		}

		// calculate mid points
		p01 = p0.average2(p1);
		p12 = p1.average2(p2);
		p02 = p0.average2(p2);

		// original normals
		n0 = normals[i0];
		n1 = normals[i1];
		n2 = normals[i2];

		// check each midpoint to see if it already exists in vertex array
		int n, pidx01, pidx12,pidx02;

		// p01
		if ( (n=smesh->check_vertex(p01)) >= 0 ) // vertex exists
		{
			pidx01 = n;
		}
		else // we need to add a vertex and a normal
		{
			// add a vertex
			pidx01 = vc;
			smesh->vertices[vc] = p01;

			// find average normal to add
			n01 = n0.average2(n1);
			smesh->normals[vc] = n01;

			++vc;
		}

		// p12
		if ( (n=smesh->check_vertex(p12)) >= 0 ) // vertex exists
		{
			pidx12 = n;
		}
		else // we need to add a vertex and a normal
		{
			// add a vertex
			pidx12 = vc;
			smesh->vertices[vc] = p12;

			// find average normal to add
			n12 = n1.average2(n2);
			smesh->normals[vc] = n12;

			++vc;
		}

		// p02
		if ( (n=smesh->check_vertex(p02)) >= 0 ) // vertex exists
		{
			pidx02 = n;
		}
		else // we need to add a vertex and a normal
		{
			// add a vertex
			pidx02 = vc;
			smesh->vertices[vc] = p02;

			// find average normal to add
			n02 = n0.average2(n2);
			smesh->normals[vc] = n02;

			++vc;
		}

		// now add four polygons to the list
		smesh->polygons[pc].vidx[0] = pidx02;
		smesh->polygons[pc].vidx[1] = i0;
		smesh->polygons[pc].vidx[2] = pidx01;
		++pc;

		smesh->polygons[pc].vidx[0] = pidx01;
		smesh->polygons[pc].vidx[1] = i1;
		smesh->polygons[pc].vidx[2] = pidx12;
		++pc;

		smesh->polygons[pc].vidx[0] = pidx12;
		smesh->polygons[pc].vidx[1] = i2;
		smesh->polygons[pc].vidx[2] = pidx02;
		++pc;

		smesh->polygons[pc].vidx[0] = pidx01;
		smesh->polygons[pc].vidx[1] = pidx12;
		smesh->polygons[pc].vidx[2] = pidx02;
		++pc;

		if ( vc >= maxverts )
		{
			cerr << "RicMesh::super_sample - number of vertices exceeded limit" << endl;
			break;
		}
	}

	// copy to properly sized mesh
	RicMesh *smesh2;
	smesh2 = new RicMesh(vc, vc, pc, p_dim);

	// copy all vertices, polygons and normals
	for ( i=0 ; i< vc ; ++i ) smesh2->vertices[i] = smesh->vertices[i];
	for ( i=0 ; i< vc ; ++i ) smesh2->normals[i] = smesh->normals[i];
	for ( i=0 ; i< vc ; ++i ) smesh2->normals[i].normalize();
	for ( i=0 ; i< pc ; ++i ) smesh2->polygons[i] = smesh->polygons[i];

	delete smesh;

	return smesh2;
}

// ---------------------------------------------------------------------------
/*!
This function super samples a mesh to produce a new mesh. Each super sampled
mesh triangle is turned into 4 smaller triangles. If a triangle is below the
minumum area, it is just copied to the new mesh. Note that this version
can create duplicate vertices which can be a problem in certain applications.

@param min_area - mininum triangle size to resample
@return new super sampled mesh
 */
RicMesh* RicMesh::super_sample2(float min_area)
{
	// needs to be a 3d mesh
	if ( p_dim != 3 ) return NULL;

	RicMesh *smesh;
	smesh = new RicMesh(v_size*8, n_size*8, p_size*4, p_dim);

	// just copy over the original vertices and normals to begin with
	for ( int i=0 ; i<v_size ; ++i )
		smesh->vertices[i] = vertices[i];
	for ( int i=0 ; i<n_size ; ++i )
		smesh->normals[i] = normals[i];

	// make four new polygons for each original polygon
	vertex p0,p1,p2,p01,p12,p02;
	vertex n0,n1,n2,n01,n12,n02;
	int	i0,i1,i2;
	int pc = 0;	// new polygon count
	int vc = v_size;	// new vertex
	float area;
	int i;
	for ( i=0 ; i<p_size ; ++i )
	{
		i0 = polygons[i].vidx[0];
		i1 = polygons[i].vidx[1];
		i2 = polygons[i].vidx[2];

		// original triangle vertices
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		// find area of polygon
		area = area_triangle(p0.pnt,p1.pnt,p2.pnt);
		if ( area < min_area ) // just copy this one and continue on
		{
			smesh->polygons[pc++] = polygons[i];
			continue;
		}

		// calculate mid points
		p01 = p0.average2(p1);
		p12 = p1.average2(p2);
		p02 = p0.average2(p2);

		// original normals
		n0 = normals[i0];
		n1 = normals[i1];
		n2 = normals[i2];

		// calculate mid points
		n01 = n0.average2(n1);
		n12 = n1.average2(n2);
		n02 = n0.average2(n2);

		// add the new vertices and normals to the list
		smesh->vertices[vc] = p01;
		smesh->normals[vc] = n01;
		smesh->vertices[vc+1] = p12;
		smesh->normals[vc+1] = n12;
		smesh->vertices[vc+2] = p02;
		smesh->normals[vc+2] = n02;


		// now add four polygons to the list
		smesh->polygons[pc].vidx[0] = vc+2;
		smesh->polygons[pc].vidx[1] = i0;
		smesh->polygons[pc].vidx[2] = vc;
		++pc;

		smesh->polygons[pc].vidx[0] = vc;
		smesh->polygons[pc].vidx[1] = i1;
		smesh->polygons[pc].vidx[2] = vc+1;
		++pc;

		smesh->polygons[pc].vidx[0] = vc+1;
		smesh->polygons[pc].vidx[1] = i2;
		smesh->polygons[pc].vidx[2] = vc+2;
		++pc;

		smesh->polygons[pc].vidx[0] = vc;
		smesh->polygons[pc].vidx[1] = vc+1;
		smesh->polygons[pc].vidx[2] = vc+2;
		++pc;
		vc += 3;

		if ( vc >= (smesh->v_size-3) )
		{
			cerr << "RicMesh::super_sample2 - number of vertices exceeded limit" << endl;
			break;
		}
	}

	// copy to properly sized mesh
	RicMesh *smesh2;
	smesh2 = new RicMesh(vc, vc, pc, p_dim);

	// copy all vertices, polygons and normals
	for ( i=0 ; i< vc ; ++i ) smesh2->vertices[i] = smesh->vertices[i];
	for ( i=0 ; i< vc ; ++i ) smesh2->normals[i] = smesh->normals[i];
	for ( i=0 ; i< vc ; ++i ) smesh2->normals[i].normalize();
	for ( i=0 ; i< pc ; ++i ) smesh2->polygons[i] = smesh->polygons[i];

	delete smesh;

	return smesh2;
}

// ---------------------------------------------------------------------------
/*!
This function looks for extremely narrow triangles in a mesh. The narrow
triangles are broken into more regular triangles in a new mesh.

@param max_ratio - maximum ratio of triangle edges to allow
@return new triangle mesh
 */
RicMesh* RicMesh::triangle_fix(float max_ratio)
{
	// needs to be a 3d mesh
	if ( p_dim != 3 ) return NULL;

	RicMesh *smesh;
	smesh = new RicMesh(v_size*10, n_size*10, p_size*20, p_dim);

	// make new sub polygons for each original polygon
	vertex p0,p1,p2,p01,p12,p02;
	vertex n0,n1,n2,n01,n12,n02;
	float l01,l12,l02;		// lengths of triangle sides
	float ratio;			// ratios of sides
	int	i0,i1,i2;
	int pc = 0;	// new polygon count
	int vc = 0;	// new vertex count
	int i,j;
	for ( i=0 ; i<p_size ; ++i )
	{
		i0 = polygons[i].vidx[0];
		i1 = polygons[i].vidx[1];
		i2 = polygons[i].vidx[2];

		// original triangle vertices
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		// original triangle normals

		// calculate side lengths
		l01 = dist(vertices[i0].pnt,vertices[i1].pnt);
		l12 = dist(vertices[i1].pnt,vertices[i2].pnt);
		l02 = dist(vertices[i0].pnt,vertices[i2].pnt);

		// figure out the shortest side
		// assign points and normals to keep in same order for rh rule
		// p1-p2 will be assigned shortest side
		if ( l01 < l12 && l01 < l02) // 0-1 shortest
		{
			ratio = max(l12/l01,l02/l01);
			p0 = vertices[i2];
			p1 = vertices[i0];
			p2 = vertices[i1];
			n0 = normals[i2];
			n1 = normals[i0];
			n2 = normals[i1];
		}
		else if ( l12 < l01 && l12 < l02) // 1-2 shortest
		{
			ratio = max(l01/l12,l02/l12);
			p0 = vertices[i0];
			p1 = vertices[i1];
			p2 = vertices[i2];
			n0 = normals[i0];
			n1 = normals[i1];
			n2 = normals[i2];
		}
		else // 0-2 shortest
		{
			ratio = max(l01/l02,l12/l02);
			p0 = vertices[i1];
			p1 = vertices[i2];
			p2 = vertices[i0];
			n0 = normals[i1];
			n1 = normals[i2];
			n2 = normals[i0];
		}

		// check to see if within limits
		if ( ratio < max_ratio )
		{
			// assign to polygon array
			smesh->polygons[pc].vidx[0] = vc;
			smesh->polygons[pc].vidx[1] = vc+1;
			smesh->polygons[pc].vidx[2] = vc+2;
			++pc;
			// assign normals and vertices
			smesh->normals[vc] = normals[i0];
			smesh->vertices[vc++] = vertices[i0];
			smesh->normals[vc] = normals[i1];
			smesh->vertices[vc++] = vertices[i1];
			smesh->normals[vc] = normals[i2];
			smesh->vertices[vc++] = vertices[i2];
			continue;
		}

		// figure out how many divisions for the long sides
		int ndiv = (int)ratio;

		// calculate increments along sides
		float inc01x,inc02x;
		float inc01y,inc02y;
		float inc01z,inc02z;
		inc01x = (p1.pnt.x-p0.pnt.x)/ndiv;
		inc02x = (p2.pnt.x-p0.pnt.x)/ndiv;
		inc01y = (p1.pnt.y-p0.pnt.y)/ndiv;
		inc02y = (p2.pnt.y-p0.pnt.y)/ndiv;
		inc01z = (p1.pnt.z-p0.pnt.z)/ndiv;
		inc02z = (p2.pnt.z-p0.pnt.z)/ndiv;

		// calculate increments along normals
		float ninc01x,ninc02x;
		float ninc01y,ninc02y;
		float ninc01z,ninc02z;
		ninc01x = (n1.pnt.x-n0.pnt.x)/ndiv;
		ninc02x = (n2.pnt.x-n0.pnt.x)/ndiv;
		ninc01y = (n1.pnt.y-n0.pnt.y)/ndiv;
		ninc02y = (n2.pnt.y-n0.pnt.y)/ndiv;
		ninc01z = (n1.pnt.z-n0.pnt.z)/ndiv;
		ninc02z = (n2.pnt.z-n0.pnt.z)/ndiv;

		vertex pp0,pp1,pp2;		// current triangle points
		vertex nn0,nn1,nn2;		// current triangle normals

		// first polygon is easy, the one away from shortest side
		pp0 = p0;
		nn0 = n0;

		pp1.pnt.x = p0.pnt.x + inc01x;
		pp1.pnt.y = p0.pnt.y + inc01y;
		pp1.pnt.z = p0.pnt.z + inc01z;

		nn1.pnt.x = n0.pnt.x + ninc01x;
		nn1.pnt.y = n0.pnt.y + ninc01y;
		nn1.pnt.z = n0.pnt.z + ninc01z;

		pp2.pnt.x = p0.pnt.x + inc02x;
		pp2.pnt.y = p0.pnt.y + inc02y;
		pp2.pnt.z = p0.pnt.z + inc02z;

		nn2.pnt.x = n0.pnt.x + ninc02x;
		nn2.pnt.y = n0.pnt.y + ninc02y;
		nn2.pnt.z = n0.pnt.z + ninc02z;

		// assign to polygon array
		smesh->polygons[pc].vidx[0] = vc;
		smesh->polygons[pc].vidx[1] = vc+1;
		smesh->polygons[pc].vidx[2] = vc+2;
		++pc;

		// assign vertices and normals to arrays
		smesh->normals[vc] = nn0;
		smesh->vertices[vc++] = pp0;
		smesh->normals[vc] = nn1;
		smesh->vertices[vc++] = pp1;
		smesh->normals[vc] = nn2;
		smesh->vertices[vc++] = pp2;

		// now repeat for each increment along the long sides
		// there will be two triangles for each increment
		for ( j=1; j<ndiv ; ++j )
		{
			// first triangle
			pp0.pnt.x = p0.pnt.x + (j*inc01x);
			pp0.pnt.y = p0.pnt.y + (j*inc01y);
			pp0.pnt.z = p0.pnt.z + (j*inc01z);

			nn0.pnt.x = n0.pnt.x + (j*ninc01x);
			nn0.pnt.y = n0.pnt.y + (j*ninc01y);
			nn0.pnt.z = n0.pnt.z + (j*ninc01z);

			pp1.pnt.x = p0.pnt.x + ((j+1)*inc01x);
			pp1.pnt.y = p0.pnt.y + ((j+1)*inc01y);
			pp1.pnt.z = p0.pnt.z + ((j+1)*inc01z);

			nn1.pnt.x = n0.pnt.x + ((j+1)*ninc01x);
			nn1.pnt.y = n0.pnt.y + ((j+1)*ninc01y);
			nn1.pnt.z = n0.pnt.z + ((j+1)*ninc01z);

			pp2.pnt.x = p0.pnt.x + ((j+1)*inc02x);
			pp2.pnt.y = p0.pnt.y + ((j+1)*inc02y);
			pp2.pnt.z = p0.pnt.z + ((j+1)*inc02z);

			nn2.pnt.x = n0.pnt.x + ((j+1)*ninc02x);
			nn2.pnt.y = n0.pnt.y + ((j+1)*ninc02y);
			nn2.pnt.z = n0.pnt.z + ((j+1)*ninc02z);

			// assign to polygon array
			smesh->polygons[pc].vidx[0] = vc;
			smesh->polygons[pc].vidx[1] = vc+1;
			smesh->polygons[pc].vidx[2] = vc+2;
			++pc;

			// assign vertices and normals to arrays
			smesh->normals[vc] = nn0;
			smesh->vertices[vc++] = pp0;
			smesh->normals[vc] = nn1;
			smesh->vertices[vc++] = pp1;
			smesh->normals[vc] = nn2;
			smesh->vertices[vc++] = pp2;

			// second triangle
			// we can cheat here by reusing pp0 and copying pp2 to pp1;
			pp1 = pp2;
			nn1 = nn2;

			pp2.pnt.x = p0.pnt.x + (j*inc02x);
			pp2.pnt.y = p0.pnt.y + (j*inc02y);
			pp2.pnt.z = p0.pnt.z + (j*inc02z);

			nn2.pnt.x = n0.pnt.x + (j*ninc02x);
			nn2.pnt.y = n0.pnt.y + (j*ninc02y);
			nn2.pnt.z = n0.pnt.z + (j*ninc02z);

			// assign to polygon array
			smesh->polygons[pc].vidx[0] = vc;
			smesh->polygons[pc].vidx[1] = vc+1;
			smesh->polygons[pc].vidx[2] = vc+2;
			++pc;

			// assign vertices and normals to arrays
			smesh->normals[vc] = nn0;
			smesh->vertices[vc++] = pp0;
			smesh->normals[vc] = nn1;
			smesh->vertices[vc++] = pp1;
			smesh->normals[vc] = nn2;
			smesh->vertices[vc++] = pp2;

		}

		// make sure we have not used too many polygons or vertices
		if ( vc > (smesh->v_size-30) || (pc > smesh->p_size-30) )
		{
			break;
		}
	}

	// copy to properly sized mesh
	RicMesh *smesh2;
	smesh2 = new RicMesh(vc, vc, pc, p_dim);

	// copy all vertices, polygons and normals
	for ( i=0 ; i< vc ; ++i ) smesh2->vertices[i] = smesh->vertices[i];
	for ( i=0 ; i< vc ; ++i ) smesh2->normals[i] = smesh->normals[i];
	for ( i=0 ; i< vc ; ++i ) smesh2->normals[i].normalize();
	for ( i=0 ; i< pc ; ++i ) smesh2->polygons[i] = smesh->polygons[i];

	delete smesh;

	return smesh2;
}

// ---------------------------------------------------------------------------
/*!
This function super samples the mesh triangles recursively until all triangles
are less than an maximum area. A mininum triangle size is set so triangles
less than that size are not super sampled.

@param min_area - mininum triangle size to resample
@param max_area - maximum triangle size allowed
@return new super sampled mesh
 */
RicMesh* RicMesh::adjust_triangle_size(float min_area, float max_area)
{
	RicMesh *m1,*m2;
	int i;

	// copy the mesh
	m1 = new RicMesh(v_size, n_size, p_size, p_dim);
	for ( i=0 ; i<v_size ; ++i ) m1->vertices[i] = vertices[i];
	for ( i=0 ; i<n_size ; ++i ) m1->normals[i] = normals[i];
	for ( i=0 ; i<p_size ; ++i ) m1->polygons[i] = polygons[i];


	float a,maxa;
	do
	{
		// calculate maximum triangle area
		maxa = 0;
		for ( i=0 ; i<m1->p_size ; ++i )
		{
			a = area_triangle(m1->vertices[m1->polygons[i].vidx[0]].pnt,
				m1->vertices[m1->polygons[i].vidx[1]].pnt,
				m1->vertices[m1->polygons[i].vidx[2]].pnt);
			if ( a > maxa ) // we are over the max so we can quit looking here
			{
				maxa = a;
				if ( maxa > max_area )
					break;
			}
		}

		if ( maxa < max_area ) // check to see if we are done
			break;

		// supersample the array
		m2 = m1->super_sample(min_area);

		// copy the mesh to the original pointer
		delete m1;
		m1 = new RicMesh(m2->v_size, m2->n_size, m2->p_size, p_dim);
		for ( i=0 ; i<m2->v_size ; ++i ) m1->vertices[i] = m2->vertices[i];
		for ( i=0 ; i<m2->n_size ; ++i ) m1->normals[i] = m2->normals[i];
		for ( i=0 ; i<m2->p_size ; ++i ) m1->polygons[i] = m2->polygons[i];
		delete m2;

	} while ( maxa > max_area );

	return m1;
}

// ---------------------------------------------------------------------------
/*!
@brief This function checks to see if a point is inside a mesh.

This function checks to see if a point is inside a mesh. It first checks to see
if the point is in the bounding box of the mesh. If not then it uses a ray from
a distant point and checks the number of intersections with mesh triangles. If
there is an odd number of intersections then the point is inside the mesh.

@param p - point to test to see if inside mesh
@return 1 if point inside mesh, 0 if not
 */
int RicMesh::point_inside_mesh(Point p)
{
	// see if point in bounding box
	if ( p.x < xmin || p.x > xmax || p.y < ymin || p.y > ymax || p.x < zmin || p.z > zmax )
		return 0;

	int nint = 0; // number of intersections
	Point pd; // point outside mesh
	pd.x = xmax+100;
	pd.y = ymax+100;
	pd.z = zmax+100;
	Point pint; // intersection point - not used

	for ( int i=0 ; i< p_size ; ++i )
	{
		Point t0,t1,t2;
		t0 = vertices[polygons[i].vidx[0]].pnt;
		t1 = vertices[polygons[i].vidx[1]].pnt;
		t2 = vertices[polygons[i].vidx[2]].pnt;

		if ( line_thru_triangle(t0, t1, t2,p, pd, &pint) )
			++nint;
	}

	// see if even or odd
	if ( nint % 2 )
		return 1; // odd
	else
		return 0; // even

}
