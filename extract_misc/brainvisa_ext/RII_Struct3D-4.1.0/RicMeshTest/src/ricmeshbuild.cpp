/***************************************************************************
 *   Copyright (C) 2007 by Bill Rogers                         *
 *   rogers@ric-rogerslinux-01                                             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//using namespace std;

#include <iostream>
#include <cstdlib>
#include "RicUtil.h"
#include "RicMeshSet.h"
#include "RicTextureSet.h"
#include "RicVolume.h"


int main(int argc, char *argv[])
{


	// read mesh
	RicMesh mesh("Bill_Lwhite.mesh");
	if ( mesh.p_size != 0 )
	{
		cout << "We just read in " << "Bill_Lwhite.mesh" <<endl;
		cout << "Number of polygons: " << mesh.p_size << endl << endl;
	}
	else
	{
		cout << "Error reading file "<< "Bill_Lwhite.mesh";
		exit(1);
	}

	// read multimesh
	RicMeshSet mBill("EH0111_Tmtktri_Mesh.mesh");
	if ( mBill.mesh[0].p_size != 0 )
	{
		cout << "We just read in "<< "EH0111_Tmtktri_Mesh.mesh" <<endl;
		cout << "Number of meshes read: "<<mBill.ntstep<<endl;
	}
	else
	{
		cout << "Error reading file "<<argv[1]<<endl;
		exit(1);
	}
	mBill.mesh[1].Write("EH0111_SC1.mesh");
	mBill.mesh[2].Write("EH0111_SC2.mesh");

	cout<<"Wrote meshes from multimesh\n\n";


	// combine two meshes
	RicMesh *CombMesh;
	CombMesh = mBill.mesh[50].Merge(&mBill.mesh[78]);

	if ( CombMesh->p_size != NULL )
	{
		cout << "Number of polygons in combined mesh: "
				<< CombMesh->p_size << endl << endl;
	}
	else
	{	return EXIT_SUCCESS;

		cout << "Error combining meshes ";
		exit(1);
	}


	// combine two meshes
	RicMesh *CombMesh2;
	CombMesh2 = CombMesh->Merge(&mBill.mesh[49]);

	if ( CombMesh2->p_size != NULL )
	{
		cout << "Number of polygons in combined mesh: "
				<< CombMesh2->p_size << endl << endl;
	}
	else
	{	return EXIT_SUCCESS;

		cout << "Error combining meshes ";
		exit(1);
	}

		// write combined mesh
	CombMesh2->Write("EH0111_S.C._left.mesh");

	return EXIT_SUCCESS;

	// now read texture
	RicTexture tex;
	tex.read_texture("BG0003_Lwhite-fem.tex");
	if ( tex.size != 0 )
	{
		cout << "We just read in " << "BG0003_Lwhite-fem.tex" <<endl;
		tex.CalcMinMaxAvg();
		cout << "Min="<<tex.min<<" Max="<<tex.max
				<<" Avg=" << tex.avg <<" Stddev="<<tex.std_dev <<endl<<endl;

	}
	else
	{
		cout << "Error reading file "<<"BG0003_Lwhite-fem.tex"<<endl;
		exit(1);
	}

	// read a multi texture
	RicTextureSet MultiRead("MultiTex.tex");
	cout << "We just read in TestMultiTex.tex" << endl;
	cout << "Number of timesteps is " << MultiRead.ntstep << endl;

	// write a texture
	MultiRead.Write("TexTest.tex");
	cout << "Wrote multitexture TexTest.tex" << endl;

	// write a single texture from the multi texture file
	MultiRead.tex[0].write_texture("SingleTex.tex");
	cout << "Wrote single texture from multitexture SingleTex.tex \n\n";

	cout << "That's all folks"<<endl;

	return EXIT_SUCCESS;
}
