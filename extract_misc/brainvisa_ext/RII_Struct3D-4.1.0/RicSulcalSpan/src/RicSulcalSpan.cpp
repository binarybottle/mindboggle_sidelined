//////////////////////////// RicSulcalSpan.cpp /////////////////////////////
// This program calculates the width of a sulcus given the
// sulcal mesh file (generally Tmtktri.mesh), the graph file, and the
// cortical mesh file.

/*! \file
Implementation file for RicSulcalSpan.cpp
Copyright (C) 2008 by Bill Rogers - Research Imaging Center - UTHSCSA
*/

/*! \mainpage
This program determines width of a sulcus based on the sulcal mesh file
(generally Tmtktri.mesh), the sulcal graph file, and the cortical mesh
file. If a single sulcal name is specified on the command line, then only
that sulcus will be measured. If no name is specified then every sulcus
in the graph file with a label name will be measured.

The program outputs average values for sulcal width, the standard deviation,
the cortical thickness along the sulcus (from cortical thickness texture file)
along with the number of measurements taken for each sulcus.
In the verbose options is selected a mesh is output for each sulcus containing
the sulcal span vectors.

If the optional transformation matrix is specified, then all the width measurements
will be transformed with the matrix inverse transform. This is used to transform the
measurements from Talairach space back to the origional space.

Command line switches:

--gm	cortical mesh file (required)

--sm	sulcal mesh file (required)

--sg	sulcal graph file (required)

--gt	cortical thickness texture file

--sn	sulcal name

--tm	transformation matrix

--mind	min distance allowed for width measurements

--maxd	max distance allowed for width measurements

--mina	min triangle area for sampling sulcal mesh

--maxa	max triangle area for sampling sulcal mesh

--sr	Maximum sulcal span segment ratio

--nsd	Number of standard deviations about avg to include

--aid	AID

--age	subject age

--sex	subject sex

-o		output file base name

-n		subject name

-v		verbose output

--bv	output for BrainVisa

--pk	special output for PeterK

 */

#include <string>
#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <cstdlib>
#include <time.h>
#include "RicMeshSet.h"
#include "RicTextureSet.h"
#include "RicUtil.h"
#include "RicGraph.h"
#include "SulcalRayTracer.h"

using namespace std;



#define MAXLABEL 100	///< max number of graph nodes allowed per sulcus
#define MINSIZE 150		///< minimum sulcus size to process

int FindTimeStep(RicMeshSet *pMesh, int tstep);

int main(int argc, char *argv[])
{
	///////// command line variables /////////////////

	string gm_name; // GM mesh name
	string sg_name; // sulcal graph name
	string sm_name; // sulcal mesh name
	string gt_name;	// cortical thickness texture name
	string outname; // base name for output files
	string sulcal_name; // sulcal name
	string matname; // name for transformation matrix file
	string sname;	// subject name
	string aid;		// subject AID name
	string age;		// subject age
	string sex;		// subject sex
	float mind = 0.0f; // minimum sulcal distance
	float maxd = 10.0f; // maximum sulcal distance
	float mina = 1.0f; // minimum sulcal mesh triangle area
	float maxa = 2.0f; // maximum sulcal mesh triangle arrea
	float spanratio = 4.0; // max allowed ratio of sulcal span segments
	float num_std = 3.0; // num of std dev about avg to allow
	bool transmat = false; // if true then use transformation matrix
	bool verbose = false; // if true output more stuff on stdout
	bool brainvisa = false; // if true output for BrainVisa
	bool peterk = false;	// if true then output for Peter K

	///////// read in the command line - tclap ///////////////

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("Acme-Presto Sulcal Width Processor", ' ', "0.66");

		// gm mesh
		TCLAP::ValueArg<string> gmname("", "gm", "Gray matter mesh file", true,
				"", "string");
		cmd.add(gmname);

		// sulcal mesh
		TCLAP::ValueArg<string> smname("", "sm", "Sulcal mesh file", true, "",
				"string");
		cmd.add(smname);

		// sulcal graph
		TCLAP::ValueArg<string> sgname("", "sg", "Sulcal graph file", true, "",
				"string");
		cmd.add(sgname);

		// gm mesh
		TCLAP::ValueArg<string> gtname("", "gt", "Gray matter thickness texture file", false,
				"", "string");
		cmd.add(gtname);

		// sulcal name
		TCLAP::ValueArg<string> sulname("", "sn", "Sulcal name", false, "",
				"string");
		cmd.add(sulname);

		// output file base name
		TCLAP::ValueArg<string> oname("o", "o", "Base output file name", false,
				"SulcalSpan", "string");
		cmd.add(oname);

		// subject name
		TCLAP::ValueArg<string> subjname("n", "name", "Subject name", false,
				"SSpanSubj", "string");
		cmd.add(subjname);

		// AID name
		TCLAP::ValueArg<string> aidname("", "aid", "AID", false,
				"SSpanAID", "string");
		cmd.add(aidname);

		// age name
		TCLAP::ValueArg<string> agename("", "age", "Subject age", false,
				"SSpanAge", "string");
		cmd.add(agename);

		// subject sex
		TCLAP::ValueArg<string> sexname("", "sex", "Subject sex", false,
				"SSpanSex", "string");
		cmd.add(sexname);

		// transformation matrix name
		TCLAP::ValueArg<string> mname("", "tm",
				"Transformation matrix file name", false, "", "string");
		cmd.add(mname);

		TCLAP::ValueArg<float> mnd("", "mind", "Minimum sulcal width (dflt 0.1)", false,
				0.1, "float");
		cmd.add(mnd);

		TCLAP::ValueArg<float> mxd("", "maxd", "Maximum sulcal width (dflt 10)", false,
				10, "float");
		cmd.add(mxd);

		TCLAP::ValueArg<float> mna("", "mina", "Minimum sulcal mesh triangle area (dflt 0.9)", false,
				0.9, "float");
		cmd.add(mna);

		TCLAP::ValueArg<float> mxa("", "maxa", "Maximum sulcal mesh triangle area (dflt 1.0)", false,
				1.0, "float");
		cmd.add(mxa);

		TCLAP::ValueArg<float> sr("", "sr", "Maximum sulcal span segment ratio (dflt 4.0)", false,
				4.0, "float");
		cmd.add(sr);

		TCLAP::ValueArg<float> nsd("", "nsd", "Number of standard deviations about avg to include (dflt 3.0)", false,
				3.0, "float");
		cmd.add(nsd);

		TCLAP::SwitchArg ver("v", "verbose", "Verbose output on stdout", false);
		cmd.add(ver);

		TCLAP::SwitchArg bv("", "bv", "Brain Visa output on stdout", false);
		cmd.add(bv);

		TCLAP::SwitchArg pk("", "pk", "Output in Peter K's special format", false);
		cmd.add(pk);

		// parse the command line
		cmd.parse(argc, argv);

		// copy command line variable to program variables
		gm_name = gmname.getValue();
		sm_name = smname.getValue();
		sg_name = sgname.getValue();
		gt_name = gtname.getValue();
		sulcal_name = sulname.getValue();
		outname = oname.getValue();
		sname = subjname.getValue();
		aid = aidname.getValue();
		age = agename.getValue();
		sex = sexname.getValue();
		matname = mname.getValue();
		mind = mnd.getValue();
		maxd = mxd.getValue();
		mina = mna.getValue();
		maxa = mxa.getValue();
		spanratio = sr.getValue();
		num_std = nsd.getValue();
		verbose = ver.getValue();
		brainvisa = bv.getValue();
		peterk = pk.getValue();

	} catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cerr << " Command Line Error" << e.error() << " for arg " << e.argId()
				<< endl;
	}

	////////////// Sanity checks on command line input ///////////

	// check for valid GM mesh file
	RicMesh gm_mesh((char*) gm_name.c_str());
	gm_mesh.calc_limits();
	if (gm_mesh.p_size != 0)
	{
		if (verbose)
			cerr << "We just read in " << gm_name << endl;
		if (verbose)
			cerr << "Number of vertices read: " << gm_mesh.v_size << endl;
	}
	else
	{
		cerr << "Error reading file" << gm_name << endl;
		exit(1);
	}

	// check for valid sulcal mesh file
	RicMeshSet s_mesh(sm_name);
	if (s_mesh.p_size != 0)
	{
		if (verbose)
			cerr << "We just read in " << sm_name << endl;
		if (verbose)
			cerr << "Number of meshes read: " << s_mesh.ntstep << endl;
	}
	else
	{
		cerr << "Error reading file" << sm_name << endl;
		exit(1);
	}

	// check for valid sulcal graph file
	RicGraph Sulcal_graph(sg_name);
	if (Sulcal_graph.nnodes != 0)
	{
		if (verbose)
			cerr << "We just read in " << sg_name << endl;
		if (verbose)
			cerr << "Number of nodes read: " << Sulcal_graph.nnodes << endl;
	}
	else
	{
		cerr << "Error reading file " << sg_name << endl;
		exit(1);
	}

	// check for valid cortical thickness texture file
	bool gt_flag = false;	// indicates a cortical thickness file
	RicTexture *gt_tex;
	gt_tex = NULL;
	if ( gt_name.length() != 0 )
	{
		gt_tex = new RicTexture;
		gt_tex->read_texture((char*)gt_name.c_str());
		if (gt_tex->size != 0)
		{
			if (verbose)
				cerr << "We just read in " << gt_name << endl;
			if (verbose)
				cerr << "Number of nodes read: " << gt_tex->size << endl;

			if ( gt_tex->size == gm_mesh.v_size ) // texture and mesh must have same size
				gt_flag = true;
			else
				cerr << "Cortical thickness texture wrong size" << endl;
		}
		else
		{
			cerr << "Error reading file " << sg_name << endl;
			exit(1);
		}
	}
	// Read the Transformation Matrix  if there is one
	CMatrix tm; // transformation matrix
	if (matname.length() != 0)
	{
		transmat = true; // we have a transformation matrix

		ifstream ifile;
		ifile.open(matname.c_str(), ifstream::in);
		if (!ifile.is_open())
		{
			cerr << "Error reading transform matrix " << matname << endl;
			exit(1);
		}

		if (verbose)
			cerr << "We just read in " << matname << endl;

		// read in the matrix to temp structure
		CMatrix tmat;
		for (int i = 0; i < 16; ++i)
			ifile >> tmat.mf[i];
		ifile.close();

		// transpose matrix to get proper orientation
		tm = tmat.Transpose();
	}

	//////////////////// Start the clock //////////////////////////
	// start a timer to see how long this all takes

	time_t starttime;
	time_t endtime;
	if ( verbose )
	{
		time(&starttime);
	}

	////////////////////// Specified Sulcus Option ////////////////////////////
	// check for special case when there is a sulcal name
	// then we just do that sulcus

	int n = 0, index = 0;
	int node_array[MAXLABEL];
	int index_array[MAXLABEL];
	float size_array[MAXLABEL];
	int total_verts=0;
	if (sulcal_name.length() != 0)
	{
		while ((n = Sulcal_graph.FindNextNodeLabel(sulcal_name)) >= 0)
		{
			if ( verbose ) cerr << "Found " << Sulcal_graph.nodes[n].label << " "
					<< Sulcal_graph.nodes[n].Tmtktri_label << " size "
					<< Sulcal_graph.nodes[n].size;

			// see if sulcus is big enough to be a keeper
			if (Sulcal_graph.nodes[n].size > MINSIZE)
			{
				node_array[index] = n;
				size_array[index] = Sulcal_graph.nodes[n].surface_area;
				// find mesh (time step) that matches Tmtktri label
				index_array[index++]
				            = FindTimeStep(&s_mesh, Sulcal_graph.nodes[n].Tmtktri_label);
				if  (verbose ) cerr << " Accepted " << endl;
				total_verts += s_mesh.mesh[n].v_size;
			}
			else
				if ( verbose ) cerr << " Rejected as too small" << endl;
		}

		if (index == 0)
		{
			cerr << "Couldn't find large enough sulci" << endl;
			exit(1);
		}

		SulcalRayTracer *ray_tracer;
		if ( gt_flag ) // if there is a boundary mesh cortical thickness map
			ray_tracer = new SulcalRayTracer(&gm_mesh, 20.0f, 0.5f, gt_tex);
		else
			ray_tracer = new SulcalRayTracer(&gm_mesh, 20.0f, 0.5f);

		// process all the sulci with the input sulcal name
		for (int i = 0; i < index; i++)
		{
			// simple insanity check to make sure that we're not off by 1 in the indexing
			if (fabs(s_mesh.mesh[index_array[i]].area() - size_array[i]) > 1.0)
			{
				cerr << "Error, two don't match "
					<< s_mesh.mesh[index_array[i]].area() << " and "
					<< size_array[i] << endl;
				cerr << "Tmtkri label=" << Sulcal_graph.nodes[node_array[i]].Tmtktri_label
					<< " Mesh time step=" << s_mesh.mesh[index_array[i]].t_step << endl;
			}

			RicMesh *sulcus;
			sulcus = s_mesh.mesh[index_array[i]].adjust_triangle_size(mina,
					maxa);

			if ( gt_flag ) // if there is a boundary mesh cortical thickness map
				ray_tracer->CalcSulcalSpanThick(sulcus);
			else
				ray_tracer->CalcSulcalSpan(sulcus);

			delete sulcus;

		}

		// if there is a transformation matrix then transform the end points
		if ( transmat ) ray_tracer->TransformEndpoints(tm);

		// calculate raw average
		ray_tracer->CalcMinMaxAvg();

		if ( verbose)
		{
			cout << "Min poly area=" << mina << " Max poly area=" << maxa << endl;

			cout << "Unfiltered average_thickness " << ray_tracer->avg_thickness;
			cout << " +/- " << ray_tracer->std_dev;
			cout << " number of vectors " << ray_tracer->number_steps << endl;
			if ( gt_flag)
				cout << "Average boundary mesh cortical thickness " << ray_tracer->boundary_thick << endl;
		}

		// filter endpoints by ratio
		ray_tracer->FilterEndpointsByRatio(spanratio);
		ray_tracer->CalcMinMaxAvg();

		if ( verbose )
		{
			cout << "Ratio filtered average_thickness " << ray_tracer->avg_thickness;
			cout << " +/- " << ray_tracer->std_dev;
			cout << " number of vectors " << ray_tracer->number_steps << endl;
			if ( gt_flag)
				cout << "Average boundary mesh cortical thickness " << ray_tracer->boundary_thick << endl;
		}

		// filter by std deviation
		ray_tracer->FilterEndpointsByStddev(num_std);
		ray_tracer->CalcMinMaxAvg();

		if ( verbose )
		{
			cout << "Stddev filtered average_thickness " << ray_tracer->avg_thickness;
			cout << " +/- " << ray_tracer->std_dev;
			cout << " number of vectors " << ray_tracer->number_steps << endl;
			if ( gt_flag)
				cout << "Average boundary cortical mesh thickness " << ray_tracer->boundary_thick << endl;

			// create an output file name
			string meshname;
			meshname = outname + "_" + sulcal_name + ".mesh";
			ray_tracer->WriteEndpointMesh(meshname);
		}

		// standard program output
		string outtext2;
		outtext2 = outname + "_" + sulcal_name + ".txt";
		FILE *fp2 = fopen((char*)outtext2.c_str(),"w");
		char buff2[512];
		if ( peterk )
		{
			sprintf(buff2,"%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%d\n",(char*)sulcal_name.c_str(),
					(char*)sname.c_str(),(char*)aid.c_str(),(char*)age.c_str(),(char*)sex.c_str(),
					ray_tracer->avg_thickness,ray_tracer->std_dev,ray_tracer->boundary_thick,
					ray_tracer->number_steps);
		}
		else
		{
			sprintf(buff2,"%s\t%f\t%f\t%f\t%d\n",(char*)sulcal_name.c_str(),
					ray_tracer->avg_thickness,ray_tracer->std_dev,ray_tracer->boundary_thick,
					ray_tracer->number_steps);
		}

		cout << buff2;
		fprintf(fp2,"%s\n",buff2);

		fclose(fp2);

		// lets see how long this took us to do
		if ( verbose )
		{
			time(&endtime);
			double elapsed = difftime(endtime, starttime);
			cerr << "Elapsed time: " << elapsed << endl;
		}
		delete ray_tracer;

		return EXIT_SUCCESS;
	}

	//////////////// end of single sulcus special case


	/////////////////// General Case, Do every sulcus ///////////////////////
	// The general case, go for every sulcus in the file measuring all those
	// that have a label name.

	// keep track of average sulcal span for all sulci
	double average,avg_tot=0;
	double std_dev,std_tot=0;
	double avg_boundary,boundary_tot=0;
	int tot_meas=0;	// total number of width measurments for all sulci
	int nsulci=0;

	// make a list of flags to indicate which nodes have been used
	int *flags;
	flags = new int[Sulcal_graph.nnodes];
	int i;
	for (i=0 ; i<Sulcal_graph.nnodes ; ++i ) flags[i]=0;

	// open the output file
	string outtext;
	outtext = outname + ".txt";
	FILE *fp = fopen((char*)outtext.c_str(),"w");

	// standard program output text
	char buff[512];

	// try every possible node
	for ( i=0 ; i<Sulcal_graph.nnodes ; ++i )
	{
		// skip emply, unknown or already used labels
		if ( Sulcal_graph.nodes[i].label == "unknown") continue;
		if ( Sulcal_graph.nodes[i].label.length() == 0 ) continue;
		if ( flags[i] == 1 ) continue;

		sulcal_name = Sulcal_graph.nodes[i].label;
		if ( verbose ) cerr << "Processing " << sulcal_name << endl;
		n = Sulcal_graph.FindNodeLabel(sulcal_name);
		flags[n] = 1;
		index = 0;
		total_verts = 0;

		// see if we have one that is big enough to be a keeper
		if (Sulcal_graph.nodes[n].size > MINSIZE)
		{
			node_array[index] = n;
			size_array[index] = Sulcal_graph.nodes[n].surface_area;
			// find mesh (time step) that matches Tmtktri label
			index_array[index++]
			            = FindTimeStep(&s_mesh, Sulcal_graph.nodes[n].Tmtktri_label);
			if ( verbose ) cerr << " Accepted " << endl;
			total_verts += s_mesh.mesh[n].v_size;
		}

		// now look for all the other sulci nodes with the same name
		while ((n = Sulcal_graph.FindNextNodeLabel(sulcal_name)) >= 0)
		{
			flags[n] = 1; // mark as used
			if ( verbose ) cerr << "Found " << Sulcal_graph.nodes[n].label << " "
					<< Sulcal_graph.nodes[n].Tmtktri_label << " size "
					<< Sulcal_graph.nodes[n].size;

			if (Sulcal_graph.nodes[n].size > MINSIZE)
			{
				node_array[index] = n;
				size_array[index] = Sulcal_graph.nodes[n].surface_area;
				// find mesh (time step) that matches Tmtktri label
				index_array[index++]
				            = FindTimeStep(&s_mesh, Sulcal_graph.nodes[n].Tmtktri_label);
				if ( verbose ) cerr << " Accepted " << endl;
				total_verts += s_mesh.mesh[n].v_size;
			}
			else
				if ( verbose ) cerr << " Rejected as too small" << endl;
		}

		if (index == 0)
		{
			if ( verbose ) cerr << "Couldn't find large enough sulci" << endl;
			continue;
		}

		SulcalRayTracer *ray_tracer;
		if ( gt_flag ) // if there is a boundary mesh cortical thickness map
			ray_tracer = new SulcalRayTracer(&gm_mesh, 20.0f, 0.5f, gt_tex);
		else
			ray_tracer = new SulcalRayTracer(&gm_mesh, 20.0f, 0.5f);

		// process all the sulci with the same name and big enough size
		for (int i = 0; i < index; i++)
		{
			// simple insanity check to make sure that we're not off by 1 in the indexing
			if (fabs(s_mesh.mesh[index_array[i]].area() - size_array[i]) > 1.0)
			{
				cerr << "Error, two don't match "
					<< s_mesh.mesh[index_array[i]].area() << " and "
					<< size_array[i] << endl;
				cerr << "Tmtkri label=" << Sulcal_graph.nodes[node_array[i]].Tmtktri_label
					<< " Mesh time step=" << s_mesh.mesh[index_array[i]].t_step << endl;
			}

			RicMesh *sulcus;
			sulcus = s_mesh.mesh[index_array[i]].adjust_triangle_size(mina,maxa);

			if ( gt_flag ) // if there is a boundary mesh cortical thickness map
				ray_tracer->CalcSulcalSpanThick(sulcus);
			else
				ray_tracer->CalcSulcalSpan(sulcus);


			delete sulcus;
		}

		// if there is a transformation matrix then transform the end points
		if ( transmat ) ray_tracer->TransformEndpoints(tm);

		// calculate raw average
		ray_tracer->CalcMinMaxAvg();

		// output a whole bunch of stuff if verbose output
		if ( verbose )
		{
			cout << "Min poly area=" << mina << " Max poly area=" << maxa << endl;

			cout << "Unfiltered average_thickness " << ray_tracer->avg_thickness;
			cout << " +/- " << ray_tracer->std_dev;
			cout << " number of vectors " << ray_tracer->number_steps << endl;
			if ( gt_flag)
				cout << "Average boundary mesh cortical thickness " << ray_tracer->boundary_thick << endl;
		}

		// filter by span ratio
		ray_tracer->FilterEndpointsByRatio(spanratio);
		ray_tracer->CalcMinMaxAvg();

		if ( verbose )
		{
			cout << "Ratio filtered average_thickness " << ray_tracer->avg_thickness;
			cout << " +/- " << ray_tracer->std_dev;
			cout << " number of vectors " << ray_tracer->number_steps << endl;
			if ( gt_flag)
				cout << "Average boundary mesh cortical thickness " << ray_tracer->boundary_thick << endl;
		}

		// filter by std deviation
		ray_tracer->FilterEndpointsByStddev(num_std);
		ray_tracer->CalcMinMaxAvg();

		if ( verbose )
		{
			cout << "Stddev filtered average_thickness " << ray_tracer->avg_thickness;
			cout << " +/- " << ray_tracer->std_dev;
			cout << " number of vectors " << ray_tracer->number_steps << endl;
			if ( gt_flag)
				cout << "Average boundary mesh cortical thickness " << ray_tracer->boundary_thick << endl;

			// write a vector mesh file for sulcal span end points if verbose
			// create an output file name
			string meshname;
			meshname = outname + "_" + sulcal_name + ".mesh";
			ray_tracer->WriteEndpointMesh(meshname);
		}

		// keep track of overall average
		avg_tot += ray_tracer->avg_tot;
		std_tot += ray_tracer->std_tot;
		boundary_tot += ray_tracer->boundary_thick;
		tot_meas += ray_tracer->number_steps;
		++nsulci;

		if ( peterk )
		{
			sprintf(buff,"%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%d\n",(char*)sulcal_name.c_str(),
					(char*)sname.c_str(),(char*)aid.c_str(),(char*)age.c_str(),(char*)sex.c_str(),
					ray_tracer->avg_thickness,ray_tracer->std_dev,ray_tracer->boundary_thick,
					ray_tracer->number_steps);
		}
		else
		{
			sprintf(buff,"%s\t%f\t%f\t%f\t%d\n",(char*)sulcal_name.c_str(),
					ray_tracer->avg_thickness,ray_tracer->std_dev,ray_tracer->boundary_thick,
					ray_tracer->number_steps);
			cout << buff;
		}

		fprintf(fp,"%s",buff);

		delete ray_tracer;

	}

	// output average for all sulci
	average = avg_tot/tot_meas;
	avg_boundary = boundary_tot/double(nsulci);

	// standard deviation for all sulci
	std_dev = sqrt(fabs(std_tot-avg_tot*avg_tot/(double)tot_meas)/((double)tot_meas-1.0));

	fprintf(fp,"Summary\t%f\t%f\t%f\t%d\n",average,std_dev,avg_boundary,tot_meas);

	fclose(fp);

	if ( peterk )
	{
		cout << sname<<"\t"<<aid<<"\t"<<age<<"\t"<<sex<<"\t"<<average<<"\t"
			<<std_dev<<"\t"<<avg_boundary<<"\t"<<tot_meas<<endl;
	}
	else
	{
		cout << sname<<"_Total\t"<<average<<"\t"<<std_dev<<"\t"<<avg_boundary
			<<"\t"<<tot_meas<<endl;
	}
	// lets see how long this took us to do
	if ( verbose )
	{
		time(&endtime);
		double elapsed = difftime(endtime, starttime);
		cerr << "Elapsed time: " << elapsed << endl;
	}
	return EXIT_SUCCESS;

}

/*!
This function returns the index of the mesh in a mesh set with
the passed time step value

@param pMesh - pointer to mesh set
@param tstep - tstep value to look for in array of meshes
@returns - array index of mesh with the time step value
*/
int FindTimeStep(RicMeshSet *pMesh, int tstep)
{
	for ( int i=0 ; i<pMesh->ntstep ; ++i )
	{
		if ( pMesh->mesh[i].t_step == tstep )
			return i;	// we found it
	}

	// we did not find it so return -1
	return -1;
}
