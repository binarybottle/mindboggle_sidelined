// --------------------------- ricvolumeconvert.cpp ------------------------
/*! \mainpage
Implementation file for RicVolumeConvert.
Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA

This program takes a volume file of one format and converts it to another
while making a valiant attempt to preserve voxel size, orientation, etc.
It will read and write GIS, ANALYZE, NEMA, and NIFTI files. (Friends don't
let friends uses ANALYZE!)

Command line switches

-i		Input volume name (required)

-t		Output format (GIS, ANALYZE, NEMA, NIFTI) (GIS default)

-o		output file name

-v		verbose output

--bv	output for BrainVisa

 */

using namespace std;
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <tclap/config.h>
#include <tclap/CmdLine.h>
#include <RicVolumeSet.h>
#include <RicUtil.h>

int main(int argc, char *argv[])
{

	///////// command line variables /////////////////

	string input_name;			// input volume name
	string output_name;			// output name
	string format;				// output format
	bool verbose=false;			// if true output more stuff on stdout
	bool brainvisa=false;		// if true output stuff for BrainVisa
	
	/// Read in the command line - tclap.

	// Every thing is wrapped in a try block
	try
	{
		// Program description
		TCLAP::CmdLine cmd("RIC Volume Transmogrification Program", ' ', "1.0");
		
		// sulcal mesh file
		TCLAP::ValueArg<string> ivol( "i","i", "Input volume file", true, "", "string");
		cmd.add(ivol);
		
		// output file base name
		TCLAP::ValueArg<string> oname( "o","o", "Base output file name", false, "", "string");
		cmd.add(oname);
		
		// output file type - overridden by output file name
		TCLAP::ValueArg<string> fname( "t","t", "Output format - GIS, NIFTI, NEMA, ANALYZE", false, "GIS", "string");
		cmd.add(fname);
		
		TCLAP::SwitchArg bv( "","bv", "BrainVisa output", false);
		cmd.add(bv);
		
		TCLAP::SwitchArg ver( "v","verbose", "Verbose output on stdout", false);
		cmd.add(ver);
		
		// parse the command line
		cmd.parse( argc, argv );
		
		// copy command line variable to program variables
		input_name = ivol.getValue();
		output_name = oname.getValue();
		format = fname.getValue();
		verbose = ver.getValue();
		brainvisa = bv.getValue();
	
	}catch (TCLAP::ArgException &e) // catch exceptions - command line mistakes
	{
		cout << " Command Line Error " << e.error() << " for arg " << e.argId() << endl; 
	}
	
/// Check format for correctness.
	int ftype;
	if ( format == "GIS" || format == "gis" )
		ftype = RIC_GIS_FILE;
	else if ( format == "NEMA" || format == "nema" )
		ftype = RIC_NEMA_FILE;
	else if ( format == "NIFTI" || format == "nifti" )
		ftype = RIC_NIFTI_FILE;
	else if ( format == "ANALYZE" || format == "analyze" )
		ftype = RIC_ANALYZE_FILE;
	else
		ftype = RIC_NO_FILE;
	
	if ( ftype == RIC_NO_FILE )
	{
		cout << "RicVolumeConvert - invalid output format: " << format << endl;
		exit(1);
	}
		
	
/// Read the input file.
	RicVolumeSet *ivol;
	ivol = new RicVolumeSet(input_name);
	if ( ivol->filetype == RIC_NO_FILE )
	{
		cout << "RicVolumeConvert - Unknown file type: " << input_name << endl;
		exit(1);
	}
	
	ivol->VolSet[0].CalcMinMaxAvg();
	
	// disclaimer for Analyze
	if ( ivol->filetype == RIC_ANALYZE_FILE )
	{
		cout << "Analyze file so all bets are off on orientation!" << endl;
	}

/// Create our output file name.
	string output_file_name;
	
	// if no name specified then just input file for base name
	if ( output_name.length() == 0 )
	{
		// look for the extension
	  	size_t found;

  		found=input_name.rfind(".");
  		output_file_name = input_name.substr(0,found);
	}
	else // just use what they give you and extension overrides format switch
	{
		output_file_name = output_name;
		
		// figure out what the extensions is
	  	size_t found;
	  	string ext;

  		found=output_file_name.rfind(".");
  		
  		ext = output_file_name.substr(found+1);
		if ( ext == "DIM" || ext == "dim" || ext == "IMA" || ext == "ima" )
			ftype = RIC_GIS_FILE;
		else if ( ext == "DES" || ext == "des" )
			ftype = RIC_NEMA_FILE;
		else if ( ext == "NII" || ext == "nii" )
			ftype = RIC_NIFTI_FILE;
		else if ( ext == "HDR" || ext == "hdr" )
			ftype = RIC_ANALYZE_FILE;
		else
			ftype = RIC_NO_FILE;
		
		if ( ftype == RIC_NO_FILE )
		{
			cout << "RicVolumeConvert - invalid output format: " << format << endl;
			exit(1);
		}
	}

/// Write the output volume.
	if ( ftype == RIC_GIS_FILE )
	{
		ivol->filetype = RIC_GIS_FILE;
		ivol->Write_GIS(output_file_name);
	}
	else if ( ftype == RIC_NEMA_FILE )
	{
		ivol->filetype = RIC_NEMA_FILE;
		ivol->Write_NEMA(output_file_name);
	}
	else if ( ftype == RIC_NIFTI_FILE )
	{
		ivol->filetype = RIC_NIFTI_FILE;
		ivol->Write_NIFTI(output_file_name);
	}
	else if ( ftype == RIC_ANALYZE_FILE )
	{	
		cout << "Analyze file so all bets are off on orientation!" << endl;
		ivol->filetype = RIC_ANALYZE_FILE;
		ivol->Write_NIFTI(output_file_name);
	}
	
	// clean up memory
	delete ivol;
	
	if ( verbose || brainvisa ) cout << "File " << input_name 
		<< " successfully transmogrified to " << output_file_name << endl;
	
	exit (0);
}
