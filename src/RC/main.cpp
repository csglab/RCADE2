// Copyright 2014 Hamed S. Najafabadi

/********************************************************************

This file is part of RCOpt.

RCOpt is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RCOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RCOpt.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

extern const char *__rf_file;
extern const char *__filter;
extern const char *__weight_file;
extern const double *__enrichment_fold;
extern const char *__output_file;


void welcome_message()
{
	cout << endl
		<< "******************************** RC version 1.01 ******************************" << endl
		<< "*.............................................................................*" << endl
		<< "* Copyright 2014 Hamed S. Najafabadi .........................................*" << endl
		<< "*.............................................................................*" << endl
		<< "*.............................................................................*" << endl
		<< endl;
}


int main( int argc, char* argv[] )
{	
	welcome_message();
	
	srand ( time(NULL) );

	if( argc <= 1 )
	// no argument is profided
	{
		print_commandline_format();
		return 0;
	}
	
	if( !read_arguments( argc, argv ) )
		return 1;
		
	print_arguments();

	//******************* open output

	ofstream ofs_PFM;
	if( !open_output( ofs_PFM, __output_file, ".PFM.txt" ) )
		return 1;

	ofstream ofs_ps;
	if( !open_output( ofs_ps, __output_file, ".ps" ) )
		return 1;

	//******************* open the result file

	s_motif *motifs[ MAX_MOTIFS ];
	int num_motifs = 0;

	cout << "Opening the randomForest output file..." << endl;
	if( !open_rndforest( __rf_file, motifs, &num_motifs, __filter ) )
		return 1;
		
	//******************* open the weight file

	if( strcmp( __weight_file, "None" ) != 0 ) // a weight file has been provided
	{
		cout << "Opening the weight file..." << endl;
		if( !open_weights( __weight_file, motifs, num_motifs ) )
			return 1;
	}
	
	//******************* calculations
	cout << "Generating PWMs..." << endl;
	if( !generate_PWMs( motifs, num_motifs, *__enrichment_fold ) )
		return 1;
	
	//******************* write the output files
	cout << "Writing to output..." << endl;
	
	// write all PFMs
	if( !write_PFMs( ofs_PFM, motifs, num_motifs ) )
		return 1;

	// write the postscript
	write_graphics( ofs_ps, motifs, num_motifs );

	cout << endl << "Job finished successfully." << endl;

	return 0;
}
