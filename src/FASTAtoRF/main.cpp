// Copyright 2014 Hamed S. Najafabadi

/********************************************************************

This file is part of FASTAtoRF.

FASTAtoRF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FASTAtoRF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FASTAtoRF.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

extern const char *__FASTA_file;
extern const int *__rnd;
extern const int *__span;
extern const int *__minl;
extern const int *__maxl;
extern const char *__output_file;


void welcome_message()
{
	cout << endl
		<< "**************************** FASTAtoRF version 1.01 ***************************" << endl
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

	ofstream ofs;
	if( !open_output( ofs, __output_file, "" ) )
		return 1;

	//******************* open the FASTA file

	s_gene *genes[ MAX_GENES ];
	int num_genes = 0;

	cout << "Opening the FASTA file..." << endl;
	if( !open_FASTA( __FASTA_file, genes, &num_genes ) )
		return 1;
		
	cout << "Finding C2H2 zinc fingers..." << endl;
	find_zfs( genes, num_genes );
	
	if( *__rnd )
	{
		cout << "Randomizing C2H2 sequences..." << endl;
		randomize_zfs( genes, num_genes );
	}
	
	//******************* write the output
	cout << "Writing the output..." << endl;
	write_ZFs( ofs, genes, num_genes, *__span, *__minl, *__maxl );
	
	cout << endl << "Job finished successfully." << endl;

	return 0;
}
