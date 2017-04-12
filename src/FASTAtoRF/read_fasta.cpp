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


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

static int _line_number; // this local variable will hold the number of lines that are read so far
	

/////////////////////////////////////////////////////////////////////////////////
bool _add_gene(
	char gene_name[],
	char gene_seq[],
	s_gene *genes[],
	int *num_genes )
// adds an entry to the gene list
//
// gene_name: the name of the new gene
{
	if(	!gene_seq[ 0 ] ) // there is no associated sequence, and thus no need to do anything
		return true; // return with no error
		
	if( !gene_name[ 0 ] ) // there is an associated sequence, but no name
	{
		cout << "ERROR: No sequence name provided (error at line " << _line_number << ")"
			<< endl;
		return false;
	}
	
	if( *num_genes >= MAX_GENES ) // too many sequences
	{
		cout << "ERROR: Maximum number of sequences (" << MAX_GENES << ") reached "
			<< "(error at line " << _line_number << ")" << endl;
		return false;
	}
	
	genes[ *num_genes ] = new s_gene;
	_COPY_STR( genes[ *num_genes ] ->name, gene_name );
	_COPY_STR( genes[ *num_genes ] ->seq, gene_seq );

	genes[ *num_genes ] ->seq_length = strlen( gene_seq );
	
	(*num_genes) ++;
	
	return true;
}


/////////////////////////////////////////////////////////////////////////////////
bool read_FASTA(
	ifstream &ifs,
	s_gene *genes[],
	int *num_genes,
	char this_EOL )
// reads sequences from &ifs
// returns true if successful, false otherwise
//
// ifs: input file stream
// genes: the list of genes whose sequences are to be read
// num_genes: will contain the number of genes that are in the list
// this_EOL: the end of line character. this character is determined separately in another function for this file stream
{
	// A line cannot be longer than MAX_LINE_LENGTH characters
	char string[ MAX_LINE_LENGTH + 1 ];
	
	char gene_name[ MAX_LINE_LENGTH + 1 ] = "";
	char gene_seq[ MAX_SEQ_LENGTH + 1 ] = "";

	// read the file line by line
	_line_number = 0; // reset the line number
	while( true )
	{
		// read the line
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );

		if( !ifs || ifs.eof() ) // the end of the file has reached
		{
			// add the last read gene to the list
			_CALL(
				_add_gene( gene_name, gene_seq, genes, num_genes )
				);

			break;
		}
			
		_line_number ++; // update the line number that was just read

		if( string[ 0 ] == '>' ) // this is a header
		{
			// add the last read gene to the list
			_CALL(
				_add_gene( gene_name, gene_seq, genes, num_genes )
				);
		
			// store the name of the new gene
			strcpy( gene_name, string + 1 );
			gene_seq[ 0 ] = 0; // reset the gene sequence
		}
		else // this line contains sequence information
			strcat( gene_seq, string );
	}
	
	return true;
}
