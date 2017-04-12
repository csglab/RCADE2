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


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

#include "structures.h"
#include "declarations.h"

#define _EOFCHECK(x)	if( !x || x.eof() ){ cout << "ERROR: Unexpected end of file." << endl << endl; return false; }
#define _PTRCHECK(x)	if( !x ){ cout << "ERROR: Unexpected end of line at line " << _line_number << "." << endl; return false; }

static int _line_number; // this local variable will hold the number of lines that are read so far

static char _tab[] = { char(9), 0 };
static char _semicolon[] = { ';', 0 };
static char _or[] = { '|', 0 };

///////////////////////////////////////////////////////////////////////////////////////////
bool read_rndforest(
	ifstream &ifs,
	s_motif *motifs[],
	int *num_motifs,
	char const *filter,
	char this_EOL )
// returns true if successful, false otherwise
{
	*num_motifs = 0;
	
	// A line cannot be longer than MAX_LINE_LENGTH characters
	char string[ MAX_LINE_LENGTH + 1 ];
	// read and ignore the first line
	ifs.getline( string, MAX_LINE_LENGTH, this_EOL );
	
	// read the file line by line
	_line_number = 0; // reset the line number
	while( true )
	{
		// read the line
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );
		if( !ifs || ifs.eof() ) // the end of the file has reached
			break;
		_line_number ++; // update the line number that was just read
		
		if( strncmp( string, "letter", 6 ) == 0 )
			continue; // this is a dummy line
			
		// get the motif/zf description
		char name[ MAX_STRING_LENGTH ];
		_PTRCHECK(
			extract_phrase( string, 0, name, _tab )
			);
			
		// match with filter, if provided
		if( filter && filter[ 0 ] &&								// a filter has been provided
			( strncmp( name, filter, strlen( filter ) ) != 0 ||		// either the filter doesn't match the start of the name
			name[ strlen( filter ) ] != '|' ) )						// or the character after the matchting string is not '|'
			continue; // ignore this entry

		// extract the motif name			
		char motif_name[ MAX_STRING_LENGTH ];
		_PTRCHECK(
			extract_phrase( name, 0, motif_name, _semicolon )
			);
			
		// find the motif that corresponds to this name
		int motif_index = find_motif( motif_name, motifs, *num_motifs );
		if( motif_index >= *num_motifs ) // this is a new motif
		{
			motifs[ motif_index ] = new s_motif;
			_COPY_STR( motifs[ motif_index ] ->name, motif_name );
			
			// extract the protein name			
			char protein_name[ MAX_STRING_LENGTH ];
			_PTRCHECK(
				extract_phrase( motif_name, 0, protein_name, _or )
				);
			_COPY_STR( motifs[ motif_index ] ->protein_name, protein_name );
			
			*num_motifs = motif_index + 1;
		}
			
		// extract the zf info
		char zf[ MAX_STRING_LENGTH ];
		_PTRCHECK(
			extract_phrase( name, 1, zf, _semicolon )
			);
		int zf_index = atoi( zf ) - 1;
		
		// extract the zf index within the original protein
		_PTRCHECK(
			extract_phrase( name, 2, zf, _semicolon )
			);
		int zf_in_protein = atoi( zf );
		if( motifs[ motif_index ] ->first_zf_in_protein < 0 ||
			motifs[ motif_index ] ->first_zf_in_protein > zf_in_protein )
			motifs[ motif_index ] ->first_zf_in_protein = zf_in_protein;
		if( motifs[ motif_index ] ->last_zf_in_protein < 0 ||
			motifs[ motif_index ] ->last_zf_in_protein < zf_in_protein )
			motifs[ motif_index ] ->last_zf_in_protein = zf_in_protein;
		
		
		// create the zf object
		motifs[ motif_index ] ->zfs[ zf_index ] = new s_c2h2;
		if( zf_index >= motifs[ motif_index ] ->num_zfs )
			motifs[ motif_index ] ->num_zfs = zf_index + 1;
			
		// read the zf PWM
		int x, n;
		int col = 13;
		for( x = 0; x < 4; x ++ )
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				char value[ MAX_STRING_LENGTH ];
				_PTRCHECK(
					extract_phrase( string, col, value, _tab )
					);
				motifs[ motif_index ] ->zfs[ zf_index ] ->PWM[ n ][ x ] = (x==3)? 0 : log10( atof( value ) * 4 + 0.0001 );
				col ++;
			}
	}
	
	// if there are any missing ZFs for each motif, create empty ZFs
	int i, j;
	for( i = 0; i < *num_motifs; i ++ )
		for( j = 0; j < motifs[ i ] ->num_zfs; j ++ )
			if( !motifs[ i ] ->zfs[ j ] )
				motifs[ i ] ->zfs[ j ] = new s_c2h2;

	// replace the motif names with user-friendly names
	for( i = 0; i < *num_motifs; i ++ )
	{
		char new_name[ MAX_STRING_LENGTH ];
		sprintf( new_name, "%s:%i-%i", motifs[ i ] ->protein_name,
			motifs[ i ] ->first_zf_in_protein, motifs[ i ] ->last_zf_in_protein );
		RELEASE( motifs[ i ] ->name );
		_COPY_STR( motifs[ i ] ->name, new_name );
	}
	
	return true;
}
