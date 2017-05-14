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
#define _PTRCHECK(x)	if( !x ){ cout << "ERROR: Unexpected end of line." << endl; return false; }

static int _line_number; // this local variable will hold the number of lines that are read so far

static char _delimiters[] = { char(9), 0 };

int _aa_index( char aa );

///////////////////////////////////////////////////////////////////////////////////////////
bool read_weights(
	ifstream &ifs,
	s_motif *motifs[],
	int num_motifs,
	char this_EOL )
// returns true if successful, false otherwise
{
	// initialize all weights as zero
	int i, j;
	for( i = 0; i < num_motifs; i ++ )
		for( j = 0; j < motifs[ i ] ->num_zfs; j ++ )
			motifs[ i ] ->zfs[ j ] ->weight = 0;
			
	// A line cannot be longer than MAX_LINE_LENGTH characters
	char string[ MAX_LINE_LENGTH + 1 ];
	
	// read the file line by line
	_line_number = 0; // reset the line number
	while( true )
	{
		// read the line
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );
		if( !ifs || ifs.eof() ) // the end of the file has reached
			break;

		_line_number ++; // update the line number that was just read
			
		if( string[ 0 ] == '#' || string[ 0 ] == 0 )
		// this is either a comment line or an empty line
			continue; // ignore the line
		
		// get the motif name	
		char motif_name[ MAX_STRING_LENGTH ];
		_PTRCHECK( extract_phrase( string, 0, motif_name, _delimiters ) );
		int motif_index = find_motif( motif_name, motifs, num_motifs );
		if( motif_index >= num_motifs )
		{
			cout << "WARNING at line " << _line_number << ": Gene " << motif_name
				<< " was not found in the sequence list." << endl;
			continue;
		}
		
		// get the ZF index
		char str_zf_index[ MAX_STRING_LENGTH ];
		_PTRCHECK( extract_phrase( string, 1, str_zf_index, _delimiters ) );
		int zf_index = atoi( str_zf_index ) - 1;
		if( zf_index >= motifs[ motif_index ] ->num_zfs )
		{
			cout << "WARNING at line " << _line_number << ": Gene " << motif_name
				<< " does not have a " << str_zf_index << "th zinc finger." << endl;
			continue;
		}
		
		// get the weight
		char weight[ MAX_STRING_LENGTH ];
		_PTRCHECK( extract_phrase( string, 2, weight, _delimiters ) );
		motifs[ motif_index ] ->zfs[ zf_index ] ->weight = atof( weight );
		if( motifs[ motif_index ] ->zfs[ zf_index ] ->weight < 0 )
			motifs[ motif_index ] ->zfs[ zf_index ] ->weight = 0;			
	}
	
	return true;
}
