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

#include "declarations.h"

static int _line_number; // this local variable will hold the number of lines that are read so far
	

/////////////////////////////////////////////////////////////////////////////////
bool _add_seq(
	char seq_name[],
	char seq[],
	s_seq *seqs[],
	int *num_seqs )
// adds an entry to the seq list
//
// seq_name: the name of the new seq
{
	if(	!seq[ 0 ] ) // there is no associated sequence, and thus no need to do anything
		return true; // return with no error
		
	if( !seq_name[ 0 ] ) // there is an associated sequence, but no name
	{
		cout << "ERROR: No sequence name provided (error at line " << _line_number << ")"
			<< endl;
		return false;
	}
	
	if( *num_seqs >= MAX_SEQS ) // too many sequences
	{
		cout << "ERROR: Maximum number of sequences (" << MAX_SEQS << ") reached "
			<< "(error at line " << _line_number << ")" << endl;
		return false;
	}
	
	seqs[ *num_seqs ] = new s_seq;
	_COPY_STR( seqs[ *num_seqs ] ->name, seq_name );
	_COPY_STR( seqs[ *num_seqs ] ->seq, seq );

	seqs[ *num_seqs ] ->seq_length = strlen( seq );
	
	// decide whether this belongs to the positive or negative set, based on the name
	if( strlen( seq_name ) >= 6 && strcmp( seq_name + strlen( seq_name ) - 6, "-dinuc" ) == 0 ) // this is a dinuc-shuffled sequence
		seqs[ *num_seqs ] ->positive = 0;
	else
		seqs[ *num_seqs ] ->positive = 1;
	
	(*num_seqs) ++;
	
	return true;
}


/////////////////////////////////////////////////////////////////////////////////
bool read_FASTA(
	ifstream &ifs,
	s_seq *seqs[],
	int *num_seqs,
	char this_EOL )
// reads sequences from &ifs
// returns true if successful, false otherwise
//
// ifs: input file stream
// seqs: the list of seqs whose sequences are to be read
// num_seqs: will contain the number of seqs that are in the list
// this_EOL: the end of line character. this character is determined separately in another function for this file stream
{
	// A line cannot be longer than MAX_LINE_LENGTH characters
	char string[ MAX_LINE_LENGTH + 1 ];
	
	char seq_name[ MAX_LINE_LENGTH + 1 ] = "";
	char seq[ MAX_SEQ_LENGTH + 1 ] = "";

	// read the file line by line
	_line_number = 0; // reset the line number
	
	while( true )
	{
		// read the line
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );

		if( !ifs || ifs.eof() ) // the end of the file has reached
		{
			// add the last read seq to the list
			_CALL(
				_add_seq( seq_name, seq, seqs, num_seqs )
				);

			break;
		}
			
		_line_number ++; // update the line number that was just read

		if( string[ 0 ] == '>' ) // this is a header
		{
			// add the last read seq to the list
			_CALL(
				_add_seq( seq_name, seq, seqs, num_seqs )
				);
		
			// store the name of the new seq
			strcpy( seq_name, string + 1 );
			seq[ 0 ] = 0; // reset the seq sequence
		}
		else // this line contains sequence information
			strcat( seq, string );
	}
	
	return true;
}
