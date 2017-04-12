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
#include <time.h>

#include "declarations.h"


///////////////////////////////////////////////////////////////////////////////////////////
bool open_output(
	ofstream &ofs,
	const char *path,
	const char *extension )
// opens an output file stream
// returns true if successful, false otherwise
//
// ofs: the file stream that will point to the output file
// path: the path and the name of the file that will be openned
// extension: the extension of the output file. if path does not contain this extension, it will be added
// example: path = "output", extension = ".txt" --> output = "output.txt"
{
	// get the length of the path and extension strings
	unsigned int path_len = strlen( path );
	unsigned int ext_len = 0;
	if( extension ) // extension is request
		ext_len = strlen( extension );

	// process the name of the output file
	char output[ MAX_STRING_LENGTH ];
	strcpy( output, path ); // copy the path into the output file name
	if( ext_len &&
		( path_len < ext_len ||
		strcmp( output + path_len - ext_len, extension ) != 0 ) )
	// a file extension is requested, and is not already appended to the end of the output file name
		strcat( output, extension );

	ofs.open( output ); // open the output file
	if( !ofs || !ofs.is_open() )
	// the file cannot be openned
	{
		cout << "ERROR: " << output << " cannot be opened." << endl;

		return false;
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
bool _open_input( ifstream &ifs, const char *path, const char *extension )
// opens an input file stream
// returns true if successful, false otherwise
// extension will be used only if 'path' itself does not exist
//
// ifs: the file stream that will point to the input file
// path: the path and the name of the input file
// extension: the default extension of the input file. if the input file itself cannot be openned, this will be appended to the end of the input file name, and will be tried again
{
	char input[ MAX_STRING_LENGTH ];
	strcpy( input, path ); // copy the path to the input file name

	// first attempt
	ifs.open( input, ios::in );
	if( !ifs || !ifs.is_open() )
	// the first attempt was unsuccessful
	{
		// second attempt
		ifs.clear(); // clear the error flag

		cout << "WARNING: " << input << " was not found. ";
		strcat( input, extension ); // append the extension to the end of input file name
		cout << "Trying " << input << "..." << endl;
		
		ifs.open( input, ios::in );
		if( !ifs || !ifs.is_open() )
		// the second attempt also failed
		{
			cout << "ERROR: " << input << " was not found." << endl;

			return false;
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
char _find_EOL( ifstream &ifs )
// finds and returns the EOL that is used in the specified file stream
// returns 0 if no EOL is found
//
// ifs: the input file stream
{
	// the list of candidate EOL characters that will be examined
	char EOL_candidates[] = { '\n', '\r' };
	int num_candidates = 2;
	
	// go to the begnning of the file
	ifs.seekg( 0, ios::beg );

	// find the EOL candidate that happens fist in the file
	while( true )
	{
		// read one character
		char ch = ifs.get();
		if( !ifs || ifs.eof() )
			break; // the end of file is reached
			
		// see whether this character matches an EOL candidate
		int i;
		for( i = 0; i < num_candidates; i ++ )
			if( ch == EOL_candidates[ i ] ) // a matching EOL candidate is found
			{
				cout << "EOL candidate no. " << i+1 << " was approved." << endl;
				ifs.seekg( 0, ios::beg ); // go to the begnning of the file
				return ch; // return this character as the EOL that is used in this file
			}
	}
	
	ifs.clear(); // clear the error flag
	ifs.seekg( 0, ios::beg ); // go to the begnning of the file
	return 0; // no EOL character was found in the file
}


///////////////////////////////////////////////////////////////////////////////////////////
bool open_rndforest(
	const char *file,
	s_motif *motifs[],
	int *num_motifs,
	char const *filter )
// opens a file and reads its content
// returns true if successful, false otherwise
{
	// open the file
	ifstream ifs;
	_CALL(
		_open_input( ifs, file, ".txt" )
		);

	// find the EOL character
	char this_EOL = _find_EOL( ifs );
	if( !this_EOL )
	// no EOL character was identified
	{
		cout << "ERROR: The file has no EOL character." << endl;
		return false;
	}

	// read the file
	cout << "Reading RF results..." << endl;
	_CALL(
		read_rndforest( ifs, motifs, num_motifs, filter, this_EOL )
		);
		
	if( *num_motifs )
		cout << *num_motifs << " motifs were read." << endl << endl;
	else
	{
		cout << "ERROR: No motifs were read. The input protein possibly does not have canonical C2H2-ZFs." << endl;
		return false;
	}

	// close the input file stream
	ifs.close();

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
bool open_weights(
	const char *weight_file,
	s_motif *motifs[],
	int num_motifs )
{
	// open the file
	ifstream ifs;
	_CALL(
		_open_input( ifs, weight_file, ".txt" )
		);

	// find the EOL character
	char this_EOL = _find_EOL( ifs );
	if( !this_EOL )
	// no EOL character was identified
	{
		cout << "ERROR: The file has no EOL character." << endl;
		return false;
	}

	// read the file
	cout << "Reading the weights..." << endl;
	_CALL(
		read_weights( ifs, motifs, num_motifs, this_EOL )
		);

	ifs.close();

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
bool open_FASTA(
	const char *file,
	s_seq *seqs[],
	int *num_seqs )
// opens a file and reads its content
// returns true if successful, false otherwise
//
{
	// open the file
	ifstream ifs;
	_CALL(
		_open_input( ifs, file, ".fasta" )
		);

	// find the EOL character
	char this_EOL = _find_EOL( ifs );
	if( !this_EOL )
	// no EOL character was identified
	{
		cout << "ERROR: The file has no EOL character." << endl;
		return false;
	}

	// read the file
	cout << "Reading the sequences..." << endl;
	_CALL(
		read_FASTA( ifs, seqs, num_seqs, this_EOL )
		);
		
	// count the number of positive instances
	int num_pos = 0;
	int i;
	for( i = 0; i < *num_seqs; i ++ )
		num_pos += seqs[ i ] ->positive;
		
	cout << *num_seqs << " sequences were read, including " << num_pos << " positive and "
		<< *num_seqs - num_pos << " negatives." << endl << endl;
		
	if( num_pos < 10 || *num_seqs - num_pos < 10 )
	{
		cout << "ERROR: At least 10 positive and 10 negative sequences are required." << endl;
		return false;
	}

	// close the input file stream
	ifs.close();

	return true;
}
