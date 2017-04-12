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

///////////////////////////////////////////////////////////////////////////////////////////
int find_gene(
	const char *name,
	s_gene *genes[],
	int num_genes )
// returns num_genes if the requested gene is not found
{
	int i;
	for( i = 0; i < num_genes; i ++ )
		if( strcmp( name, genes[ i ] ->name ) == 0 )
			break;

	return i;
}

///////////////////////////////////////////////////////////////////////////////////////////
bool _is_delimiter(
	char ch,
	const char *delimiters,
	int num_delimiters )
// determines whether a character is a delimiter or not
// returns true if ch is a delimiter, false otherwise
//
// ch: the query character
// delimiters: the list of delimiters
// num_delimiters: the size of the delimiter list
{
	// examine the delimiters one by one to find a match
	for( int i = 0; i < num_delimiters; i ++ )
		if( ch == delimiters[ i ] )
			return true; // a matching delimiter is found

	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////
char* extract_phrase(
	const char *line,
	int phrase_index,
	char *phrase,
	const char *delimiters )
// extracts a phrase from a delimitted text
// returns NULL if unsuccessful, 'phrase' otherwise
//
// line: the line from which the phrase will be extracted
// phrase_index: the index of the phrase that will be extracted; 0 means that first phrase, and so on
// phrase: the string in which the extracted phrase will be stored
// delimiters: the list of valid delimiters, separating different phrases within the text
{
	int len = strlen( line );
	int num_delimiters = strlen( delimiters );

	// find the first position after 'phrase_index'th occurrance of a delimiter
	int curr_index = 0;
	int start_pos;
	for( start_pos = 0; start_pos < len && curr_index < phrase_index; start_pos ++ )
		if( _is_delimiter( line[ start_pos ], delimiters, num_delimiters ) )
			curr_index ++;

	// return NULL if the requested phrase is not found
	if( start_pos >= len )
		return NULL;

	// extract the phrase
	int pos;
	for( pos = start_pos; pos < len; pos ++ )
		if( _is_delimiter( line[ pos ], delimiters, num_delimiters ) )
			break;
		else
			phrase[ pos - start_pos ] = line[ pos ];

	phrase[ pos - start_pos ] = 0;

	// return the extracted phrase
	return phrase;
}
