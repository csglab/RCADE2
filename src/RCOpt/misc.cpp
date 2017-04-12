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

///////////////////////////////////////////////////////////////////////////////////////////
double normal_cum_p( double z )
// returns P( Z > z )
// this code is taken from http://finance.bi.no/~bernt/gcc_prog/recipes/recipes/node23.html
{
	if( z > 200 )
		z = 200;
	if( z < -200 )
		z = -200;
		
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	
	double a = fabs( z );
	double t = 1.0 / ( 1.0 + a*p );
	double b = c2 * exp( - z*z / 2.0 );
	double n = ( ( ( ( b5 * t + b4 ) * t + b3 ) * t + b2 ) * t + b1 ) * t;
	n *= b;
	
	if( z < 0.0 )
		n = 1.0 - n;
		
	return n;
}

///////////////////////////////////////////////////////////////////////////////////////////
int find_motif(
	const char *name,
	s_motif *motifs[],
	int num_motifs )
// returns num_motifs if the requested motif is not found
{
	int i;
	for( i = 0; i < num_motifs; i ++ )
		if( strcmp( name, motifs[ i ] ->name ) == 0 )
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
