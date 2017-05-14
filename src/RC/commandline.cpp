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


#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"
#include "commandline.h"


/******* This is the order of variables for structure s_command_parameter:

	const char *notes;
	const char *prefix;
	int type;
	const bool mandatory;
	bool set;
	int v_int;
	double v_double;
	char v_string[ MAX_STRING_LENGTH ];

*******/


/////////////////// Here is where the commandline parameters must be defined /////////////
////////////////// This part is all that needs to be modified in a new program ///////////

static s_command_parameter _command_parameters[] = {
	{	"Input file containing the randomForest results",
		"-rf", _TYPE_STRING, true, false, 0, 0.0, "None" },
	{	"The filter for ZFPs",
		"-filter", _TYPE_STRING, true, true, 0, 0.0, "" },
	{	"Input file containing the ZF weights",
		"-w", _TYPE_STRING, true, true, 0, 0.0, "None" },
	{	"Enrichment fold expected for each ZF (negative means PWMs should not be normalized)",
		"-enf", _TYPE_DOUBLE, true, true, 0, 50.0, "None" },
	{	"Output file",
		"-out"  , _TYPE_STRING, true, false, 0, 0.0, "None" } };

int _num_parameters = 5;

const char *__rf_file = _command_parameters[ 0 ].v_string;
const char *__filter = _command_parameters[ 1 ].v_string;
const char *__weight_file = _command_parameters[ 2 ].v_string;
const double *__enrichment_fold = & _command_parameters[ 3 ].v_double;
const char *__output_file = _command_parameters[ 4 ].v_string;

///////////////////////////////////////////////////////////////////////////////////////////



#define _CHECK_ARG	{ if( i >= argc-1 ){ \
	cout << "ERROR: Command line parameter \"" << argv[i] \
	<< "\" requires an argument." << endl << endl; \
	return false; } }



///////////////////////////////////////////////////////////////////////////////////////////
int _find_parameter( char *prefix )
// finds the parameter that has the specified prefix
// returns the parameter index if successful, otherwise returns _num_paramters
//
// prefix: the prefix of the parameter that is being searched
{
	int i;
	for( i = 0; i < _num_parameters; i ++ )
	{
		if( strcmp( _command_parameters[ i ].prefix, prefix ) == 0 )
			break;
	}
	
	return i;
}



///////////////////////////////////////////////////////////////////////////////////////////
bool read_arguments(
	int argc,
	char* argv[] )
// reads the arguments from command line input
// returns true if successful, false otherwise
//
// argc: the number of arguments
// argv: the argument list
{
	int i;
	for( i = 1; i < argc; i ++ )
	{
		int parameter_index = _find_parameter( argv[i] );
		if( parameter_index >= _num_parameters )
		{
			cout << "ERROR: Unknown command line parameter: " << argv[i] << endl << endl;
			return false;
		}
		
		switch( _command_parameters[ parameter_index ].type )
		{
			case _TYPE_INT:

				_CHECK_ARG;
		
				// convert the next argument to an integer number
				{
					std::string test_string( argv[i+1] ); // create a string object
					std::istringstream iss( test_string ); // create a string stream

					if( !( iss >> _command_parameters[ parameter_index ].v_int ) )
					// the conversion was not successful
					{
						cout << "ERROR: " << argv[i+1] << " is not an integer number." << endl;
						return false;
					}
				}	
				
				// ignore the next argv, as it was already read as the argument for this parameter
				i ++;

				break;

			case _TYPE_DOUBLE:

				_CHECK_ARG;
				
				// convert the next argument to an integer number
				{
					std::string test_string( argv[i+1] ); // create a string object
					std::istringstream iss( test_string ); // create a string stream

					if( !( iss >> _command_parameters[ parameter_index ].v_double ) )
					// the conversion was not successful
					{
						cout << "ERROR: " << argv[i+1] << " is not a double-precision number." << endl;
						return false;
					}
				}	
				
				// ignore the next argv, as it was already read as the argument for this parameter
				i ++;

				break;

			case _TYPE_STRING:

				_CHECK_ARG;
				
				strcpy( _command_parameters[ parameter_index ].v_string, argv[i+1] );
				
				// ignore the next argv, as it was already read as the argument for this parameter
				i ++;

				break;			
		}

		_command_parameters[ parameter_index ].set = true;
	}

	for( i = 0; i < _num_parameters; i ++ )
		if( _command_parameters[ i ].mandatory && !_command_parameters[ i ].set )
		// this parameter is mandatory, but is not determined
		{
			cout << "ERROR: Command line parameter \"" << _command_parameters[ i ].prefix
				<< "\" is mandatory." << endl << endl;
				
			return false;
		}
	
	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////
void print_arguments()
// prints the arguments that are read from the command line input
{
	cout << "Received parameters:" << endl;
	
	int i;
	for( i = 0; i < _num_parameters; i ++ )
	{
		cout << _command_parameters[ i ].notes << ": ";
		switch( _command_parameters[ i ].type )
		{
			case _TYPE_INT: cout << _command_parameters[ i ].v_int << endl; break;
			case _TYPE_DOUBLE: cout << _command_parameters[ i ].v_double << endl; break;
			case _TYPE_STRING: cout << _command_parameters[ i ].v_string << endl; break;
		}
	}
	
	cout << endl;
}



///////////////////////////////////////////////////////////////////////////////////////////
void print_commandline_format()
// prints the list of command line paramters and arguments
{
	cout << "::::::::::::::::::::::: List of command line parameters :::::::::::::::::::::::"
		<< endl << endl;
	
	int i;
	for( i = 0; i < _num_parameters; i ++ )
	{
		cout << _command_parameters[ i ].prefix << ": "
			<< _command_parameters[ i ].notes << "." << endl;
		
		switch( _command_parameters[ i ].type )
		{
			case _TYPE_INT:
				cout << "Type: Integer" << endl;
				cout << "Default: " << _command_parameters[ i ].v_int << endl;
				break;
				
			case _TYPE_DOUBLE:
				cout << "Type: Double-precision floating point" << endl;
				cout << "Default: " << _command_parameters[ i ].v_double << endl;
				break;
				
			case _TYPE_STRING:
				cout << "Type: String" << endl;
				cout << "Default: " << _command_parameters[ i ].v_string << endl;
				break;
		}
		
		cout << endl;
	}
	
	cout << endl;
}
