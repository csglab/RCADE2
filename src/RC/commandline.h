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

#ifndef _H_COMMANDLINE // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_COMMANDLINE // this macro indicates that this file is now included


#include "structures.h"


#define _TYPE_NULL		0  // this type of parameters does not need any argument
#define _TYPE_INT		1  // this type of parameters need an integer argument
#define _TYPE_DOUBLE	2  // this type of parameters need a double argument
#define _TYPE_STRING	3  // this type of parameters need a string argument

struct s_command_parameter
{
	const char *notes;
	const char *prefix;

	int type;
	
	const bool mandatory;
	bool set;
	
	int v_int;
	double v_double;
	char v_string[ MAX_STRING_LENGTH ];
};

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
