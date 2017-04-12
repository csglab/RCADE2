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

void ps_scale( ofstream &ofs, double x_scale, double y_scale )
{
	ofs << x_scale << " " << y_scale << " scale" << endl;

}

void ps_write_string_rotated( ofstream &ofs, double x, double y, const char *string, double red, double green, double blue, double angle )
{
	ofs << x << " " << y << " translate" << endl
		<< angle << " rotate" << endl
		<< "newpath" << endl
		<< 0 << " " << 0 << " moveto" << endl
		<< "(" << string << ") true charpath" << endl
		<< "closepath" << endl
		<< red << " " << green << " " << blue << " setrgbcolor" << endl
		<< "fill" << endl
		<< -angle << " rotate" << endl
		<< -x << " " << -y << " translate" << endl;
}

void ps_write_string( ofstream &ofs, const char *string, double red, double green, double blue )
{
	if( !string )
		return;
		
	ofs << "newpath" << endl
		<< "0 0 moveto" << endl
		<< "(" << string << ") true charpath" << endl
		<< red << " " << green << " " << blue << " setrgbcolor" << endl
		<< "fill" << endl;
}

void ps_write_char( ofstream &ofs, char ch, double red, double green, double blue )
{
	ofs << "newpath" << endl
		<< "0 0 moveto" << endl
		<< "(";

	if( ch == ')' || ch == '(' )
		ofs << "\\";
		
	ofs << ch << ") true charpath" << endl
		<< red << " " << green << " " << blue << " setrgbcolor" << endl
		<< "fill" << endl;
}

void ps_translate( ofstream &ofs, double x, double y )
{
	ofs << x << " " << y << " " << "translate" << endl;
}

void ps_rotate( ofstream &ofs, double angle )
{
	ofs << angle << " rotate" << endl;
}


void ps_line( ofstream &ofs, double x1, double y1, double x2, double y2, double width, double red, double green, double blue )
{
	ofs << "newpath" << endl;
	ofs << x1 << " " << y1 << " moveto" << endl;
	ofs << x2 << " " << y2 << " lineto" << endl;
	ofs << width << " setlinewidth" << endl;
	ofs << red << " " << green << " " << blue << " setrgbcolor" << endl;
	ofs << "stroke" << endl;
}


void ps_box_path( ofstream &ofs, double left, double bottom, double width, double height, double red, double green, double blue )
{
	ofs << "newpath" << endl
		<< left << " " << bottom << " moveto" << endl
		<< width << " " << 0 << " rlineto" << endl
		<< 0 << " " << height << " rlineto" << endl
		<< -width << " " << 0 << " rlineto" << endl
		<< 0 << " " << -height << " rlineto" << endl
		<< "closepath" << endl
		<< 0.2 << " setlinewidth" << endl
		<< red << " " << green << " " << blue << " setrgbcolor" << endl
		<< "stroke" << endl;
}

void ps_font( ofstream &ofs, const char *font, double fontsize )
{
	ofs << "/" << font << " findfont" << endl
	<< fontsize*0.2 << " scalefont" << endl
	<< "setfont" << endl;
}
