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

#define FONT_UNIT_HEIGHT	double(1.0/7.0)
#define FONT_UNIT_WIDTH		double(1.0/8.0)


#define PAGE_WIDTH	(612.0)
#define PAGE_HEIGHT	(792.0)
#define TOP_MARGIN				(PAGE_HEIGHT/10.0)
#define BOTTOM_MARGIN			TOP_MARGIN
#define LEFT_MARGIN				(PAGE_WIDTH/20.0)
#define RIGHT_MARGIN			LEFT_MARGIN

#define MIN_BASE_HEIGHT		(0.1)
#define COLUMN_WIDTH		(PAGE_WIDTH / 80)
#define LOGO_HEIGHT			(COLUMN_WIDTH * 5)
#define BASE_FONT_SIZE		( COLUMN_WIDTH * 0.9 / FONT_UNIT_WIDTH )

struct s_rgb
{
	double r, g, b;
};

// nucleotide order: "ACGT";
static s_rgb __n_letter_colors[] = 
	{
		{ 8.0, 0.0, 0.0 },
		{ 0.0, 0.0, 8.0 },
		{ 1.0, 0.7, 0.0 },
		{ 0.0, 0.5, 0.0 }
	};
	
void ps_scale( ofstream &ofs, double x_scale, double y_scale );
void ps_write_string_rotated( ofstream &ofs, double x, double y, const char *string, double red, double green, double blue, double angle );
void ps_write_string( ofstream &ofs, const char *string, double red, double green, double blue );
void ps_write_char( ofstream &ofs, char ch, double red, double green, double blue );
void ps_translate( ofstream &ofs, double x, double y );
void ps_rotate( ofstream &ofs, double angle );
void ps_line( ofstream &ofs, double x1, double y1, double x2, double y2, double width, double red, double green, double blue );
void ps_box_path( ofstream &ofs, double left, double bottom, double width, double height, double red, double green, double blue );
void ps_font( ofstream &ofs, const char *font, double fontsize );
void write_motif( ofstream &ofs, double *PWM[], int PWM_width, double axis_fontsize, int reverse );
