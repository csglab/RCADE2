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
#include "graphics.h"

struct s_base
{
	int index;
	double value;
};

////////////////////////////////////////////////////////////////
int _compare_bases( const void *a, const void *b )
{
	s_base *base1 = (s_base *) a;
	s_base *base2 = (s_base *) b;
	double value1 = base1 ->value;
	double value2 = base2 ->value;
	
	if( value1 > value2 )
		return 1;
	if( value1 < value2 )
		return -1;
	
	return 0;
}


////////////////////////////////////////////////////////////////
void _write_x_axis(
	ofstream &ofs,
	int PWM_width,
	int reverse,
	double fontsize )
// write the horizontal axis values
{
	double x = 0;

	ps_translate( ofs, COLUMN_WIDTH / 2 - fontsize * FONT_UNIT_WIDTH / 2, 0 );
	x += ( COLUMN_WIDTH / 2 - fontsize * FONT_UNIT_WIDTH / 2 );

	ps_font( ofs, "Helvetica", fontsize );

	for( int i = 0; i < PWM_width; i ++ )
	{
		int x_value = i+1;
		if( reverse )
			x_value = PWM_width - i;
			
		ps_write_char( ofs, '0' + x_value%10, 0, 0, 0 );
		ps_translate( ofs, COLUMN_WIDTH, 0 );
		x += COLUMN_WIDTH;
	}
	ps_translate( ofs, -x, 0 );
}

	
////////////////////////////////////////////////////////////////
void _write_y_axis(
	ofstream &ofs,
	double max_sum_positive,
	double max_sum_negative,
	double fontsize,
	int PWM_width )
// write the horizontal axis values
{
	char string[ MAX_STRING_LENGTH ];
	ps_font( ofs, "Helvetica", fontsize );

	ps_translate( ofs, COLUMN_WIDTH * ( PWM_width + 1 ), - fontsize * FONT_UNIT_HEIGHT / 2.5 );
	sprintf( string, "0" );
	ps_write_string( ofs, string, 0, 0, 0 );
	
	ps_translate( ofs, 0, LOGO_HEIGHT / 2 );
	sprintf( string, "+%0.2f", max_sum_positive );
	ps_write_string( ofs, string, 0, 0, 0 );

	ps_translate( ofs, 0, - LOGO_HEIGHT );
	sprintf( string, "%0.2f", -max_sum_negative );
	ps_write_string( ofs, string, 0, 0, 0 );

	ps_translate( ofs, - COLUMN_WIDTH * ( PWM_width + 1 ), LOGO_HEIGHT / 2 );
	ps_translate( ofs, 0, fontsize * FONT_UNIT_HEIGHT / 2.5 );
}

////////////////////////////////////////////////////////////////
void write_motif(
	ofstream &ofs,
	double *PWM[],
	int PWM_width,
	double axis_fontsize,
	int reverse )
{
	double x = 0;
	double y = 0;

	_write_x_axis( ofs, PWM_width, reverse, axis_fontsize );
	ps_translate( ofs, 0, axis_fontsize * FONT_UNIT_HEIGHT );
	y += ( axis_fontsize * FONT_UNIT_HEIGHT );
	ps_translate( ofs, 0, MIN_BASE_HEIGHT + LOGO_HEIGHT / 2 );
	y += ( MIN_BASE_HEIGHT + LOGO_HEIGHT / 2 );

	double max_sum_positive = 0;
	double max_sum_negative = 0;
	int i;
	for( i = 0; i < PWM_width; i ++ )
	{
		// calculate the sum of elements in this column of the motif
		double sum_positive = 0;
		double sum_negative = 0;
		int j;
		for( j = 0; j < NUM_N_LETTERS; j ++ )
			if( PWM[j][i] >= 0 )
				sum_positive += PWM[j][i];
			else
				sum_negative += fabs( PWM[j][i] );
				
		if( max_sum_positive < sum_positive )
			max_sum_positive = sum_positive;
		if( max_sum_negative < sum_negative )
			max_sum_negative = sum_negative;
	}
	
	max_sum_positive += 0.0000001;
	max_sum_negative += 0.0000001;
	
	_write_y_axis( ofs, max_sum_positive, max_sum_negative, axis_fontsize,
		PWM_width );

	ps_font( ofs, "Helvetica-Bold", BASE_FONT_SIZE );
	
	int col;
	for( col = 0; col < PWM_width; col ++ )
	{
		i = col;
		if( reverse )
			i = PWM_width - col - 1;
			
		double positive_y = MIN_BASE_HEIGHT;
		double negative_y = MIN_BASE_HEIGHT;
		
		s_base ranked_bases[ NUM_N_LETTERS ];
		int row;
		for( row = 0; row < NUM_N_LETTERS; row ++ )
		{
			ranked_bases[ row ].index = row;
			ranked_bases[ row ].value = fabs( PWM[ row ][ i ] );
		}
	
		// now sort by abs of PWM value
		qsort( ranked_bases, NUM_N_LETTERS, sizeof(s_base), _compare_bases );

		for( row = 0; row < NUM_N_LETTERS; row ++ )
		{
			int j = ranked_bases[ row ].index;
			
			if( PWM[j][i] >= 0 )
			{
				double height = PWM[j][i] / max_sum_positive * LOGO_HEIGHT / 2;

				// set the scale according to calculated height
				double newscale = height / (BASE_FONT_SIZE*FONT_UNIT_HEIGHT) - MIN_BASE_HEIGHT;
				if( newscale <= 0 )
					continue;
					
				// move up
				ps_translate( ofs, 0, positive_y );

				ps_scale( ofs, 1, newscale );
	
				// write the letter
				if( reverse )
					ps_write_char( ofs, __n_letters[ 3-j ],
						__n_letter_colors[ 3-j ].r,
						__n_letter_colors[ 3-j ].g,
						__n_letter_colors[ 3-j ].b );
				else
					ps_write_char( ofs, __n_letters[ j ],
						__n_letter_colors[ j ].r,
						__n_letter_colors[ j ].g,
						__n_letter_colors[ j ].b );
	
				// reset the scale to default
				ps_scale( ofs, 1, 1.0/(newscale) );
	
				// move back
				ps_translate( ofs, 0, -positive_y );

				// update the current position along y axis
				positive_y += height;
			}
			else
			{
				double height = fabs( PWM[j][i] ) / max_sum_negative * LOGO_HEIGHT / 2;

				// set the scale according to calculated height
				double newscale = height / (BASE_FONT_SIZE*FONT_UNIT_HEIGHT) - MIN_BASE_HEIGHT;
				if( newscale <= 0 )
					continue;
					
				// move down
				ps_translate( ofs, 0, -negative_y );

				ps_scale( ofs, 1, -newscale );
	
				// write the letter
				if( reverse )
					ps_write_char( ofs, __n_letters[ 3-j ],
						( __n_letter_colors[ 3-j ].r + 4.0 ) / 5,
						( __n_letter_colors[ 3-j ].g + 4.0 ) / 5,
						( __n_letter_colors[ 3-j ].b + 4.0 ) / 5 );
				else
					ps_write_char( ofs, __n_letters[ j ],
						( __n_letter_colors[ j ].r + 4.0 ) / 5,
						( __n_letter_colors[ j ].g + 4.0 ) / 5,
						( __n_letter_colors[ j ].b + 4.0 ) / 5 );
	
				// reset the scale to default
				ps_scale( ofs, 1, -1.0/(newscale) );
	
				// move back
				ps_translate( ofs, 0, negative_y );

				// update the current position along y axis
				negative_y += height;
			}
		}

		// move right
		ps_translate( ofs, COLUMN_WIDTH, 0 );
		x += COLUMN_WIDTH;
	}

	// restore the position along x axis
	ps_translate( ofs, -x, -y );
}
