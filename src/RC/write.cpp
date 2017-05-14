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

#define _PTRCHECK(x)	if( !x ){ cout << "ERROR: Unexpected end of string." << endl; return false; }
static char _orsign[] = { '|', 0 };

///////////////////////////////////////////////////////////////////////////////////////////
bool write_PFMs(
	ofstream &ofs,
	s_motif *motifs[],
	int num_motifs )
{
	// write the PWM for each motif
	int i;
	for( i = 0; i < num_motifs; i ++ )
	{
		ofs << "TF" << char(9) << "Unknown" << endl
			<< "TF Name" << char(9) << "Unknown" << endl
			<< "Gene" << char(9) << motifs[ i ] ->protein_name << endl
			<< "Motif" << char(9) << motifs[ i ] ->name << endl
			<< "Family" << char(9) << "C2H2 ZF" << endl
			<< "Species" << char(9) << "Unknown" << endl;

		ofs << "Pos";
		int n;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
			ofs << char(9) << __n_letters[ n ];
		ofs << endl;

		int j;
		for( j = 0; j < motifs[ i ] ->PWM_width; j ++ )
		{
			ofs << j+1;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
				ofs << char(9) << motifs[ i ] ->PFM[ n ][ j ];
			ofs << endl;
		}

		ofs << endl << endl;
	}
	
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
void write_graphics(
	ofstream &ofs,
	s_motif *motifs[],
	int num_motifs )
{
	// find the maximum number of zinc fingers per motif
	int max_PWM_width = -1;
	int i;
	for( i = 0; i < num_motifs; i ++ )
		if( max_PWM_width < motifs[ i ] ->PWM_width )
			max_PWM_width = motifs[ i ] ->PWM_width;
			
	// initialize the PWM
	double *PWM[ NUM_N_LETTERS ];
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		PWM[ n ] = new double[ max_PWM_width ];

	// draw the logos
	char string[ MAX_STRING_LENGTH + 1 ];
	for( i = 0; i < num_motifs; i ++ )
	{
		// copy the information for the relevant PWM
		int motif_index = i;
		int PWM_width = motifs[ motif_index ] ->PWM_width;
		int n, x;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
			for( x = 0; x < PWM_width; x ++ )
				PWM[ n ][ x ] = motifs[ motif_index ] ->PWM[ n ][ x ];

		// if this is the beginning of the page, add margin					
		if( i%8 == 0 )
			//ps_translate( ofs, LEFT_MARGIN, PAGE_HEIGHT - TOP_MARGIN - LOGO_HEIGHT );
			ps_translate( ofs, PAGE_WIDTH / 5, PAGE_HEIGHT - TOP_MARGIN - LOGO_HEIGHT );

		// write the motif name
		ps_translate( ofs, 0, LOGO_HEIGHT / 2 + 10 );
		ps_font( ofs, "Helvetica-Bold", BASE_FONT_SIZE * 0.75 );
		ps_write_string( ofs, motifs[ motif_index ] ->name, 0, 0, 0 );
		ps_translate( ofs, 0, -LOGO_HEIGHT / 2 - 10 );

		// write the motif
		ps_translate( ofs, 0, -LOGO_HEIGHT * 0.6 );
		write_motif( ofs, PWM, PWM_width, BASE_FONT_SIZE / 2, 0 );

		// draw a box
		ps_box_path( ofs,
			-COLUMN_WIDTH / 2, -LOGO_HEIGHT * 0.05,
			PWM_width * COLUMN_WIDTH + COLUMN_WIDTH,
			LOGO_HEIGHT * 1.2,
			0.7, 0.7, 0.7 );
				
		// make room for the next logo
		ps_translate( ofs, 0, -BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5 );

		// if this is the optimized PWM, make room for the title of the next motif
		ps_translate( ofs, 0, -15 );
	
		if( i%8 == 7 )
			ofs << "showpage" << endl;
	}

	if( i%8 != 0 )
		ofs << "showpage" << endl;
}
