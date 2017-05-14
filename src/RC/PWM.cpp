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

int _l = 4; // the length of the recognition sequence for each ZF

//////////////////////////////////////////////////////////
double _get_max_range(
	double *PWM[],
	int PWM_width )
// determines the maximum difference between the max and min values of any column in the PWM
{
	double max_range = 0;	
	// determine the min and max scores in each position, based on which the overall min/max will be deteremined
	int n, x;
	for( x = 0; x < PWM_width; x ++ )
	{
		// get the min and max in position x
		double this_min, this_max;	
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			if( n == 0 ||
				this_min > PWM[ n ][ x ] )
				this_min = PWM[ n ][ x ];
			if( n == 0 ||
				this_max < PWM[ n ][ x ] )
				this_max = PWM[ n ][ x ];
		}
		
		double this_range = this_max - this_min;
		if( x == 0 ||
			max_range < this_range )
			max_range = this_range;
	}
	
	return max_range;
}

//////////////////////////////////////////////////////////
void _convert_PWM_to_likelihood(
	double *PWM[],
	int PWM_width,
	double enrichment_fold )
{	
	double c = _get_max_range( PWM, PWM_width ) / log( enrichment_fold );
	c += 0.00000001; // this pseudocount ensures that there will be no DIV/0
	
	//////////////
	//c = 0.746986509;
	//////////////
	
	// convert to likelihoods
	int x, n;
	for( x = 0; x < PWM_width; x ++ )
	{
		// convert this column values
		
		// pass 1: convert to frequencies
		// initialize the frequencies
		double sum_exp = 0;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			PWM[ n ][ x ] = exp( PWM[ n ][ x ] / c );
			sum_exp += PWM[ n ][ x ];
		}
		// pass 2: convert to log-likelihoods
		// normalize the frequencies, and then covert to log-likelihood
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			PWM[ n ][ x ] /= sum_exp; // normalize the frequencies
			PWM[ n ][ x ] = log10( PWM[ n ][ x ] * NUM_N_LETTERS ); // convert to log-likelihood
		}
	}
}


//////////////////////////////////////////////////////////
void _get_PWM(
	s_motif *motif,
	double *PWM[],
	int *PWM_width,
	double enrichment_fold )
{
	// determine the PWM width based on the number of zinc fingers and the bait length in the model
	*PWM_width = motif ->num_zfs * 3 + _l - 3;
	// set the PWM to zero
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		// reset all the numbers of this row to zero
		memset( PWM[ n ], 0, sizeof(double) * (*PWM_width) );
	}

	int offset = *PWM_width - _l; // the offset of the zinc finger PWM relative to the overall PWM
	
	// concatenate the PWMs of the ZFs
	int i, x;
	for( i = 0; i < motif ->num_zfs; i ++ )
	{
		for( n = 0; n < NUM_N_LETTERS; n ++ )	// examine each nucleotide
			for( x = 0; x < _l; x ++ ) // examine each bait position
			//for( x = 0; x < 3; x ++ ) // examine only the first three positions
				PWM[ n ][ x + offset ] +=
					motif ->zfs[ i ] ->weight * motif ->zfs[ i ] ->PWM[ n ][ x ];
				
		offset -= 3; // update the offset
	}
	
	// the PWM should be converted to likelihood-based PWM
	if( enrichment_fold > 0 )
		_convert_PWM_to_likelihood( PWM, *PWM_width, enrichment_fold );
}


///////////////////////////////////////////////////////////////////////////////////////////
bool generate_PWMs(
	s_motif *motifs[],
	int num_motifs,
	double enrichment_fold )
{
	// find the maximum number of zinc fingers per motif
	int max_num_zfs = 0;
	int i;
	for( i = 0; i < num_motifs; i ++ )
		if( max_num_zfs < motifs[ i ] ->num_zfs )
			max_num_zfs = motifs[ i ] ->num_zfs;
			
	// initialize the PWM
	double *PWM[ NUM_N_LETTERS ];
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		PWM[ n ] = new double[ max_num_zfs * 3 + _l - 3 ];
		
	// store the PWM for each motif
	for( i = 0; i < num_motifs; i ++ )
	{
		// get the PWM for this motif
		int PWM_width;
		_get_PWM( motifs[ i ], PWM, &PWM_width, enrichment_fold );
		
		// store the PWM and the PFM for this motif
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			motifs[ i ] ->PWM[ n ] = new double[ PWM_width ];
			motifs[ i ] ->PFM[ n ] = new double[ PWM_width ];

			int j;
			for( j = 0; j < PWM_width; j ++ )
			{
				motifs[ i ] ->PWM[ n ][ j ] = PWM[ n ][ j ];
				motifs[ i ] ->PFM[ n ][ j ] = pow( 10, PWM[ n ][ j ] ) / 4.0;
			}
		}

		motifs[ i ] ->PWM_width = PWM_width - 1;
	}
	
	return true;
}
