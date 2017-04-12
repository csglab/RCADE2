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

int _calcp = 1;

////////////////////////////////////////////////////////////////////////////////////////
void _calculate_dist_CovXY(
	double *PWM1[],
	double *PWM2[],
	int x,
	double *mean,
	double *var )
{
	const char *shuffling_indices[] = { "0123", "0132", "0213", "0231", "0312", "0321", "1023", "1032", "1203", "1230", "1302", "1320", "2013", "2031", "2103", "2130", "2301", "2310", "3012", "3021", "3102", "3120", "3201", "3210" };
	
	*mean = 0;
	*var = 0;
	int count = 0;
	
	int i, j;
	for( i = 0; i < 24; i ++ )
		for( j = 0; j < 24; j ++ )
		// this is one possible combination from 4!x4! combinations
		{
			// calculate EXY, EX, and EY for position x, given this combination
			double EXY = 0;
			double EX = 0;
			double EY = 0;
		
			int n;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				int n1 = shuffling_indices[ i ][ n ] - '0';
				int n2 = shuffling_indices[ j ][ n ] - '0';
				
				EXY += PWM1[ n1 ][ x ] * PWM2[ n2 ][ x ];
				EX += PWM1[ n1 ][ x ];
				EY += PWM2[ n2 ][ x ];
			}
		
			EXY /= NUM_N_LETTERS;
			EX /= NUM_N_LETTERS;
			EY /= NUM_N_LETTERS;
			double CovXY = EXY - EX * EY;
			
			// for now, calculate sum_CovXY and sum_(CovXy^2)
			*mean += CovXY;
			*var += CovXY * CovXY;
			count ++;
		}
		
	// now, calculate E_CovXY and E_(CovXY^2)
	(*mean) /= count;
	(*var) /= count;
	
	// now, calculate Var_CovXY = E_(CovXY^2) - (E_CovXY)^2
	(*var) -= (*mean) * (*mean);
}			

////////////////////////////////////////////////////////////////////////////////////////
double get_correlation(
	double *PWM1[],
	double *PWM2[],
	int PWM_width,
	double *mean,
	double *stdev )
// returns the correlation of the two PWMs
{
	// the statistics for the two PWMs
	double CovXY = 0;
	double VarX = 0;
	double VarY = 0;
		
	// the statistics for the shuffled versions of the two PWMs
	*mean = 0;
	*stdev = 0;
	double mean_CovXY = 0;
	double var_CovXY = 0;
	
	int x;
	for( x = 0; x < PWM_width; x ++ )
	{
		double EXY = 0;
		double EX = 0;
		double EY = 0;
		double EX2 = 0;
		double EY2 = 0;
		
		int n;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			EXY += PWM1[ n ][ x ] * PWM2[ n ][ x ];
			EX += PWM1[ n ][ x ];
			EY += PWM2[ n ][ x ];
			EX2 += PWM1[ n ][ x ] * PWM1[ n ][ x ];
			EY2 += PWM2[ n ][ x ] * PWM2[ n ][ x ];
		}
		
		EXY /= NUM_N_LETTERS;
		EX /= NUM_N_LETTERS;
		EY /= NUM_N_LETTERS;
		EX2 /= NUM_N_LETTERS;
		EY2 /= NUM_N_LETTERS;
		
		CovXY += ( EXY - EX * EY );
		VarX += ( EX2 - EX * EX );
		VarY += ( EY2 - EY * EY );
		
		
		// calculate the statistics for shuffled versions of these two PWMs at this position
		double this_mean_CovXY = 0;
		double this_var_CovXY = 0;
		if( _calcp )
			_calculate_dist_CovXY( PWM1, PWM2, x, &this_mean_CovXY, &this_var_CovXY );
		mean_CovXY += this_mean_CovXY;
		var_CovXY += this_var_CovXY;
	}
	
	if( VarX <= 0 || VarY <= 0 )
		return 0;
		
	*mean = mean_CovXY / sqrt( VarX * VarY );
	*stdev = sqrt( var_CovXY / ( VarX * VarY ) );
		
	return CovXY / sqrt( VarX * VarY );
}
