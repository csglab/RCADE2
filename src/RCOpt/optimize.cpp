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

extern const int *__mode;

#define MAX_ITER			1000
#define MAX_UNCHANGED		20
#define CONVERGENCE_RATE	0.8

////////////////////////////////////////////////////////////////////////////////
int _compare_motif_scores( const void *a, const void *b )
{
	s_motif *motif1 = *( (s_motif **) a );
	s_motif *motif2 = *( (s_motif **) b );

	if( motif1 ->score > motif2 ->score )
		return -1;
	else if( motif1 ->score < motif2 ->score )
		return 1;
	
	return 0;
}

//////////////////////////////////////////////////////////
int _n_index( char nucleotide );

//////////////////////////////////////////////////////////
void _update_PFM(
	s_seq *seqs[],
	int num_seqs,
	double *PFM[],
	int PFM_width,
	int max_PFM_width )
{
	// store the current PFM
	static double *prev_PFM[ NUM_N_LETTERS ];
	int n, x;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		if( !prev_PFM[ n ] )
			prev_PFM[ n ] = new double[ max_PFM_width ];
		memcpy( prev_PFM[ n ], PFM[ n ], sizeof(double) * PFM_width  );
	}

	// reset the PFM to zero
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		memset( PFM[ n ], 0, sizeof(double) * PFM_width );
		
	// the sum of columns at each position of the PFM
	int counts[ MAX_SEQ_LENGTH ];
	memset( counts, 0, sizeof(int) * PFM_width );
	
	int i;
	for( i = 0; i < num_seqs; i ++ )
		if( seqs[ i ] ->leading_edge )
		{
			int j;
			for( j = 0; j < PFM_width; j ++ )
			{
				int n_index = _n_index( seqs[ i ] ->seq[ seqs[ i ] ->max_pos + j ] );
				if( n_index < 0 ) // this base is unknown
					continue;
					
				if( seqs[ i ] ->max_dir == FORWARD_ONLY )
				{
					PFM[ n_index ][ j ] ++;
					counts[ j ] ++;
				}
				else
				{
					PFM[ NUM_N_LETTERS - n_index - 1 ][ PFM_width - j - 1 ] ++;
					counts[ PFM_width - j - 1 ] ++;
				}
			}
		}
		
	// normalize the PFM
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		for( i = 0; i < PFM_width; i ++ )
			PFM[ n ][ i ] = ( PFM[ n ][ i ] + 0.25 ) / ( counts[ i ] + 1.0 );

	// average the new PFM with the previous PFM
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		for( i = 0; i < PFM_width; i ++ )
			PFM[ n ][ i ] = ( PFM[ n ][ i ] * CONVERGENCE_RATE + prev_PFM[ n ][ i ] * (1.0-CONVERGENCE_RATE) );

}

//////////////////////////////////////////////////////////
bool _calc_AUC(
	s_seq *seqs[],
	int num_seqs,
	int num_pos,
	int num_neg,
	double *auc,
	double *p,
	int *num_leading_edge )
// returns true if the leading-edge set has changed
{
	// find the leading edge set
	// also, calculate the AUC
	double ES = 0;
	double fc_pos = 0;
	double fc_neg = 0;
	double max_ES = 0;
	double sum_rank = 0;
	int leading_edge = 0;
	int j;
	for( j = 0; j < num_seqs; j ++ )
	{
		if( seqs[ j ] ->positive )
		{
			fc_pos ++;
			sum_rank += (j+1);
		}
		else
			fc_neg ++;
			
		ES = fc_pos / num_pos - fc_neg / num_neg;
		
		if( j == 0 || max_ES < ES )
		{
			max_ES = ES;
			leading_edge = j;
		}
	}
	double U = num_pos * num_neg + num_pos * ( num_pos + 1.0 ) / 2.0 - sum_rank;
	*auc = U / ( num_pos * num_neg );
	
	// determine whether the AUC is significant
	double mU = num_pos * num_neg / 2.0;
	double sU = sqrt( num_pos * num_neg * ( num_pos + num_neg + 1.0 ) / 12.0 );
	double zU = ( U - mU ) / sU;
	*p = normal_cum_p( zU );

	// mark the leading edge
	bool changed = false;
	*num_leading_edge = 0;
	for( j = 0; j < num_seqs; j ++ )
		if( seqs[ j ] ->positive )
		{
			if( j <= leading_edge && seqs[ j ] ->leading_edge == 0 ) // this is a new leading edge sequence
			{
				changed = true;
				seqs[ j ] ->leading_edge = 1;
			}
			else if( j > leading_edge && seqs[ j ] ->leading_edge == 1 ) // this is a new non-leading edge sequence
			{
				changed = true;
				seqs[ j ] ->leading_edge = 0;
			}
			
			(*num_leading_edge) += seqs[ j ] ->leading_edge;
		}
	
	return changed;
}

//////////////////////////////////////////////////////////
int _optimize_PFM(
	s_seq *seqs[],
	int num_seqs,
	int num_pos,
	int num_neg,
	double *ini_PWM[],
	double *PFM[],
	int PFM_width,
	int max_PWM_width,
	double *ini_AUC,
	double *ini_p,
	double *opt_AUC,
	double *opt_p )
// returns 1 if optimized, 0 otherwise
{
	*ini_AUC = -1;
	*ini_p = 1;
	*opt_AUC = -1;
	*opt_p = 1;

	// initialize internal variables and buffers
	double best_AUC = -1;
	double best_AUC_p = -1;
	double best_score = -1;
	static double *best_PFM[ NUM_N_LETTERS ];
	static double *PWM[ NUM_N_LETTERS ];
	int n, x;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		if( !best_PFM[ n ] )
			best_PFM[ n ] = new double[ max_PWM_width ];
		memset( best_PFM[ n ], 0, sizeof(double) * max_PWM_width  );
		if( !PWM[ n ] )
			PWM[ n ] = new double[ max_PWM_width ];
		memset( PWM[ n ], 0, sizeof(double) * max_PWM_width  );
	}
		
	
	// initialize the leading edge set to null
	int i;
	for( i = 0; i < num_seqs; i ++ )
		seqs[ i ] ->leading_edge = 0;
		
	int counter = 0;
	for( i = 0; i < MAX_ITER; i ++ )
	{
		// scan and sort the sequences by their affinity
		scan_sequences( seqs, num_seqs, PFM, PFM_width, BOTH_STRANDS );
		
		// calculate AUC and related statistics
		double auc = 0;
		double p_auc = 1;
		int num_leading_edge = 0;
		bool changed = _calc_AUC( seqs, num_seqs, num_pos, num_neg, &auc, &p_auc, &num_leading_edge );

		// calculate the correlation of this motif with the initial motif
		for( n = 0; n < NUM_N_LETTERS; n ++ )
			for( x = 0; x < PFM_width; x ++ )
				PWM[ n ][ x ] = log10( PFM[ n ][ x ] * 4 );			
		double mean_correl = 0, stdev_correl = 1;
		double correl = get_correlation( ini_PWM, PWM, PFM_width, &mean_correl, &stdev_correl );
		double z_correl = ( correl - mean_correl ) / stdev_correl;
		double p_correl = normal_cum_p( z_correl );
		
		double score;
		if( *__mode == 1 )
			score = -log10( p_auc );
		else if( *__mode == 2 )
			score = log10( p_auc ) * log10( p_correl );
		else if( *__mode == 3 )
			score = best_score + 2;
		else
		{
			cout << "ERROR: Invalid optimization mode requested." << endl;
			return 0;
		}

		// store the AUC and its p-value
		if( i == 0 )
		{
			*ini_AUC = auc;
			*ini_p = p_auc;
		}
		
		double prev_best_score = best_score;
		if( i == 0 || best_score < score ) // this is the best PFM so far
		{
			best_score = score;
			best_AUC = auc;
			best_AUC_p = p_auc;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
				memcpy( best_PFM[ n ], PFM[ n ], sizeof(double)*PFM_width );
		}
			
		//cout << "Leading edge size: " << num_leading_edge
		//	<< " ... AUC: " << auc
		//	<< " ... p-AUC: " << p_auc
		//	<< " ... correl: " << correl
		//	<< " ... p-correl: " << p_correl << endl;
		if( p_auc > 1e-4 )
		{
			//cout << "AUC is insignificant. Optimization halted." << endl;
			return 0;
		}
		
		// check whether substantial improvement is achieved
		if( *__mode != 3 && score < prev_best_score + 1 ) // if the objective is not to converge, there must be at least one order of magnitude improvement in the log(p)
			{
				//cout << "Optimization converged." << endl;
				break;
			}
		
		// check whether the algorithm converged
		if( !changed )
		{
			counter ++;
			if( counter >= MAX_UNCHANGED ) // for MAX_UNCHANGED cycles, nothing changed, and so the algorithm converged
			{
				//cout << "Optimization converged." << endl;
				break;
			}
		}
		else
			counter = 0;

		if( i < MAX_ITER - 1 ) // there will be at least one more iteration
		// update the PFM based on the leading edge sequences
			_update_PFM( seqs, num_seqs, PFM, PFM_width, max_PWM_width );
	}
	
	// restore the best PFM
	*opt_AUC = best_AUC;
	*opt_p = best_AUC_p;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		memcpy( PFM[ n ], best_PFM[ n ], sizeof(double)*PFM_width );

	return 1;
}


//////////////////////////////////////////////////////////
void optimize_PFMs(
	s_seq *seqs[],
	int num_seqs,
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
	double *PFM[ NUM_N_LETTERS ];
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		PFM[ n ] = new double[ max_PWM_width ];
		
	// find the number of positive and negative sequences
	int num_pos = 0;
	int num_neg = 0;
	for( i = 0; i < num_seqs; i ++ )
		if( seqs[ i ] ->positive )
			num_pos ++;
		else
			num_neg ++;
		
	// optimize each PWM
	for( i = 0; i < num_motifs; i ++ )
	{
		//cout << "Optimizing " << motifs[ i ] ->name << " ..." << endl;
		// initialize the PFM as the current PFM of this motif
		for( n = 0; n < NUM_N_LETTERS; n ++ )
			memcpy( PFM[ n ] , motifs[ i ] ->PFM[ n ], sizeof(double)*motifs[ i ] ->PWM_width );

		// optimize this PFM
		motifs[ i ] ->PFM_optimized =
			_optimize_PFM( seqs, num_seqs, num_pos, num_neg,
				motifs[ i ] ->PWM, PFM, motifs[ i ] ->PWM_width,
				max_PWM_width,
				&motifs[ i ] ->ini_AUC, &motifs[ i ] ->ini_p,
				&motifs[ i ] ->opt_AUC, &motifs[ i ] ->opt_p );
		
		// store the optimized PFM, convert to PWM, and calculate correlation with the original PWM
		if( motifs[ i ] ->PFM_optimized )
		{
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				motifs[ i ] ->opt_PFM[ n ] = new double[ motifs[ i ] ->PWM_width ];
				memcpy( motifs[ i ] ->opt_PFM[ n ] , PFM[ n ], sizeof(double)*motifs[ i ] ->PWM_width );
				
				int x;
				for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
					PFM[ n ][ x ] = log10( PFM[ n ][ x ] * 4 );
			}
			
			double mean_correl = 0;
			double stdev_correl = 1;
			motifs[ i ] ->correl = get_correlation( motifs[ i ] ->PWM, PFM,
				motifs[ i ] ->PWM_width, &mean_correl, &stdev_correl );
			double z_correl = ( motifs[ i ] ->correl - mean_correl ) / stdev_correl;
			motifs[ i ] ->p_correl = normal_cum_p( z_correl );
			
			//motifs[ i ] ->score =
			//	- log10( motifs[ i ] ->opt_p ) * log10( motifs[ i ] ->ini_p ) * log10( motifs[ i ] ->p_correl );

			if( *__mode == 2 ) // if the objective is to optimize motif similarity in addition to AUC, use both for sorting motifs
				motifs[ i ] ->score =
					log10( motifs[ i ] ->opt_p ) * log10( motifs[ i ] ->p_correl );
			else // otherwise, sort motifs by AUC
				motifs[ i ] ->score =
					-log10( motifs[ i ] ->opt_p );

		}
		
		//cout << endl;
	}
	
	// sort the motifs by descending order of their scores
	qsort( motifs, num_motifs, sizeof(s_motif*), _compare_motif_scores );
}