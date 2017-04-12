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

//////////////////////////////////////////////////////////
int _n_index( char nucleotide )
// returns the index of the specified nucleotide
{
	int i;
	for( i = 0; i < NUM_N_LETTERS; i ++ )
		if( nucleotide == __n_letters[ i ] )
			return i;
			
	return -1;
}

//////////////////////////////////////////////////////////
double _get_binding_p(
	char *seq,
	double *PFM[],
	int PFM_width,
	int direction,
	int *max_dir )
// This function returns the score of a particular sequence stretch relative to the given PWM
{
	double score_f = 1;
	double score_r = 1;
	
	int i;
	for( i = 0; i < PFM_width; i ++ )
	{
		int pos_f = i;
		int pos_r = PFM_width - i - 1;
		
		double this_score_f = 0;
		double this_score_r = 0;
		
		int n_index = _n_index( seq[ i ] );

		if( n_index < 0 ) // the nucleotide is unknown
		// thus, the scores will be averaged over this position of the PWM
			for( n_index = 0; n_index < NUM_N_LETTERS; n_index ++ )
			{
				this_score_f += PFM[ n_index ][ pos_f ] / double(NUM_N_LETTERS);
				this_score_r += PFM[ n_index ][ pos_r ] / double(NUM_N_LETTERS);
			}
		else // the nucleotide is known
		{
			this_score_f = PFM[ n_index ][ pos_f ];
			this_score_r = PFM[ NUM_N_LETTERS - n_index - 1 ][ pos_r ];
		}
		
		// since we're dealing with PFMs, the probabilities are multiplied
		score_f *= this_score_f;
		score_r *= this_score_r;
	}

	if( direction == FORWARD_ONLY ) // only the forward direction is requested
	{
		*max_dir = FORWARD_ONLY;
		return score_f;
	}
	else if( direction == REVERSE_ONLY ) // only the reverse direction is requested
	{
		*max_dir = REVERSE_ONLY;
		return score_r;
	}

	// any other value for direction means that both directions should be considered
	if( score_r > score_f )
		*max_dir = REVERSE_ONLY;
	else
		*max_dir = FORWARD_ONLY;
	return score_f + score_r;
}

///////////////////////////////////////////////////////////////////////////////////////////
void _scan_sequence(
	s_seq *seq,
	double *PFM[],
	int PFM_width,
	int direction )
// This function scans the given sequence for the instances of the motif
{
	// calculate the sum of the scores for this motif over the specified sequence
	// also, store the position of the window that has the maximum score
	double sum_score = 0;
	double max_score = 0;
	int max_pos = -1;
	int max_dir = -1;
	int i;
	for( i = 0; i <= seq ->seq_length - PFM_width; i ++ )
	{
		// get the score of this position
		int dir = -1;
		double score = _get_binding_p( seq ->seq + i, PFM, PFM_width, direction, &dir );
		sum_score += score; // the final affinity score is the sum of affinity scores across the sequence
		
		// update the max score
		if( i == 0 || // either this is the first score
			max_score < score ) // or this is the best score so far
		{
			max_score = score;
			max_pos = i;
			max_dir = dir;
		}
	}
	
	seq ->score = sum_score;
	seq ->max_pos = max_pos;
	seq ->max_dir = max_dir;
}

////////////////////////////////////////////////////////////////////////////////
int _compare_affinities( const void *a, const void *b )
{
	s_seq *seq1 = *( (s_seq **) a );
	s_seq *seq2 = *( (s_seq **) b );

	if( seq1 ->score > seq2 ->score )
		return -1;
	else if( seq1 ->score < seq2 ->score )
		return 1;
	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
void scan_sequences(
	s_seq *seqs[],
	int num_seqs,
	double *PFM[],
	int PFM_width,
	int direction )
{
	// scan the sequences
	int i;
	for( i = 0; i < num_seqs; i ++ )
		_scan_sequence( seqs[ i ], PFM, PFM_width, direction );
		
	// sort the sequences by descending order of affinities
	qsort( seqs, num_seqs, sizeof(s_seq*), _compare_affinities );
}
