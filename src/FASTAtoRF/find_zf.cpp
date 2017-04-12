// Copyright 2014 Hamed S. Najafabadi

/********************************************************************

This file is part of FASTAtoRF.

FASTAtoRF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FASTAtoRF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FASTAtoRF.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

#define VALIDATE_POS(x,l)		((x)<l)

/////////////////////////////////////////////////////////////////////////////////
bool _is_c2h2(
	char *sequence,
	int start,
	int length,
	int *helix_start,
	int *last_pos )
// Uses the pattern #-X-C-X(1-5)-C-X3-#-X5-#-X2-H-X(3-6)-[H/C] to find C2H2s
// See http://pfam.sanger.ac.uk/family/PF00096
//  The pattern in regular expression: ..C.{1,5}C.{12}H.{3,6}[HC], with ..C.{1,5}C.{12}H.{3,6}H having priority
//
// returns true if the specified start position on the sequence corresponds to the
//  beginning of a C2H2, otherwise false
{
	int pos;
	if( !VALIDATE_POS( pos = start + 2, length ) || // the sequence is too short to be a ZF
		sequence[ pos ] != 'C' ) // the amino acid #2 (i.e. the third amino acid) must be C
		return false;
		
	// first, search for ..C.{1,5}C.{12}H.{3,6}H
	int gap1, gap2; // the length of the gaps between the pair of C's, and the pair of H's, respectively
	for( gap1 = 1; gap1 <= 5; gap1 ++ )
	{
		*helix_start = start + 3 + gap1;
		
		if( !VALIDATE_POS( pos = *helix_start, length ) || // this position is not valid
			sequence[ pos ] != 'C' ) // this position does not have the correct amino acid
			continue;
			
		if( !VALIDATE_POS( pos = *helix_start + 13, length ) || // this position is not valid
			sequence[ pos ] != 'H' ) // this position does not have the correct amino acid
			continue;

		for( gap2 = 3; gap2 <= 6; gap2 ++ )
		{
			*last_pos = *helix_start + 14 + gap2;

			if( !VALIDATE_POS( pos = *last_pos, length ) || // this position is not valid
				sequence[ pos ] != 'H' ) // this position does not have the correct amino acid
				continue;
				
			// all conditions are met
			return true;
		}
	}

	// now, search for ..C.{1,5}C.{12}H.{3,6}C
	for( gap1 = 1; gap1 <= 5; gap1 ++ )
	{
		*helix_start = start + 3 + gap1;
		
		if( !VALIDATE_POS( pos = *helix_start, length ) || // this position is not valid
			sequence[ pos ] != 'C' ) // this position does not have the correct amino acid
			continue;
			
		if( !VALIDATE_POS( pos = *helix_start + 13, length ) || // this position is not valid
			sequence[ pos ] != 'H' ) // this position does not have the correct amino acid
			continue;

		for( gap2 = 3; gap2 <= 6; gap2 ++ )
		{
			*last_pos = *helix_start + 14 + gap2;

			if( !VALIDATE_POS( pos = *last_pos, length ) || // this position is not valid
				sequence[ pos ] != 'C' ) // this position does not have the correct amino acid
				continue;
				
			// all conditions are met
			return true;
		}
	}
	
	return false;		
}			

/////////////////////////////////////////////////////////////////////////////////
void _copy_zf_sequence(
	s_c2h2 *zf,
	char *gene_seq,
	int zf_start,
	int helix_start,
	int last_pos )
// Copies the info of the found zf
// Reminder: The C2H2 pattern in regular expression: ..C.{1,5}C.{12}H.{3,6}[HC]
//
// This function converts the sequence into:
//
//  ..C-----C............H------[HC]
//  01234567890123456789012345678
//
//  Therefore, the information about the amino acids between the pair of C's, and between
//  the pair of H's is discarded
{
	zf ->residues = new char[ 30 ]; // allocate memory
	
	// copy the info
	zf ->seq_length = 29;
	zf ->zf_start = zf_start;
	zf ->helix_start = helix_start;
	zf ->last_pos = last_pos;
	
	// copy the first three positions
	memcpy( zf ->residues, gene_seq + zf_start, sizeof(char) * 3 );
	
	// set the next 5 positions to - (gap)
	memset( zf ->residues + 3, '-', sizeof(char) * 5 );
	
	// copy the helix region
	memcpy( zf ->residues + 8, gene_seq + helix_start, sizeof(char) * 14 );
	
	// set the next 6 positions to - (gap)
	memset( zf ->residues + 22, '-', sizeof(char) * 6 );
	
	// copy the last position
	zf ->residues[ 28 ] = gene_seq[ last_pos ];

	zf ->residues[ 29 ] = 0; // terminate the string	
}

/////////////////////////////////////////////////////////////////////////////////
void find_zfs(
	s_gene *genes[],
	int num_genes )
{
	//cout << "GENE" << char(9) << "ZF_INDEX"
	//	<< char(9) << "ZF_START"
	//	<< char(9) << "ZF_END"
	//	<< char(9) << "ZF_SEQ" << endl;

	// examine all genes
	int i;
	for( i = 0; i < num_genes; i ++ )
	{
		genes[ i ] ->num_zfs = 0;
		
		int j;
		for( j = 0; j < genes[ i ] ->seq_length; j ++ )
		{
			int zf_start = j;
			int helix_start;
			int last_pos;
			
			if( _is_c2h2( genes[ i ] ->seq, zf_start, genes[ i ] ->seq_length,
				&helix_start, &last_pos ) )
			// this is a C2H2
			{
				if( genes[ i ] ->num_zfs && // there are previously found ZFs
					zf_start <= genes[ i ] ->zfs[ genes[ i ] ->num_zfs - 1 ] ->last_pos ) // this zf overlaps with the previous one
				{
					cout << "WARNING: Overlapping ZFs are found in " << genes[ i ] ->name
						<< ". Only the first ZF instance is considered." << endl;
				}
				else if( genes[ i ] ->num_zfs < MAX_ZF_PER_GENE ) // there is still room for more ZFs
				// add the found c2h2 to the list
				{
					int index = genes[ i ] ->num_zfs;
					
					genes[ i ] ->zfs[ index ] = new s_c2h2; // create the new ZF
					_copy_zf_sequence(
						genes[ i ] ->zfs[ index ],
						genes[ i ] ->seq,
						zf_start,
						helix_start,
						last_pos ); // copy the sequence, according to the canonical format

					genes[ i ] ->num_zfs ++; // update the number of ZFs
				
					//cout << genes[ i ] ->name << char(9) << index+1
					//	<< char(9) << genes[ i ] ->zfs[ index ] ->zf_start+1
					//	<< char(9) << genes[ i ] ->zfs[ index ] ->last_pos+1
					//	<< char(9) << genes[ i ] ->zfs[ index ] ->residues << endl;
				}
				else // too many ZFs
				{
					cout << "WARNING: Too many ZFs for " << genes[ i ] ->name << ". "
						<< "Some ZFs are ignored." << endl;
					break; // ignore the rest of the sequence
				}
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////
void randomize_zfs(
	s_gene *genes[],
	int num_genes )
{
	s_c2h2 *zfs[ MAX_ZFS ];
	int num_zfs = 0;
	
	// create a list of all zfs
	int i, j;
	for( i = 0; i < num_genes; i ++ )
		for( j = 0; j < genes[ i ] ->num_zfs; j ++ )
		{
			zfs[ num_zfs ] = genes[ i ] ->zfs[ j ];
			num_zfs ++;
		}
	
	// shuffle each position of the ZFs independently	
	for( i = 0; i < _L; i ++ )
		for( j = 0; j < num_zfs - 1; j ++ )
		{
			int rnd_index = ( rand() % (num_zfs-j) ) + j;
	
			char swap = zfs[ j ] ->residues[ i ];
			zfs[ j ] ->residues[ i ] = zfs[ rnd_index ] ->residues[ i ];
			zfs[ rnd_index ] ->residues[ i ] = swap;
		}
}
