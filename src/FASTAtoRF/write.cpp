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

//////////////////////////////////////////////////////////
bool _write_contig(
	ofstream &ofs,
	s_gene *gene,
	int contig_index,
	int start_index_within_contig,
	int zf_span,
	int min_linker_length,
	int max_linker_length,
	int *contig_found )
// returns true if the contig_index was found, false otherwise
{
	// determine the first and last zinc fingers that are in the requested contig
	int first_zf_index = -1;
	int last_zf_index = -1;
	int i;
	int this_contig_index = 0;
	for( i = 0; i < gene ->num_zfs; i ++ )
	{
		if( this_contig_index == contig_index ) // this is the requested contig
		{
			if(	first_zf_index < 0 ) // this is the first zf of this contig
				first_zf_index = i;
				
			last_zf_index = i; // so far, this is the last zf of this contig
		}
		
		if( i < gene ->num_zfs - 1 ) // this is not the last zf
		{
			// the linker length is the distance between this zf and the next zf
			int linker_length =
				gene ->zfs[ i + 1 ] ->zf_start - gene ->zfs[ i ] ->last_pos - 1;
			
			if( linker_length < min_linker_length ||
				linker_length > max_linker_length ) // the linker length is out of the scope, and therefore the next zf belongs to a different contig
				this_contig_index ++;
		}
	}
	
	if( first_zf_index < 0 ) // if the fist_zf_index is not set, it means that the requested contig was not found
	{
		*contig_found = 0; // this tells the caller that the contig was not found
		return false;
	}
	*contig_found = 1;
	
	if( start_index_within_contig < 0 ) // the start index cannot be negative
		return false;
		
	if( zf_span <= 0 &&
		start_index_within_contig > 0 ) // if no particular span is requested, the start index can only be zero
		return false;
		
	// adjust the first and last ZFs based on the requested span
	first_zf_index += start_index_within_contig; // start from the requested zinc finger within the contig
	if( zf_span > 0 ) // only a few ZFs from this contig should be included
	{
		if( first_zf_index + zf_span - 1 > last_zf_index )
			return false;
			
		last_zf_index = first_zf_index + zf_span - 1;
	}
			

	// write the info for each zf
	for( i = 0; i < gene ->num_zfs; i ++ )
	{
		// is this zf in the requested contig?
		if( i >= first_zf_index && i <= last_zf_index )
		{

			ofs << gene ->name << "|" << contig_index+1 << "|" << start_index_within_contig+1 << "x" << zf_span
				<< ";" << i-first_zf_index+1 << ";" << i+1;

				int j;
				for( j = 0; j < 12; j ++ )
					if( strchr( __aa_letters, gene ->zfs[ i ] ->residues[ j + 9 ] ) )
						ofs << char(9) << gene ->zfs[ i ] ->residues[ j + 9 ];
					else
						ofs << char(9) << "-";
				for( j = 0; j < 18; j ++ )
					ofs << char(9) << 0;
					
			ofs << endl;
		}
	}
		
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
void write_ZFs(
	ofstream &ofs,
	s_gene *genes[],
	int num_genes,
	int zf_span,
	int min_linker_length,
	int max_linker_length )
{
	ofs << "Name"
		<< char(9) << "p1" << char(9) << "p2" << char(9) << "p3" << char(9) << "p4" << char(9) << "p5" << char(9) << "p6" << char(9) << "p7" << char(9) << "p8" << char(9) << "p9" << char(9) << "p10" << char(9) << "p11" << char(9) << "p12"
		<< char(9) << "A1" << char(9) << "C1" << char(9) << "G1" << char(9) << "T1" << char(9) << "A2" << char(9) << "C2" << char(9) << "G2" << char(9) << "T2" << char(9) << "A3" << char(9) << "C3" << char(9) << "G3" << char(9) << "T3" << char(9) << "A4" << char(9) << "C4" << char(9) << "G4" << char(9) << "T4"
		<< char(9) << "Bind" << char(9) << "Accuracy" << endl;

	// write the ZFs for each gene
	int i, j, k;
	for( i = 0; i < num_genes; i ++ )
		if( genes[ i ] ->num_zfs ) // this gene has at least one C2H2 zinc finger
		{
			int contig_found = 1;
			
			for( j = 0; contig_found == 1; j ++ ) // try all possible contigs
				for( k = 0; ; k ++ ) // try all possible start indices
					if( !_write_contig( ofs,
						genes[ i ], j, k, zf_span,
						min_linker_length, max_linker_length,
						&contig_found ) )
						break; // this contig was not found
		}
	
	// write a set of dummy ZFs so that all the predictors have matching factors as the randomForest
	for( i = 0; i < NUM_AA_LETTERS; i ++ )
	{
		ofs << "letter" << __aa_letters[ i ];
		for( j = 0; j < 12; j ++ )
			ofs << char(9) << __aa_letters[ i ];
		for( j = 0; j < 18; j ++ )
			ofs << char(9) << 0;
		ofs << endl;
	}			
}
