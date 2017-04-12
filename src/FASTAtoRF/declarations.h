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

#ifndef _H_DECLARATIONS // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_DECLARATIONS // this macro indicates that this file is now included


#include "structures.h" // this header file is required for definition of structures


bool read_arguments( int argc, char* argv[] );
void print_arguments();
void print_commandline_format();

void find_zfs( s_gene *genes[], int num_genes );
void randomize_zfs( s_gene *genes[], int num_genes );

int find_gene( const char *name, s_gene *genes[], int num_genes );
char* extract_phrase( const char *line, int phrase_index, char *phrase, const char *delimiters );

bool open_output( ofstream &ofs, const char *path, const char *extension );
bool open_FASTA( const char *FASTA_file, s_gene *genes[], int *num_genes );

bool read_FASTA( ifstream &ifs, s_gene *genes[], int *num_genes, char this_EOL );

void write_ZFs( ofstream &ofs, s_gene *genes[], int num_genes, int zf_span, int min_linker_length, int max_linker_length );

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
