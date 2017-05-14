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

#ifndef _H_DECLARATIONS // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_DECLARATIONS // this macro indicates that this file is now included


#include "structures.h" // this header file is required for definition of structures


bool read_arguments( int argc, char* argv[] );
void print_arguments();
void print_commandline_format();

int find_motif( const char *name, s_motif *motifs[], int num_motifs );
char* extract_phrase( const char *line, int phrase_index, char *phrase, const char *delimiters );

bool open_output( ofstream &ofs, const char *path, const char *extension );
bool open_rndforest( const char *FASTA_file, s_motif *motifs[], int *num_motifs, char const *filter );
bool open_weights( const char *weight_file, s_motif *motifs[], int num_motifs );

bool read_rndforest( ifstream &ifs, s_motif *motifs[], int *num_motifs, char const *filter, char this_EOL );
bool read_weights( ifstream &ifs, s_motif *motifs[], int num_motifs, char this_EOL );

bool generate_PWMs( s_motif *motifs[], int num_motifs, double enrichment_fold );

bool write_PFMs( ofstream &ofs, s_motif *motifs[], int num_motifs );
void write_graphics( ofstream &ofs, s_motif *motifs[], int num_motifs );


//*******************************************************************
#endif // this is to make sure that this file will not be included twice
