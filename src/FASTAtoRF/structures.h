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

#ifndef _H_STRUCTURES // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_STRUCTURES // this macro indicates that this file is now included

#ifndef _NAMESPACE
using namespace std;
#define _NAMSPACE
#endif


#include <string.h>
#include <stdlib.h>

// ************************************************ macro definitions

// constants
#define MAX_GENES				200000
#define MAX_SEQ_LENGTH			100000
#define MAX_FEATURES			100000
#define MAX_LINE_LENGTH			100000
#define MAX_STRING_LENGTH		1000
#define MAX_ZF_PER_GENE			100
#define MAX_ZFS					100000
#define MAX_BAIT_LENGTH			10

#define _L	29 // the length of each ZF sequence


// macros
#define RELEASE(x)		do{ if(x) delete (x); }while(false) // this do-while loop is called only once
#define _CALL(x)		do{ if(!(x)) return false; }while(false)
#define _COPY_STR(x,y)	do{ x=new char[strlen(y)+1]; strcpy((x),(y)); }while(false)
#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))

// ********************************************** typedef definitions

typedef unsigned char BYTE;

// ********************************************** local variables
static char __aa_letters[] = "-ACDEFGHIKLMNPQRSTVWY";
#define NUM_AA_LETTERS	21


// ******************************************** structure definitions

// this structure will hold the information for the c2h2 zinc fingers
struct s_c2h2
{
	s_c2h2()
	{
		residues = NULL; // no memory allocated yet
		seq_length = 0;
		
		zf_start = -1;
		helix_start = -1;
		last_pos = -1;
	}
	~s_c2h2()
	{
		RELEASE( residues );
	}
	
	char *residues; // the sequence of the zinc finger
	int seq_length; // the length of the zinc finger

	int zf_start; // the start position of the zf in the original sequence
	int helix_start; // the start position of the helix region (including the second C) in the original sequence
	int last_pos; // the last position of the zf
};

// this structure will hold the information for the input genes
struct s_gene
{
	s_gene()
	{
		// initialize all the variables, setting them to zero
		
		name = NULL; // no memory allocated yet
		index = -1;
		seq = NULL; // no memory allocated yet
		seq_length = 0; // the initial length of the sequence is zero
		
		memset( zfs, 0, sizeof(s_c2h2*) * MAX_ZF_PER_GENE ); // no memory allocated yet
		num_zfs = 0;
	}
	~s_gene()
	{
		// release all the allocated memory
	
		RELEASE( name );
		RELEASE( seq );		
	}

	char *name; // the name of this gene, as read from the input FASTA file
	int index;

	char *seq; // the gene, as read from the input FASTA file
	int seq_length; // the length of the nucleotide sequence
	
	s_c2h2 *zfs[ MAX_ZF_PER_GENE ]; // the zinc fingers within this gene
	int num_zfs; // the number of zinc fingers
};

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
