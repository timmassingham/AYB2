/** 
 * \file nuc.h
 * Public parts of Nucleotide Base and Phred Quality Score Class.
 *//* 
 *  Created : 2010
 *  Author : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the AYB base-calling software.
 *
 *  AYB is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AYB is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AYB.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef NUC_H_
#define NUC_H_

#include <stdint.h>
#include "xio.h"

typedef char NUC;
#define NBASE       4

#define NUC_AMBIG   4
#define NUC_A       0
#define NUC_C       1
#define NUC_G       2
#define NUC_T       3

typedef char PHREDCHAR;
#define NULL_PHRED  32                          // ascii space
#define MIN_PHRED   33                          // start of printable chars
#define MAX_PHRED   126                         // end of printable chars
#define ERR_PHRED   0

/* standard functions */
void show_NUC(XFILE * fp, const NUC nuc);
void show_PHREDCHAR(XFILE * fp, const PHREDCHAR nuc);
NUC read_NUC(XFILE * fp);
PHREDCHAR read_PHREDCHAR(XFILE * fp);

/* use type safe ARRAY construct to define arrays of NUC and PHREDCHAR */
#define X(A) A ## NUC
#include "array.def"
#undef X
#define X(A) A ## PHREDCHAR
#include "array.def"
#undef X


/* function prototypes */
NUC nuc_from_char( const char c) __attribute__((const));  // Strictly not constant but is constant or abort
char char_from_nuc(const NUC nuc) __attribute__((const));
ARRAY(NUC) nucs_from_string( const char * nucstr );
NUC complement(const NUC nuc) __attribute__((const));
ARRAY(NUC) reverse_complement(const ARRAY(NUC) nucs);
PHREDCHAR phredchar_from_char( const char c)  __attribute__((const)); 
PHREDCHAR phredchar_from_prob( const real_t p)  __attribute__((const));
real_t quality_from_prob(const real_t p) __attribute__((const)) __attribute__((const));
PHREDCHAR phredchar_from_quality( real_t qual) __attribute__((const)) __attribute__((const));

#endif /* NUC_H_ */
