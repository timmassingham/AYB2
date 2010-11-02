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

#include <stdbool.h>
#include <stdint.h>
#include "xio.h"

typedef char NUC;
#define NBASE       4                           ///< Number of nucleotide bases.

#define NUC_AMBIG   4                           ///< Indicates unable to determine base.
#define NUC_A       0                           ///< Position for base A.
#define NUC_C       1                           ///< Position for base C.
#define NUC_G       2                           ///< Position for base G.
#define NUC_T       3                           ///< Position for base T.

typedef char PHREDCHAR;
#define NULL_PHRED  32                          ///< Ascii space character.
#define MIN_PHRED   33                          ///< Start of printable characters.
#define MAX_PHRED   126                         ///< End of printable characters.
#define ERR_PHRED   0

#define MIN_QUALITY 0.0                         ///< Minimum quality indicates unable to compute.

/* standard functions */
void show_NUC(XFILE * fp, const NUC nuc);
void show_PHREDCHAR(XFILE * fp, const PHREDCHAR nuc);
NUC read_NUC(XFILE * fp);
PHREDCHAR read_PHREDCHAR(XFILE * fp);

/** Returns true if NUC set to ambiguous. */
static inline bool isambig(const NUC nuc) { return ((nuc - NUC_AMBIG) == 0); }

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
bool has_ambiguous_base(const NUC * restrict nucs, const uint_fast32_t n);

PHREDCHAR phredchar_from_char( const char c)  __attribute__((const)); 
PHREDCHAR phredchar_from_prob( const real_t p)  __attribute__((const));
real_t quality_from_prob(const real_t p) __attribute__((const));
PHREDCHAR phredchar_from_quality( real_t qual) __attribute__((const));

#endif /* NUC_H_ */
