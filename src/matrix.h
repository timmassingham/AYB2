/**
 * \file matrix.h
 * Public parts of Matrix Class.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
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

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdint.h>
#include <stdbool.h>
#include "xio.h"
#include "utility.h"

/*  Specification did not specify what type of input should be taken. Use a
 * typedef so it can be easily changed, although care must be taken with the
 * input and output routines (printf etc) to make sure types match.
 */

typedef int16_t int_t;                          ///< Define size of integer data type to use.
#define INT_FORMAT " %6d"                       ///< Integer data type output format.
#define INT_MIN INT16_MIN                       ///< Integer data type minimum value.
#define INT_MAX INT16_MAX                       ///< Integer data type maximum value.

/**
 * Matrix structure includes size and whether data stored as real or integer.
 * Use of integer reduces memory requirements.
 */
struct _matrix_str {
    int    nrow, ncol;
    real_t * x;
    int_t  * xint;
    bool   useint;
    };

// Make future abstraction easier
typedef struct _matrix_str * MAT;

// standard functions
MAT new_MAT( const int nrow, const int ncol );
MAT new_MAT_int( const int nrow, const int ncol, const bool useint );
MAT free_MAT( MAT mat );
MAT copy_MAT( const MAT mat);
void show_MAT( XFILE * fp, const MAT mat, const uint32_t mrow, const uint32_t mcol);
void show_MAT_rownum( XFILE * fp, const MAT mat, const uint32_t mrow, const uint32_t mcol, bool rownum);

// standard variations
MAT new_MAT_from_array( const uint32_t nrow, const uint32_t ncol, const real_t * x);
MAT coerce_MAT_from_array(const uint32_t nrow, const uint32_t ncol, real_t * x);
MAT coerce_MAT_from_intarray(const uint32_t nrow, const uint32_t ncol, int_t * x);
MAT identity_MAT( const int nrow);
MAT copyinto_MAT( MAT matout, const MAT matin);
MAT append_columns(MAT matout, const MAT matin, int colstart, int colend);
MAT set_MAT( MAT mat, const real_t x);

// stream i/o
int_t clipint(long int val);
int count_line_columns(const int nrow, char *ptr);
MAT new_MAT_from_line(const int nrow, int *ncol, char *ptr);
void write_MAT_to_line (XFILE * fp, const MAT mat);
MAT read_MAT_from_column_file(XFILE * fp);
void write_MAT_to_column_file(XFILE * fp, const MAT mat, bool freeformat);

// Identities
bool is_square(const MAT mat);

// Operations
MAT vectranspose( const MAT mat, const unsigned int p );
MAT reshape_MAT( MAT mat, const int nrow);
MAT cholesky( MAT mat);
MAT invert_cholesky( MAT mat);
MAT trim_MAT( MAT mat, const int mrow, const int mcol, const bool forwards);
MAT * block_diagonal_MAT( const MAT mat, const int n);
MAT scale_MAT(MAT mat, const real_t f);
MAT transpose_inplace( MAT mat);
MAT transpose( const MAT mat);
MAT invert(const MAT mat);
MAT invert_symmetric(const MAT mat);

real_t xMy( const real_t * x, const MAT M, const real_t * y);
real_t normalise_MAT(MAT mat, const real_t delta_diag);

#endif /* MATRIX_H_ */
