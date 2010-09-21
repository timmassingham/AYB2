/** 
 * \file call_bases.h
 * Public parts of Call Bases.
 *//* 
 *  Created : 20 Apr 2010
 *  Author  : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
 *
 *  This file is part of the AYB base calling software.
 *
 *  AYB is free software; you can redistribute it and/or modify
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

#ifndef CALL_BASES_H_
#define CALL_BASES_H_

#include "matrix.h"
#include "nuc.h"
#include "utility.h"

/** Pair type allows return of base and quality. */
struct basequal { NUC base; PHREDCHAR qual;};

/* function prototypes */
NUC call_base_simple( const real_t * restrict p);
struct basequal call_base_null(void);
struct basequal call_base( const real_t * restrict p, const real_t lambda, const MAT omega);
bool set_mu(const char *mu_str);

#endif /* CALL_BASES_H_ */
