/** 
 * \file mixnormal.h
 * Public parts of routines for fitting a mixed normal distribution.
 *//* 
 *  Created : 11 May 2011
 *  Author : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
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

#ifndef MIXNORMAL_H_
#define MIXNORMAL_H_

#include "utility.h"
#include <stdbool.h>


/** Data for a mixed normal fitted distribution. */
typedef struct {
        unsigned int nmix;
        real_t * prob;
        real_t * mean;
        real_t * sd;
} * NormMixParam;


/* function prototypes */

// standard functions
NormMixParam new_NormMixParam( const unsigned int nmix);
void free_NormMixParam(NormMixParam param);
void show_NormMixParam( FILE * fp, const NormMixParam param);

NormMixParam fit_mixnormal(const real_t * x, const bool * allowed, const unsigned int n, const unsigned int nmix, const unsigned niter);

#endif /* MIXNORMAL_H_ */

