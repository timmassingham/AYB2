/** 
 * \file lambda.c
 * Lambda Calculation.
 *//* 
 *  Created : 23 Apr 2010
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

#include <math.h>
#include "lambda.h"

/* constants */
/* None      */


/* members */
/* None    */


/* private functions */
/* None              */


/* public functions */

/**
 * Estimate brightness of a cluster using a Ordinary Least Squares approach.
 * Unlike the routine in the original AYB code, this uses processed intensities.
 * \n Used for initial estimate of brightness.
 *   -  p            Matrix of processed intensities
 *   -  base         Array of base calls, one per cycle
 */
real_t estimate_lambdaOLS( const MAT p, const NUC * base){
    validate(NULL!=p,NAN);
    validate(NBASE==p->nrow,NAN);
    validate(NULL!=base,NAN);
    const uint32_t ncycle = p->ncol;

    real_t numerator = 0.0;
    for (uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        const int cybase = base[cycle];
        numerator += p->x[cycle*NBASE+cybase];
    }

    real_t lambda = numerator / ncycle;
    return (lambda>0.)?lambda:0.;
}
