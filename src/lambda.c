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

/** Cauchy weighting function for IWLS. */
static inline real_t cauchy(const real_t xsqr, const real_t v){
    return 1.0/(1.0+ xsqr/v);
}


/* public functions */

/**
 * Estimate brightness of a cluster using a Ordinary Least Squares approach.
 * Unlike the routine in the original AYB code, this uses processed intensities.
 * \n Used for initial estimate of brightness.
 * - p:          Matrix of processed intensities
 * - base:       Array of base calls, one per cycle
 */
real_t estimate_lambdaOLS( const MAT p, const NUC * base){
    validate(NULL!=p,NAN);
    validate(NBASE==p->nrow,NAN);
    validate(NULL!=base,NAN);
    const uint_fast32_t ncycle = p->ncol;

    real_t numerator = 0.0, denominator = 0.0;
    for (uint_fast32_t cycle=0 ; cycle<ncycle ; cycle++){
        const int cybase = base[cycle];
        if (!isambig(cybase)){
            numerator += p->x[cycle*NBASE+cybase];
            denominator += 1.0;
        }
    }

    real_t lambda = (denominator>0.0)?(numerator / denominator):0.0;
    return (lambda>0.)?lambda:0.;
}

/**
 * Estimate brightness of a cluster using a Weighted Least Squares approach.
 * (Actually a single step of an Iteratively reWeighted Least Squares.)
 * Unlike the routine in the original AYB code, this uses processed intensities.
 * - p            Matrix of processed intensities
 * - base         Array of base calls, one per cycle
 * - oldlambda    Previous estimate of lambda
 * - v            Array of cycle specific scales, length ncycle
 */
real_t estimate_lambdaWLS( const MAT p, const NUC * base, const real_t oldlambda, const real_t * v){
    validate(NULL!=p,NAN);
    validate(NBASE==p->nrow,NAN);
    validate(NULL!=base,NAN);
    validate(oldlambda>=0.,NAN);
    validate(NULL!=v,NAN);
    const uint_fast32_t ncycle = p->ncol;

    real_t numerator = 0.0, denominator=0.0;
    for (uint_fast32_t cycle=0 ; cycle<ncycle ; cycle++){
        const int cybase = base[cycle];
        if(!isambig(cybase)){
            // Calculate Sum Squared Error, then weight
            real_t sse = 0.0;
            for ( int j=0 ; j<NBASE ; j++){
                sse += p->x[cycle*NBASE+j] * p->x[cycle*NBASE+j];
            }
            sse -= 2.0 * oldlambda * p->x[cycle*NBASE+cybase];
            sse += oldlambda*oldlambda;
            real_t w = cauchy(sse,v[cycle]);
            // Accumulate
            numerator += p->x[cycle*NBASE+cybase] * w;
            denominator += w;
        }
    }

    real_t lambda = (denominator>0.)?(numerator / denominator):0.;
    return (lambda>0.)?lambda:0.;
}

/**
 * Estimate lambda by least square.
 * vec(I_i - N) = lambda_i A vec(S_i)
 * \n Solution is y^t A s / s^tA^tAs where y = Vec(I_i - N) and s = Vec(S_i)
 */
real_t estimate_lambda_A ( const MAT intensity, const MAT N, const MAT At, const NUC * base){
    if(NULL==intensity || NULL==N || NULL==At || NULL==base){ return NAN; }
    const uint_fast32_t ncycle = intensity->ncol;

    // Calculate A vec(S_i)
    const int lda = NBASE*ncycle;
    real_t As[lda];
    for ( int i=0 ; i<lda ; i++){
        As[i] = 0.;
        for ( int j=0 ; j<ncycle ; j++){
            const int idx = j*NBASE+base[j];
            As[i] += At->x[i*lda+idx];
        }
    }
    // Numerator and denominator of solution
    real_t sAAs = 0.0;
    real_t sAy = 0.0;
    for ( int i=0 ; i<lda ; i++){
        sAAs += As[i]*As[i];
        sAy += (intensity->xint[i]-N->x[i]) * As[i];
    }

    // Ensure that lambda is sufficiently positive.
    real_t lambda = sAy/sAAs;
    return (lambda>0.)?lambda:0.;
}
