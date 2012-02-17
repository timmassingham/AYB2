/** 
 * \file call_bases.c
 * Call Bases.
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

#include <assert.h>
#include <math.h>
#include "call_bases.h"


/* constants */
/* None      */

/* members */

static real_t Mu = 1e-5;                        ///< Adjusts range of quality scores (for old call_base).
static real_t PolyQual = 38.0;                  ///< Generalised error value, adjusts level of quality scores.


/* private functions */

/**
 * Subcalculation for call_bases.
 * Returns -2*I_b^t om p + lambda * I_b^t om I^b, where I_b is indicator vector of base b.
 * - p:        Intensities of i^th cycle
 * - lambda:   Brightness of cluster
 * - b:        Putative base
 * - om:       Pointer to first entry of i^th block on the diagonal of Omega; lda of Omega is NBASE*ncycle
 */
static real_t baselike (const real_t * restrict p, const real_t lambda, const NUC b, const real_t * om, const int ncycle){
    if (NULL==p || NULL==om || !isfinite(lambda) || NUC_AMBIG==b) { return NAN;}

    const int lda = NBASE*ncycle;
    real_t like = lambda*om[b*lda+b];
    for ( int i=0 ; i<NBASE ; i++){
        like -= 2.0* p[i]*om[b*lda+i];
    }
    return like;
}

/**
 * Subcalculation for call_bases.
 * Returns lambda I_b^t om I_a - p_{i+1}^t om p_{i} - p_{i}^t om p_{i+1}^t, where I_b/a is indicator vector of base b/a.
 * - p:        Intensities of i^th cycle
 * - lambda:   Brightness of cluster
 * - a:        Putative base for cycle i
 * - b:        Putative base for cycle i+1
 * - om:       Pointer to first entry of i^th block on the lower diagonal of Omega; lda of Omega is NBASE*ncycle
 */
static real_t crosslike(const real_t * restrict p, const real_t lambda, const NUC a, const NUC b, const real_t * om, const int ncycle){
    if (NULL==p || NULL==om || !isfinite(lambda) || NUC_AMBIG==a || NUC_AMBIG==b) { return NAN;}

    const int lda = NBASE*ncycle;
    real_t like = lambda*om[a*lda+b];
    for ( int i=0 ; i<NBASE ; i++){
        like -= p[i]*om[i*lda+b];
        like -= p[NBASE+i]*om[a*lda+i];
    }
    return like;
}

/** Maximum of n real_ts. Should possibly be moved into utility.h library. */
static inline int max_real_t(const real_t * restrict p, const uint32_t n){
    validate(NULL!=p,-1);
    real_t m = p[0];
    int idx = 0;
    for ( uint32_t i=1 ; i<n ; i++){
        if(p[i]>m){ m=p[i]; idx=i;}
    }
    return idx;
}


/* public functions */

/**
 * Call base from processed intensities using maximum intensity.
 * Used for initial base call.
 */
NUC call_base_simple( const real_t * restrict p){
    return max_real_t(p,NBASE);
}

/** Return a nodata base call, used when missing data. */
NUC call_base_nodata(void){
    return NUC_AMBIG;
}

/** Return a null base call, used when insufficient data available. */
struct basequal call_base_null(void){
    struct basequal b = {0, MIN_QUALITY};
    return b;
}

/**
 * Call base from processed intensities using minimum Least Squares.
 * Also returns a quality score.
 * - p:       Processed intensities for given cycle
 * - lambda:  Cluster brightness
 * - omega:   Cycle specific inverse covariance matrix
 */
struct basequal call_base( const real_t * restrict p, const real_t lambda, const real_t * restrict penalty, const MAT omega){
    assert(NULL!=p);
    assert(NULL!=omega);
    assert(NBASE==omega->nrow);
    assert(NBASE==omega->ncol);

    if(0==lambda){
        return call_base_null();
    }

    int call = 0;
    real_t minstat = HUGE_VAL;
    real_t stat[NBASE] = { 0.,0.,0.,0.};
    for ( int i=0 ; i<NBASE ; i++){
        stat[i] = lambda * omega->x[i*NBASE+i];
        for ( int j=0 ; j<NBASE ; j++){
            stat[i] -= 2.0 * p[j] * omega->x[i*NBASE+j];
        }
        stat[i] *= lambda;
        stat[i] += penalty[i];
        if(stat[i]<minstat){
            minstat = stat[i];
            call = i;
        }
    }

    /* Summation of probabilities for normalisation,
     * having removed factor exp(-0.5*(K+minstat))
     */
    real_t tot = 0.;
    for ( int i=0 ; i<NBASE ; i++){
        tot += exp(-0.5*(stat[i]-minstat));
    }

    real_t K = xMy(p,omega,p);
    real_t maxprob = exp(-0.5*(K+minstat));

    /* Calculate posterior probability in numerically stable fashion
     * Note that maxp can be extremely small.
     */
    const real_t exp_pen = exp(-0.5*penalty[call]);
    real_t post_prob = (maxprob<Mu) ?
        // Case probabilities small compared to mu
        (exp_pen*Mu + maxprob ) / (4.0*Mu + maxprob*tot) :
        // Case probabilities large compared to mu
        (exp_pen*Mu/maxprob + 1.0) / (4.0*Mu/maxprob + tot);

    struct basequal b = {call, quality_from_prob(post_prob)};
    return b;
}

/**
 * Call bases from processed intensities.
 * Uses a dynamic programming algorithm (Viterbi) to find the bases that minimise the
 * squared error:
\verbatim
        (p-lambda*s)^t Omega (p-lambda*s).
      = p^t Om p + lambda ( 2 * s^t Omega p + lambda^2 s^t Omega s )
                          |------------------ A -------------------|
\endverbatim
 * Sufficient to minimise A.
 * Omega is block-tridiagonal, that is composed from 4x4 blocks, e.g.
\verbatim
        Om_{11} Om_{21}^t 0
        Om_{21} Om_{22}   Om_{32}^t
        0       Om_{32}   Om_{33}
\endverbatim
 * Return value is p^t Om p + lambda * min A.
 */  
real_t call_bases( const MAT p, const real_t lambda, const MAT omega, NUC * base){
    if (NULL==base || NULL==p || NULL==omega) { return NAN; }

    const int ncycle = p->ncol;
    const int lda = omega->ncol;

    // array contains accumulation information and trace-back directions
    real_t array[NBASE*ncycle];
    // Initialise. First elements of array contain contribution from first cycle.
    for ( int b=0 ; b<NBASE ; b++){
        array[b] = baselike(p->x,lambda,b,omega->x,ncycle);
    }

    // Forwards piece of algorithm. For each cycle, for each base, find the base in the
    // previous cycle that minimises
    for ( int cy=1; cy<ncycle ; cy++){
        NUC precall[NBASE];
        for ( int b=0 ; b<NBASE ; b++){
            // Contribution from calling b at cycle cy, independent from other cycles
            real_t diag = baselike(p->x+cy*NBASE,lambda,b,omega->x+cy*NBASE*lda+cy*NBASE,ncycle);
            // Find call prev at previous cycle that minimises cost of calling b and prev
            real_t minstat = HUGE_VAL;
            int minidx = 0;
            real_t stat[NBASE];
            for ( int prev=0 ; prev<NBASE ; prev++){
                stat[prev] = diag +                     // Cost of calling b at cycle cy
                        array[(cy-1)*NBASE+prev] +      // Cost of calling prev at cycle (cy-1)
                        // Adjustment for calling both.
                        2.0*crosslike(p->x+(cy-1)*NBASE,lambda,prev,b,omega->x+(cy-1)*NBASE*lda+cy*NBASE,ncycle);
                // Keep track of best previous call
                if(stat[prev]<minstat){
                    minidx = prev;
                    minstat = stat[prev];
                }
            }
            // Save best call for previous cycle given call b at this cycle.
            array[cy*NBASE+b] = minstat;
            precall[b] = minidx;
        }
        // No longer need previous entries of array. Use memory to save trace-back
        // information.
        for ( int b=0 ; b<NBASE ; b++){
            array[(cy-1)*NBASE+b] = precall[b];
        }
    }

    // Backwards piece of algorithm.
    // Find best call on last cycle
    real_t minstat = HUGE_VAL;
    int minidx = 0;
    for ( int b=0 ; b<NBASE ; b++){
        if(array[(ncycle-1)*NBASE+b]<minstat){ minidx = b; minstat = array[(ncycle-1)*NBASE+b];}
    }
    // Trace calls backward using previously stored information
    base[ncycle-1] = minidx;
    for ( int cy=(ncycle-2) ; cy>=0 ; cy--){
        base[cy] = (NUC)array[cy*NBASE+base[cy+1]];
    }

    return xMy(p->x,omega,p->x) + lambda * minstat;
}

/** 
 * Calculate the call_bases return value for a supplied sequence. 
 * See call_bases for algorithm. 
 */
real_t calculate_lss(const MAT p, const real_t lambda, const MAT omega, const NUC * base) {
    if (NULL==base || NULL==p || NULL==omega) { return NAN; }

    const int ncycle = p->ncol;
    const int lda = omega->ncol;

	real_t res = 0.0;
    for ( int cy=0; cy<ncycle ; cy++){
		res += baselike(p->x+cy*NBASE,lambda,base[cy],omega->x+cy*NBASE*lda+cy*NBASE,ncycle);
	}
    for ( int cy=1; cy<ncycle ; cy++){
		res += 2.0*crosslike(p->x+(cy-1)*NBASE,lambda,base[cy-1],base[cy],omega->x+(cy-1)*NBASE*lda+cy*NBASE,ncycle);
	}
	
	return xMy(p->x,omega,p->x) + lambda * res;
}

/** 
 * Posterior probabilities via fwds/bwds.
 * Repeats much of call_bases.
 */
void call_qualities_post(const MAT p, const real_t lambda, const MAT omega, const real_t effDF, NUC * base, real_t * qual){
	if(NULL==base || NULL==p || NULL==omega){ return; }

    const int ncycle = p->ncol;
    const int lda = omega->ncol;
    const real_t polyErr = -expm1(-PolyQual/10.0 * log(10.0));

    // arrays contain accumulation information
    real_t farray[4*ncycle], barray[4*ncycle];
    // Calculate costs
    real_t basecost[4*ncycle], crosscost[16*ncycle];
    for ( int b=0 ; b<NBASE ; b++){
        basecost[b] = baselike(p->x,lambda,b,omega->x,ncycle);
    }
    for ( int cy=1 ; cy<ncycle ; cy++){
        for ( int b=0 ; b<NBASE ; b++){
            basecost[4*cy+b] = baselike(p->x+cy*NBASE,lambda,b,omega->x+cy*NBASE*lda+cy*NBASE,ncycle);
            for ( int b2=0 ; b2<NBASE ; b2++){
                crosscost[16*cy+4*b+b2] = crosslike(p->x+(cy-1)*NBASE,lambda,b,b2,omega->x+(cy-1)*NBASE*lda+cy*NBASE,ncycle);
            }
        }
    }

    // Forwards piece of algorithm.
    // Initalise
    for ( int b=0 ; b<4 ; b++){ farray[b] = 0.0; }
    // Iteration
    for ( int cy=1; cy<ncycle ; cy++){
        for ( int b=0 ; b<NBASE ; b++){
            farray[4*cy+b] = HUGE_VAL;
            // Find minimum for past calls give current call
            for ( int prev=0 ; prev<4 ; prev++){
                real_t callLS = basecost[4*(cy-1)+prev] + farray[4*(cy-1)+prev] + 2.0*crosscost[16*cy+4*prev+b];
                if(callLS<farray[4*cy+b]){ farray[4*cy+b] = callLS; }
            }
        }
    }
    real_t fwdslike = HUGE_VAL;
    for ( int b=0 ; b<4 ; b++){
        if(basecost[(ncycle-1)*4+b]+farray[(ncycle-1)*4+b]<fwdslike){
            fwdslike = basecost[(ncycle-1)*4+b]+farray[(ncycle-1)*4+b];
        }
    }

    // Backwards piece of algorithm.
    // Initialise
    for ( int b=0 ; b<4 ; b++){ barray[(ncycle-1)*4+b] = 0.0; }
    // Iteration
    for ( int cy=(ncycle-2) ; cy>=0 ; cy--){
        for ( int b=0 ; b<NBASE ; b++){
            barray[4*cy+b] = HUGE_VAL;
            // Find minimum for future calls given current call
            for ( int nxt=0 ; nxt<4 ; nxt++){
                real_t callLS = basecost[4*(cy+1)+nxt] + barray[4*(cy+1)+nxt] + 2.0*crosscost[16*(cy+1)+4*b+nxt];
                if(callLS<barray[4*cy+b]){ barray[4*cy+b] = callLS; }
            }
        }
    }
    real_t bkdslike = HUGE_VAL;
    for ( int b=0 ; b<4 ; b++){
        if(basecost[b]+barray[b]<bkdslike){
            bkdslike = basecost[b]+barray[b];
        }
    }
	//fprintf(stdout,"%e %e\n",fwdslike,bkdslike);

    real_t pOp = xMy(p->x,omega,p->x);
    //fprintf(stdout,"%e\n",pOp+lambda*(farray[base[0]] + barray[base[0]] + basecost[base[0]]));
    for ( int cy=0 ; cy<ncycle ; cy++){
        int b = base[cy];
        real_t xmax = farray[4*cy+b] + barray[4*cy+b] + basecost[4*cy+b];
        xmax = pOp + lambda * xmax;
        real_t sum = 0.0;
        for ( int b=0 ; b<4 ; b++){
            real_t stat = farray[4*cy+b] + barray[4*cy+b] + basecost[4*cy+b];
            stat *= lambda;
            stat += pOp;
            sum += exp(-0.5*(1.0+effDF)*(log1p(stat)-log1p(xmax)));
        }
        real_t prob = 1.0 / sum;
        prob *= polyErr;
        qual[cy] = quality_from_prob(prob);
    }
}

/** Return value of generalised error. */
real_t get_generr(void) {
    return PolyQual;
}

/** Set value for generalised error. */
bool set_generr(const char *generr_str) {

    char *endptr;
    PolyQual = strtor(generr_str, &endptr);

    return (PolyQual > 0);
}

/** Return value of Mu. */
real_t get_mu(void) {
    return Mu;
}

/** Set value for Mu. */
bool set_mu(const char *mu_str) {

    char *endptr;
    Mu = strtor(mu_str, &endptr);

    return (Mu > 0);
}

