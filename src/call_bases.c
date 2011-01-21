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
#include "dirio.h"
#include "message.h"
/* Include set calibration table to be used when none given */
#include "../tables/newcalibrationS2.tab"

/* constants */

/** Quality calibration conversion items. */
typedef enum QCCIdT {E_INTERCEPT, E_SLOPE, E_PRIOR_BASE, E_BASE_NEXT, E_PRIOR_BASE_NEXT} QCCID;
/** Name text for quality calibration conversion table messages and output. */
static const char *MESS_TEXT[] = {"Intercept", "Slope", "Prior-Base", "Base-Next", "Prior-Base-Next"};

/* members */

static real_t Mu = 1e-5;                        ///< Adjusts range of quality scores.


/* private functions */

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

/**
 * Read a single value from a column matrix file.
 * Returns false if fails to read.
 */
static bool get_matrix_single(XFILE * fp, real_t *value) {

    MAT mat = read_MAT_from_column_file(fp);
    if (mat == NULL) {return false;}

    /* matrix read returns at least one element */
    *value = mat->x[0];
    free_MAT(mat);
    return true;
}

/**
 * Read an array of values from a column matrix file, and apply as conversion.
 * Returns false if fails to read or matrix wrong size.
 */
static bool get_and_convert_array(XFILE * fp, const int nrow, const int ncol, const real_t scale, real_t *values) {

    MAT mat = read_MAT_from_column_file(fp);
    if (mat == NULL) {return false;}
    if ((mat->nrow != nrow) || (mat->ncol != ncol)) {
        free_MAT(mat);
        return false;
    }

    /* convert values and store back in same place */
    const int nelt = nrow *ncol;
    for (int idx = 0; idx < nelt; idx++) {
        values[idx] = values[idx] * scale + mat->x[idx];
    }

    free_MAT(mat);
    return true;
}

/** Output a set of array values as a list of columns. */
static void show_array (XFILE * fp, real_t *values, const int nrow, const int ncol, const char *name) {

    xfprintf(fp, "# %s\n", name);
    MAT mat = coerce_MAT_from_array(nrow, ncol, values);
    write_MAT_to_column_file (fp, mat, true);

    /* only free top structure as values refer to array */
    xfree(mat);
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
 * Adjust quality score for base using a linear calibration and neighbours.
 * Vectors used are defined in calibration table (e.g. tables/newcalibrationS2.tab)
 * or a path given from command line.
 * First and last bases of a read are special cases, dealt with by setting the
 * prior or next base (respectively) to be NUC_AMBIG.
 */
real_t adjust_quality(const real_t qual, const NUC prior, const NUC base, const NUC next){
   if(isambig(base)){ return MIN_QUALITY; }
   real_t new_qual = calibration_intercept + calibration_scale * qual;
   if(!isambig(prior)){
      new_qual += calibration_baseprior_adj[prior*NBASE+base];
   }
   if(!isambig(next)){
      new_qual += calibration_basenext_adj[next*NBASE+base];
   }
   if(!isambig(next) && !isambig(prior)){
      new_qual += calibration_priorbasenext_adj[(next*NBASE+prior)*NBASE+base];
   }
   return new_qual;
}

/** Return value of Mu */
real_t get_mu(void) {
    return Mu;
}

/** Set value for Mu. */
bool set_mu(const char *mu_str) {

    char *endptr;
    Mu = strtor(mu_str, &endptr);

    return (Mu > 0);
}

/**
 * Read in quality calibration conversion table if supplied and use to convert the default values.
 * Conversion formula is generally:
 * - value = original value * new scale + new value
 * - exception: no addition term for new scale
 *
 * Output resultant values if message level at least Debug.
 * Returns false if file supplied but failed to read.
 */
bool read_quality_table(void) {

    if (!matrix_from_file(E_QUALTAB)) {return true;}

    XFILE *fp = open_matrix(E_QUALTAB);
    if (fp == NULL) {return false;}

    real_t new_intercept, new_scale;
    int index = 0;

    /* read intercept and scale */
    index = E_INTERCEPT;
    if (!get_matrix_single(fp, &new_intercept)) {goto cleanup;}
    index = E_SLOPE;
    if (!get_matrix_single(fp, &new_scale)) {goto cleanup;}

    /* convert intercept and scale */
    calibration_intercept = calibration_intercept * new_scale + new_intercept;
    calibration_scale = calibration_scale * new_scale;

    /* read and convert adjustments */
    index = E_PRIOR_BASE;
    if (!get_and_convert_array(fp, NBASE, NBASE, new_scale, calibration_baseprior_adj)) {goto cleanup;}
    index = E_BASE_NEXT;
    if (!get_and_convert_array(fp, NBASE, NBASE, new_scale, calibration_basenext_adj)) {goto cleanup;}
    index = E_PRIOR_BASE_NEXT;
    if (!get_and_convert_array(fp, NBASE, NBASE * NBASE, new_scale, calibration_priorbasenext_adj)) {goto cleanup;}

    /* output resultant values */
    if (get_message_level() >= MSG_DEBUG) {
        XFILE *fpout = open_output_blk("calibration.tab", BLK_SINGLE);
        if (!xfisnull(fpout)) {
            /* show_array needs an array */
            real_t value[1];
            value[0] = calibration_intercept;
            show_array (fpout, value, 1, 1, MESS_TEXT[E_INTERCEPT]);
            value[0] = calibration_scale;
            show_array (fpout, value, 1, 1, MESS_TEXT[E_SLOPE]);
            show_array (fpout, calibration_baseprior_adj, NBASE, NBASE, MESS_TEXT[E_PRIOR_BASE]);
            show_array (fpout, calibration_basenext_adj, NBASE, NBASE, MESS_TEXT[E_BASE_NEXT]);
            show_array (fpout, calibration_priorbasenext_adj, NBASE, NBASE * NBASE, MESS_TEXT[E_PRIOR_BASE_NEXT]);
        }
        xfclose(fpout);
    }

    xfclose(fp);
    return true;

 cleanup:
    message(E_BAD_INPUT_SS, MSG_FATAL, "quality calibration conversion table", MESS_TEXT[index]);
    xfclose(fp);
    return false;
}
