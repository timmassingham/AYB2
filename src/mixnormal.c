/** 
 * \file mixnormal.c
 * Routines for fitting a mixed normal distribution.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mixnormal.h"

/* constants */
/* None      */

/* members */
/* None */


/* private functions */

/** 
 * Initialise the distribution data structure for the supplied data. 
 * Parameters are number of values supplied and number of distributions in the mix.
 * Returns a new mixed normal data structure.
 */
static NormMixParam initialise_mixture(const real_t * x, const unsigned int n, const unsigned int nmix){
	if(NULL==x){ return NULL;}
	if ((n==0) || (nmix==0)) { return NULL;}
	NormMixParam param = new_NormMixParam(nmix);
	if(NULL==param){ return NULL; }

	real_t mean = 0.0;
	real_t var = 0.0;
	for ( int i=0 ; i<n ; i++){
		mean += x[i];
		var += x[i]*x[i];
	}
	mean /= n;
	var = var/n - mean*mean;
	real_t sd = sqrt(var);

	for ( int j=0 ; j<nmix ; j++){
		param->prob[j] = 1.0/nmix;
		param->mean[j] = mean + (2.0-4.0*((j+1.0)/(nmix+1.0))) * sd;
		param->sd[j] = sd/nmix;
	}

	return param;
}

/**
 * Perform one optimisation loop for the supplied data.
 * Parameters are number of values supplied and mixed normal data structure.
 */
static real_t emstep_mixnormal(const real_t * x, const unsigned int n, NormMixParam param){
	if(NULL==param){ return NAN; }
	if(NULL==x){ return NAN; }
	const unsigned int nmix = param->nmix;
	real_t m_mean[nmix];
	real_t m_var[nmix];
	real_t m_we[nmix];
	real_t logprob[nmix];
	real_t logsd[nmix];
	for ( int j=0 ; j<nmix ; j++){
		m_mean[j] = m_var[j] = m_we[j] = 0.0;
		logprob[j] = log(param->prob[j]);
		logsd[j] = log(param->sd[j]);
	}
	real_t loglike = 0.0;
	for ( int i=0 ; i<n ; i++){
		// E-step, conditional probability
		real_t cprob[nmix];
		real_t tot = 0.0;
		real_t maxp = -HUGE_VAL;
		for ( int j=0 ; j<nmix ; j++){
			real_t tstat = (x[i]-param->mean[j])/param->sd[j];
			cprob[j] = logprob[j] -0.5*tstat*tstat - logsd[j];
			if(cprob[j]>maxp){ maxp = cprob[j];}
		}
		for ( int j=0 ; j<nmix ; j++){
			cprob[j] = exp(cprob[j]-maxp);
			tot += cprob[j];
		}
		loglike += log(tot) + maxp;
		for ( int j=0 ; j<nmix ; j++){
			cprob[j] /= tot;
		}

		// M-step, accumulate for estimates of mixture mean and variance
		for ( int j=0 ; j<nmix ; j++){
			m_mean[j] += cprob[j] * x[i];
			m_var[j] += cprob[j] * x[i] * x[i];
			m_we[j] += cprob[j];
		}
	}
	
	for ( int j=0 ; j<nmix ; j++){
		param->mean[j] = m_mean[j] / m_we[j];
		param->sd[j] = sqrt(m_var[j]/m_we[j] - param->mean[j]*param->mean[j]);
		param->prob[j] = m_we[j]/n;
	}

	return loglike - 0.5 * log(2.0*M_PI) * n;
}


/* public functions */

NormMixParam new_NormMixParam( const unsigned int nmix){
	NormMixParam param = calloc(1,sizeof(*param));
	if(NULL==param){ return NULL;}
	param->nmix = nmix;
	param->prob = calloc(nmix,sizeof(real_t));
	param->mean = calloc(nmix,sizeof(real_t));
	param->sd = calloc(nmix,sizeof(real_t));

	if (NULL==param->prob || NULL==param->mean || NULL==param->sd ){
		free_NormMixParam(param);
		param = NULL;
	}
	return param;
}

void free_NormMixParam(NormMixParam param){
	if(NULL==param){ return; }
	free(param->prob);
	free(param->mean);
	free(param->sd);
	free(param);
}

void show_NormMixParam( FILE * fp, const NormMixParam param){
	fprintf(fp,"Normal mixture, %u mixtures.\n",param->nmix);
	for ( int i=0 ; i<param->nmix ; i++){
		fprintf(fp,"%3d: %6.4f %6.4f %6.4f\n",i+1,param->prob[i],param->mean[i],param->sd[i]);
	}
}

/** 
 * Fit a mixed normal distribution to the supplied data.
 * Parameters are number of values supplied, 
 * number of distributions in the mix and number of iterations to use.
 * Returns a new mixed normal data structure.
 */
NormMixParam fit_mixnormal(const real_t * x, const unsigned int n, const unsigned int nmix, const unsigned niter){
	if(NULL==x){ return NULL;}
	NormMixParam param = initialise_mixture(x,n,nmix);
	if(NULL==param){ return NULL; }
	for ( int i=0 ; i<niter ; i++){
		emstep_mixnormal(x,n,param);
	}
	return param;
}


#ifdef TEST
#include <err.h>

int main ( int argc, char * argv[]){
    if(argc<4){
	    /* arguments are number of distributions in mix, number of iterations and input file of reals */
		errx(EXIT_FAILURE, "Usage: test-mixnormal nmix niter filename");
	}
	
	const unsigned nmix = atoi(argv[1]);
	const unsigned niter = atoi(argv[2]);
	if(nmix==0){ errx(EXIT_FAILURE,"zero nmix parameter supplied\n");}
	FILE * fp = fopen(argv[3],"r");
	if(NULL==fp){ errx(EXIT_FAILURE,"Failed to open %s\n",argv[3]);}

	unsigned int nelt = 0;
	fscanf(fp,"%u",&nelt);
	real_t * x = calloc(nelt,sizeof(real_t));
	for ( int i=0 ; i<nelt ; i++){
		fscanf(fp,REAL_FORMAT_IN,&x[i]);
	}
	fclose(fp);

    if (nelt==0) {
        fputs("No values in input file\n", stdout);
        return EXIT_FAILURE;
    }
            
    NormMixParam param = NULL;
    fputs("Initialise distribution null values\n", stdout);
	param = initialise_mixture(NULL,nelt,nmix);
	if (param==NULL) {
        fputs("Return value null, ok\n", stdout);
    }
    else {
        fputs("Return value not null, not ok\n", stdout);
        free_NormMixParam(param);
	}
	
    fputs("Initialise distribution no values\n", stdout);
	param = initialise_mixture(x,0,nmix);
	if (param==NULL) {
        fputs("Return value null, ok\n", stdout);
    }
    else {
        fputs("Return value not null, not ok\n", stdout);
        free_NormMixParam(param);
	}
	
    fputs("Initialise distribution nmix zero\n", stdout);
	param = initialise_mixture(x,nelt,0);
	if (param==NULL) {
        fputs("Return value null, ok\n", stdout);
    }
    else {
        fputs("Return value not null, not ok\n", stdout);
        free_NormMixParam(param);
	}
	
    fputs("Fit mixed normal distribution\n", stdout);
	param = initialise_mixture(x,nelt,nmix);
	for ( int i=0 ; i<niter ; i++){
		show_NormMixParam(stdout,param);
		real_t l = emstep_mixnormal(x,nelt,param);
		fprintf(stdout,"loglike = %f\n",l);
	}
	show_NormMixParam(stdout,param);
    free_NormMixParam(param);

    return EXIT_SUCCESS;
}
#endif /* TEST */

