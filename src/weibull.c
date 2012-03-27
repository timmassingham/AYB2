/** 
 * \file weibull.c
 * Routines for the Weibull distribution
 *//* 
 *  Created : 2010
 *  Author : Tim Massingham
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

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <err.h>
#include "statistics.h"
#include "weibull.h"

/** Cumulative density function of Weibull distribution
 * Four possible cases depending on tail or log, deal with
 * all of them accurately.
 * x is non-negative.
 * Shape and scale are strictly positive parameters of Weibull.
 * tail = true: calulate tail probability.
 * logp = true: return logarithm of cumulative density.
 */
real_t pweibull (  real_t x,  real_t shape,  real_t scale,  bool tail,  bool logp){
    validate(x>=0,NAN);
    validate(shape>0,NAN);
    validate(scale>0,NAN);
        
    real_t res = -pow(x/scale,shape);
    if( (true==tail) && (true==logp) ){ return res; }
    
    if ( (false==tail) && (false==logp) ){ return -expm1(res); } 
    
    res = exp(res);
    if ( true==tail ){ return res;}
    
    return log1p(-res);
}

/** Inverse CDF of Weibull
 * Dealing with log and tails accurately.
 * p is a probability (belongs to closed interval [0,1]).
 * shape, scale are strictly positive parameters of Weibull.
 * tail = true: p is the tail probability of a Weibull.
 * logp = true: p is given its logarithm.
 */
real_t qweibull ( real_t p,  real_t shape,  real_t scale,  bool tail,  bool logp){
    validate(shape>0,NAN);
    validate(shape>0,NAN);
    validate(logp || (p>=0 && p<=1),NAN);
    validate(!logp || p<=0., NAN);
    
    // Default cases
    if(NAN==p){ return NAN; }
    if(p<=0. && !logp){ return 0.;}
    if(p>=1. && !logp){ return HUGE_VALR;}
    if(p>=0. && logp){ return HUGE_VALR;}
  
    real_t res = NAN;
    if(logp){
       res = (true==tail) ? -p : -log(-expm1(p));
    } else {
       res = (true==tail) ? -log(p) : -log1p(-p);
    }
    
    return scale*pow(res,1.0/shape);
}

/** Weibull density distribution
 * x non-negative.
 * Shape and scale are strictly positive parameters of Weibull.
 * logd = true: return logarithm of density.
 */
real_t dweibull (  real_t x, real_t shape,  real_t scale,  bool logd ){
    validate(x>=0.0,NAN);
    validate(shape>0,NAN);
    validate(shape>0,NAN);

    if ( NAN==x){ return NAN;}
    real_t scaledx = x/scale;
    if(true==logd){ return log(shape/scale) + (shape-1) * log(scaledx) - pow(scaledx,shape);}
    
    real_t res = pow(scaledx,shape-1.);
    return (shape/scale) * res * exp( -res*scaledx );
}

/** Compare to reals
 */
int cmp_real_t ( const void * a, const void * b){
    real_t ad = *(real_t *)a;
    real_t bd = *(real_t *)b;
    
    if(ad>bd) return 1;
    if(ad<bd) return -1;
    return 0;
}

/** Fit Weibull distribution to data using a least-squares method
 * Since CDF = 1 - exp( -(x/l)^k )
 * Then log(-log(1-CDF)) = k log x - k log l
 * is a linear relationship; fit to data using OLS.
 * Returns a pair_real structure containing the fitted shape and 
 * scale.
 */
pair_real fit_weibull ( const real_t * x_orig, const bool * allowed, const uint_fast32_t n ){
    real_t * x_sorted = calloc(n,sizeof(real_t));
    real_t * y_sorted = calloc(n,sizeof(real_t));
    if(NULL==x_sorted || NULL==y_sorted){
        errx(EXIT_FAILURE,"Failed to allocate memory in %s (%s:%d)",__func__,__FILE__,__LINE__);
    }
    uint_fast32_t nallowed = 0;
    for( uint_fast32_t i=0 ; i<n ; i++){
        if(!allowed[i]){ continue; }
        x_sorted[nallowed] = log(x_orig[i]);
	nallowed++;
    }
    for( uint_fast32_t i=0 ; i<nallowed ; i++){
        y_sorted[i] = log(-log( (nallowed-i)/(nallowed+1.0) ));
    }
    qsort(x_sorted,nallowed,sizeof(*x_sorted),cmp_real_t);

    real_t * res = linearRegression(x_sorted,y_sorted,NULL,nallowed,NULL);
    pair_real resf = {NAN,NAN};
    if(NULL!=res){
        resf.e1 = res[0];
        resf.e2 = exp(-res[1]/res[0]);
        free(res);
    }
    
    free(x_sorted);
    free(y_sorted);
    return resf;
}

/** Fit Weibull distribution to data using weighted least-squares method
 * Tests suggest that the fit achieved using this method agrees more
 * closely with the maximum likelihood fitting than the OLS method used 
 * by fit_weibull although it is possible that the fit is less robust
 * if the data in not from a Weibull distribution.
 * 
 * Since CDF = 1 - exp( -(x/l)^k )
 * Then log(-log(1-CDF)) = k log x - k log l
 * is a linear relationship; fit to data using WLS.
 * 
 * The weights for the regression are calculated by taking the variance
 * of each order statistic and transforming using the delta method.
 *
 * Returns a pair_real structure containing the fitted shape and 
 * scale.
 */
pair_real wfit_weibull ( const real_t * x_orig, const bool * allowed, const uint_fast32_t n ){
    real_t * x_sorted = calloc(n,sizeof(real_t));
    real_t * y_sorted = calloc(n,sizeof(real_t));
    real_t * w_sorted = calloc(n,sizeof(real_t));
    if(NULL==x_sorted || NULL==y_sorted || NULL==w_sorted){
        errx(EXIT_FAILURE,"Failed to allocate memory in %s (%s:%d)",__func__,__FILE__,__LINE__);
    }
    uint_fast32_t nallowed = 0;
    for( uint_fast32_t i=0 ; i<n ; i++){
	if(!allowed[i]){ continue; }
        x_sorted[nallowed] = log(x_orig[i]);
	nallowed++;
    }
    for( uint_fast32_t i=0 ; i<nallowed ; i++){
        y_sorted[i] = log(-log( (nallowed-i)/(nallowed+1.0) ));
        /*  Calculate weight.
         * k'th order statistic is distributed Beta(k,n+1-k),
         * which has variance V = k * (n+1-k) / ( (n+2)(n+1)^2 )
         *  The variance of log(-log(1-p)) is therefore
         * V/((1-p)log(1-p))^2 by the delta method.
         *  Finally, the weight for the linear regression is
         * inversely proportional to the variance.
         */
        real_t p = (i+1.0)/(nallowed+1.0);
        w_sorted[i] = (i+1.0)*(nallowed-i)/((nallowed+1.0)*(nallowed+1.0)*(nallowed+2.0));
	w_sorted[i] /= (1-p)*(1-p)*log1p(-p)*log1p(-p);
        w_sorted[i] = 1.0/w_sorted[i];
    }
    qsort(x_sorted,nallowed,sizeof(*x_sorted),cmp_real_t);

    real_t * res = wLinearRegression(w_sorted,x_sorted,y_sorted,NULL,nallowed,NULL);
    pair_real resf = {NAN,NAN};
    if(NULL!=res){
        resf.e1 = res[0];
        resf.e2 = exp(-res[1]/res[0]);
        free(res);
    }

    free(w_sorted);
    free(y_sorted);
    free(x_sorted);
    return resf;
}

    
    

#ifdef TEST
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/** Test whether suffix of string matches given value.
 */
bool hassuffix( const char * str, const char * suf ){
    const uint_fast32_t suf_len = strlen(suf);
    const uint_fast32_t str_len = strlen(str);
    if ( suf_len>str_len){ return false;}
    return (0==strcmp(str+(str_len-suf_len),suf));
}
    
int main ( int argc, char * argv[] ){
    if(1==argc){
        fputs("Usage:\n"
"qweibull x shape scale tail logp\n"
"pweibull p shape scale tail logp\n"
"dweibull x shape scale log\n"
"fweibull ndata file\n"
,stderr);
        return EXIT_FAILURE;
    }
    
    real_t x,p,shape,scale;
    unsigned int tail,logp;
    printf("%s",argv[0]);

    if ( hassuffix(argv[0],"qweibull") ){
        sscanf(argv[1],REAL_FORMAT_IN,&x);
        sscanf(argv[2],REAL_FORMAT_IN,&shape);
        sscanf(argv[3],REAL_FORMAT_IN,&scale);
        sscanf(argv[4],"%u",&tail);
        sscanf(argv[5],"%u",&logp);
        printf("(%e,%e,%e,%u,%u)=",x,shape,scale,tail,logp);
        printf("%e\n",qweibull(x,shape,scale,tail,logp));
        return EXIT_SUCCESS;
    }

    if ( hassuffix(argv[0],"pweibull") ){
        sscanf(argv[1],REAL_FORMAT_IN,&p);
        sscanf(argv[2],REAL_FORMAT_IN,&shape);
        sscanf(argv[3],REAL_FORMAT_IN,&scale);
        sscanf(argv[4],"%u",&tail);
        sscanf(argv[5],"%u",&logp);
        printf("(%e,%e,%e,%u,%u)=",p,shape,scale,tail,logp);
        printf("%e\n",pweibull(p,shape,scale,tail,logp));
        return EXIT_SUCCESS;
    }
    
    if ( hassuffix(argv[0],"dweibull") ){
        sscanf(argv[1],REAL_FORMAT_IN,&x);
        sscanf(argv[2],REAL_FORMAT_IN,&shape);
        sscanf(argv[3],REAL_FORMAT_IN,&scale);
        sscanf(argv[4],"%u",&logp);
        printf("(%e,%e,%e,%u)=",x,shape,scale,logp);
        printf("%e\n",dweibull(x,shape,scale,logp));
        return EXIT_SUCCESS;
    }

    if ( hassuffix(argv[0],"fweibull") ){
        unsigned int n;
        sscanf(argv[1],"%u",&n);
        real_t * x = calloc(n,sizeof(real_t));
        FILE * fp = fopen(argv[2],"r");
        for ( int i=0 ; i<n ; i++){
            fscanf(fp,REAL_FORMAT_IN,&x[i]);
        }

	// OLS fit
	real_t loglike = 0.;
        pair_real param = fit_weibull(x,n);
        for ( uint_fast32_t i=0 ; i<n ; i++){
		loglike -= dweibull(x[i],param.e1,param.e2,true);
	}
        printf("\nOLS: shape = %f  scale = %f  loglike = %f\n",param.e1,param.e2,loglike);
	// WLS fit
        loglike = 0.;
	param = wfit_weibull(x,n);
        for ( uint_fast32_t i=0 ; i<n ; i++){
                loglike -= dweibull(x[i],param.e1,param.e2,true);
        }
        printf("WLS: shape = %f  scale = %f  loglike = %f\n",param.e1,param.e2,loglike);
	return EXIT_SUCCESS;
    }
    
    return EXIT_FAILURE;
}
#endif

