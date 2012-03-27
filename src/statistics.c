/**
 * \file statistics.c
 * Statistical Functions.
 *//*
 *  Created : 2010
 *  Author : Tim Massingham
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
 
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "statistics.h"


real_t sum( const real_t * x, const bool * allowed, const uint_fast32_t n){
    if(NULL==x){return NAN;}
    
    real_t sum = 0.0;
    real_t c = 0.;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t y = x[i] - c;
        volatile real_t t = sum + y;
        c = (t-sum) - y;
        sum = t;
    }
    return sum;
}

int cmpReal(const void * x, const void * y){
   if( *(real_t*)(x) == *(real_t*)(y) ){ return 0;}
   return *(real_t*)(x)>*(real_t*)(y)?1:-1;
}

bool isodd(const uint_fast32_t n){
   return n%2;
}

// Median by sorting.
// Divide and conquer would be quicker (linear vs n log(n) )
real_t median( const real_t * x, const bool * allowed, const uint_fast32_t n){
   if(NULL==x){ return NAN;}
   if(0==n){ return NAN;}

   real_t * xc = calloc(n,sizeof(double));
   if(NULL==xc){ return NAN;}

   uint_fast32_t nallowed = 0;
   for ( int i=0 ; i<n ; i++){
       if(NULL!=allowed && !allowed[i]){ continue; }
       xc[nallowed] = x[i];
       nallowed++;
   }
   qsort(xc,nallowed,sizeof(real_t),cmpReal);

   int minIdx = (nallowed-1)/2;
   real_t ret = isodd(nallowed)?xc[minIdx]:0.5*(xc[minIdx]+xc[minIdx+1]);
   free(xc);

   return ret;
}

    
/*  Mean by compensated summation
 * See: http://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
real_t mean( const real_t * x, const bool * allowed, const uint_fast32_t n){
    if(NULL==x){return NAN;}
   
    uint_fast32_t nallowed = 0; 
    real_t sum = 0.0;
    real_t c = 0.;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t y = x[i] - c;
        volatile real_t t = sum + y;
        c = (t-sum) - y;
        sum = t;
	nallowed++;
    }
    return sum/nallowed;
}

real_t meanxx( const real_t * x, const bool * allowed, const uint_fast32_t n){
    if(NULL==x){return NAN;}
   
    uint_fast32_t nallowed = 0; 
    real_t sum = 0.0;
    real_t c = 0.;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t y = x[i]*x[i] - c;
        volatile real_t t = sum + y;
        c = (t-sum) - y;
        sum = t;
	nallowed++;
    }
    return sum/nallowed;
}

real_t meanxy( const real_t * x, const real_t * y, const bool * allowed, const uint_fast32_t n){
    if(NULL==x){return NAN;}
    if(NULL==y){return NAN;}
    
    uint_fast32_t nallowed = 0;
    real_t sum = 0.0;
    real_t c = 0.;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t r = x[i]*y[i] - c;
        volatile real_t t = sum + r;
        c = (t-sum) - r;
        sum = t;
	nallowed++;
    }
    return sum/nallowed;
}


/*  Weighted mean by compensated summation
 * See: http://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
real_t wmean( const real_t * w, const real_t * x, const bool * allowed, const uint_fast32_t n){
    if(NULL==w || NULL==x ){return NAN;}
    
    real_t sum = 0.0;
    real_t wsum = 0.0;
    real_t c = 0.;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t y = w[i]*x[i] - c;
        volatile real_t t = sum + y;
        c = (t-sum) - y;
        sum = t;
        // Sum of weights not compensated
        wsum += w[i];
    }
    return sum/wsum;
}

real_t wmeanxx( const real_t * w, const real_t * x, const bool * allowed, const uint_fast32_t n){
    if( NULL==w || NULL==x ){return NAN;}
    
    real_t sum = 0.0;
    real_t c = 0.;
    real_t wsum = 0.0;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t y = w[i]*x[i]*x[i] - c;
        volatile real_t t = sum + y;
        c = (t-sum) - y;
        sum = t;
        // Sum of weights not compensated
        wsum += w[i];
    }
    return sum/wsum;
}

real_t wmeanxy( const real_t * w, const real_t * x, const real_t * y, const bool * allowed, const uint_fast32_t n){
    if(NULL==x){return NAN;}
    if(NULL==y){return NAN;}
    
    real_t sum = 0.0;
    real_t c = 0.;
    real_t wsum = 0.0;
    for ( uint_fast32_t i=0; i<n ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        volatile real_t r = w[i]*x[i]*y[i] - c;
        volatile real_t t = sum + r;
        c = (t-sum) - r;
        sum = t;
        // Sum of weights not compensated
        wsum += w[i];
    }
    return sum/wsum;
}


/* Variance by compensated summation
 */
real_t variance( const real_t * x, const bool * allowed, const uint_fast32_t n){
        if(NULL==x){return NAN;}
        
        const real_t m = mean(x,allowed,n);
	uint_fast32_t nallowed = 0;
        real_t v = 0.0;
        real_t c = 0.;
        for ( uint_fast32_t i=0; i<n ; i++){
            if(NULL!=allowed && !allowed[i]){ continue; }
            volatile real_t y = (x[i] - m)*(x[i]-m) - c;
            volatile real_t t = v + y;
            c = (t-v) - y;
            v = t;
	    nallowed++;
        }
        return v/(nallowed-1);
}

/* Weighted variance by compensated summation
 */
real_t wvariance( const real_t * w, const real_t * x, const bool * allowed, const uint_fast32_t n){
        if(NULL==w || NULL==x){return NAN;}
        
        const real_t m = wmean(x,w,allowed,n);
        real_t v = 0.;
        real_t c = 0.;
        real_t wsum = w[0];
	uint_fast32_t nallowed = 0;
        for ( uint_fast32_t i=0; i<n ; i++){
            if(NULL!=allowed && !allowed[i]){ continue; }
            volatile real_t y = w[i]*(x[i] - m)*(x[i]-m) - c;
            volatile real_t t = v + y;
            c = (t-v) - y;
            v = t;
            // Sum of weights not compensated
            wsum += w[i];
	    nallowed++;
        }
        return v/(nallowed-1.) * nallowed/wsum;
}

/* Returns [ slope, constant, minLS ]
 * where minLS = sum_i w_i res_i^2
 */
real_t * wLinearRegression( const real_t * w, const real_t * x, const real_t * y, const bool * allowed, const uint_fast32_t n, real_t *res){
    if(NULL==w){ return NULL;}
    if(NULL==x){ return NULL;}
    if(NULL==y){ return NULL;}
    if(n==0){return NULL;} // Matter of definition

    if(NULL==res){ res = calloc(3,sizeof(real_t));}
    if(NULL==res){ return NULL;}
    
    real_t wmx  = wmean(w,x,allowed,n);
    real_t wmy  = wmean(w,y,allowed,n);
    real_t wmxx = wmeanxx(w,x,allowed,n); // Could be made more generic, to wmean(w,f(x),n);
    real_t wmxy = wmeanxy(w,x,y,allowed,n);
    real_t wmyy = wmeanxx(w,y,allowed,n);
    real_t ws   = sum(w,allowed,n);
    
    res[0] = (wmxy-wmx*wmy)/(wmxx-wmx*wmx); // Solve LS problem for slope
    res[1] = wmy - res[0]*wmx; // Error remaining
    res[2] = ws*(wmyy - res[0] * wmxy - res[1] * wmy);
    
    return res;
}

real_t * linearRegression( const real_t * x, const real_t * y, const bool * allowed, const uint_fast32_t n, real_t *res){
    if(NULL==x){ return NULL;}
    if(NULL==y){ return NULL;}
    if(n==0){return NULL;} // Matter of definition

    if(NULL==res){ res = calloc(3,sizeof(real_t));}
    if(NULL==res){ return NULL;}
    
    real_t mx  = mean(x,allowed,n);
    real_t my  = mean(y,allowed,n);
    real_t mxx = meanxx(x,allowed,n); // Could be made more generic, to wmean(w,f(x),n);
    real_t mxy = meanxy(x,y,allowed,n);
    real_t myy = meanxx(y,allowed,n);
    uint_fast32_t nallowed = 0;
    if(NULL==allowed){
	    nallowed = n;
    } else {
        for ( uint_fast32_t i=0 ; i<n ; i++){
	    nallowed += allowed[i];
	}
    }
    
    res[0] = (mxy-mx*my)/(mxx-mx*mx); // Solve LS problem for slope
    res[1] = my - res[0]*mx; // Error remaining
    res[2] = nallowed*(myy - res[0] * mxy - res[1] * my);
    
    return res;
}

real_t * linearResiduals( const real_t * x, const real_t * y, const bool * allowed, const real_t * param, const uint_fast32_t nobs, real_t * resid){
    if(NULL==x){return NULL;}
    if(NULL==y){return NULL;}
    if(NULL==param){return NULL;}
    
    if(NULL==resid){
        resid = calloc(nobs,sizeof(real_t));
        if(NULL==resid){ return NULL;}
    }
    for ( uint_fast32_t i=0 ; i<nobs ; i++){
	if(NULL!=allowed && !allowed[i]){ continue; }
        resid[i] = y[i] - param[0]*x[i] - param[1];
    }
    
    return resid;
}
            
    

/* Code for iterated least squares estimation */
real_t __attribute__((const)) cauchy( const real_t xsqr, const real_t v){
   const real_t ratio = xsqr/v;
   return 1./(1.+0.5*ratio);
}

real_t __attribute__((const)) tukey_biweight( const real_t xsqr, const real_t v){
   const real_t ratio = xsqr/v;
   return (ratio < 36.)?(1.-ratio/36.)*(1.-ratio/36.):0.;
}

real_t * iwlsLinearRegression(real_t (*f)(const real_t,const real_t), const real_t * x, const real_t * y, const bool * allowed, const uint_fast32_t niter, const uint_fast32_t n, real_t * res){
    if(NULL==x){return NULL;}
    if(NULL==y){return NULL;}

    if(NULL==res){
        res = calloc(3,sizeof(real_t));
        if(NULL==res){return NULL;}
    }
    real_t * w = calloc(n,sizeof(real_t));
    if(NULL==w){return NULL;}
    for ( uint_fast32_t i=0 ; i<n ; i++){ w[i] = 1.;}
    
    
    // Initial linear regression
    res = wLinearRegression(w,x,y,allowed,n,res);
    real_t * resid = NULL;
    for ( uint_fast32_t iter=0 ; iter<niter ; iter++){
        // First, robust estimate of scale
        resid = linearResiduals(x,y,allowed,res,n,resid);
        real_t v = wvariance(w,resid,allowed,n);
        // New weights
        for( uint_fast32_t i=0 ; i<n ; i++){
            w[i] = f(resid[i]*resid[i],v);
        }
        // Weighted regression
        res = wLinearRegression(w,x,y,allowed,n,res);
    }
    xfree(resid);
    return res;        
}

#ifdef TEST
#include <stdio.h>
#include <stdlib.h>

real_t * readmatrix( FILE * fp, int nrow, int ncol ){
    if(NULL==fp){ return NULL;}
    real_t * mat = calloc(nrow*ncol,sizeof(real_t));
    if(NULL==mat){ return NULL;}
    for ( int row=0 ; row<nrow ; row++){
        for ( int col=0 ; col<ncol ; col++){
            fscanf(fp,"%lf",&mat[col*nrow+row]);
        }
        fprintf(stdout,"Row %d",row);
        for ( int col=0 ; col<ncol ; col++){
            fprintf(stdout,"\t%f",mat[col*nrow+row]);
        }
        fputc('\n',stdout);
    }
    
    return mat;
}
        
int main( int argc, char *argv[]){
    int nobs;
    if ( 3!=argc){ fputs("Usage: wLR nobs filename\n",stderr); return EXIT_FAILURE;}
 
    FILE * fp = fopen(argv[2],"r");
    sscanf(argv[1],"%d",&nobs);
    real_t * vars = readmatrix(fp,nobs,3);

    fputs("* Standard linear regression.\n",stdout);
    real_t * res = linearRegression(vars+nobs,vars+2*nobs,nobs,NULL);
    fprintf(stdout,"Slope = %f\nConstant = %f\ndiffLS = %f\n",res[0],res[1],res[2]);

    fputs("* Weighted linear regression, using weights from file.\n",stdout);
    res = wLinearRegression(vars,vars+nobs,vars+2*nobs,nobs,res);
    fprintf(stdout,"Slope = %f\nConstant = %f\ndiffLS = %f\n",res[0],res[1],res[2]);
    
    fputs("* Iteratively reWeighted Least Squares, Cauchy.\n",stdout); 
    res = iwlsLinearRegression(cauchy,vars+nobs,vars+2*nobs,100,nobs,res);
    fprintf(stdout,"Slope = %f\nConstant = %f\ndiffLS = %f\n",res[0],res[1],res[2]);
    fputs("* Iteratively reWeighted Least Squares, Tukey biweight.\n",stdout); 
    res = iwlsLinearRegression(tukey_biweight,vars+nobs,vars+2*nobs,100,nobs,res);
    fprintf(stdout,"Slope = %f\nConstant = %f\ndiffLS = %f\n",res[0],res[1],res[2]);
    
    return EXIT_SUCCESS;
}

#endif
