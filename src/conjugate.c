/**
 * \file conjugate.c
 * Calculate restricted fitted omega.
 *//*
 *  Created : 15 Jun 2011
 *  Author  : Hazel Marsden
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
#include "conjugate.h"
#include "lapack.h"
#include "utility.h"

/* constants */

static const int BLOCK_SIZE = 4;             ///< Same as NBASE.

/** Multiplier that is almost one to ensure that lambda can never actually be equal to min_cut. */ 
static const real_t ALMOST_ONE = 1.0 - 3.0e-8;

/* members */
/* None    */


/* private functions */

/**
 *  Update gamma_i using Polak and Ribiere method.
 *  gamma_i = (gradNew - gradOld) . gradNew / gradOld.gradOld
 */
static real_t gammaI( const real_t * gradOld, const real_t * gradNew, const unsigned int n){
    if(NULL==gradOld || NULL==gradNew){ return NAN;}

    real_t dggn=0,gogo=0;
    for ( int i=0 ; i<n ; i++){
        dggn += (gradNew[i]-gradOld[i])*gradNew[i];
        gogo += gradOld[i]*gradOld[i];
    }
    return dggn/gogo;
}

/**
 *  Calculate the objective function.
 *  Objective is tr(V U^tU ) - 2 log det (U).
 *  - u is the (upper) Choleky factorisation of Omega
 *  - info is the matrix V
 */
static real_t objective(const real_t *u, const unsigned int np , void * info){
    if(NULL==u|| NULL==info){ return NAN; }
    const MAT V = (MAT) info;

    const int n = V->nrow;
    real_t res = 0.;
    
    // Form trace(VU^tU) = 1^t(V o U^tU)1
    // Copy u into tmp and use LAPACK routine to calculate tmp := u^t tmp
    const real_t alpha = 1.0;
    real_t tmp[np];
    for( int i=0 ; i<np ; i++){ tmp[i] = u[i];}
    trmm(LAPACK_LEFT,LAPACK_UPPER,LAPACK_TRANS,LAPACK_NONUNITTRI,&n,&n,&alpha,u,&n,tmp,&n);
    for ( int i=0 ; i<np ; i++){
        res += V->x[i] * tmp[i];
    }
    // - 2 log det (U). U upper triangular => product of diagional entries is det
    for ( int i=0 ; i<n ; i++){
        res -= 2.0*log(u[i*n+i]);
    }
    return res;
}

/**
 *  Make a matrix x encoded as an n*n array into a block tridiagonal matrix.
 *  Note: uses BLOCK_SIZE which is same as NBASE.
 *        Memory access pattern is poor and could be improved.
 */
static void make_block( real_t * x, const int n){
    if(NULL==x){ return; }

    // Number of cycles (blocks)
    const int ncy = n/BLOCK_SIZE;
    for ( int i=0 ; i<ncy ; i++){
        for ( int i2=0 ; i2<BLOCK_SIZE; i2++){
            for ( int j=0 ; j<(i-1) ; j++){
                for ( int j2=0 ; j2<BLOCK_SIZE ; j2++){
                    x[(i*BLOCK_SIZE+i2)*n+j*BLOCK_SIZE+j2] = 0.0;
                    x[(j*BLOCK_SIZE+j2)*n+i*BLOCK_SIZE+i2] = 0.0;
                }
            }
        }
    }
}


/**
 *  Calculate the gradient of the objective function.
 *  Objective is tr(V U^tU ) - 2 log det (U).
 *  \n Gradient is 2 UV^t - 2 U^{-t} then forced to have block-tridiagonal
 *  upper triangular structure since these elements must stay zero 
 *  (gradient zero).
 *  - u is the Choleky factorisation of Omega
 *  - info is the matrix V
 */
static void gradObj(const real_t *u,const unsigned int np, real_t * g, void * info){
    if(NULL==u||NULL==info||NULL==g){ return; }
    const MAT V = (MAT)info;
    const int n = V->nrow;

    // Calculate full derivative
    // g := UV
    for( int i=0 ; i<np ; i++){ g[i] = V->x[i];}
    const real_t alpha = 2.0;
    trmm(LAPACK_LEFT,LAPACK_UPPER,LAPACK_NOTRANS,LAPACK_NONUNITTRI,&n,&n,&alpha,u,&n,g,&n);

    // Add together g := (tmp-g) = 2*(VU - U^{-t})
    // Note: U^{-t} is lower triangular, so only diagonal elements are needed.
    // Diagonal elements of inverse of triangular matrix are inverse of elements of diagonal
    for( int i=0 ; i<n ; i++){
        g[i*n+i] -=  2.0/u[i*n+i];
    }

    // Make matrix upper triangular
    for ( int i=0 ; i<n ; i++){
        for (int j=i+1 ; j<n ; j++){
            g[i*n+j] = 0.0;
        }
    }
    // Make matrix block tridiagonal
    make_block(g,n);
}


/**
 *  Change in the value of the objective function, starting from lambda=0.
 *  Note the starting point and direction are implicitly encoded in the
 *  scalars a, b and the vector ratio[]. See linemin_obj for details.
 */
static real_t deltaObj(real_t lambda, real_t a, real_t b, real_t * ratio, int n){
    if(NULL==ratio || !isfinite(lambda) || !isfinite(a) || !isfinite(b)){ return NAN; }

    real_t res = lambda*(a + b * lambda);
    for ( int i=0 ; i<n ; i++){
        res -= 2.0 *log1p(lambda*ratio[i]);
    }
    return res;
}

/**
 *  Gradient of change in the objective function, starting from lambda=0.
 *  Note the starting point and direction are implicitly encoded in the
 *  scalars a, b and the vector ratio[]. See linemin_obj for details.
 */
static real_t ddeltaObj(real_t lambda, real_t a, real_t b, real_t * ratio, int n){
    if(NULL==ratio || !isfinite(lambda) || !isfinite(a) || !isfinite(b)){ return NAN; }

    real_t res = a + 2.0 * b * lambda;
    for ( int i=0 ; i<n ; i++){
        res -= 2.0*ratio[i]/(1.0+lambda*ratio[i]);
    }
    return res;
}

/**
 *  Curvature of change in the objective function, starting from lambda=0.
 *  Note the starting point and direction are implicitly encoded in the
 *  scalars a, b and the vector ratio[]. See linemin_obj for details.
 */
static real_t d2deltaObj(real_t lambda, real_t a, real_t b, real_t * ratio, int n){
    if(NULL==ratio || !isfinite(lambda) || !isfinite(a) || !isfinite(b)){ return NAN; }

    real_t res = 2.0 * b;
    for ( int i=0 ; i<n ; i++){
        real_t r = ratio[i] / (1.0+lambda*ratio[i]);
        res += 2.0*r*r;
    }
    return res;
}


/**
 *  Ratio for Newton-Raphson iteration.
 *  Gradient / Curvature for the change in the objective function, starting from lambda=0.
 *  Slightly more efficient than calculating both separately.
 */
static real_t newtonObj(real_t lambda, real_t a, real_t b, real_t * ratio, int n){
    if(NULL==ratio || !isfinite(lambda) || !isfinite(a) || !isfinite(b)){ return NAN; }

    real_t top=0.0,bot=0.0;
    top = a + 2.0 * b * lambda;
    bot = 2.0 * b;
    for ( int i=0 ; i<n ; i++){
        real_t r = ratio[i]/(1.0+lambda*ratio[i]);
        top -= 2.0*r;
        bot += 2.0*r*r;
    }
    return top/bot;
}

/**
 *  Coordinate transform to ensure that lambda never crosses min_cut.
 *  (-inf,0) -> (min_cut*ALMOST_ONE,0)
 */
static real_t transform( real_t lam, real_t min_cut){
    return isfinite(min_cut) ? ALMOST_ONE*(-expm1(lam)*min_cut) : lam;
}

/**
 *  Derviative of coordinate transform to ensure that lambda never crosses min_cut.
 *  (-inf,0) -> (min_cut*ALMOST_ONE,0)
 */
static real_t dtransform( real_t lamt, real_t min_cut){
    return isfinite(min_cut) ? ALMOST_ONE*(lamt-min_cut) : 1.0;
}

/**
 *  Minimise objective function along a line.
 *  Start at u[] and minimise in direction d[].
 *  Method is Newton-Raphson, protected to ensure that object always decreases.
 *  The coordinates are transformed to ensure a cut-point in the surface
 *  (determinant of u is zero) is never crossed.
\verbatim
   Objective is tr(V(U+lam*D)^t(U+lam*D) - 2.0 * log det (U+lam*D)
              = tr(VU^tU) + lam tr(VD^tU) + lam tr(VU^tD) + lam^2 tr(VD^tD) - 2 \sum_i log(U_{ii}+lam*D_{ii})
              = tr(VU^tU) + lam a + lam^2 b - 2 \sum_i log(U_{ii}+lam*D_{ii})
   Change from lam=0
              = lam a + lam^2 b - 2 \sum_i ( log(U_{ii}+lam*D_{ii}) - log(U_{ii}) )
              = lam a + lam^2 b - 2 \sum_i log1p(lam*r_i) where r_i = D_{ii}/U_{ii}
\endverbatim
 *  Returns improvement from lam=0.
 */
static real_t linemin_obj(real_t *u,const unsigned int np, const real_t * d, void * info){
    
    MAT V = (MAT)info;
    const int n = V->nrow;
    // Find cuts in surface. Value of lambda such that any diagonal element
    // of u + lambda *d becomes zero
    real_t min_cut = -HUGE_VAL;
    for ( int i=0 ; i<n ; i++){
        real_t cut = -u[i*n+i]/d[i*n+i];
        // Only interested in negative cuts, since that is direction being optimised
        if(cut<0 && cut>min_cut){ min_cut=cut;}
    }

    // Initialise coefficients
    // First a = tr(V(U^tD+D^tU)). Note: tr(A^tB) = 1^t(A o B) 1 and V is symmetric
    real_t a = 0.0;
    real_t tmp[np];
    for( int i=0 ; i<np ; i++){ tmp[i] = d[i];}
    const real_t alpha = 1.0;
    // Calculate tmp = U^tD using LAPACK routine for triangular matrix
    trmm(LAPACK_LEFT,LAPACK_UPPER,LAPACK_TRANS,LAPACK_NONUNITTRI,&n,&n,&alpha,u,&n,tmp,&n);
    // tmp = U^tD + D^tU
    for( int i=0 ; i<n ; i++){
       for ( int j=0 ; j<=i ; j++){
           tmp[i*n+j] += tmp[j*n+i];
           tmp[j*n+i] = tmp[i*n+j];
       }
    }
    // tr( V tmp) = 1^t(V o tmp ) 1
    for ( int i=0 ; i<n ; i++){
        for ( int j=0 ; j<n ; j++){
            a += V->x[i*n+j] * tmp[i*n+j];
        }
    }

    // Second b = tr(V(D^tD)). 
    real_t b = 0.0;
    for( int i=0 ; i<np ; i++){ tmp[i] = d[i];}
    // Calculate tmp = D^tD using LAPACK routine for triangular matrix
    trmm(LAPACK_LEFT,LAPACK_UPPER,LAPACK_TRANS,LAPACK_NONUNITTRI,&n,&n,&alpha,d,&n,tmp,&n);
    // tr( V tmp) = 1^t(V o tmp ) 1
    for ( int i=0 ; i<n ; i++){
        for ( int j=0 ; j<n ; j++){
            b += V->x[i*n+j] * tmp[i*n+j];
        }
    }
    // Ratios r_i
    real_t ratio[n];
    for ( int i=0 ; i<n ; i++){
        ratio[i] = d[i*n+i]/u[i*n+i];
    }

    // Newton-Raphson iteration 
    real_t lambda = 0.0;
    int max_it = 10;
    real_t oldVal = 0.0;
    for ( int i=0 ; i<max_it ; i++){
        // Transformed lambda
        real_t lambda_tran = transform(lambda,min_cut);
        // Newton adjustment to Lambda
        // (bound by 2.0 to reduce possibility of diverging to infinity).
        real_t adj = newtonObj(lambda_tran,a,b,ratio,n)/dtransform(lambda_tran,min_cut);
        if(fabs(adj)>2.0){
            adj = (adj>0.0)?2.0:-2.0;
        }

        // Ensure solution has improved
        // Repeatedly halve adjustment until definitely have improvement.
        real_t len = 2.0;
        real_t newVal = 0.0;
        real_t lambda_new = 0.0;
        do {
            len *= 0.5;
            if(len<1e-12){ len=0.0; break; } // If len too small, give up

            // Adjust lambda and transform
            lambda_new = lambda - len * adj;
            real_t lambda_new_tran = transform(lambda_new,min_cut);
            // Calculate change in object and see if improved
            newVal = deltaObj(lambda_new_tran,a,b,ratio,n);
        } while (newVal>oldVal);

        lambda = lambda_new;
        if(oldVal-newVal<3e-8){ break; }

        oldVal = newVal;
    }

    // Update u to new point.
    const real_t lambda_tran = transform(lambda,min_cut);
    for( int i=0 ; i<np ; i++){
        u[i] += lambda_tran * d[i];
    }

    return deltaObj(lambda_tran,a,b,ratio,n);
}


/* public functions */

/**
 *  Fit restricted omega to V using conjugate gradient method.
 *  Objective function is tr(V Omega) - log det (Omega).
 *  Need to ensure that Om is a valid restricted information matrix
 *  -- symmetric positive definite and block tridiagonal.
 *  Omega is block-tridiagonal, that is composed from 4x4 blocks, e.g.
\verbatim
        Om_{11} Om_{21}^t 0
        Om_{21} Om_{22}   Om_{32}^t
        0       Om_{32}   Om_{33}
\endverbatim
 *  To ensure this, fit_omega works with the Cholesky decomposition
 *  Omega = U^tU.
 *  It can be shown (considering the block LU decomposition formulas
 *  and noting that the inverse of U is also upper triangular) that
 *  the Cholesky decomposition of Omega must retain the block tridiagonal
 *  form (although all entries in one or the other triangle must be zero).
 *
 *  \n Objective is tr(V U^tU ) - 2 log det (U).
 */  
MAT fit_omega(const MAT V, MAT omega){
    if(NULL==V){ return NULL;}
    const int n = V->nrow;

    // Initial guess is either initialOmega (previous solution)
    // or a diagonal matrix consisting of entries which are the
    // inverse of the diagonal entries of V.
    //   Note: initialising to identity matrix doesn't work very
    // well due to scaling issues.
    if(omega==NULL){
        omega = new_MAT(n,n);
        if(NULL==omega){ return NULL; }
        for ( int i=0 ; i<n ; i++){ omega->x[i*n+i] = 1.e-5 + 1.0/(1.0e-5 + V->x[i*n+i]); }
    }
    // Get Cholesky factorisation
    cholesky(omega);

    // Initial gradient 
    MAT grad = new_MAT(n,n);
    gradObj(omega->x,n*n,grad->x,(void *) V);
    // Memory for gradient update and current direction
    MAT grad2 = new_MAT(n,n);
    MAT d = copy_MAT(grad);

    // Conjugate gradient optimisation
    int max_it = 400;
    real_t initObj = objective(omega->x,n*n,(void *)V);
    for ( int i=0 ; i<max_it ; i++){
        // Minimise along direction.
        real_t delta = linemin_obj(omega->x,n*n,d->x,(void *)V);
        // Exit if change is small enough.
        if(fabs(delta)<1e-5*fabs(initObj)){ break; }

        // Update direction using Polak-Ribiere
        gradObj(omega->x,n*n,grad2->x,(void *) V);
        real_t gi = gammaI(grad->x,grad2->x,n*n);
        for ( int j=0 ; j<n*n ; j++){
            d->x[j] = grad2->x[j] + gi*d->x[j];
            grad->x[j] = grad2->x[j];
        }
    }

    free_MAT(grad);
    free_MAT(grad2);

    // Calculate omega := u^t u, where u is Cholesky factorisation
    // Using d as temporary memory. Contains Cholesky factorisation of Omega
    for ( int i=0 ; i<n*n ; i++){ d->x[i] = omega->x[i];}
    const real_t alpha = 1.0;
    // Calculate omega := u^t u using LAPACK routine. d and omega both initially contain u
    trmm(LAPACK_LEFT,LAPACK_UPPER,LAPACK_TRANS,LAPACK_NONUNITTRI,&n,&n,&alpha,d->x,&n,omega->x,&n);
    free_MAT(d);

    return omega;
}

    
#ifdef TEST
#include <stdlib.h>
#include <err.h>

const real_t vArr[] = {
    2.063893564, 1.160316934, -0.869556069, -1.108887191, 2.258075616,
    0.068911827, -0.714806040, 0.998297497, 1.091978610, 2.086450099,
    -0.095339025, -0.717392638, -1.446990545, 0.511447253, 1.150830116,
    2.311288953, 1.160316934, 1.006077872, -0.680621929, -0.564876754,
    1.377256286, -0.106242066, -0.479989751, 0.594936668, 0.684660423,
    1.334965581, -0.012061631, -0.284491352, -0.825357140, 0.288930543,
    0.578279500, 1.502367433, -0.869556069, -0.680621929, 0.797009212,
    0.369124001, -0.857357530, 0.279774777, 0.278455136, -0.591874401,
    -0.345236895, -1.180812049, -0.143870104, 0.077849617, 0.359907094,
    -0.417614047, -0.455974822, -1.031163024, -1.108887191, -0.564876754,
    0.369124001, 0.904927298, -1.396968922, -0.059661284, 0.392275881,
    -0.604721450, -0.714600299, -1.234621287, 0.156520671, 0.531024226,
    0.952267193, -0.298238779, -0.598263588, -1.447630998, 2.258075616,
    1.377256286, -0.857357530, -1.396968922, 2.892487390, 0.113370407,
    -0.918879096, 1.194234508, 1.440743317, 2.434619647, -0.235961202,
    -0.953580659, -1.859784872, 0.565171705, 1.227122299, 2.926987751,
    0.068911827, -0.106242066, 0.279774777, -0.059661284, 0.113370407,
    0.405398841, 0.002476389, -0.163844786, 0.117638671, -0.302057303,
    -0.182950280, -0.123197821, -0.226029888, -0.183819424, 0.083337908,
    0.007748548, -0.714806040, -0.479989751, 0.278455136, 0.392275881,
    -0.918879096, 0.002476389, 0.458693064, -0.380952381, -0.458718498,
    -0.755866368, 0.034722658, 0.374333678, 0.585577179, -0.111285031,
    -0.379853154, -0.832860226, 0.998297497, 0.594936668, -0.591874401,
    -0.604721450, 1.194234508, -0.163844786, -0.380952381, 0.755570735,
    0.524770905, 1.325908215, 0.024372059, -0.372717814, -0.656624057,
    0.418392807, 0.548090266, 1.293403422, 1.091978610, 0.684660423,
    -0.345236895, -0.714600299, 1.440743317, 0.117638671, -0.458718498,
    0.524770905, 0.810797776, 1.151555748, -0.181145502, -0.485862103,
    -0.972303698, 0.196065345, 0.561057468, 1.477799503, 2.086450099,
    1.334965581, -1.180812049, -1.234621287, 2.434619647, -0.302057303,
    -0.755866368, 1.325908215, 1.151555748, 2.705137974, 0.025013778,
    -0.707207748, -1.434649420, 0.752528700, 1.101922849, 2.651984103,
    -0.095339025, -0.012061631, -0.143870104, 0.156520671, -0.235961202,
    -0.182950280, 0.034722658, 0.024372059, -0.181145502, 0.025013778,
    0.208459434, 0.085203299, 0.231473131, 0.078947790, -0.033341825,
    -0.213835305, -0.717392638, -0.284491352, 0.077849617, 0.531024226,
    -0.953580659, -0.123197821, 0.374333678, -0.372717814, -0.485862103,
    -0.707207748, 0.085203299, 0.611815976, 0.727670567, -0.100674764,
    -0.441197798, -0.817965456, -1.446990545, -0.825357140, 0.359907094,
    0.952267193, -1.859784872, -0.226029888, 0.585577179, -0.656624057,
    -0.972303698, -1.434649420, 0.231473131, 0.727670567, 1.425114366,
    -0.229561851, -0.765665285, -1.843120492, 0.511447253, 0.288930543,
    -0.417614047, -0.298238779, 0.565171705, -0.183819424, -0.111285031,
    0.418392807, 0.196065345, 0.752528700, 0.078947790, -0.100674764,
    -0.229561851, 0.399133302, 0.325160699, 0.716417884, 1.150830116,
    0.578279500, -0.455974822, -0.598263588, 1.227122299, 0.083337908,
    -0.379853154, 0.548090266, 0.561057468, 1.101922849, -0.033341825,
    -0.441197798, -0.765665285, 0.325160699, 0.797520921, 1.234101682,
    2.311288953, 1.502367433, -1.031163024, -1.447630998, 2.926987751,
    0.007748548, -0.832860226, 1.293403422, 1.477799503, 2.651984103,
    -0.213835305, -0.817965456, -1.843120492, 0.716417884, 1.234101682,
    3.242693662
};

real_t omArr[] = {
    11.4968758009, 0.0182690275, 3.8550846192, 1.4096675416, -2.4696355299,
    -3.3922273755, 0.7882389718, 1.1136483476, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0182690275, 13.8926190070, 7.6212607813, -2.2938435623,
    -4.4284025520, 5.6444517283, 0.9055813717, 4.4308224903, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 3.8550846192, 7.6212607813, 16.5801535174,
    -1.1408584745, -3.0400974709, 2.7464495336, -2.6359135286, 3.8705849913,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 1.4096675416, -2.2938435623,
    -1.1408584745, 7.6649890079, 0.5177368268, -1.5919179555, 1.7168056426,
    -0.9240682006, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -2.4696355299,
    -4.4284025520, -3.0400974709, 0.5177368268, 18.0596324000, -3.2882785798,
    6.5904983987, -4.0120845577, -6.4959705229, -0.5286384483, 3.5464304683,
    1.5405487881, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    -3.3922273755, 5.6444517283, 2.7464495336, -1.5919179555, -3.2882785798,
    16.7634355861, -5.1060489664, -1.5244953762, -2.0631188052, 8.5545795324,
    1.3891479396, -0.5231096395, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.7882389718, 0.9055813717, -2.6359135286, 1.7168056426,
    6.5904983987, -5.1060489664, 13.9082039692, 2.9456425930, 2.8387511237,
    -3.9608667989, 1.3219484749, -3.8331951065, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 1.1136483476, 4.4308224903, 3.8705849913,
    -0.9240682006, -4.0120845577, -1.5244953762, 2.9456425930, 18.0930242167,
    9.1670942998, -7.8062211002, 2.7483603023, 4.3781397326, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, -6.4959705229, -2.0631188052, 2.8387511237,
    9.1670942998, 28.3041673559, -6.3185345715, 5.1891298494, 4.3088464016,
    -4.9085602142, 9.5889394929, 2.3575773352, -9.8064868772, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, -0.5286384483, 8.5545795324,
    -3.9608667989, -7.8062211002, -6.3185345715, 13.2697860780, -3.0354579183,
    -0.2981326680, 3.6816284053, 1.6313900639, -0.7211693403, -0.1357010780,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 3.5464304683,
    1.3891479396, 1.3219484749, 2.7483603023, 5.1891298494, -3.0354579183,
    19.0936578361, 6.5140658964, -1.3872799658, -3.6840630263, 0.2445515141,
    0.2489848243, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    1.5405487881, -0.5231096395, -3.8331951065, 4.3781397326, 4.3088464016,
    -0.2981326680, 6.5140658964, 13.5433250618, -5.0260812134, -0.0007715702,
    3.6474993863, -3.7380628661, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    -4.9085602142, 3.6816284053, -1.3872799658, -5.0260812134, 11.7963722286,
    -2.2269460787, -2.9421146681, 3.3171112384, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 9.5889394929, 1.6313900639, -3.6840630263, -0.0007715702,
    -2.2269460787, 25.0374602411, -5.0691072841, -9.7727852921, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 2.3575773352, -0.7211693403, 0.2445515141,
    3.6474993863, -2.9421146681, -5.0691072841, 10.6287970672, -0.0408203854,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, -9.8064868772, -0.1357010780,
    0.2489848243, -3.7380628661, 3.3171112384, -9.7727852921, -0.0408203854,
    13.2194288752
};

// Objective is 72.68812

int main(int argc, char * argv[]){
    if(argc==1){
        MAT V = new_MAT_from_array(16,16,vArr);
        /*MAT Om0 = new_MAT_from_array(16,16,omArr);
        cholesky(Om0);
        show_MAT(xstdout,Om0,0,0);
        fprintf(stdout,"Objective = %e\n",objective(Om0->x,16*16,(void *)V));

        MAT grad = new_MAT(16,16);
        gradObj(Om0->x,16*16,grad->x,(void *) V);
        show_MAT(xstdout,grad,0,0);

        MAT grad2 = new_MAT(16,16);
        MAT d = copy_MAT(grad);
        for ( int i=0 ; i<30 ; i++){
            real_t delta = linemin_obj(Om0->x,16*16,d->x,(void *)V);
            fprintf(stdout,"Improvement is %e\n",delta);
            fprintf(stdout,"Objective = %e\n",objective(Om0->x,16*16,(void *)V));
            // Update direction
            gradObj(Om0->x,16*16,grad2->x,(void *) V);
            real_t gi = gammaI(grad->x,grad2->x,16*16);
            for ( int i=0 ; i<16*16 ; i++){
                d->x[i] = grad2->x[i] + gi*d->x[i];
                grad->x[i] = grad2->x[i];
            }
        }*/

        MAT omega = fit_omega(V,NULL);
        MAT omegaInv = invert_symmetric(omega);
        show_MAT(xstdout,V,0,0);
        show_MAT(xstdout,omegaInv,0,0);
        show_MAT(xstdout,omega,0,0);
    } 

    else {  
        //read from file
        FILE * fp = fopen(argv[1],"r");
        int nr,nc;
        fscanf(fp,"%d%d",&nr,&nc);
        if(nr!=nc){ errx(EXIT_FAILURE,"nrow!=ncol");}
        real_t * x = calloc(nr*nc,sizeof(real_t));
        for ( int i=0 ; i<nr*nc ; i++){
               fscanf(fp,REAL_FORMAT_IN,&x[i]);
        }
        MAT m = new_MAT_from_array(nr,nc,x);
        MAT omega2 = fit_omega(m,NULL);
        MAT omegaInv2 = invert_symmetric(omega2);
        //show_MAT(xstdout,m,0,0);
        //show_MAT(xstdout,omegaInv2,0,0);
        show_MAT(xstdout,omega2,0,0);
    }

    return EXIT_SUCCESS;
}
#endif
