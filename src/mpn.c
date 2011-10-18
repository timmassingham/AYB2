/** 
 * \file mpn.c
 * Calculations of terms for Parameter Estimation.
 *//* 
 *  Created : 9 Jun 2010
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
#include <strings.h>
#include "lapack.h"
#include "mpn.h"
#include "nuc.h"
#include "statistics.h"
#include "utility.h"


/* constants */
/* None      */


/* members */
/* None    */


/* private functions */

static inline real_t sum_squares(const real_t * x, const uint32_t n){
    validate(NULL!=x,NAN);
    real_t res = 0.;
    for ( uint32_t i=0 ; i<n ; i++){
        res += x[i] * x[i];
    }
    return res;
}


/* public functions */

// Crude method of obtaining weights
MAT calculateWe( const MAT lssi, MAT we){
    validate(NULL!=lssi,NULL);
    const uint32_t ncluster = lssi->nrow;
    if(NULL==we){
        we = new_MAT(ncluster,1);
        validate(NULL!=we,NULL);
    }
    memset(we->x, 0, ncluster*sizeof(real_t));

    real_t meanLSSi = mean(lssi->x,ncluster);
    real_t varLSSi = variance(lssi->x,ncluster);
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        const real_t d = lssi->x[cl]-meanLSSi;
        we->x[cl] = cauchy(d*d,varLSSi);
    }
    return we;
}


//MAT calculateIbar( const ARRAY(int16_t) intmat, const MAT we, MAT Ibar){
MAT calculateIbar( const TILE tile, const MAT we, MAT Ibar){
    validate(NULL!=tile,NULL);
    validate(NULL!=we,NULL);
    const uint32_t ncluster = tile->ncluster;
    const uint32_t ncycle = tile->ncycle;

    if(NULL==Ibar){
        Ibar = new_MAT(NBASE,ncycle);
        validate(NULL!=Ibar,NULL);
    }
    validate(Ibar->nrow==NBASE,NULL);
    validate(Ibar->ncol==ncycle,NULL);
    memset(Ibar->x, 0, Ibar->nrow*Ibar->ncol*sizeof(real_t));

    unsigned int cl = 0;
    LIST(CLUSTER) node = tile->clusterlist;
    while (NULL != node && cl < ncluster){
        for( uint32_t idx=0 ; idx<NBASE*ncycle ; idx++){
            Ibar->x[idx] += node->elt->signals->xint[idx] * we->x[cl];
        }

        /* next cluster */
        node = node->nxt;
        cl++;
    }

    return Ibar;
}

MAT calculateSbar( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const uint32_t ncycle, MAT Sbar){
    validate(NULL!=lambda,NULL);
    validate(NULL!=we,NULL);
    validate(NULL!=bases.elt,NULL);
    const uint32_t ncluster = lambda->nrow;
    validate(ncluster==we->nrow,NULL);
    if(NULL==Sbar){
        Sbar = new_MAT(NBASE,ncycle);
        validate(NULL!=Sbar,NULL);
    }
    validate(Sbar->nrow==NBASE,NULL);
    validate(Sbar->ncol==ncycle,NULL);
    memset(Sbar->x, 0, Sbar->nrow*Sbar->ncol*sizeof(real_t));

    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            int base = bases.elt[cl*ncycle+cy];
            if(!isambig(base)){
                Sbar->x[cy*NBASE+base] += we->x[cl] * lambda->x[cl];
            }
        }
    }

    return Sbar;
}
real_t calculateWbar( const MAT we){
    validate(NULL!=we,NAN);
    const uint32_t ncluster = we->nrow;

    real_t wbar = 0.;
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        wbar += we->x[cl];
    }
    return wbar;
}


MAT calculateJ( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const uint_fast32_t ncycle, MAT J){
    validate(NULL!=lambda,NULL);
    validate(NULL!=we,NULL);
    validate(NULL!=bases.elt,NULL);
    const uint_fast32_t ncluster = lambda->nrow;

    if(NULL==J){
        J = new_MAT(NBASE*NBASE,ncycle*ncycle);
        validate(NULL!=J,NULL);
    }
    memset(J->x, 0, J->nrow*J->ncol*sizeof(real_t));

    const uint_fast32_t lda = NBASE*NBASE;
    for ( uint_fast32_t cl=0 ; cl<ncluster ; cl++){
        const real_t welam = we->x[cl] * lambda->x[cl] * lambda->x[cl];
        for ( uint_fast32_t cy=0 ; cy<ncycle ; cy++){
            const uint_fast32_t base = bases.elt[cl*ncycle+cy];
            if(!isambig(base)){
                const uint_fast32_t offset = cy*ncycle*NBASE*NBASE + base*NBASE;
                for ( uint_fast32_t cy2=0 ; cy2<ncycle ; cy2++){
                    const uint_fast32_t base2 = bases.elt[cl*ncycle+cy2];
                    if(!isambig(base2)){
                        J->x[offset+cy2*lda+base2] += welam;
                    }
                }
            }
        }
    }

    return J;
}

// tmp should be NBASE*NBASE+ncycle*ncycle
real_t calculateDeltaLSE(const MAT Mt, const MAT P, const MAT N, const MAT J, const MAT K, real_t * tmp){
    // BLAS definitions
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;
    // Sanity checks
    validate(NULL!=Mt,NAN);
    validate(NULL!=P,NAN);
    validate(NULL!=N,NAN);
    validate(NULL!=J,NAN);
    validate(NULL!=K,NAN);
    validate(NULL!=tmp,NAN);
    const int ncycle = P->nrow;
    const int nbase = NBASE;

    real_t delta = -sum_squares(N->x,ncycle*NBASE); // tr(N^tN) = 1^t (N o N) 1 = sum_squares
    delta -= 2.0 * xMy(Mt->x,K,P->x);
    // PP^t
    gemm(LAPACK_NOTRANS,LAPACK_TRANS,&ncycle,&ncycle,&ncycle,&alpha,P->x,&ncycle,P->x,&ncycle,&beta,tmp+NBASE*NBASE,&ncycle);
    // MtM
    gemm(LAPACK_NOTRANS,LAPACK_TRANS,&nbase, &nbase, &nbase, &alpha,Mt->x,&nbase, Mt->x,&nbase, &beta,tmp,&nbase);
    delta += xMy(tmp,J,tmp+NBASE*NBASE);

    return delta;
}

MAT calculateK( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const TILE tile, const uint_fast32_t ncycle, MAT K){
    LIST(CLUSTER) node = NULL;
    uint_fast32_t cl = 0;
    real_t * tmp = NULL;

    validate(NULL!=lambda,NULL);
    validate(NULL!=we,NULL);
    validate(NULL!=tile,NULL);
    const uint_fast32_t ncluster = tile->ncluster;
    const uint_fast32_t lda = NBASE*NBASE;

    if(NULL==K){
        K = new_MAT(NBASE*NBASE,ncycle*ncycle);
        if(NULL==K){ goto cleanup; }
    }
    memset(K->x, 0, K->nrow*K->ncol*sizeof(real_t));

    tmp = calloc(NBASE*ncycle,sizeof(real_t));
    if(NULL==tmp){ goto cleanup; }

    node = tile->clusterlist;
    while (NULL != node && cl < ncluster){
        const bool has_ambig = has_ambiguous_base(bases.elt+cl*ncycle, ncycle);
        const real_t welam = we->x[cl] * lambda->x[cl];
        for ( uint_fast32_t i=0 ; i<(NBASE*ncycle) ; i++){
           tmp[i] = welam * node->elt->signals->xint[i];
        }
        if (has_ambig) {
            for ( uint_fast32_t cy=0 ; cy<ncycle ; cy++){
                const uint_fast32_t ioffset = cy*NBASE;
                for ( uint_fast32_t cy2=0 ; cy2<ncycle ; cy2++){
                    const uint_fast32_t base = bases.elt[cl*ncycle+cy2];
                    if(!isambig(base)){
                        const uint_fast32_t koffset = cy*lda*ncycle + cy2*lda + base;
                        for ( uint_fast32_t ch=0 ; ch<NBASE ; ch++){
                           K->x[ koffset + ch*NBASE] += tmp[ioffset + ch];
                        }
                    }
                }
            }
        }
        else {
            for ( uint_fast32_t cy=0 ; cy<ncycle ; cy++){
                const uint_fast32_t ioffset = cy*NBASE;
                for ( uint_fast32_t cy2=0 ; cy2<ncycle ; cy2++){
                    const uint_fast32_t base = bases.elt[cl*ncycle+cy2];
                    const uint_fast32_t koffset = cy*lda*ncycle + cy2*lda + base;
                    for ( uint_fast32_t ch=0 ; ch<NBASE ; ch++){
                        K->x[ koffset + ch*NBASE] += tmp[ioffset + ch];
                    }
                }
            }
        }

        /* next cluster */
        node = node->nxt;
        cl++;
    }
    free(tmp);

    return K;

cleanup:
    free(tmp);
    free_MAT(K);
    return NULL;
}


// tmp should be NBASE*NBASE+ncycle*ncycle
MAT calculateMlhs( const MAT var, const real_t wbar, const MAT SbarT, const MAT P, const MAT Jt, real_t * tmp, MAT lhs){
    validate(NULL!=var,NULL);
    validate(SbarT!=NULL,NULL);
    validate(P!=NULL,NULL);
    validate(Jt!=NULL,NULL);
    validate(tmp!=NULL,NULL);
    const int unit = 1;
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;
    const int ncycle = P->nrow;

    if(NULL==lhs){
        lhs = new_MAT(NBASE+ncycle,NBASE+ncycle);
        validate(NULL!=lhs,NULL);
    }
    memset(lhs->x, 0, lhs->nrow*lhs->ncol*sizeof(real_t));

    const uint32_t lda = NBASE+ncycle;
    // Reshape Jvec(P diag(var) Pt) via P diag(sqrt(v)) %*% ( P diag(sqrt(v)) )^t
    // Scale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] *= sqrt(var->x[cy]);
        }
    }
    gemm(LAPACK_NOTRANS,LAPACK_TRANS,&ncycle,&ncycle,&ncycle,&alpha,P->x,&ncycle,P->x,&ncycle,&beta,tmp+NBASE*NBASE,&ncycle);
    gemv(LAPACK_TRANS,&Jt->nrow,&Jt->ncol,&alpha,Jt->x,&Jt->nrow,tmp+NBASE*NBASE,&unit,&beta,tmp,&unit);
    for ( uint32_t cy=0 ; cy<NBASE ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[cy*lda+base] = tmp[cy*NBASE+base];
        }
    }
    // Unscale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] /= sqrt(var->x[cy]);
        }
    }

    // Sbar %*% P %*% diag(var)
    gemm(LAPACK_TRANS,LAPACK_NOTRANS,&SbarT->ncol,&P->ncol,&P->nrow,&alpha,SbarT->x,&SbarT->nrow,P->x,&P->nrow,&beta,lhs->x+NBASE*lda,&lhs->nrow);
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[NBASE*lda+cy*lda+base] *= var->x[cy];
        }
    }
    // Copy Sbar %*% P into new bit of array
{
    const uint32_t offset1 = lda * NBASE;
    const uint32_t offset2 = NBASE;
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[offset2+base*lda+cy] = lhs->x[offset1+cy*lda+base];
        }
    }
}
    // Id matrix
{
    const uint32_t offset = NBASE*lda+NBASE;
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        lhs->x[offset+cy*lda+cy] = wbar*var->x[cy];
    }
}

    return lhs;
}


// tmp should be NBASE*NBASE long
MAT calculateMrhs( const MAT var, const MAT IbarT, const MAT P, const MAT Kt, real_t * tmp, MAT rhs){
    validate(NULL!=var,NULL);
    validate(NULL!=IbarT,NULL);
    validate(NULL!=P,NULL);
    validate(NULL!=Kt,NULL);
    validate(NULL!=tmp,NULL);
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;

    const uint32_t ncycle = P->nrow;
    const int lda = NBASE + ncycle;
    if(NULL==rhs){
        rhs = new_MAT(lda,NBASE);
        validate(NULL!=rhs,NULL);
    }
    memset(rhs->x, 0, rhs->nrow*rhs->ncol*sizeof(real_t));

    // Reshape KVec(Pdiag(v))
    // Scale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] *= var->x[cy];
        }
    }
    gemv(LAPACK_TRANS,&Kt->nrow,&Kt->ncol,&alpha,Kt->x,&Kt->nrow,P->x,LAPACK_UNIT,&beta,tmp,LAPACK_UNIT);
    for ( uint32_t base1=0 ; base1<NBASE ; base1++){
        for ( uint32_t base2=0 ; base2<NBASE ; base2++){
            rhs->x[base1*lda+base2] = tmp[base1*NBASE+base2];
        }
    }
    // Unscale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] /= var->x[cy];
        }
    }

    // Copy in diag(v)*IbarT
    for ( uint32_t base=0 ; base<NBASE ; base++){
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            rhs->x[NBASE+base*lda+cy] = IbarT->x[base*ncycle+cy] * var->x[cy];
        }
    }
    return rhs;
}


// tmp should be NBASE*NBASE+ncycle*ncycle
MAT calculatePlhs( const real_t wbar, const MAT Sbar, const MAT Mt, const MAT J, real_t * tmp, MAT lhs){
    validate(Sbar!=NULL,NULL);
    validate(Mt!=NULL,NULL);
    validate(J!=NULL,NULL);
    validate(tmp!=NULL,NULL);
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;
    const int ncycle = Sbar->ncol;
    const int nbase = NBASE;

    if(NULL==lhs){
        lhs = new_MAT(ncycle,ncycle);
        validate(NULL!=lhs,NULL);
    }
    memset(lhs->x, 0, lhs->nrow*lhs->ncol*sizeof(real_t));

    const int lda = ncycle;
    // Reshape Jtvec(MtM)
    gemm(LAPACK_NOTRANS, LAPACK_TRANS,&nbase, &nbase, &nbase, &alpha,Mt->x,&nbase, Mt->x,  &nbase, &beta, tmp+ncycle*ncycle,&nbase);
    gemv(LAPACK_TRANS,&J->nrow,&J->ncol,&alpha,J->x,&J->nrow,tmp+ncycle*ncycle,LAPACK_UNIT,&beta,tmp,LAPACK_UNIT);
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            lhs->x[cy*lda+cy2] = tmp[cy*ncycle+cy2];
        }
    }

    return lhs;
}


// tmp should be ncycle*ncycle+NBASE*ncycle long
MAT calculatePrhs( const MAT Ibar, const MAT Mt, const MAT Sbar, const MAT N, const MAT K, real_t * tmp, MAT rhs){
    validate(NULL!=Ibar,NULL);
    validate(NULL!=Mt,NULL);
    validate(NULL!=Sbar,NULL);
    validate(NULL!=N,NULL);
    validate(NULL!=K,NULL);
    validate(NULL!=tmp,NULL);
    const real_t alpha = 1.0;
    const real_t negalpha = -1.0;
    const real_t  beta = 0.0;

    const int ncycle = Ibar->ncol;
    const int lda = ncycle;
    const int nbase = Mt->nrow;
    if(NULL==rhs){
        rhs = new_MAT(lda,ncycle);
        validate(NULL!=rhs,NULL);
    }
    memset(rhs->x, 0, rhs->nrow*rhs->ncol*sizeof(real_t));

    // Reshape KtVec(M)
    gemv(LAPACK_TRANS,&K->nrow,&K->ncol,&alpha,K->x,&K->nrow,Mt->x,LAPACK_UNIT,&beta,tmp,LAPACK_UNIT);
    //  - Sbar^t M^t N
    gemm(LAPACK_NOTRANS,LAPACK_NOTRANS,&nbase,&ncycle,&nbase,&alpha,Mt->x,&nbase,N->x,&nbase,&beta,tmp+ncycle*ncycle,&nbase);
    gemm(LAPACK_TRANS,LAPACK_NOTRANS,&ncycle,&ncycle,&nbase,&negalpha,Sbar->x,&nbase,tmp+ncycle*ncycle,&nbase,&alpha,tmp,&ncycle);

    for ( uint32_t cy1=0 ; cy1<ncycle ; cy1++){
        for ( uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            rhs->x[cy1*lda+cy2] = tmp[cy1*ncycle+cy2];
        }
    }

    return rhs;
}

/**
 * Calculates matrix J used in calculateLhs.
 * J is the matrix \\sum_i we_i lambda_i lambda_i Vec(S_i) Vec(S_i)^t.
 */
MAT calculateNewJ(const MAT lambda, const ARRAY(NUC) bases, const MAT we, const int ncycle, MAT newJ){
   if(NULL==lambda || NULL==bases.elt || NULL==we){ return NULL;}

   // Allocate memory if necessary and initialise to zero
   if(NULL==newJ){
      newJ = new_MAT(ncycle*NBASE,ncycle*NBASE);
      if(NULL==newJ){ return NULL; }
   }
   memset(newJ->x, 0, newJ->nrow*newJ->ncol*sizeof(real_t));

   const uint_fast32_t lda = ncycle * NBASE;
   const uint_fast32_t ncluster = we->nrow;
   for ( uint_fast32_t cl=0 ; cl<ncluster ; cl++){
      const real_t eltmult = we->x[cl] * lambda->x[cl] * lambda->x[cl];
      for ( uint_fast32_t i=0 ; i<ncycle ; i++){
         const uint_fast32_t base = bases.elt[cl*ncycle+i];
         if(!isambig(base)){
             const uint_fast32_t idx1 = i*NBASE+base;
             for ( uint_fast32_t j=0 ; j<ncycle ; j++){
                const uint_fast32_t base2 = bases.elt[cl*ncycle+j];
                if(!isambig(base2)){
                    const uint_fast32_t idx2 = j*NBASE+base2;
                    newJ->x[idx1*lda+idx2] += eltmult;
                }
             }
         }
      }
   }

   return newJ;
}

/**
 * Calculates matrix K used in calculateRhs.
 * K is the matrix \\sum_i we_i lambda_i Vec(S_i) Vec(I_i)^t.
 * First calculate its transpose (better memory layout).
 */
MAT calculateNewK(const MAT lambda, const ARRAY(NUC) bases, const TILE tile, const MAT we, const int ncycle, MAT newK){
    if(NULL==lambda || NULL==bases.elt || NULL==tile || NULL==we){ return NULL;}

    // Allocate memory if necessary and initialise to zero
    if(NULL==newK){
        newK = new_MAT(ncycle*NBASE,ncycle*NBASE);
        if(NULL==newK){ return NULL; }
    }
    memset(newK->x, 0, newK->nrow*newK->ncol*sizeof(real_t));

	// Calculate transpose
    const uint_fast32_t ncluster = tile->ncluster;
    const uint_fast32_t lda = ncycle * NBASE;
    uint_fast32_t cl = 0;
    LIST(CLUSTER) node = tile->clusterlist;
    while (NULL != node && cl < ncluster){
        for ( uint_fast32_t i=0 ; i<ncycle ; i++){
            const uint_fast32_t base = bases.elt[cl*ncycle+i];
            if(!isambig(base)){
                const uint_fast32_t col = i*NBASE + base;
                const real_t colmult = we->x[cl] * lambda->x[cl];
                for ( uint_fast32_t j=0 ; j<lda ; j++){
                    newK->x[col*lda+j] += node->elt->signals->xint[j] * colmult;
                }
            }
        }

        /* next cluster */
        node = node->nxt;
        cl++;
    }
    // Transpose (square) matrix newK
    transpose_inplace(newK);

    return newK;
}

/**
 * Calculates Lhs of equation to solve for \\hat{A}^t
 * where Lhs \\hat{A}^t = Rhs.
\verbatim
                 _                          _
Lhs is a matrix |      J        Vec(S_bar)   |
                |_  Vec(S_bar)^t      wbar  _|
\endverbatim
 */
MAT calculateLhs( const real_t wbar,const MAT J, const MAT Sbar, MAT lhs){
    if( NULL==J || NULL==Sbar){ return NULL; }

    // Allocate memory for lhs if necessary. One more row and column than J
    if(NULL==lhs){
        lhs = new_MAT(J->nrow+1,J->nrow+1);
        if(NULL==lhs){ return NULL;}
    }
    const int lda = J->nrow+1;
    const int nelt = J->nrow;

    for ( int i=0 ; i<nelt ; i++){
        // Copy in J
        for ( int j=0 ; j<nelt ; j++){
            lhs->x[i*lda+j] = J->x[i*nelt+j];
        }
        // Copy Vec(Sbar) and Vec(Sbar)^t
        lhs->x[i*lda+nelt] = Sbar->x[i];
        lhs->x[nelt*lda+i] = Sbar->x[i];
    }

    // wbar
    lhs->x[lda*lda-1] = wbar;
    return lhs;
}

/**
 * Calculates Rhs of equation to solve for \\hat{A}^t
 * where Lhs \\hat{A}^t = Rhs.
\verbatim
                 _               _
Rhs is a matrix |      K          |
                |_  Vec(I_bar)^t _|
\endverbatim
 */
MAT calculateRhs( const MAT K, const MAT Ibar, MAT rhs){
    if( NULL==K || NULL==Ibar ){ return NULL; }

    if( NULL==rhs){
        rhs = new_MAT(K->nrow+1,K->nrow);
    }
    const int nelt = K->nrow;
    const int lda = K->nrow+1;

    for ( int i=0 ; i<nelt ; i++){
        for ( int j=0 ; j<nelt ; j++){
            rhs->x[i*lda+j] = K->x[i*nelt+j];
        }
        rhs->x[i*lda+nelt] = Ibar->x[i];
    }

    return rhs;
}

/** 
 * Solve system of linear equations using Cholesky decomposition.
 * Wrapper for LAPACK routine.
 * Assumes that lhs is positive-definite, which problems being considered are.
 * Both lhs and rhs are destructively updated.
 * Result is stored in rhs.
 */
int solverChol( MAT lhs, MAT rhs, real_t * null){
    validate(NULL!=lhs,-4);
    validate(NULL!=rhs,-6);
    validate(lhs->nrow==lhs->ncol,-2);
    validate(lhs->ncol==rhs->nrow,-3);

    int info = 0;
    posv(LAPACK_UPPER, &lhs->ncol, &rhs->ncol, lhs->x, &lhs->nrow, rhs->x, &rhs->nrow, &info);
    if(0!=info){ fprintf(stderr,"Error in solver, info = %d\n",info);}
    return info;
}

/** 
 * Solve system of linear equations using SVD.
 * Wrapper for LAPACK routine.
 * Both lhs and rhs are destructively updated.
 * Result is stored in rhs.
 * tmp should be 6*N.
 */
int solverSVD(MAT lhs, MAT rhs, real_t * tmp, const real_t delta_diag) {
    validate(NULL!=lhs,-4);
    validate(NULL!=rhs,-6);
    validate(NULL!=tmp,-11);

    const int N = lhs->nrow;
    int INFO=0,RANK=0,IWORK=5*N;
    real_t RCOND = 3e-8;

    /* add a diagonal offset to avoid failure to solve */
    if(0.0 != delta_diag){
        for (int i = 0; i < N; i++){
            lhs->x[i * N + i] += delta_diag;
        }
    }

    gelss(&lhs->nrow,&lhs->ncol,&rhs->ncol,lhs->x,&lhs->nrow,
                      rhs->x,&rhs->nrow,tmp,&RCOND,&RANK,tmp+N,&IWORK,&INFO);
    return INFO;
}

/** 
 * Solve system of linear equations using SVD, setting negative results to zero.
 * Wrapper for LAPACK routine.
 * Both lhs and rhs are destructively updated.
 * Result is stored in rhs.
 * tmp should be 6*N.
 */
int solverZeroSVD(MAT lhs, MAT rhs, real_t * tmp, const real_t delta_diag){
    int info = solverSVD(lhs,rhs,tmp, delta_diag);
    if(info==0){
        // Success
        const int N = lhs->nrow;
        for ( uint32_t i=0 ; i<N*N ; i++){
            if(rhs->x[i]<0){ rhs->x[i] = 0.;}
        }
    }
    return info;
}


#ifdef FORTRAN
/** 
 * Solve system of linear equations using non-negative least squares.
 * Wrapper for Fortran routine in s/dnnls.f.
 * Both lhs and rhs are destructively updated.
 * Result is stored in rhs.
 * lhs must be square, tmp should be N*rhs->ncol.
 */
int solverNNLS(MAT lhs, MAT rhs, real_t *tmp, const real_t delta_diag){
#ifdef NFORTRAN
    /* Eclipse IDE does not handle fortran */
    return 0;
#else
    validate(NULL!=lhs,-4);
    validate(NULL!=rhs,-6);
    validate(NULL!=tmp,-11);

    const int N = lhs->nrow;
    const int rcol = rhs->ncol;
    /* lhs square and rhs row match */
    validate(N==lhs->ncol,-2);
    validate(N==rhs->nrow,-3);

    real_t RNORM;
    real_t W[N];
    real_t ZZ[N];
    int INDEX[N],MODE=0;

    /* add a diagonal offset to avoid failure to solve */
    if(0.0 != delta_diag){
        for (int i = 0; i < N; i++){
            lhs->x[i * N + i] += delta_diag;
        }
    }

    /* nnls only calculates for a vector so loop for each column of rhs
       results collected sequentially in tmp
       lhs and rhs destructively updated so use working copies, 
       all of lhs and column at a time of rhs 
    */
    real_t * lhs_tmp = malloc(N*N*sizeof(real_t));
    real_t * rhs_tmp = malloc(N*sizeof(real_t));
    for( int cy=0 ; cy<rcol ; cy++){
        memcpy(lhs_tmp,lhs->x,N*N*sizeof(real_t));
        memcpy(rhs_tmp,rhs->x+cy*N,N*sizeof(real_t));
        nnls(lhs_tmp,&N,&N,&N,rhs_tmp,tmp+cy*N,&RNORM,W,ZZ,INDEX,&MODE);
    }
    /* result is expected in rhs */
    memcpy(rhs->x, tmp, N*rcol*sizeof(real_t));

    free(lhs_tmp);
    free(rhs_tmp);
    return MODE;
#endif
}
#endif


#ifdef TEST
#include <stdio.h>
/* Inputs */
static  NUC base_array[] = {3, 2, 3, 3, 3, 1, 3, 1, 3, 0, 1, 3, 3, 0, 2};
static  real_t we[] = {0.33081400, 0.29916646, 0.79687810, 0.04376649, 0.04898477};
static  real_t lambda[] = {2.852168, 11.890647, 10.675880, 5.542096, 3.245923};
static  real_t ints_t[] = {
    0.32822142, 0.3611550, 0.4056395, 0.6217474, 0.3701723, 0.54075999, 0.9004344, 0.8360398, 0.7729722, 0.9753776, 0.06152802, 0.5637200,
    0.03651698, 0.8598354, 0.4986751, 0.7010750, 0.3627614, 0.06679484, 0.9104111, 0.8990147, 0.4424424, 0.2002516, 0.11924224, 0.2969899,
    0.54164990, 0.9133029, 0.5054887, 0.5377883, 0.2405506, 0.58718250, 0.3922185, 0.1585344, 0.8682236, 0.3153632, 0.86569760, 0.8310858,
    0.61345680, 0.9237610, 0.9891983, 0.1324787, 0.3127209, 0.69892800, 0.2603790, 0.2256760, 0.8400769, 0.9111117, 0.27118720, 0.5472079,
    0.94282040, 0.4265982, 0.1394416, 0.7856237, 0.2671884, 0.11008570, 0.4843887, 0.1912447, 0.7967254, 0.5434273, 0.83477990, 0.7742892
};
static real_t P[] = {
    0.9445303038, 0.70435005, 0.74991005,
    0.0006557733, 0.05396262, 0.05993442,
    0.0018982966, 0.08522693, 0.54507703
};
static real_t M[] = {
     0.8390111, 0.8550323, 0.6927066, 0.7199632,
     0.4642314, 0.8833080, 0.1549395, 0.1361410,
     0.1612448, 0.3195291, 0.5993373, 0.4147634,
     0.8614796, 0.3521542, 0.2954812, 0.7250016
};

/* Expected Results, 2dp output

J matrix:
1:     1.34     0.00     0.00     0.00     0.52     0.00     0.00     0.00     0.00
2:     0.00     1.34     0.00     0.00     0.00     0.00     0.00     0.00     0.00
3:     0.00     0.00     0.00     0.00     0.00     0.52     0.00     0.00     0.00
4:     0.00     0.00     1.34     0.52     0.00     0.00     0.00     0.00     0.00
5:     0.00     0.00     0.00     1.34     0.00     0.00     0.00     0.00     0.00
6:     0.00     0.00     0.00     0.00    92.17     0.00     0.00     0.00    42.30
7:     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00
8:     0.00     0.00     0.00    90.82     0.00    92.17    42.30    42.30     0.00
9:     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.52     0.00
10:     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00
11:     0.00     0.00     0.00     0.00     2.69     0.00     0.00     0.00     0.52
12:     0.00     0.00     0.00     2.69     0.00     2.69     0.52     0.00     0.00
13:     0.00     0.52     0.00     0.00     0.00     0.00     1.34     0.00     0.00
14:     0.00    90.82    42.30     0.00     0.00    42.30     0.00    92.17     0.00
15:     0.00     2.69     0.52     0.00     0.00     0.00     0.00     2.69     0.00
16:   136.33    42.30    93.51    42.30    42.30     0.00    93.51     0.00    94.86
K matrix:
1:     0.15     0.15     0.00     0.08     0.04     0.00     0.20     0.13     0.00
2:     0.00     4.76     0.13     0.00     2.12     1.29     0.00     7.59     1.57
3:     0.00     0.31     0.15     0.00     0.35     0.04     0.00     0.73     0.13
4:     5.20     0.13     5.07     3.73     1.29     2.47     9.82     1.57     8.32
5:     0.22     0.07     0.00     0.17     0.02     0.00     0.22     0.09     0.00
6:     0.00     7.99     3.06     0.00     5.16     0.24     0.00     2.90     0.71
7:     0.00     0.34     0.07     0.00     0.51     0.02     0.00     0.92     0.09
8:    11.24     3.06     8.33     5.76     0.24     5.68     4.40     0.71     3.82
9:     0.24     0.02     0.00     0.06     0.08     0.00     0.07     0.13     0.00
10:     0.00     4.54     1.77     0.00     3.40     3.24     0.00     7.43     0.42
11:     0.00     0.38     0.02     0.00     0.85     0.08     0.00     0.06     0.13
12:     6.48     1.77     4.92     7.50     3.24     4.25     7.98     0.42     7.49
13:     0.03     0.12     0.00     0.05     0.03     0.00     0.13     0.12     0.00
14:     0.00     4.61     2.49     0.00     1.40     3.20     0.00     7.20     1.06
15:     0.00     0.59     0.12     0.00     0.79     0.03     0.00     0.53     0.12
16:     7.78     2.49     5.19     5.37     3.20     2.19     8.78     1.06     7.73

Sbar matrix:
1:     0.24     0.16     0.00
2:     0.00     8.75     3.56
3:     0.00     0.94     0.16
4:    13.17     3.56     9.69
Ibar matrix:
1:     0.62     0.45     1.16
2:     1.17     0.70     0.70
3:     0.74     0.92     0.80
4:     0.89     0.69     1.00
Wbar:    1.52

M lhs matrix:
1:     1.46     0.89     0.30     1.30     0.34     0.01     0.01
2:     0.89    83.17     0.00   168.15     8.83     0.69     2.68
3:     0.30     0.00     1.81     3.71     0.78     0.06     0.17
4:     1.30   168.15     3.71   413.88    22.21     0.78     5.61
5:     0.34     8.83     0.78    22.21     1.52     0.00     0.00
6:     0.01     0.69     0.06     0.78     0.00     1.52     0.00
7:     0.01     2.68     0.17     5.61     0.00     0.00     1.52
M rhs matrix:
1:     0.26     0.27     0.26     0.13
2:     5.14     8.85     5.77     6.57
3:     0.48     0.44     0.41     0.66
4:    13.71    21.53    15.63    17.63
5:     0.62     1.17     0.74     0.89
6:     0.45     0.70     0.92     0.69
7:     1.16     0.70     0.80     1.00
M solution matrix:
1:     0.10    -0.03     0.25     0.07
2:    -0.06    -0.04    -0.04    -0.03
3:    -0.05    -0.40     0.34     0.35
4:     0.03    -0.01     0.09     0.07
5:     0.38     1.39    -0.88    -0.46
6:     0.31     0.50     0.56     0.42
7:     0.78     0.62     0.21     0.42

P lhs matrix:
1:   204.92   144.94   177.23
2:   144.94   161.29   117.57
3:   177.23   117.57   184.56
P rhs matrix:
1:    16.50    11.64    19.23
2:    15.04    11.92    11.75
3:    16.01     8.79    18.08
P ls solution matrix:
1:    -0.03     0.04     0.08
2:     0.07     0.06    -0.03
3:     0.07    -0.03     0.04
P zero solution matrix:
1:     0.00     0.04     0.08
2:     0.07     0.06     0.00
3:     0.07     0.00     0.04
P nnls solution matrix (success=1):
1:     0.00     0.01     0.03
2:     0.06     0.07     0.00
3:     0.05     0.00     0.07
*/
int main ( void){
    const uint32_t ncluster = 5;
    const uint32_t ncycle = 3;

    /* use null variance for now */
    MAT var = new_MAT(ncycle, 1);
    set_MAT(var, 1.0);

    /* put input values into required data structures */
    fputs("Inputs:\n", stdout);
    ARRAY(NUC) bases = coerce_ARRAY(NUC)(ncluster * ncycle, base_array);
    fputs("Bases array:\n", stdout);
    show_ARRAY(NUC)(xstdout, bases, "", ncluster * ncycle);
    MAT matWe = coerce_MAT_from_array(ncluster,1,we);
    fputs("Weight matrix:\n",stdout);
    show_MAT(xstdout,matWe,0,0);
    MAT matLambda = coerce_MAT_from_array(ncluster,1,lambda);
    fputs("Lambda matrix:\n",stdout);
    show_MAT(xstdout,matLambda,0,0);

    /* values originally supplied as small reals; increase size and change to integer */
    int_t * tmp2 = malloc(ncluster * ncycle * NBASE * sizeof(int_t));
    for (unsigned int i = 0; i < ncluster * ncycle * NBASE; i++) {
        tmp2[i] = (int_t)(roundr(ints_t[i] * 1000));
    }
    TILE tileInts = coerce_TILE_from_array(ncluster, ncycle, tmp2);
    fputs("Intensity tile:\n", stdout);
    show_TILE(xstdout, tileInts, 10);

    MAT matM = coerce_MAT_from_array(NBASE,NBASE,M);
    MAT matMt = transpose(matM);
    fputs("Initial M matrix:\n",stdout);
    show_MAT(xstdout,matM,0,0);
    MAT matP = coerce_MAT_from_array(ncycle,ncycle,P);
    fputs("Initial P matrix:\n",stdout);
    show_MAT(xstdout,matP,0,0);
    MAT N = new_MAT(NBASE,ncycle);
    fputs("Initial N matrix:\n",stdout);
    show_MAT(xstdout,N,0,0);

    /* pre-calculation terms */
    fputs("\nResults:\n", stdout);
    MAT J = calculateJ(matLambda,matWe,bases,ncycle,NULL);
    MAT Jt = transpose(J);
    fputs("J matrix:\n",stdout);
    show_MAT(xstdout,J,0,0);

    MAT K = calculateK(matLambda,matWe,bases,tileInts,ncycle,NULL);
    MAT Kt = transpose(K);
    fputs("K matrix:\n",stdout);
    show_MAT(xstdout,K,0,0);

    MAT matSbar = calculateSbar(matLambda, matWe, bases, ncycle, NULL);
    MAT matSbarT = transpose(matSbar);
    fputs("Sbar matrix:\n",stdout);
    show_MAT(xstdout,matSbar,0,0);

    MAT matIbar = calculateIbar(tileInts, matWe, NULL);
    MAT matIbarT = transpose(matIbar);
    fputs("Ibar matrix:\n",stdout);
    show_MAT(xstdout,matIbar,0,0);

    real_t Wbar = calculateWbar(matWe);
#ifdef NDEBUG
    xfprintf(xstdout, "Wbar:%#8.2f\n", Wbar);
#else
    xfprintf(xstdout, "Wbar:%#12.6f\n", Wbar);
#endif

    /* calculate M, P */
    real_t * tmp = calloc(NBASE*NBASE*ncycle*ncycle, sizeof(real_t));

    MAT lhs = calculateMlhs(var,Wbar,matSbarT,matP,Jt,tmp,NULL);
    fputs("M lhs matrix:\n",stdout);
    show_MAT(xstdout,lhs,0,0);
    MAT rhs = calculateMrhs(var,matIbarT,matP,Kt,tmp,NULL);
    fputs("M rhs matrix:\n",stdout);
    show_MAT(xstdout,rhs,0,0);

    int retM = solverSVD(lhs, rhs, tmp, 0.0);
    xfprintf(xstdout,"M solution matrix (info=%d):\n",retM);
    show_MAT(xstdout,rhs,0,0);

    MAT Plhs = calculatePlhs(Wbar,matSbar,matMt,J,tmp,NULL);
    fputs("P lhs matrix:\n",stdout);
    show_MAT(xstdout,Plhs,0,0);
    MAT Prhs = calculatePrhs(matIbar,matMt,matSbar,N,K,tmp,NULL);
    fputs("P rhs matrix:\n",stdout);
    show_MAT(xstdout,Prhs,0,0);

    /* copy to allow multiple P solvers */
    MAT Plhs_copy = copy_MAT(Plhs);
    MAT Prhs_copy = copy_MAT(Prhs);

    /* SVD */
    int retP = solverSVD(Plhs_copy, Prhs_copy, tmp, 0.0);
    xfprintf(xstdout,"P ls solution matrix (info=%d):\n",retP);
    show_MAT(xstdout,Prhs_copy,0,0);

    /* ZeroSVD */
    Plhs_copy = copyinto_MAT(Plhs_copy, Plhs);
    Prhs_copy = copyinto_MAT(Prhs_copy, Prhs);

    retP = solverZeroSVD(Plhs_copy, Prhs_copy, tmp, 0.0);
    xfprintf(xstdout,"P zero solution matrix (info=%d):\n",retP);
    show_MAT(xstdout,Prhs_copy,0,0);

#ifdef FORTRAN
    /* NNLS */
    Plhs_copy = copyinto_MAT(Plhs_copy, Plhs);
    Prhs_copy = copyinto_MAT(Prhs_copy, Prhs);

    retP = solverNNLS(Plhs_copy, Prhs_copy, tmp, 0.0);
    xfprintf(xstdout,"P nnls solution matrix (info=%d):\n",retP);
    show_MAT(xstdout,Prhs_copy,0,0);
#endif

    /* nothing freed but exiting anyway */
    return EXIT_SUCCESS;
}
#endif
