/**
 * \file lapack.h
 * Description.
 *//* 
 *  Created : 2010
 *  Author : Tim Massingham
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
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

#ifndef LAPACK_H_
#define LAPACK_H_

#ifdef __linux
#define F77_NAME(A) A ## _
#else
#define F77_NAME(A) A
#endif

static const char LAPACK_UPPER[] = {'U'};
static const char LAPACK_LOWER[] = {'L'};
static const char LAPACK_UNITTRI[] = {'U'};
static const char LAPACK_NONUNITTRI[] = {'N'};
static const char LAPACK_TRANS[] = {'T'};
static const char LAPACK_NOTRANS[] = {'N'};
static const char LAPACK_LEFT[] = {'L'};
static const char LAPACK_RIGHT[] = {'R'};
static const int  LAPACK_UNIT[] = {1};


// Float functions
void F77_NAME(spotrf)( const char * uplo, const int * n, float * A, const int * lda, int * info );
void F77_NAME(spotri)( const char * uplo, const int * n, float * A, const int * lda, int * info );
void F77_NAME(strtri)( const char * uplo, const char * diag, const int * n, float * a, const int * lda, int * info);
void F77_NAME(strmm)(  const char * direct, const char * uplo, const char * trans, const char *diag, const int * m , const int * n, 
                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb);

void F77_NAME(sposv)(  const char * uplo, const int * N, const int * NRHS, float * A, const int * lda, float * B, const int * ldb, const int * info);
void F77_NAME(sgelss)( const int * M, const int * N, const int * NRHS, float * A,
                       const int * lda, float * B, const int * ldb, float * S,
                       float * RCOND, int * rank, float * WORK, int * LWORK, int * INFO);
void F77_NAME(sgemm)(  const char * transa, const char * transb,
                       const int * M, const int * N, const int *K,
                       const float * alpha, const float * A, const int * lda,
                       const float * B, const int * ldb, const float * beta,
                       float * C, const int * ldc);
void F77_NAME(sgemv)(  const char * trans, const int * M, const int * N,
                       const float * alpha, const float * A, const int * lda,
                       const float * X, const int * incx, const float * beta,
                       float * y, const int * incy);

void F77_NAME(sgetrf)( const int * M, const int * N, float * A, const int * lda,
                       int * ipiv, int * info);
void F77_NAME(sgetri)( const int * N, float * A, const int * lda, const int * ipiv,
                       float * work, const int * lwork, int * info);
void F77_NAME(sgetrs)( const char * trans, const int * N, const int * NRHS, const float * A, const int * lda, 
                       const int * ipiv, float * B, const int * ldb, int * info);

void F77_NAME(ssyr)(   const char * uplo, const int * N, const float * alpha, const float * x,
                       const int * incx, float * A, const int * lda);

// Non-LAPACK routine for non-negative least squares
void F77_NAME(snnls)(  float *A, const int * MDA, const int * M, const int* N, float * B, float * X, 
                       float * RNORM, float * W, float * ZZ, int * INDEX, int * MODE);

// Double functions
void F77_NAME(dpotrf)( const char * uplo, const int * n, double * A, const int * lda, int * info );
void F77_NAME(dpotri)( const char * uplo, const int * n, double * A, const int * lda, int * info );
void F77_NAME(dtrtri)( const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);
void F77_NAME(dtrmm)(  const char * direct, const char * uplo, const char * trans, const char *diag, const int * m , const int * n, 
                       const double * alpha, const double * A, const int * lda, const double * B, const int * ldb);

void F77_NAME(dposv)(  const char * uplo, const int * N, const int * NRHS, double * A, const int * lda, double * B, const int * ldb, const int * info);
void F77_NAME(dgelss)( const int * M, const int * N, const int * NRHS, double * A,
                       const int * lda, double * B, const int * ldb, double * S,
                       double * RCOND, int * rank, double * WORK, int * LWORK, int * INFO);
void F77_NAME(dgemm)(  const char * transa, const char * transb,
                       const int * M, const int * N, const int *K,
                       const double * alpha, const double * A, const int * lda,
                       const double * B, const int * ldb, const double * beta,
                       double * C, const int * ldc);
void F77_NAME(dgemv)(  const char * trans, const int * M, const int * N,
                       const double * alpha, const double * A, const int * lda,
                       const double * X, const int * incx, const double * beta,
                       double * y, const int * incy);

void F77_NAME(dgetrf)( const int * M, const int * N, double * A, const int * lda,
                       int * ipiv, int * info);
void F77_NAME(dgetri)( const int * N, double * A, const int * lda, const int * ipiv,
                       double * work, const int * lwork, int * info);
void F77_NAME(dgetrs)( const char * trans, const int * N, const int * NRHS, const double * A, const int * lda, 
                       const int * ipiv, double * B, const int * ldb, int * info);

void F77_NAME(dsyr)(   const char * uplo, const int * N, const double * alpha, const double * x,
                       const int * incx, double * A, const int * lda);

// Non-LAPACK routine for non-negative least squares
void F77_NAME(dnnls)(  double *A, const int * MDA, const int * M, const int* N, double * B, double * X, 
                       double * RNORM, double * W, double * ZZ, int * INDEX, int * MODE);


// Generic definitions
#ifdef USEFLOAT
    #define potrf   F77_NAME(spotrf)
    #define potri   F77_NAME(spotri)
    #define trtri   F77_NAME(strtri)
    #define trmm    F77_NAME(strmm)
    #define posv    F77_NAME(sposv)
    #define gelss   F77_NAME(sgelss)
    #define gemm    F77_NAME(sgemm)
    #define gemv    F77_NAME(sgemv)
    #define getrf   F77_NAME(sgetrf)
    #define getri   F77_NAME(sgetri)
    #define getrs   F77_NAME(sgetrs)
    #define syr     F77_NAME(ssyr)
    #define nnls    F77_NAME(snnls)
#else
    #define potrf   F77_NAME(dpotrf)
    #define potri   F77_NAME(dpotri)
    #define trtri   F77_NAME(dtrtri)
    #define trmm    F77_NAME(dtrmm)
    #define posv    F77_NAME(dposv)
    #define gelss   F77_NAME(dgelss)
    #define gemm    F77_NAME(dgemm)
    #define gemv    F77_NAME(dgemv)
    #define getrf   F77_NAME(dgetrf)
    #define getri   F77_NAME(dgetri)
    #define getrs   F77_NAME(dgetrs)
    #define syr     F77_NAME(dsyr)
    #define nnls    F77_NAME(dnnls)
#endif



#endif /* LAPACK_H_ */
