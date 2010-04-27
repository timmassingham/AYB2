/**
 * \file lapack.h
 * Description.
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

#ifndef LAPACK_H_
#define LAPACK_H_

#ifdef __linux
#define F77_NAME(A) A ## _
#else
#define F77_NAME(A) A
#endif

const char LAPACK_UPPER[] = {'U'};
const char LAPACK_LOWER[] = {'L'};
const char LAPACK_UNITTRI[] = {'U'};
const char LAPACK_NONUNITTRI[] = {'N'};

// Float functions
void F77_NAME(spotrf)( const char * uplo, const int * n, float * A, const int * lda, int * info );
void F77_NAME(strtri)( const char * uplo, const char * diag, const int * n, float * a, const int * lda, int * info);

void F77_NAME(sgetrf)(const int * M, const int * N, float * A, const int * lda,
                      int * ipiv, int * info);

void F77_NAME(sgetri)(const int * N, float * A, const int * lda, const int * ipiv,
                      float * work, const int * lwork, int * info);

// Double functions
void F77_NAME(dpotrf)( const char * uplo, const int * n, double * A, const int * lda, int * info );
void F77_NAME(dtrtri)( const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);

void F77_NAME(dgetrf)(const int * M, const int * N, double * A, const int * lda,
                      int * ipiv, int * info);

void F77_NAME(dgetri)(const int * N, double * A, const int * lda, const int * ipiv,
                      double * work, const int * lwork, int * info);



// Generic definitions
#ifdef USEFLOAT
    #define potrf   F77_NAME(spotrf)
    #define trtri   F77_NAME(strtri)
    #define getrf   F77_NAME(sgetrf)
    #define getri   F77_NAME(sgetri)
#else
    #define potrf   F77_NAME(dpotrf)
    #define trtri   F77_NAME(dtrtri)
    #define getrf   F77_NAME(dgetrf)
    #define getri   F77_NAME(dgetri)
#endif

#endif /* LAPACK_H_ */
