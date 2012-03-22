/**
 * \file matrix.c
 * Matrix Class.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
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

/* Standard copyright header */
/* Description
 * Include statement about row/col major format
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <err.h>
#include <string.h>
#include <math.h>
#include <ctype.h>          // for isblank
#include "matrix.h"
#include "lapack.h"


#define WARN_MEM(A) warn("Failed to allocation memory for %s at %s:%d.\n",(A),__FILE__,__LINE__)

/* constants */
/* None      */

/* members */
/* none */


/* private functions */
/* undetermined */


/* public functions */
/* undetermined */


/* First deal with allocation and deallocation of matrices  */
/* Allocate memory for matrix of a specified size */
MAT new_MAT( const int nrow, const int ncol ){
    /* default is real type */
    MAT mat = new_MAT_int(nrow, ncol, false);
    return mat;
}

/**
 * Create a new matrix structure and allocate memory for specified size.
 * Create integer array if requested else real.
 */
MAT new_MAT_int( const int nrow, const int ncol, const bool useint ) {
    MAT mat = malloc(sizeof(*mat));
    if ( NULL==mat){
        WARN_MEM("matrix");
        return NULL;
    }
    mat->ncol=ncol;
    mat->nrow=nrow;
    mat->useint = useint;
    /* Number of rows or columns might be zero but probably an error, so warn
     * as such. Want to avoid malloc(0) since this is "implementation defined"
     * in the C standard, may be a real address that should not be used or NULL
     */
     if ( 0==ncol || 0==nrow ){
         warn("One of dimensions of matrix equal to zero. Setting memory to NULL.\n");
         mat->x = NULL;
         mat->xint = NULL;
     }
     else {
         /* Usual case. Use calloc rather than malloc so elements are initialised */
         bool memfail = false;
         if (useint) {
             mat->xint = calloc(nrow*ncol,sizeof(int_t));
             memfail = (NULL==mat->xint);
         }
         else {
             mat->x = calloc(nrow*ncol,sizeof(real_t));
             memfail = (NULL==mat->x);
         }

         if (memfail) {
             WARN_MEM("matrix elements");
             xfree(mat);
             mat = NULL;
         }
     }

     return mat;
}

/* Free memory allocated for matrix */
MAT free_MAT ( MAT mat ){
    if(NULL==mat){ return NULL; }
    /* free real or int array */
    if (mat->useint) {
        if ( NULL!=mat->xint){
            xfree(mat->xint);
        }
    }
    else {
        /* Memory for elements may be NULL if nrow or ncol equals zero */
        if ( NULL!=mat->x){
            xfree(mat->x);
        }
    }
    xfree(mat);
    return NULL;
}

MAT copy_MAT( const MAT mat){
    if(NULL==mat){ return NULL;}

    MAT newmat = new_MAT_int(mat->nrow, mat->ncol, mat->useint);
    if(NULL==newmat){ return NULL;}

    /* copy real or int array */
    if (mat->useint) {
        memcpy(newmat->xint,mat->xint,mat->nrow*mat->ncol*sizeof(int_t));
    }
    else {
        memcpy(newmat->x,mat->x,mat->nrow*mat->ncol*sizeof(real_t));
    }
    return newmat;
}

void show_MAT ( XFILE * fp, const MAT mat, const uint_fast32_t mrow, const uint_fast32_t mcol){
    /* original default was with rownum */
    show_MAT_rownum( fp, mat, mrow, mcol, true);
}

void show_MAT_rownum( XFILE * fp, const MAT mat, const uint_fast32_t mrow, const uint_fast32_t mcol, bool rownum) {
    if(NULL==fp){ return;}
    if(NULL==mat){ return;}

#ifdef NDEBUG
    char fmt[] = " %#8.2f";
#else
    char fmt[] = " %#12.6f";
#endif
    const uint_fast32_t nrow = mat->nrow;
    const uint_fast32_t ncol = mat->ncol;
#ifdef NDEBUG
    const uint_fast32_t maxrow = (mrow!=0 && mrow<nrow)?mrow:nrow;
    const uint_fast32_t maxcol = (mcol!=0 && mcol<ncol)?mcol:ncol;
#else
    const uint_fast32_t maxrow = nrow;
    const uint_fast32_t maxcol = ncol;
#endif
    for( int row=0 ; row<maxrow ; row++){
        if (rownum) {
            xfprintf(fp,"%d:",row+1);
        }
        for ( int col=0 ; col<maxcol ; col++){
            if (mat->useint) {
                xfprintf(fp, INT_FORMAT, mat->xint[col * nrow + row]);
            }
            else {
                xfprintf(fp, fmt, mat->x[col*nrow+row]);
            }
        }
        if(maxcol<ncol){ xfprintf(fp,"\t... (%u others)",ncol-maxcol); }
        xfputc('\n',fp);
    }
    if( maxrow<nrow){ xfprintf(fp,"... (%u others)\n",nrow-maxrow); }
}

/**
 * Create a new real value matrix from a supplied real array.
 */
MAT new_MAT_from_array( const uint_fast32_t nrow, const uint_fast32_t ncol, const real_t * x){
    if(NULL==x){ return NULL;}
    MAT mat = new_MAT(nrow,ncol);
    if(NULL==mat){return NULL;}
    memcpy(mat->x,x,nrow*ncol*sizeof(real_t));
    return mat;
}

/**
 * Create a matrix from a supplied real array.
 * Resultant matrix references supplied array values directly so only free top structure.
 */
MAT coerce_MAT_from_array(const uint_fast32_t nrow, const uint_fast32_t ncol, real_t * x){
    assert(NULL!=x);
    MAT mat = malloc(sizeof(*mat));
    if(NULL==mat){ return NULL; }
    mat->nrow = nrow;
    mat->ncol = ncol;
    mat->useint = false;
    mat->x = x;
    return mat;
}

/**
 * Create a matrix from a supplied integer array.
 * Resultant matrix references supplied array values directly so only free top structure.
 */
MAT coerce_MAT_from_intarray(const uint_fast32_t nrow, const uint_fast32_t ncol, int_t * x){
    assert(NULL!=x);
    MAT mat = malloc(sizeof(*mat));
    if(NULL==mat){ return NULL; }
    mat->nrow = nrow;
    mat->ncol = ncol;
    mat->useint = true;
    mat->xint = x;
    return mat;
}

/** Create a new real value identity matrix of the specified size. */
MAT identity_MAT( const int nrow){
    MAT mat = new_MAT(nrow,nrow);
    validate(NULL!=mat,NULL);
    for ( int i=0 ; i<nrow ; i++){
        mat->x[i*nrow+i] = 1.0;
    }
    return mat;
}

/** Copy values from a real value matrix into another of the same size. */
MAT copyinto_MAT( MAT newmat, const MAT mat){
    if(NULL==newmat || NULL==mat){ return NULL;}
    if(newmat->nrow!=mat->nrow){ return NULL;}
    if(newmat->ncol!=mat->ncol){ return NULL;}
    memcpy(newmat->x,mat->x,mat->nrow*mat->ncol*sizeof(real_t));
    return newmat;
}

/**
 * Append matrix columns, from colstart to colend inclusive, from matin onto matout.
 * matout may be empty, in which case it is created with same number of rows as matin.
 * Returns matout unchanged if:
 * - matin is null
 * - colstart is greater than colend (warning issued)
 * - matin and matout do not have the same number of rows (warning issued)
 * - matin and matout do no have the same value type (warning issued)
 * - matout is not empty and extend fails
 *
 * Returns NULL if:
 *  - matout is empty and create fails
 *
 * Parameter validation:
 * - colstart is coerced to zero if negative and a warning issued.
 * - colend is coerced to last column if greater and a warning issued.
 */
MAT append_columns(MAT matout, const MAT matin, int colstart, int colend) {
    int colout = 0;

    /* validate parameters */
    validate(NULL != matin, matout);
    int nrow = matin->nrow;
    int ncol = matin->ncol;

    if (colstart < 0) {
        warnx("Matrix append_columns: column start increased from %d to 0.", colstart);
        colstart = 0;
    }
    if (colstart > colend) {
        warnx("Matrix append_columns: column start (%d) greater than column end (%d).", colstart, colend);
        return matout;
    }
    if (colend >= ncol) {
        warnx("Matrix append_columns: column end reduced from %d to maximum (%d).", colend, ncol - 1);
        colend = ncol - 1;
    }

    if (matout != NULL) {
        if (matout->nrow != nrow) {
            warnx("Matrix append_columns: different number of rows in source (%d) and target (%d).", nrow, matout->nrow );
            return matout;
        }
        if (matout->useint != matin->useint) {
            warnx("Matrix append_columns: different value types in source and target.");
            return matout;
        }
        /* target column to start append */
        colout = matout->ncol;
    }

    /* calculate size of appended matrix and create or extend */
    int newcol = colout + colend - colstart + 1;
    if (matout == NULL) {
        matout = new_MAT_int(nrow, newcol, matin->useint);
        validate(NULL != matout, NULL);
    }
    else {
        if (matout->useint) {
            int_t * tmp = realloc(matout->xint, nrow * newcol * sizeof(int_t));
            if (NULL == tmp) {
                WARN_MEM("matrix append");
                return matout;
            }
            matout->xint = tmp;
        }
        else {
            real_t * tmp = realloc(matout->x, nrow * newcol * sizeof(real_t));
            if (NULL == tmp) {
                WARN_MEM("matrix append");
                return matout;
            }
            matout->x = tmp;
        }
        matout->ncol = newcol;
    }

    /* copy selected columns to append */
    if (matout->useint) {
        memcpy(matout->xint + colout * nrow, matin->xint + colstart * nrow,
               nrow * (colend - colstart + 1) * sizeof(int_t));
    }
    else {
        memcpy(matout->x + colout * nrow, matin->x + colstart * nrow,
               nrow * (colend - colstart + 1) * sizeof(real_t));
    }

    return matout;
}

/** Set all elements in a supplied real value matrix to the specified value. */
MAT set_MAT( MAT mat, const real_t x){
    if(NULL==mat){ return NULL;}
    const uint_fast32_t nelt = mat->nrow * mat->ncol;
    for ( uint_fast32_t i=0 ; i<nelt ; i++){
        mat->x[i] = x;
    }
    return mat;
}

/** Ensure a value is within integer value range. */
int_t clipint(long int val) {
    if      (val < INT_MIN) { val = INT_MIN; }
    else if (val > INT_MAX) { val = INT_MAX; }

    return (int_t)val;
}

/** Count the number of sets of columns in a tab separated line of char (standard illumina format). */
int count_line_columns(const int nrow, char *ptr) {

    if (ptr == NULL) {return 0;}
    if (nrow <= 0) {return 0;}

    int nc = 0;
    bool found = true;
    /* read until run out */
    while (found) {
        if(ptr[0] != '\t'){
            found = false;
        }
        else {
            int cnt = 0;
            for (int nr = 0; nr < nrow; nr++) {
                if (ptr[0] == 0) {
                    found = false;
                    break;
                }
                else {
                    strtor(ptr, &ptr);
                    cnt++;
                }
            }
            /* only count complete columns */
            if (cnt == nrow) {++nc;}
        }
    }

    return nc;
}

/**
 * Create a new integer value matrix from tab separated sets of columns in a line of char.
 * Return actual number of columns found as reference parameter if not enough.
 */
MAT new_MAT_from_line(const int nrow, int *ncol, char *ptr){

    if (ptr == NULL) {return NULL;}
    MAT mat = new_MAT_int(nrow, *ncol, true);
    if(mat == NULL) {return NULL;}

    int nc = -1;
    bool found = true;
    /* read until number required or run out */
    while (found && (++nc < *ncol)) {

        /* should start with a tab */
        if(ptr[0] != '\t'){
            found = false;
        }
        else {
            for (int nr = 0; nr < nrow; nr++) {
                if (ptr[0] == 0) {
                    found = false;
                    break;
                }
                else {
                    /*
                     * read real value, round to integer and ensure within range
                     * Probably safe to assume value is less than max longint
                     */
                    mat->xint[nc * nrow + nr] = clipint((long int)(roundr(strtor(ptr, &ptr))));
                }
            }
        }
    }

    if (!found) {
        /* resize the matrix to match number of columns found */
        mat = trim_MAT(mat, nrow, nc, true);
        *ncol = nc;
    }

    return mat;
}

/** Write a real value matrix to file in a single line as tab separated sets of columns (standard illumina format). */
void write_MAT_to_line (XFILE * fp, const MAT mat) {
    if (NULL == fp) {return;}
    if (NULL == mat) {return;}

    const int nrow = mat->nrow;
    const int ncol = mat->ncol;

    for (int col = 0; col < ncol; col++) {
        /* each set begins with a tab and then separated with space */
        xfprintf(fp, "\t%0.1f", mat->x[col * nrow]);
        for (int row = 1; row < nrow; row++) {
            xfprintf(fp, " %0.1f", mat->x[col * nrow + row]);
        }
    }
    /* end with a new line */
    xfprintf(fp, "\n");
}

/**
 * Get the next non-comment line from the supplied file.
 * Leading whitespace is ignored, but blank and empty lines are returned.
 * Also returns length in len.
 */
static char * getnextline(XFILE * fp, size_t * len ){

    static const char COMMENTCHAR = '#';
    if (xfisnull(fp)) {return NULL;}

    char * line = NULL;
    bool found = false;

    while (!found) {
        line = xfgetln(fp, len);
        if (line == NULL) {
            found = true;
        }
        else {
            /* skip whitespace */
            int i = -1;
            while ((++i < *len) && (isblank(line[i])))
                ;
            /* check for comment indicator at first non-whitespace char */
            if (line[i] == COMMENTCHAR) {
                xfree(line);
            }
            else {
                found = true;
            }
        }
    }

    return line;
}

/**
 * Create a new matrix from a list of columns in a file, one column per row.
 * First row contains number of rows and columns.
 * Possible errors: number of rows and columns not properly specified
 *                  not enough elements
 *                  if too many elements then remainder ignored
 */
MAT read_MAT_from_column_file(XFILE * fp) {

    if (xfisnull(fp)) {return NULL;}

    /* get first line from file */
    char * line = NULL;
    size_t len = 0;
    line = getnextline(fp, &len);
    if (line == NULL) {
        warnx("General matrix read error");
        return NULL;
    }

    /* first line is number of rows and columns */
    char * ptr = line;
    int nrow = 0, ncol = 0;
    nrow = strtol(ptr, &ptr, 0);
    ncol = strtol(ptr, &ptr, 0);
    xfree(line);
    line = NULL;
    if ((nrow <= 0) || (ncol <= 0)) {
        warnx("Matrix size incorrectly specified: read in as %d by %d", nrow, ncol);
        return NULL;
    }

    MAT mat = new_MAT(nrow, ncol);
    if(mat == NULL) {return NULL;}

    bool found = true;
    int nc = -1;
    real_t val;
    char * prevptr = NULL;
    while (found && (++nc < ncol)) {
        /* get each column line */
        line = getnextline(fp, &len);
        if (line == NULL) {
            warnx("General matrix read error");
            found = false;
        }
        else {
            ptr = line;

            /* extract values for this column and store */
            for (int nr = 0; nr < nrow; nr++) {
                if (ptr[0] == 0) {
                    found = false;
                    break;
                }
                else {
                    prevptr = ptr;
                    val = strtor(ptr, &ptr);
                    /* check for invalid value */
                    if ((val == 0.0) && (ptr == prevptr)) {
                        found = false;
                        break;
                    }
                    else {
                        mat->x[nc * nrow + nr] = val;
                    }
                }
            }
            xfree(line);
        }
    }

    if (!found) {
        /* failed to read enough elements */
        warnx("Insufficient matrix data or incorrect file format; "
                "expected %d column rows but found only %d", ncol, nc);
        mat = free_MAT(mat);
    }
    return mat;
}

/**
 * Write a matrix to a file as a list of columns, one column per row.
 * First row contains number of rows and columns.
 * Can be formatted in fixed field size to 2 dp or free format.
 */
void write_MAT_to_column_file(XFILE * fp, const MAT mat, bool freeformat) {
    if (NULL == fp) {return;}
    if (NULL == mat) {return;}

    char fmt[] = " %8.2f";
    if (freeformat) {
        strcpy(fmt, "%f ");
    }

    const int nrow = mat->nrow;
    const int ncol = mat->ncol;

    /* first line is number of rows and columns */
    xfprintf(fp, "%d %d\n", nrow, ncol);

    for (int col = 0; col < ncol; col++) {
        for (int row = 0; row < nrow; row++) {
            xfprintf(fp, fmt, mat->x[col * nrow + row]);
        }
        /* next column line */
        xfprintf(fp, "\n");
    }
}

bool is_square(const MAT mat){
    validate(NULL!=mat,false);
    if(mat->nrow!=mat->ncol){ return false;}
    return true;
}

/*  Vec transpose operation */
MAT vectranspose ( const MAT mat, const unsigned int p ){
    /* Simple checking of arguments */
    if ( NULL==mat){
        warn("Attempting to apply vec-transpose operation to a NULL matrix");
        return NULL;
    }
    if ( 0==mat->nrow || 0==mat->ncol ){
        warn("Vec-transpose of matrix with no columns or no rows is not completely defined\n");
    }
    if ( 0!=(mat->nrow % p) ){
        warn("Invalid application of vec-transpose(%u) to a matrix with %u rows.\n",p,mat->nrow);
        return NULL;
    }


    MAT vtmat = new_MAT(p*mat->ncol,mat->nrow/p);
    if ( NULL==vtmat){
        return NULL;
    }

    /*  Iterate through columns of matrix doing necessary rearrangement.
     * Note: assumed column-major format
     * Each column is split into imaginary subcolumns of length p,
     * which will form the submatrix corresponding to that column.
     * Refer to external documentation for definition of vec-tranpose
     * (give document reference here).
     */
    /* Routine is valid for matrices with zero columns or rows since it
     * it will never attempt to access elements in this case
     */
     for ( unsigned int col=0 ; col<mat->ncol ; col++){
         /* offset represents where in the "matrix stack" that the column
          * should be formed into.*/
         unsigned int offset = col*p;
         for ( unsigned int subcol=0 ; subcol<(mat->nrow/p) ; subcol++){
             for ( unsigned int i=0 ; i<p ; i++){
                 vtmat->x[ (subcol*vtmat->nrow) + offset + i ]
                        = mat->x[col*mat->nrow + subcol*p + i];
             }
         }
     }

     return vtmat;
}


/**
 * Make a matrix in lower triangular form into a
 * symmetric matrix by copying the lower triangle
 * into the upper, transposing
 *
 * @param mat	Matrix in lower triangular form
 * return 	mat, or NULL if mat is not a square matrix.
 */
MAT symmeteriseL2U( MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);
    const uint_fast32_t n = mat->ncol;
    for ( uint_fast32_t col=0 ; col<n ; col++){
        for ( uint_fast32_t row=col ; row<n ; row++){
            mat->x[row*n+col] = mat->x[col*n+row];
        }
    }
    return mat;
}

/** 
 * Cholesky decomposition of matrix. Done in place.
 */
MAT cholesky( MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);
    
    // LAPACK routine for Cholesky decomposition
    int info=0;
    int n = mat->nrow;
    potrf(LAPACK_UPPER,&mat->nrow,mat->x,&mat->nrow,&info);
    if(info!=0){ warnx("potrf in %s returned %d\n",__func__,info);}
    
    // Zero entries in lower diagonal (not referenced in potrf, so retain initial values)
    for ( int i=0 ; i<n ; i++){
        for ( int j=i+1 ; j<n ; j++){
           mat->x[i*n+j] = 0;
        }
    } 
    return mat;
}

// Invert cholesky factorisation
MAT invert_cholesky( MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);
    int info=0;
    trtri(LAPACK_LOWER,LAPACK_NONUNITTRI,&mat->nrow,mat->x,&mat->nrow,&info);
    symmeteriseL2U(mat);
    return mat;
}

MAT reshape_MAT( MAT mat, const int nrow){
    validate(NULL!=mat,NULL);
    validate(((mat->nrow*mat->ncol)%nrow)==0,NULL);
    mat->ncol = (mat->nrow*mat->ncol)/nrow;
    mat->nrow = nrow;
    return mat;
}

/**
 * Reduce the size of a matrix.
 * Adjust the elements array to match.
 */
MAT trim_MAT( MAT mat, const int mrow, const int mcol, const bool forwards){
    validate(NULL!=mat,NULL);
    validate(mrow>=0,NULL);
    validate(mcol>=0,NULL);
    validate(mrow<=mat->nrow,NULL);
    validate(mcol<=mat->ncol,NULL);
    if(forwards==false){ errx(EXIT_FAILURE,"Forwards==false not implemented in %s (%s:%d)\n",__func__,__FILE__,__LINE__);}
    for ( uint_fast32_t col=0 ; col<mcol ; col++){
        uint_fast32_t midx = col*mrow;
        uint_fast32_t nidx = col*mat->nrow;
        if (mat->useint) {
            memmove(mat->xint+midx,mat->xint+nidx,mrow*sizeof(int_t));
        }
        else {
            memmove(mat->x+midx,mat->x+nidx,mrow*sizeof(real_t));
        }
    }
    mat->nrow = mrow;
    mat->ncol = mcol;
    return mat;
}

MAT * block_diagonal_MAT( const MAT mat, const int n){
    validate(NULL!=mat,NULL);
    validate(mat->ncol==mat->nrow,NULL);    // Ensure symmetry
    const int nelts = mat->ncol / n;        // Number of blocks on diagonal
    validate((mat->ncol % n)==0,NULL);      // Is parameter valid?

    // Create memory
    MAT * mats = calloc(nelts,sizeof(*mats));
    if(NULL==mats){ goto cleanup; }
    for ( uint_fast32_t i=0 ; i<nelts ; i++){
        mats[i] = new_MAT(n,n);
        if(NULL==mats[i]){ goto cleanup;}
    }
    // Copy into diagonals
    for ( uint_fast32_t i=0 ; i<nelts ; i++){
        for ( uint_fast32_t col=0 ; col<n ; col++){
            const uint_fast32_t oldcol = i*n+col;
            for ( uint_fast32_t row=0 ; row<n ; row++){
                const uint_fast32_t oldrow = i*n+row;
                mats[i]->x[col*n+row] = mat->x[oldcol*mat->nrow+oldrow];
            }
        }
    }
    return mats;

cleanup:
    if(NULL!=mats){
        for ( uint_fast32_t i=0 ; i<nelts ; i++){
            free_MAT(mats[i]);
        }
    }
    xfree(mats);
    return NULL;
}

MAT scale_MAT(MAT mat, const real_t f){
    validate(NULL!=mat,NULL);
    const uint_fast32_t nelt = mat->ncol * mat->nrow;
    for ( uint_fast32_t elt=0 ; elt<nelt ; elt++){
            mat->x[elt] *= f;
    }
    return mat;
}

/** Return the transpose of a supplied square matrix, replacing the original. */
MAT transpose_inplace( MAT mat){
    validate(NULL!=mat,NULL);
    const uint_fast32_t ncol = mat->ncol;
    const uint_fast32_t nrow = mat->nrow;
    validate(ncol==nrow,NULL);

    for ( uint_fast32_t col=0 ; col<ncol ; col++){
        for ( uint_fast32_t row=0 ; row<col ; row++){
            real_t x = mat->x[col*nrow + row];
            mat->x[col*nrow + row] = mat->x[row*ncol + col];
            mat->x[row*ncol + col] = x;
        }
    }
    mat->nrow = ncol;
    mat->ncol = nrow;
    return mat;
}

/** Create a new matrix which is the transpose of a supplied matrix. */
MAT transpose( const MAT mat){
    validate(NULL!=mat,NULL);
    MAT tmat = new_MAT(mat->ncol,mat->nrow);
    validate(NULL!=tmat,NULL);
    for ( uint_fast32_t col=0 ; col<mat->ncol ; col++){
        for ( uint_fast32_t row=0 ; row<mat->nrow ; row++){
            tmat->x[row*mat->ncol+col] = mat->x[col*mat->nrow+row];
        }
    }
    return tmat;
}

/** Create a new matrix which is the inverse of a supplied square matrix. */
MAT invert(const MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);

    int INFO;
    int N;
    int * IPIV=NULL;
    int LWORK=0;
    real_t * WORK=NULL,WORKSIZE=0;

    N = mat->nrow;
    // Get temporary memory for inversion
    LWORK = -1;
    getri(&N,mat->x,&N,IPIV,&WORKSIZE,&LWORK,&INFO);
    LWORK=(int)WORKSIZE;
    WORK = malloc(LWORK*sizeof(real_t));
    IPIV = malloc(N*sizeof(int));

    MAT matinv = copy_MAT(mat);
    // LU decomposition required for inversion
    getrf(&N,&N,matinv->x,&N,IPIV,&INFO);
    // Invert
    getri(&N,matinv->x,&N,IPIV,WORK,&LWORK,&INFO);

    xfree(IPIV);
    xfree(WORK);

    return matinv;
}

/** Create a new matrix which is the inverse of a supplied symmetric positive definite matrix. */
MAT invert_symmetric(const MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);

    int info = 0;
    int N = mat->nrow;

    MAT matinv = copy_MAT(mat);
    validate(NULL!=matinv,NULL);

    // Cholesky factorisation required for inversion
    potrf(LAPACK_LOWER,&N,matinv->x,&N,&info);
    if (info!=0) {
        warnx("Supplied matrix is not positive definite, order: %d", info);
        free_MAT(matinv);
        return NULL;
    }

    // Invert
    potri(LAPACK_LOWER,&N,matinv->x,&N,&info);
    if (info!=0) {
        warnx("Cholesky diagonal is zero, element: %d", info);
        free_MAT(matinv);
        return NULL;
    }

    // Only lower triangle returned
    symmeteriseL2U(matinv);
    return matinv;
}

/**
 * Return a vector * matrix * vector product.
 * Assumes the supplied vectors are the appropriate size.
 */
real_t xMy( const real_t * x, const MAT M, const real_t * y){
    validate(NULL!=x,NAN);
    validate(NULL!=M,NAN);
    validate(NULL!=y,NAN);
    const uint_fast32_t ncol = M->ncol;
    const uint_fast32_t nrow = M->nrow;

    real_t res = 0.;
    for ( uint_fast32_t col=0 ; col<ncol ; col++){
        real_t rowtot = 0.;
        for ( uint_fast32_t row=0 ; row<nrow ; row++){
            rowtot += x[row] * M->x[col*nrow+row];
        }
        res += rowtot * y[col];
    }
    return res;
}

real_t normalise_MAT(MAT mat, const real_t delta_diag){
    validate(NULL!=mat,NAN);
    validate(mat->nrow==mat->ncol,NAN);
    const int n = mat->nrow;
    int piv[n];

    if(0.0!=delta_diag){
        for ( uint_fast32_t i=0 ; i<n ; i++){
            mat->x[i*n+i] += delta_diag;
        }
    }

    MAT mcopy = copy_MAT(mat);
    int info = 0;
    getrf(&n,&n,mcopy->x,&n,piv,&info);

    real_t logdet = 0.;
    for ( uint_fast32_t i=0 ; i<n ; i++){
        logdet += log(fabs(mcopy->x[i*n+i]));
    }
    free_MAT(mcopy);

    real_t f = 1e-5 + exp(logdet/n);
    scale_MAT(mat,1./f);
    return f;
}

/**
 * LU decomposition of matrix.
 * Routine returns structure rather than a matrix since pivoting information must also be transferred.
 */
struct structLU LUdecomposition( const MAT mat){
    // NULL structure to be returned on error
    struct structLU structLUnill = {NULL,NULL};
    validate(NULL!=mat,structLUnill);
    validate(mat->nrow==mat->ncol,structLUnill);
    const int n = mat->nrow;

    int * piv = calloc(n,sizeof(int));;
    MAT mcopy = copy_MAT(mat);
    if(NULL==piv || NULL==mcopy){ goto cleanup; }
   
    // Call LAPACK routine for LU decomposition
    int info = 0;
    getrf(&n,&n,mcopy->x,&n,piv,&info);
    if(info!=0){ warnx("getrf in %s returned %d\n",__func__,info);}

    return (struct structLU) { mcopy, piv};

cleanup:
    return structLUnill;
}


#ifdef TEST

/* Create a new integer value matrix as a copy of the one supplied (real or integer).*/
MAT copy_MAT_int(const MAT mat) {
    if(NULL==mat) { return NULL; }

    MAT newmat = new_MAT_int(mat->nrow, mat->ncol, true);
    if(NULL==newmat) { return NULL; }

    if (mat->useint) {
        /* can do a memory copy */
        memcpy(newmat->xint, mat->xint, mat->nrow * mat->ncol * sizeof(int_t));
    }
    else {
        /* copy values, ensuring within range */
        int nrow = mat->nrow;
        for (int col = 0; col < mat->ncol; col++) {
            for (int row = 0; row < nrow; row++) {
                newmat->xint[col * nrow + row] = clipint((long int)(roundr(mat->x[col * nrow + row])));
            }
        }
    }
    return newmat;
}

int main (int argc, char * argv[]){
    if(argc!=3){
        fputs("Usage: test appendto appendfrom (filenames)\n",stdout);
        return EXIT_FAILURE;
    }

    /* read inputs */
    XFILE * fp = xfopen(argv[1], XFILE_UNKNOWN, "r");
    MAT matout = read_MAT_from_column_file(fp);
    fp = xfclose (fp);
    fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    MAT matin = read_MAT_from_column_file(fp);
    fp = xfclose (fp);

    /* use copy of target file to preserve original */
    MAT outcopy;

    fputs("Append all:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, matin, 0, matin->ncol - 1);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append some:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, matin, 1, matin->ncol - 2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append to null:\n", stdout);
//    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, matin, 1, matin->ncol - 2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append from null:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, NULL, 1, matin->ncol - 2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append colstart > colend:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, matin, matin->ncol/2 + 1, matin->ncol/2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append colstart negative:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, matin, -1, matin->ncol - 2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append colend larger than source:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = append_columns(outcopy, matin, 1, matin->ncol);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append row mismatch:\n", stdout);
    outcopy = copy_MAT(matout);
    outcopy = trim_MAT(outcopy, outcopy->nrow - 1, outcopy->ncol - 1, true);
//    outcopy = new_MAT(matin->nrow + 1, matin->ncol);
    outcopy = append_columns(outcopy, matin, 1, matin->ncol - 2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    fputs("Append type mismatch:\n", stdout);
    outcopy = copy_MAT_int(matout);
    outcopy = append_columns(outcopy, matin, 1, matin->ncol - 2);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);

    matout->x[0] = (real_t)(INT16_MIN *2);
    matout->x[1] = (real_t)(INT16_MAX *2);
    fprintf(stdout, "Values outside integer range: %0.2f %0.2f; see clipped values below\n", matout->x[0], matout->x[1]);

    fputs("Append integer type:\n", stdout);
    outcopy = copy_MAT_int(matout);
    MAT outcopyin = copy_MAT_int(matin);
    outcopy = append_columns(outcopy, outcopyin, 0, outcopyin->ncol - 1);
    show_MAT(xstdout, outcopy, outcopy->nrow, outcopy->ncol);
    outcopy = free_MAT(outcopy);
    outcopyin = free_MAT(outcopyin);

    free_MAT(matin);
    free_MAT(matout);
}

#endif
