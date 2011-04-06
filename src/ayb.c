/** 
 * \file ayb.c
 * Ayb Class.
 * Handles the data required for modelling and top level modelling.
 * Per run data stored as local global members but AYB data passed in as pointer to structure.
 *//* 
 *  Created : 16 Mar 2011
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
#include "ayb.h"
#include "call_bases.h"
#include "cif.h"
#include "cluster.h"
#include "dirio.h"
#include "intensities.h"
#include "lambda.h"
#include "message.h"
#include "mpn.h"
#include "nuc.h"
#include "qual_table.h"
#include "statistics.h"


/** AYB structure contains the data required for modelling. */
struct AybT {
    uint32_t ncluster;
    uint32_t ncycle;
    TILE tile;
    ARRAY(NUC) bases;
    ARRAY(PHREDCHAR) quals;
    MAT M, P, N;
    MAT lambda;
    MAT we, cycle_var;
};

/** Structure for output of final working values. */
struct WorkT {
    XFILE *fp;
    CIFDATA cif;
};
typedef struct WorkT * WORKPTR;


/* constants */

static const unsigned int AYB_NITER = 20;       ///< Number of parameter estimation loops.
static const real_t DELTA_DIAG = 1.0;           ///< Delta for solver routines.

/** Initial Crosstalk matrix if not read in, fixed values of approximately the right shape. */
static const real_t INITIAL_CROSSTALK[] = {
    2.0114300, 1.7217841, 0.06436576, 0.1126401,
    0.6919319, 1.8022413, 0.06436576, 0.0804572,
    0.2735545, 0.2252802, 1.39995531, 0.9976693,
    0.2896459, 0.2413716, 0.11264008, 1.3194981
};

/** Name text for matrix messages. Match to IOTYPE enum in dirio. */
static const char *MATRIX_TEXT[] = {"Crosstalk", "Noise", "Phasing"};
/** Possible solver routine text. Match to SOLVERS list. Used to match --solver argument. */
static const char *SOLVER_TEXT[] = {"ls", "zero", "nnls"};

/* members */

static real_t BasePenalty[NBASE] = {0.0,0.0,0.0,0.0}; ///< Penalty for least-squares base-calling.
static MAT Matrix[E_MNP];                       ///< Predetermined matrices.
static bool ShowWorking = false;                ///< Set to output final working values.

/** Possible solver routines and number. */
typedef enum SolvIdT {E_LS, E_ZERO, E_NNLS, E_SOLVER_NUM} SOLVID;
static SOLVID SolverIndex = E_ZERO;             ///< Selected solver routine index.
/** Template for P solver routine. */
typedef int (*SOLVER)(MAT , MAT, real_t *, const real_t);
/** Possible solver routines. */
SOLVER SOLVERS[] = { solverSVD, solverZeroSVD, solverNNLS };
static SOLVER SolverRoutine = solverZeroSVD;    ///< Selected solver routine, match to index.


/* private functions */

/**
 * Accumulate full covariance matrix (ncycle * NBASE x ncycle * NBASE).
 * \n Note: If V is NULL, the required memory is allocated.
 * - p:        Matrix of processed intensities
 * - lambda:   Brightness of cluster
 * - base:     Current base call
 * - V:        Array of covariance matrices of size 1; single element used for accumulation
 */
static MAT * accumulate_all_covariance( const real_t we, const MAT p, const real_t lambda, const NUC * base, MAT * V) {
    validate(NULL!=p, NULL);
    validate(NBASE==p->nrow, NULL);
    validate(lambda>=0.0, NULL);

    const uint32_t ncycle = p->ncol;
    const uint32_t ncyclebase = ncycle * NBASE;

    // Allocate memory for V if necessary
    // for compatibility of return value with normal accumulate use an array of size 1
    if ( NULL==V ){
        V = calloc(1, sizeof(*V));
        if(NULL==V) {return NULL;}
        V[0] = new_MAT(ncyclebase, ncyclebase);
        if(NULL==V[0]) {
            xfree(V);
            return NULL;
        }
    }

    // Perform accululation. V += R R^t
    // Note: R = P - \lambda I_b, where I_b is unit vector with b'th elt = 1
    //  so R R^t = P P^t - \lambda I_b P^t
    //                   - \lambda P I_b^t + \lambda^2 I_b I_b^2
    // P P^t
    for (uint32_t i = 0; i < ncyclebase; i++) {
        for (uint32_t j = 0; j < ncyclebase; j++) {
            (V[0])->x[j * ncyclebase + i] += we * p->x[i] * p->x[j];
        }
    }

    // \lambda I_b P^t and \lambda P I_b^t
    for (uint32_t cy1 = 0; cy1 < ncycle; cy1++){

        if (isambig(base[cy1])) { continue;}

        uint32_t ibase1 = cy1 * NBASE + base[cy1];
        for (uint32_t j = 0; j < ncyclebase; j++) {
            (V[0])->x[j * ncyclebase + ibase1] -= we * lambda * p->x[j];
            (V[0])->x[ibase1 * ncyclebase + j] -= we * lambda * p->x[j];
        }

        // \lambda^2 I_b I_b^t
        for (uint32_t cy2 = 0; cy2 < ncycle; cy2++) {

            if (isambig(base[cy2])) { continue;}

            uint32_t ibase2 = cy2 * NBASE + base[cy2];
            (V[0])->x[ibase2 * ncyclebase + ibase1] += we * lambda * lambda;
        }
    }

    return V;
}

/**
 * Accumulate required variances (inner summation of variance calculation).
 * \n Note: If V is NULL, the required memory is allocated.
 * - p:        Matrix of processed intensities
 * - lambda:   Brightness of cluster
 * - base:     Current base call
 * - V:        Array of covariance matrices into which accumulation occurs
 */
static MAT * accumulate_covariance( const real_t we, const MAT p, const real_t lambda, const NUC * base, MAT * V){
    validate(NULL!=p,NULL);
    validate(NBASE==p->nrow,NULL);
    validate(lambda>=0.,NULL);

    const uint32_t ncycle = p->ncol;
    // Allocate memory for V if necessary
    if ( NULL==V ){
        V = calloc(ncycle+1,sizeof(*V));
        if(NULL==V){ return NULL;}
        for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
            V[cycle] = new_MAT(NBASE,NBASE);
            if(NULL==V[cycle]){ goto cleanup;}
        }
    }
    // Perform accululation. V += R R^t
    // Note: R = P - \lambda I_b, where I_b is unit vector with b'th elt = 1
    //  so R R^t = P P^t - \lambda I_b P^t
    //                   - \lambda P I_b^t + \lambda^2 I_b I_b^2
    for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        const int cybase = base[cycle];
        // P P^t
        for ( uint32_t i=0 ; i<NBASE ; i++){
            for ( uint32_t j=0 ; j<NBASE ; j++){
                V[cycle]->x[i*NBASE+j] += we * p->x[cycle*NBASE+i] * p->x[cycle*NBASE+j];
            }
        }

        if (isambig(cybase)) { continue;}

        // \lambda I_b P^t and \lambda P I_b^t
        for ( uint32_t i=0 ; i<NBASE ; i++){
            V[cycle]->x[cybase*NBASE+i] -= we * lambda * p->x[cycle*NBASE+i];
            V[cycle]->x[i*NBASE+cybase] -= we * lambda * p->x[cycle*NBASE+i];
        }
        // \lambda^2 I_b I_b^t
        V[cycle]->x[cybase*NBASE+cybase] += we * lambda*lambda;
    }
    // Create residuals
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        if (!isambig(base[cy])) { p->x[cy*NBASE+base[cy]] -= lambda;}
    }

    return V;

cleanup:
    for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        free_MAT(V[cycle]);
    }
    xfree(V);
    return NULL;
}

/**
 * Initialise crosstalk (M), Noise (N) and Phasing (P) matrices.
 * May use read in values or initialise using an internal method.
 */
static MAT init_matrix(MAT mat, const IOTYPE idx) {

    if (mat == NULL) {return NULL;}

    if (Matrix[idx] == NULL) {
        /* initialise internally */
        unsigned int nrow = mat->nrow;
        switch (idx) {
            case E_CROSSTALK:
                /* shouldn't reach this as crosstalk always set up in read_matrices */
                /* initial crosstalk from default array */
                free_MAT(mat);
                mat = new_MAT_from_array(NBASE, NBASE, INITIAL_CROSSTALK);
                break;

            case E_NOISE:
                /* initial noise is zero */
                set_MAT(mat, 0.0);
                break;

            case E_PHASING:
                /* initial phasing */
                #warning "Phasing not properly initialised yet"
                free_MAT(mat);
                mat = identity_MAT(nrow);
                break;

            default: ;
        }
    }

    else {
        /* use read in matrix */
        unsigned int ncol = mat->ncol;
        mat = copyinto_MAT(mat, Matrix[idx]);
        if (mat == NULL) {
            message (E_MATRIXINIT_SDD, MSG_ERR, MATRIX_TEXT[idx], ncol, Matrix[idx]->ncol);
        }
    }
    return mat;
}

/** Returns true if array of data is missing (all zero). */
static inline bool nodata(const real_t * sig, const int num) {

    for (unsigned int i = 0; i < num; i++) {
        if (sig[i] != 0) {return false;}
    }
    return true;
}

/**
 * Read in any external crosstalk (M), Noise (N) and Phasing (P) matrices.
 * Returns false if failed to read a supplied matrix file.
 */
static bool read_matrices(void) {
    XFILE *fpmat = NULL;
    bool found = true;

    for (IOTYPE idx = (IOTYPE)0; idx < E_MNP; idx++) {
        if (matrix_from_file(idx)) {
            fpmat = open_matrix(idx);
            if (fpmat == NULL) {
                found = false;
            }
            else {
                /* read in matrix */
                Matrix[idx] = read_MAT_from_column_file(fpmat);
                if (Matrix[idx] == NULL) {
                    found = false;
                }
            }
            xfclose(fpmat);
        }
        else {
            /* no input file specified, initialise using internal method later */
            Matrix[idx] = NULL;
        }

        /* exit if this iteration failed */
        if (!found) {
            break;
        }
    }

    return found;
}

/** Calculate new weights. */
static real_t update_cluster_weights(AYB ayb){
    validate(NULL!=ayb,NAN);
    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle   = ayb->ncycle;
    real_t sumLSS = 0.;

    MAT e = NULL;
    /*  Calculate least squares error, using ayb->we as temporary storage */
    unsigned int cl = 0;
    LIST(CLUSTER) node = ayb->tile->clusterlist;
    while (NULL != node && cl < ncluster){
        ayb->we->x[cl] = 0.;
        MAT cycle_ints = node->elt->signals;
        NUC * cycle_bases = ayb->bases.elt + cl*ncycle;
        e = expected_intensities(ayb->lambda->x[cl],cycle_bases,ayb->M,ayb->P,ayb->N,e);
        for( uint32_t idx=0 ; idx<NBASE*ncycle ; idx++){
            real_t tmp = cycle_ints->x[idx] - e->x[idx];
            ayb->we->x[cl] += tmp*tmp;
        }
        sumLSS += ayb->we->x[cl];

        /* next cluster */
        node = node->nxt;
        cl++;
    }
    free_MAT(e);

    /* Calculate weight for each cluster */
    real_t meanLSSi = mean(ayb->we->x,ncluster);
    real_t varLSSi = variance(ayb->we->x,ncluster);
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        const real_t d = ayb->we->x[cl]-meanLSSi;
        ayb->we->x[cl] = cauchy(d*d,varLSSi);
    }
    //xfputs("Cluster weights:\n",xstderr);
    //show_MAT(xstderr,ayb->we,8,1);
    return sumLSS;
}

/* Functions for final processed intensities output */

/**
 * Initialise for final processed intensities output.
 * Open output file to be in same format as intensities input file.
 * If cif then also create cif structure.
 */
static WORKPTR open_processed(const AYB ayb, const int blk) {

    WORKPTR work = calloc(1, sizeof(*work));
    work->fp = open_output_blk("pif", blk);

    switch (get_input_format()) {
        case E_TXT:
            /* nothing to do */
            break;

        case E_CIF:
            work->cif = create_cif(ayb->ncycle, ayb->ncluster);
            break;

        default: ;
    }

    return work;
}

/** Output/store a line of processed intensities in same format as intensities input file. */
static void write_processed(WORKPTR work, const AYB ayb, const CLUSTER cluster, const uint32_t cl, MAT pcl_int) {

    const uint32_t ncycle = ayb->ncycle;

    switch (get_input_format()) {
        case E_TXT:
            if (!xfisnull(work->fp)) {
                write_lane_tile (work->fp, ayb->tile);
                write_coordinates (work->fp, cluster);
                write_MAT_to_line(work->fp, pcl_int);
            }
            break;

        case E_CIF:
            for (uint32_t cy = 0; cy < ncycle; cy++) {
                for (uint32_t base = 0; base < NBASE; base++){
                    cif_set_from_real (work->cif, cl, base, cy, pcl_int->x[cy * NBASE + base]);
                }
            }
            break;

        default: ;
    }
}

/** Output final model values. */
static void output_final(const AYB ayb, const int blk) {
    XFILE *fpfin = NULL;

    /* final model values */
    fpfin = open_output_blk("final", blk);
    if (!xfisnull(fpfin)) {
        show_AYB(fpfin, ayb, false);
    }
    xfclose(fpfin);

    /* final M, N, P in input format */
    fpfin = open_output_blk("M", blk);
    if (!xfisnull(fpfin)) {
        write_MAT_to_column_file (fpfin, ayb->M, false);
    }
    xfclose(fpfin);

    fpfin = open_output_blk("N", blk);
    if (!xfisnull(fpfin)) {
        write_MAT_to_column_file (fpfin, ayb->N, false);
    }
    xfclose(fpfin);

    fpfin = open_output_blk("P", blk);
    if (!xfisnull(fpfin)) {
        write_MAT_to_column_file (fpfin, ayb->P, false);
    }
    xfclose(fpfin);
}

/**
 * Finish off and close up final processed intensities output.
 * Output is done here for cif and final model.
 * If status is error then just free resources.
 */
static WORKPTR close_processed(WORKPTR work, const AYB ayb, const int blk, const int status) {

    if (status != DATA_ERR) {
        switch (get_input_format()) {
            case E_TXT:
                /* nothing to do */
                break;

            case E_CIF:
                if (!xfisnull(work->fp)) {
                    writeCIFtoStream(work->cif, work->fp);
                }
                break;

            default: ;
        }

        /* final model values */
        output_final(ayb, blk);
    }

    free_cif(work->cif);
    work->fp = xfclose(work->fp);
    return work;
}


/* public functions */

/* standard functions */

AYB new_AYB(const uint32_t ncycle, const uint32_t ncluster){
    AYB ayb = malloc(sizeof(*ayb));
    if(NULL==ayb){ return NULL;}
    ayb->ncycle = ncycle;
    ayb->ncluster = ncluster;
    ayb->tile = new_TILE();
    ayb->bases = new_ARRAY(NUC)(ncluster*ncycle);
    ayb->quals = new_ARRAY(PHREDCHAR)(ncluster*ncycle);
    ayb->M = new_MAT(NBASE,NBASE);
    ayb->P = new_MAT(ncycle,ncycle);
    ayb->N = new_MAT(NBASE,ncycle);
    ayb->lambda = new_MAT(ncluster,1);
    ayb->we = new_MAT(ncluster,1);
    ayb->cycle_var = new_MAT(ncycle,1);
    if( NULL==ayb->tile || NULL==ayb->bases.elt || NULL==ayb->quals.elt
            || NULL==ayb->M || NULL==ayb->P || NULL==ayb->N
            || NULL==ayb->lambda || NULL==ayb->we || NULL==ayb->cycle_var){
        goto cleanup;
    }

    return ayb;

cleanup:
    free_AYB(ayb);
    return NULL;
}

AYB free_AYB(AYB ayb){
    if(NULL==ayb){ return NULL;}
    free_TILE(ayb->tile);
    free_ARRAY(NUC)(ayb->bases);
    free_ARRAY(PHREDCHAR)(ayb->quals);
    free_MAT(ayb->M);
    free_MAT(ayb->P);
    free_MAT(ayb->N);
    free_MAT(ayb->we);
    free_MAT(ayb->cycle_var);
    free_MAT(ayb->lambda);
    xfree(ayb);
    return NULL;
}

AYB copy_AYB(const AYB ayb){
    if(NULL==ayb){return NULL;}
    AYB ayb_copy = malloc(sizeof(*ayb));
    if(NULL==ayb_copy){ return NULL;}

    ayb_copy->ncycle = ayb->ncycle;
    ayb_copy->ncluster = ayb->ncluster;

    ayb_copy->tile = copy_TILE(ayb->tile);
    if(NULL!=ayb->tile && NULL==ayb_copy->tile){ goto cleanup;}

    ayb_copy->bases = copy_ARRAY(NUC)(ayb->bases);
    if(NULL!=ayb->bases.elt && NULL==ayb_copy->bases.elt){ goto cleanup;}

    ayb_copy->quals = copy_ARRAY(PHREDCHAR)(ayb->quals);
    if(NULL!=ayb->quals.elt && NULL==ayb_copy->quals.elt){ goto cleanup;}

    ayb_copy->M = copy_MAT(ayb->M);
    if(NULL!=ayb->M && NULL==ayb_copy->M){ goto cleanup;}

    ayb_copy->P = copy_MAT(ayb->P);
    if(NULL!=ayb->P && NULL==ayb_copy->P){ goto cleanup;}

    ayb_copy->N = copy_MAT(ayb->N);
    if(NULL!=ayb->N && NULL==ayb_copy->N){ goto cleanup;}

    ayb_copy->we = copy_MAT(ayb->we);
    if(NULL!=ayb->we && NULL==ayb_copy->we){ goto cleanup;}

    ayb_copy->cycle_var = copy_MAT(ayb->cycle_var);
    if(NULL!=ayb->cycle_var && NULL==ayb_copy->cycle_var){ goto cleanup;}

    ayb_copy->lambda = copy_MAT(ayb->lambda);
    if(NULL!=ayb->lambda && NULL==ayb_copy->lambda){ goto cleanup;}

    return ayb_copy;

cleanup:
    free_AYB(ayb_copy);
    return NULL;
}

void show_AYB(XFILE * fp, const AYB ayb, bool showall){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    xfprintf(fp,"%u cycles from %u clusters\n",ayb->ncycle,ayb->ncluster);
    xfputs("M:\n",fp); show_MAT(fp,ayb->M,NBASE,NBASE);
    xfputs("P:\n",fp); show_MAT(fp,ayb->P,ayb->ncycle,ayb->ncycle);
    xfputs("N:\n",fp); show_MAT(fp,ayb->N,NBASE,ayb->ncycle);
    xfputs("we:\n",fp); show_MAT(fp,ayb->we,ayb->ncluster,1);
    xfputs("cycle_var:\n",fp); show_MAT(fp,ayb->cycle_var,ayb->ncycle,1);
    xfputs("lambda:\n",fp); show_MAT(fp,ayb->lambda,ayb->ncluster,1);
    if (showall) {
        xfputs("Bases:\n",fp); show_ARRAY(NUC)(fp,ayb->bases,"",ayb->ncycle);
        xfputc('\n',fp);
        xfputs("Quality:\n",fp); show_ARRAY(PHREDCHAR)(fp,ayb->quals,"",ayb->ncycle*10);
        xfputc('\n',fp);
#ifdef NDEBUG
        xfputs("Tile:\n",fp); show_TILE(fp,ayb->tile,10);
#else
        xfputs("Tile:\n",fp); show_TILE(fp,ayb->tile,ayb->ncluster);
#endif
    }
    xfputc('\n',fp);
}

/* access functions */

/** Return array of non-zero lambdas and how many in num. */
real_t * get_AYB_lambdas(const AYB ayb, uint32_t *num) {

    real_t * lambdas = calloc(ayb->lambda->nrow, sizeof(real_t));
    *num = 0;

    for (uint32_t i = 0; i < ayb->lambda->nrow; i++) {
        if (ayb->lambda->x[i] != 0.0) {
            lambdas[(*num)++] = ayb->lambda->x[i];
        }
    }
    return lambdas;
}

/** Return the number of clusters. */
uint32_t get_AYB_ncluster(const AYB ayb) {
    return ayb->ncluster;
}

/** Return the number of cycles. */
uint32_t get_AYB_ncycle(const AYB ayb) {
    return ayb->ncycle;
}

/** Replace any existing tile with the supplied one. */
AYB replace_AYB_tile(AYB ayb, const TILE tile) {

    /* free any tile memory, is allocated in new_AYB */
    free_TILE(ayb->tile);
    ayb->tile = copy_TILE(tile);
    return ayb;
}

/** Show bases in a single line. */
void show_AYB_bases(XFILE * fp, const AYB ayb, const uint32_t cl) {

    const uint32_t ncycle = ayb->ncycle;

    for (uint32_t cy = 0; cy < ncycle; cy++){
        show_NUC(fp, ayb->bases.elt[cl * ncycle + cy]);
    }
}

/** Show qualities in a single line. */
void show_AYB_quals(XFILE * fp, const AYB ayb, const uint32_t cl) {

    const uint32_t ncycle = ayb->ncycle;

    for (uint32_t cy = 0; cy < ncycle; cy++){
        show_PHREDCHAR(fp, ayb->quals.elt[cl * ncycle + cy]);
    }
}


/*
 * Routines to calculate covariance of errors, using the "fast approach".
 * Working with processed intensities.
 */

/**
 * Calculate covariance of (processed) residuals.
 * Returns a pointer to array of matrices, one per cycle.
 */
MAT * calculate_covariance(AYB ayb, bool all){

    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle = ayb->ncycle;

    MAT * V = NULL;                     // memory allocated in accumulate
    MAT pcl_int = NULL;                 // Shell for processed intensities
    MAT Minv_t = NULL, Pinv_t = NULL;   // Temporary matrices to speed calculation
    /* Create t(inverse(M)) and t(inverse(P)) */
    Minv_t = transpose_inplace(invert(ayb->M));
    if (NULL == Minv_t) { goto cleanup; }
    Pinv_t = transpose_inplace(invert(ayb->P));
    if (NULL == Pinv_t) { goto cleanup; }

    real_t wesum = 0.;

    unsigned int cl = 0;
    LIST(CLUSTER) node = ayb->tile->clusterlist;
    while (NULL != node && cl < ncluster){
        const NUC * cl_bases = ayb->bases.elt + cl * ncycle;
        pcl_int = process_intensities(node->elt->signals, Minv_t, Pinv_t, ayb->N, pcl_int);
        if (NULL == pcl_int) { goto cleanup; }

        /* add this cluster values */
        if (all) {
            V = accumulate_all_covariance(ayb->we->x[cl], pcl_int, ayb->lambda->x[cl], cl_bases, V);
        }
        else {
            V = accumulate_covariance(ayb->we->x[cl], pcl_int, ayb->lambda->x[cl], cl_bases, V);
        }
        if (NULL == V) { goto cleanup; }

        /* sum denominator */
        wesum += ayb->we->x[cl];

        /* next cluster */
        node = node->nxt;
        cl++;
    }

    pcl_int = free_MAT(pcl_int);
    Pinv_t = free_MAT(Pinv_t);
    Minv_t = free_MAT(Minv_t);

    /* scale sum of squares to make covariance */
    if (all) {
        if (NULL == V || NULL == V[0]) { goto cleanup; }
        scale_MAT(V[0], 1.0/wesum);
    }
    else {
        for ( uint32_t cy = 0; cy < ncycle; cy++){
            if (NULL == V || NULL == V[cy]) { goto cleanup; }
            /* add a diagonal offset to help with bad data */
            for ( uint32_t i = 0; i < NBASE; i++){
                V[cy]->x[i * NBASE + i] += 1.0;
            }
            scale_MAT(V[cy], 1.0/wesum);
        }
    }
    return V;

cleanup:
    if (NULL != V) {
        const uint32_t nc =  all ? 1 : ncycle;
        for (uint32_t cy = 0; cy < nc; cy++) {
            free_MAT(V[cy]);
        }
        xfree(V);
    }
    free_MAT(pcl_int);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
    return NULL;
}

/**
 * Call bases. Includes calculate covariance and estimate lambda.
 * Last iteration parameter for final working.
 * Returns number of zero lambdas or data error indication.
 */
int estimate_bases(AYB ayb, const int blk, const bool lastiter, const bool showdebug) {

    validate(NULL != ayb, 0);
    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle = ayb->ncycle;
    MAT pcl_int = NULL;                 // Shell for processed intensities
    MAT Pinv_t=NULL, Minv_t=NULL;       // Temporary matrices to speed computation
    MAT * V = NULL;                     // Array of covariance matrices
    real_t * qual = NULL;
    int ret_count = 0;
    bool show_processed = (ShowWorking && lastiter);
    WORKPTR work = NULL;

#ifndef NDEBUG
    XFILE *fpi2 = NULL;
    XFILE *fpout = NULL;
    if (showdebug) {
        fpout = open_output("ayb2");
        if (!xfisnull(fpout)) {
            show_AYB(fpout, ayb, false);
        }
        fpout = xfclose(fpout);
    }
#endif

    /* calculate covariance */
    V = calculate_covariance(ayb, false);

    if (V == NULL) {
        /* set calls to null and terminate processing */
        unsigned int cl = 0;

        LIST(CLUSTER) node = ayb->tile->clusterlist;
        while (NULL != node && cl < ncluster){
            NUC * cl_bases = ayb->bases.elt + cl * ncycle;
            PHREDCHAR * cl_quals = ayb->quals.elt + cl * ncycle;

            /* call null base for each cycle */
            for (uint32_t cy = 0; cy < ncycle; cy++){
                struct basequal bq = call_base_null();
                cl_bases[cy] = bq.base;
                cl_quals[cy] = MIN_PHRED;// bq.qual is in quality-score space, want PHRED
            }

            /* next cluster */
            node = node->nxt;
            cl++;
        }

        ret_count = DATA_ERR;
        goto cleanup;
    }

#ifndef NDEBUG
    if (showdebug) {
        fpout = open_output("cov");
        if (!xfisnull(fpout)) {
            xfputs("covariance:\n", fpout);
            for (uint32_t cy = 0; cy < ncycle; cy++){
                show_MAT(fpout, V[cy], NBASE, NBASE);
            }
        }
        fpout = xfclose(fpout);
    }
#endif

    /* scale is variance of residuals; get from V matrices */
    for (uint32_t cy = 0; cy < ncycle; cy++){
        ayb->cycle_var->x[cy] = 0.;
        for (uint32_t b = 0; b < NBASE; b++){
            ayb->cycle_var->x[cy] += V[cy]->x[b * NBASE + b];
        }
        /* hmhm I have no idea why a reciprocal is taken here */
//        ayb->cycle_var->x[cy] = 1.0/ayb->cycle_var->x[cy];
    }

    /* invert variance matrices to get omega */
    for (uint32_t cy = 0; cy < ncycle; cy++){
        MAT a = V[cy];
        V[cy] = invert(a);
        free_MAT(a);
    }

    /* process intensities then estimate lambda and call bases for each cluster */
    Minv_t = transpose_inplace(invert(ayb->M));
    Pinv_t = transpose_inplace(invert(ayb->P));

#ifndef NDEBUG
    if (showdebug) {
        fpi2 = open_output("pi2");
        fpout = open_output("lam2");
        if (!xfisnull(fpout)) {
            xfputs("lambda:\n", fpout);
        }
    }
#endif

    /* output processed intensities if requested and final iteration */
    if (show_processed) {
        work = open_processed(ayb, blk);
    }

    unsigned int cl = 0;

    /* Temporary storage for (real_t) quality values */
    qual = calloc(ncycle, sizeof(real_t));
    if (NULL == qual) {
        ret_count = DATA_ERR;
        goto cleanup;
    }

    LIST(CLUSTER) node = ayb->tile->clusterlist;
    while (NULL != node && cl < ncluster){
        NUC * cl_bases = ayb->bases.elt + cl * ncycle;
        PHREDCHAR * cl_quals = ayb->quals.elt + cl * ncycle;
        pcl_int = process_intensities(node->elt->signals, Minv_t, Pinv_t, ayb->N, pcl_int);
        if (NULL == pcl_int) {
            ret_count = DATA_ERR;
            goto cleanup;
        }

        if (show_processed) {
            write_processed(work, ayb, node->elt, cl, pcl_int);
        }

        /* estimate lambda using Weighted Least Squares */
//        ayb->lambda->x[cl] = estimate_lambdaGWLS(pcl_int, cl_bases, ayb->lambda->x[cl], ayb->cycle_var->x, V);
        ayb->lambda->x[cl] = estimate_lambdaWLS(pcl_int, cl_bases, ayb->lambda->x[cl], ayb->cycle_var->x);
        if (ayb->lambda->x[cl] == 0.0) {
            ret_count++;
        }

#ifndef NDEBUG
    if (showdebug) {
        if (!xfisnull(fpi2)) {
            show_MAT(fpi2, pcl_int, pcl_int->nrow, pcl_int->ncol);
        }
        if (!xfisnull(fpout)) {
            xfprintf(fpout, "%d: %#12.6f\n", cl + 1, ayb->lambda->x[cl]);
        }
    }
#endif

        /* call bases for each cycle */
        for (uint32_t cy = 0; cy < ncycle; cy++){
            /* deal differently with missing data */
            if (nodata(node->elt->signals->x + cy * NBASE, NBASE)) {
                cl_bases[cy] = call_base_nodata();
                qual[cy] = MIN_QUALITY;
            }
            else {
                struct basequal bq = call_base(pcl_int->x+cy * NBASE, ayb->lambda->x[cl], BasePenalty, V[cy]);
                cl_bases[cy] = bq.base;
                qual[cy] = bq.qual;
            }
        }

        if (ayb->lambda->x[cl] != 0.0) {
            /* adjust quality values; first and Last bases of read are special cases */
            qual[0] = adjust_first_quality(qual[0], cl_bases[0], cl_bases[1]);
            for (uint32_t cy=1; cy<(ncycle-1); cy++) {
                qual[cy] = adjust_quality(qual[cy], cl_bases[cy-1], cl_bases[cy], cl_bases[cy+1]);
            }
            qual[ncycle-1] = adjust_last_quality(qual[ncycle-1], cl_bases[ncycle-2], cl_bases[ncycle-1]);

        }

        /* convert qualities to phred char */
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            cl_quals[cy] = phredchar_from_quality(qual[cy]);
        }

        /* repeat estimate lambda with the new bases */
        ayb->lambda->x[cl] = estimate_lambdaWLS(pcl_int, cl_bases, ayb->lambda->x[cl], ayb->cycle_var->x);

        /* next cluster */
        node = node->nxt;
        cl++;
    }

#ifndef NDEBUG
    if (showdebug) {
        xfclose(fpi2);
        xfclose(fpout);
    }
#endif

/* cleanup for success and error states */
cleanup:
    if (show_processed) {
        work = close_processed(work, ayb, blk, ret_count);
        xfree(work);
    }

    xfree(qual);
    free_MAT(pcl_int);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
    if (NULL != V) {
        for (uint32_t cy = 0; cy < ncycle; cy++) {
            free_MAT(V[cy]);
        }
        xfree(V);
    }
    return ret_count;
}

/**
 * Code for parameter estimation loop.
 * Equivalent to parameter_estimation_loop in old AYB.
 * - Updates: M, P, N
 * - Recalculates weights
 * - Scales all lambda by factor
 */
real_t estimate_MPN(AYB ayb){
    validate(NULL!=ayb,NAN);
    const uint32_t ncycle = ayb->ncycle;
    /*  Calculate new weights */
    //timestamp("Updating weights\n",stderr);
    real_t sumLSS = update_cluster_weights(ayb);

    /*  Precalculate terms for iteration */
    //timestamp("Calculating matrices\n",stderr);
    //timestamp("J\t",stderr);
    MAT J = calculateJ(ayb->lambda,ayb->we,ayb->bases,ncycle,NULL);
    MAT Jt = transpose(J);
    //timestamp("K\t",stderr);
    MAT K = calculateK(ayb->lambda,ayb->we,ayb->bases,ayb->tile,ncycle,NULL);
    //timestamp("Others\n",stderr);
    MAT Kt = transpose(K);
    MAT Sbar = calculateSbar(ayb->lambda,ayb->we,ayb->bases,ncycle,NULL);
    MAT SbarT = transpose(Sbar);
    MAT Ibar = calculateIbar(ayb->tile,ayb->we,NULL);
    MAT IbarT = transpose(Ibar);
    real_t Wbar = calculateWbar(ayb->we);
    real_t lambdaf = 1.;
    real_t * tmp = calloc(ncycle*ncycle*NBASE*NBASE,sizeof(real_t));
    /* Convenience terms for M, P and N */
    MAT matMt = transpose_inplace(ayb->M);
    MAT matP = ayb->P;
    MAT matN = ayb->N;

    //xfputs("Staring main loop",xstderr);
    /*  Main iteration */
    MAT mlhs=NULL,mrhs=NULL,plhs=NULL,prhs=NULL;
    const uint32_t lda = NBASE + ncycle;
    real_t det=0.;
    //xfprintf(xstderr,"Starting\tdelta = %e\n",calculateDeltaLSE( matMt, matP, matN, J, K, tmp));
    for( uint32_t i=0 ; i<AYB_NITER ; i++){
        /*
         *  Solution for phasing
         */
        plhs = calculatePlhs(Wbar,Sbar,matMt,J,tmp,plhs);
        prhs = calculatePrhs(Ibar,matMt,Sbar,matN,K,tmp,prhs);
        SolverRoutine(plhs,prhs,tmp, DELTA_DIAG);
        for ( uint32_t i=0 ; i<ncycle ; i++){
            for(uint32_t j=0 ; j<ncycle ; j++){
                matP->x[i*ncycle+j] = prhs->x[i*ncycle+j];
            }
        }

        // Scaling so det(P) = 1
        det = normalise_MAT(matP,3e-8);
        scale_MAT(J,det*det); scale_MAT(Jt,det*det);
        scale_MAT(K,det);     scale_MAT(Kt,det);
        scale_MAT(Sbar,det);  scale_MAT(SbarT,det);
        lambdaf *= det;

        /*
         *  Solution for cross-talk and constant noise
         */
        mlhs = calculateMlhs(ayb->cycle_var,Wbar,SbarT,matP,Jt,tmp,mlhs);
        mrhs = calculateMrhs(ayb->cycle_var,IbarT,matP,Kt,tmp,mrhs);
        solverSVD(mlhs,mrhs,tmp,DELTA_DIAG);

        for(uint32_t i=0 ; i<NBASE ; i++){
            for(uint32_t j=0 ; j<NBASE ; j++){
                matMt->x[i*NBASE+j] = mrhs->x[i*lda+j];
            }
        }
        for(uint32_t i=0 ; i<NBASE ; i++){
            for(uint32_t j=0 ; j<ncycle ; j++){
                matN->x[j*NBASE+i] = mrhs->x[i*lda+j+NBASE];
            }
        }
        // Scaling so det(M) = 1
        det = normalise_MAT(matMt,3e-8);
        scale_MAT(J,det*det); scale_MAT(Jt,det*det);
        scale_MAT(K,det);     scale_MAT(Kt,det);
        scale_MAT(Sbar,det);  scale_MAT(SbarT,det);
        lambdaf *= det;
        //xfprintf(xstderr,"... done %u\tdelta = %e\n",i,calculateDeltaLSE( matMt, matP, matN, J, K, tmp));
    }
    real_t delta = calculateDeltaLSE( matMt, matP, matN, J, K, tmp);

    // Transpose Mt back to normal form
    matMt = transpose_inplace(matMt);
    if (NULL==matMt) { message(E_MATRIX_FAIL_S, MSG_ERR, "transposed"); }

    // Scale lambdas by factor
    scale_MAT(ayb->lambda,lambdaf);

    // Clean-up memory
    free_MAT(prhs);
    free_MAT(plhs);
    free_MAT(mrhs);
    free_MAT(mlhs);
    xfree(tmp);
    free_MAT(IbarT);
    free_MAT(Ibar);
    free_MAT(SbarT);
    free_MAT(Sbar);
    free_MAT(Kt);
    free_MAT(K);
    free_MAT(Jt);
    free_MAT(J);

    //xfprintf(xstderr,"Initial %e\tImprovement %e\t = %e\n",sumLSS,delta,sumLSS-delta);
    //xfprintf(xstderr,"Updated weights %e\n", update_cluster_weights(ayb));
    return sumLSS-delta;
}

/**
 * Set initial values for the model.
 * Returns false if one of the initial matrices is wrong dimension or process intensities fails.
 */
bool initialise_model(AYB ayb, const bool showdebug) {

    /* initial M, P, N */
    ayb->M = init_matrix(ayb->M, E_CROSSTALK);
    if (ayb->M == NULL) {
        message (E_MATRIX_FAIL_S, MSG_ERR, MATRIX_TEXT[E_CROSSTALK]);
        return false;
    }
    ayb->N = init_matrix(ayb->N, E_NOISE);
    if (ayb->N == NULL) {
        message (E_MATRIX_FAIL_S, MSG_ERR, MATRIX_TEXT[E_NOISE]);
        return false;
    }
    ayb->P = init_matrix(ayb->P, E_PHASING);
    if (ayb->P == NULL) {
        message (E_MATRIX_FAIL_S, MSG_ERR, MATRIX_TEXT[E_PHASING]);
        return false;
    }

    /* Initial weights and cycles are all equal. Arbitrarily one */
    set_MAT(ayb->we, 1.0);
    set_MAT(ayb->cycle_var, 1.0);

    /* process intensities then call initial bases and lambda for each cluster */
    MAT pcl_int = NULL;                     // Shell for processed intensities
    MAT Minv_t = transpose_inplace(invert(ayb->M));
    MAT Pinv_t = transpose_inplace(invert(ayb->P));
    bool ret = true;

#ifndef NDEBUG
    XFILE *fpout = NULL;
    if (showdebug) {
        fpout = open_output("pi");
    }
#endif

    unsigned int cl = 0;
    LIST(CLUSTER) node = ayb->tile->clusterlist;
    while (NULL != node && cl < ayb->ncluster){
        pcl_int = process_intensities(node->elt->signals, Minv_t, Pinv_t, ayb->N, pcl_int);
        if (NULL == pcl_int) {
            ret = false;
            goto cleanup;
        }

#ifndef NDEBUG
    if (showdebug) {
        if (!xfisnull(fpout)) {
            show_MAT(fpout, pcl_int, pcl_int->nrow, pcl_int->ncol);
        }
    }
#endif

        /* call initial bases for each cycle */
        NUC * cl_bases = ayb->bases.elt + cl * ayb->ncycle;
        PHREDCHAR * cl_quals = ayb->quals.elt + cl * ayb->ncycle;
        for ( uint32_t cy = 0; cy < ayb->ncycle; cy++){
            /* deal differently with missing data */
            if (nodata(node->elt->signals->x + cy * NBASE, NBASE)) {
                cl_bases[cy] = call_base_nodata();
            }
            else {
                cl_bases[cy] = call_base_simple(pcl_int->x + cy * NBASE);
            }
            cl_quals[cy] = MIN_PHRED;
        }
        /* initial lambda */
        ayb->lambda->x[cl] = estimate_lambdaOLS(pcl_int, cl_bases);

        /* next cluster */
        node = node->nxt;
        cl++;
    }

#ifndef NDEBUG
    if (showdebug) {
        xfclose(fpout);
    }
#endif

cleanup:
    free_MAT(pcl_int);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
    return ret;
}

/**
 * Set the GC composition of the sequence.
 * Reads proportion from string and creates a penalty vector for the
 * least square statistic such that the penalty is zero for equally likely bases.
 * Returns true if valid value (0 to 1 non-inclusive).
 */
bool set_composition(const char *comp_str){

    char *endptr;
    const double gc = strtod(comp_str, &endptr);

    /* limit GC to a proportion and avoid log overflow */
    if (gc <= 0.0 || gc >= 1.0) {return false;}

    const double pen_mid = 2.0 * log(0.5);
    BasePenalty[NUC_A] = BasePenalty[NUC_T] = -2.0 * log1p(-gc) + pen_mid;
    BasePenalty[NUC_C] = BasePenalty[NUC_G] = -2.0 * log(gc) + pen_mid;
    return true;
}

/** Set show working flag. */
void set_show_working(void) {

    ShowWorking = true;
}

/**
 * Set solver routine to use for estimation of P.
 * Text must match one of options specified in SOLVER_TEXT, case insensitive.
 * Returns true if match found.
 */
bool set_solver(const char *solver_str) {

    /* match to options */
    int matchidx = match_string(solver_str, SOLVER_TEXT, E_SOLVER_NUM);
    if (matchidx>=0) {
        SolverIndex = (SOLVID)matchidx;
        SolverRoutine = SOLVERS[matchidx];
        return true;
    }
    else {
        return false;
    }
}

/**
 * Start up; call at program start after options.
 * Issues option info messages.
 * Returns true if any M, N, P matrix and quality table reads are successful.
 */
bool startup_ayb(void) {

    message(E_OPT_SELECT_SS, MSG_INFO, "P solver", SOLVER_TEXT[SolverIndex]);
    message(E_OPT_SELECT_SE, MSG_INFO, "Mu value", get_mu());
    if (BasePenalty[0] != 0.0) {
        message(E_OPT_SELECT_SE, MSG_INFO, "Genome GC Composition", exp(log(0.5) - BasePenalty[NUC_C] / 2.0));
    }

    /* read any M, N, P */
    if (!read_matrices()) {
        return false;
    }

    /* read any quality calibration table */
    return read_quality_table();
}

/** Tidy up; call at program shutdown. */
void tidyup_ayb(void) {

    /* free memory */
    for (IOTYPE idx = (IOTYPE)0; idx < E_MNP; idx++) {
        Matrix[idx] = free_MAT(Matrix[idx]);
    }
}
