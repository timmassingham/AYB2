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
#include "aybthread.h"
#include <math.h>
#include "ayb.h"
#include "call_bases.h"
#include "cif.h"
#include "cluster.h"
#include "conjugate.h"
#include "dirio.h"
#include "intensities.h"
#include "lambda.h"
#include "lapack.h"
#include "message.h"
#include "mpn.h"
#include "nuc.h"
#include "qual_table.h"
#include "spikein.h"
#include "statistics.h"


/** AYB structure contains the data required for modelling. */
struct AybT {
    uint_fast32_t ncluster;
    uint_fast32_t ncycle;
    TILE tile;
    ARRAY(NUC) bases;
    ARRAY(PHREDCHAR) quals;
//    MAT M, P, N;
    MAT N;
    MAT At;
    MAT lambda;
    MAT lss;
    MAT we, cycle_var;
    MAT omega;
    bool *spiked;
};

/** Structure for spike-in quality counts. */
struct QSpikeT {
    uint_fast32_t obs;
    uint_fast32_t diff;
};

typedef struct QSpikeT * QSPIKEPTR;

/** Structure for output of final working values. */
struct WorkT {
    XFILE *fp;
    CIFDATA cif;
};
typedef struct WorkT * WORKPTR;


/* constants */

static const unsigned int AYB_NITER = 20;       ///< Number of parameter estimation loops.
static const real_t DELTA_DIAG = 1.0;           ///< Delta for solver routines.
static const real_t RIDGE_VAL = 100000.0;       ///< At and N solver constant.

/** Initial Crosstalk matrix if not read in, fixed values of approximately the right shape. */
static const real_t INITIAL_CROSSTALK[] = {
    2.0114300, 1.7217841, 0.06436576, 0.1126401,
    0.6919319, 1.8022413, 0.06436576, 0.0804572,
    0.2735545, 0.2252802, 1.39995531, 0.9976693,
    0.2896459, 0.2413716, 0.11264008, 1.3194981
};

/** Name text for matrix messages. Match to IOTYPE enum in dirio. */
static const char *MATRIX_TEXT[] = {"Crosstalk", "Noise", "Parameter A"};

/* members */

static MAT Matrix[E_MNP];                       ///< Predetermined matrices.
static MAT Initial_At = NULL;                   ///< Initial parameter matrix.
static bool FixedParam = false;                 ///< Use fixed supplied parameter matrices.
static bool ShowWorking = false;                ///< Set to output final working values.
static bool SpikeIn = false;                    ///< Use spike-in data.
static bool SpikeFound = false;                 ///< Spike-in data found for this tile block.
static bool SpikeCalib = false;                 ///< Calibrate qualities using spike-in data.


/* private functions */

/**
 * Accumulate covariance matrix (ncycle * NBASE x ncycle * NBASE).
 * Only lower triangular matrix is accumulated.
 * Note: If V is NULL, the required memory is allocated.
 * - p:        Matrix of processed intensities
 * - lambda:   Brightness of cluster
 * - base:     Current base call
 * - do_full:  Calculate full covariance matrix, or just band suitable for fitting Omega
 * - V:        Covariance matrix used for accumulation
 */
static MAT  accumulate_covariance( const real_t we, const MAT p, const real_t lambda, const NUC * base, const bool do_full, MAT V) {
    validate(NULL!=p, NULL);
    validate(NBASE==p->nrow, NULL);
    validate(lambda>=0.0, NULL);

    const uint_fast32_t ncycle = p->ncol;
    const int lda = ncycle * NBASE;

    // Allocate memory for V if necessary
    if (NULL==V) {
        V = new_MAT(lda, lda);
        if (NULL==V) { return NULL; }
    }

    // Perform accumulation. V += R R^t
    // Note: R = P - \lambda I_b, where I_b is unit vector with b'th elt = 1
    for (uint_fast32_t cy = 0; cy < ncycle; cy++){
        if (!isambig(base[cy])) {
            p->x[cy * NBASE + base[cy]] -= lambda;
        }
    }
    if ( do_full ){
    	syr(LAPACK_LOWER, &lda, &we, p->x, LAPACK_UNIT, V->x, &lda);
    } else {
	    // Only update base*base blocks on the diagonal or immediately above and below it
	    // Update more elements than necessary but excess are ignored later.
	    for ( int i=0 ; i<lda ; i++){
		    for ( int j=i ; j<i+2*NBASE ; j++ ){
			    if(j>=lda){ break; }
			    V->x[i*lda+j] += we * p->x[i] * p->x[j];
		    }
	    }
    }

    return V;
}

/** Calibrate qualities using spike-in data.
 *     if Qp = predicted quality, Nerrs = number different, Nobs = number total 
 * \n and Qa = actual quality (to find) then:
 * \n Pp = 10 ^ (-Qp/10)
 * \n Pa = (Nerrs + Pp)/(Nobs + 1)
 * \n Qa = -10 * log10(Pa)
 */
static void calibrate_by_spikein(AYB ayb, const int blk, QSPIKEPTR qspike) {

    validate(NULL != ayb, );
    validate(NULL != qspike, );
    
    /* calculate actual quality for each predicted */
    int qact[MAX_QUALITY + 1];
    for (int q = 0; q <= MAX_QUALITY; q++) {
        real_t pp = pow(10, -(real_t)q/10);
        real_t pa = (qspike[q].diff + pp) / (qspike[q].obs + 1);
        real_t qa = -10 * log10(pa);
        qact[q] = (int)(qa + 0.5);
    }

    if (SpikeCalib) {
        /* replace existing phred chars with calibrated quality ones */
        const uint_fast32_t ncluster = ayb->ncluster;
        const uint_fast32_t ncycle = ayb->ncycle;

        for (uint_fast32_t cl = 0; cl < ncluster; cl++) {
            PHREDCHAR * cl_quals = ayb->quals.elt + cl * ncycle;
            for (uint_fast32_t cy = 0; cy < ncycle; cy++){
                int q = qualint_from_phredchar(cl_quals[cy]);
                cl_quals[cy] = phredchar_from_quality(qact[q]);
            }
        }    
    }
    else {
        /* output counts to a file */
        XFILE *fp = NULL;
        fp = open_output_blk("qspike", blk);
        if (!xfisnull(fp)) {
            for (int q = 0; q <= MAX_QUALITY; q++) {
                xfprintf(fp, "%d\t%d\t%d\n", q, qspike[q].diff, qspike[q].obs);
            }
        }
        fp = xfclose(fp);
    }
}

/** Calibrate qualities using calibration tables. */
static void calibrate_by_table(const uint_fast32_t ncycle, const real_t lambda, NUC * cl_bases, real_t * qual) {            

    if (lambda != 0.0) {
        /* adjust quality values; first and Last bases of read are special cases */
        qual[0] = adjust_first_quality(qual[0], cl_bases[0], cl_bases[1]);
        for (uint_fast32_t cy=1; cy<(ncycle-1); cy++) {
            qual[cy] = adjust_quality(qual[cy], cl_bases[cy-1], cl_bases[cy], cl_bases[cy+1]);
        }
        qual[ncycle-1] = adjust_last_quality(qual[ncycle-1], cl_bases[ncycle-2], cl_bases[ncycle-1]);
    }
}

/**
 * Initialise crosstalk (M), Noise (N) and parameter (A) matrices.
 * May use read in values or initialise using an internal method.
 */
static MAT init_matrix(MAT mat, const IOTYPE idx, const MAT M) {

    if (mat == NULL) {return NULL;}

    if (Matrix[idx] == NULL) {
        /* initialise internally */
        switch (idx) {
            case E_CROSSTALK:
                /* initial crosstalk from default array */
                free_MAT(mat);
                mat = new_MAT_from_array(NBASE, NBASE, INITIAL_CROSSTALK);
                break;

            case E_NOISE:
                /* initial noise is zero */
                set_MAT(mat, 0.0);
                break;

            case E_PARAMA:
                if (M == NULL) {return NULL;}
                const uint_fast32_t lda = mat->nrow;
                if (mat->ncol != lda) {return NULL;}
                const uint_fast32_t ncycle = lda / NBASE;

                /* initial A has blocks of initial M down the diagonal */
                for (uint_fast32_t cy = 0; cy < ncycle; cy++) {
                    uint_fast32_t offset = cy * NBASE;
                    for (uint_fast32_t i = 0; i < NBASE; i++) {
                        for (uint_fast32_t j = 0; j < NBASE; j++) {
                            mat->x[(offset + i) * lda + (offset + j)] = M->x[i * NBASE + j];
                        }
                    }
                }
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
static inline bool nodata(const int_t * sig, const int num) {

    for (unsigned int i = 0; i < num; i++) {
        if (sig[i] != 0) {return false;}
    }
    return true;
}

/**
 * Read in any external crosstalk (M), Noise (N) and Param (A) matrices.
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

/** 
 * Read in any spike-in data. 
 * For each spike-in cluster store bases and flag to indicate a spike-in cluster.
 * Issues a warning message if a cluster number is out of range.
 */
static void read_spikein_data(AYB ayb, const int blk) {
    
    SpikeFound = false;
    if (!SpikeIn) { return; }
    validate(NULL != ayb, );
    if (!read_spikein_file(ayb->ncycle, blk)) { return; }

    SpikeFound = true;
    const uint_fast32_t ncluster = ayb->ncluster;
    const uint_fast32_t ncycle = ayb->ncycle;

    SPIKEIN spikein = get_next_spikein();
    while (spikein != NULL) {
        /* store the bases for the given cluster */
        uint_fast32_t cl = spikein->knum;
        if (cl < ncluster) {
            NUC * cl_bases = ayb->bases.elt + cl * ncycle;
            for (uint_fast32_t cy = 0; cy < ncycle; cy++){
                cl_bases[cy] = spikein->kbases.elt[cy];
            }

            /* set the spike-in cluster flag */
            ayb->spiked[cl] = true;
        }
        else {
            message(E_BAD_CLUSTER_SD, MSG_WARN, "spike-in file", cl + 1);
        }
        spikein = get_next_spikein();
    }

    /* don't reference the spike-in list again */
    tidyup_spikein();
}

/** Set all calls to null (before terminating processing on error). */
static void set_null_calls(AYB ayb) {

    const uint_fast32_t ncluster = ayb->ncluster;
    const uint_fast32_t ncycle = ayb->ncycle;

    unsigned int cl = 0;
    LIST(CLUSTER) node = ayb->tile->clusterlist;
    while (NULL != node && cl < ncluster){
        NUC * cl_bases = ayb->bases.elt + cl * ncycle;
        PHREDCHAR * cl_quals = ayb->quals.elt + cl * ncycle;

        /* call null base for each cycle */
        for (uint_fast32_t cy = 0; cy < ncycle; cy++){
            struct basequal bq = call_base_null();
            cl_bases[cy] = bq.base;
            cl_quals[cy] = MIN_PHRED;// bq.qual is in quality-score space, want PHRED
        }

        /* next cluster */
        node = node->nxt;
        cl++;
    }
}

/**
 * Calculate and store the least squares error for a cluster.
 * ip is the pre-calculated processed intensities for the given cluster.
 */
static void store_cluster_error(AYB ayb, MAT ip, unsigned int cl){
    validate(NULL!=ayb,);
    validate(NULL!=ip,);
    validate(cl < ayb->ncluster,);
    const uint_fast32_t ncycle = ayb->ncycle;

    ayb->we->x[cl] = 0.;
    NUC * cycle_bases = ayb->bases.elt + cl*ncycle;

    for ( uint_fast32_t i=0 ; i<ncycle ; i++){
        ip->x[i*NBASE+cycle_bases[i]] -= ayb->lambda->x[cl];
    }
    for( uint_fast32_t idx=0 ; idx<NBASE*ncycle ; idx++){
        ayb->we->x[cl] += ip->x[idx]*ip->x[idx];
    }
}

/** Calculate new weights assuming least squares error already stored. */
static real_t update_cluster_weights(AYB ayb){
    validate(NULL!=ayb,NAN);
    const uint_fast32_t ncluster = ayb->ncluster;
    real_t sumLSS = 0.;

    /* Calculate weight for each cluster */
    real_t meanLSSi = mean(ayb->we->x,ncluster);
    real_t varLSSi = variance(ayb->we->x,ncluster);
    if( varLSSi != 0.0){
        for ( uint_fast32_t cl=0 ; cl<ncluster ; cl++){
            sumLSS += ayb->we->x[cl];
            const real_t d = ayb->we->x[cl]-meanLSSi;
            ayb->we->x[cl] = cauchy(d*d,varLSSi);
	}
    } else {
	// If variance is zero, set all weights to one.
	for ( uint_fast32_t cl=0 ; cl<ncluster ; cl++){
            ayb->we->x[cl] = 1.0;
	}
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
static void write_processed(WORKPTR work, const AYB ayb, const CLUSTER cluster, const uint_fast32_t cl, MAT pcl_int) {

    const uint_fast32_t ncycle = ayb->ncycle;

    switch (get_input_format()) {
        case E_TXT:
            if (!xfisnull(work->fp)) {
                write_lane_tile (work->fp, ayb->tile);
                write_coordinates (work->fp, cluster);
                write_MAT_to_line(work->fp, pcl_int);
            }
            break;

        case E_CIF:
            for (uint_fast32_t cy = 0; cy < ncycle; cy++) {
                for (uint_fast32_t base = 0; base < NBASE; base++){
                    cif_set_from_real (work->cif, cl, base, cy, pcl_int->x[cy * NBASE + base]);
                }
            }
            break;

        default: ;
    }
}

/** Output final model values. */
static void output_final(const AYB ayb, const real_t effDF, const int blk) {
    XFILE *fpfin = NULL;

    /* final model values */
    fpfin = open_output_blk("final", blk);
    if (!xfisnull(fpfin)) {
        show_AYB(fpfin, ayb, false);
    }
    xfprintf(fpfin, "effDF: %#0.2f\n", effDF);
    xfclose(fpfin);

    /* final N, A in input format */
    fpfin = open_output_blk("N", blk);
    if (!xfisnull(fpfin)) {
        write_MAT_to_column_file (fpfin, ayb->N, false);
    }
    xfclose(fpfin);

    fpfin = open_output_blk("A", blk);
    if (!xfisnull(fpfin)) {
        /* input A is not transposed */
        transpose_inplace(ayb->At);
        write_MAT_to_column_file (fpfin, ayb->At, false);
    }
    xfclose(fpfin);
}

/**
 * Finish off and close up final processed intensities output.
 * Output is done here for cif and final model.
 * If status is error then just free resources.
 */
static WORKPTR close_processed(WORKPTR work, const AYB ayb, const real_t effDF, const int blk, const int status) {

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
        output_final(ayb, effDF, blk);
    }

    free_cif(work->cif);
    work->fp = xfclose(work->fp);
    return work;
}


/* public functions */

/* standard functions */

AYB new_AYB(const uint_fast32_t ncycle, const uint_fast32_t ncluster){
    AYB ayb = malloc(sizeof(*ayb));
    if(NULL==ayb){ return NULL;}
    ayb->ncycle = ncycle;
    ayb->ncluster = ncluster;
    ayb->tile = new_TILE();
    ayb->bases = new_ARRAY(NUC)(ncluster*ncycle);
    ayb->quals = new_ARRAY(PHREDCHAR)(ncluster*ncycle);
//    ayb->M = new_MAT(NBASE,NBASE);
//    ayb->P = new_MAT(ncycle,ncycle);
    ayb->N = new_MAT(NBASE,ncycle);
    ayb->At = new_MAT(NBASE*ncycle,NBASE*ncycle);
    ayb->lambda = new_MAT(ncluster,1);
    ayb->lss = new_MAT(ncluster,1);
    ayb->we = new_MAT(ncluster,1);
    ayb->cycle_var = new_MAT(ncycle,1);
    ayb->omega = NULL;
    ayb->spiked = calloc(ncluster, sizeof(bool));
    
    if( NULL==ayb->tile || NULL==ayb->bases.elt || NULL==ayb->quals.elt
//            || NULL==ayb->M || NULL==ayb->P || NULL==ayb->N
            || NULL==ayb->N || NULL==ayb->At
            || NULL==ayb->lambda || NULL==ayb->lss || NULL==ayb->we || NULL==ayb->cycle_var || NULL==ayb->spiked){
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
//    free_MAT(ayb->M);
//    free_MAT(ayb->P);
    free_MAT(ayb->N);
    free_MAT(ayb->At);
    free_MAT(ayb->lambda);
    free_MAT(ayb->lss);
    free_MAT(ayb->we);
    free_MAT(ayb->cycle_var);
    free_MAT(ayb->omega);
    xfree(ayb->spiked);
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

//    ayb_copy->M = copy_MAT(ayb->M);
//    if(NULL!=ayb->M && NULL==ayb_copy->M){ goto cleanup;}

//    ayb_copy->P = copy_MAT(ayb->P);
//    if(NULL!=ayb->P && NULL==ayb_copy->P){ goto cleanup;}

    ayb_copy->N = copy_MAT(ayb->N);
    if(NULL!=ayb->N && NULL==ayb_copy->N){ goto cleanup;}

    ayb_copy->At = copy_MAT(ayb->At);
    if(NULL!=ayb->At && NULL==ayb_copy->At){ goto cleanup;}

    ayb_copy->lambda = copy_MAT(ayb->lambda);
    if(NULL!=ayb->lambda && NULL==ayb_copy->lambda){ goto cleanup;}

    ayb_copy->lss = copy_MAT(ayb->lss);
    if(NULL!=ayb->lss && NULL==ayb_copy->lss){ goto cleanup;}

    ayb_copy->we = copy_MAT(ayb->we);
    if(NULL!=ayb->we && NULL==ayb_copy->we){ goto cleanup;}

    ayb_copy->cycle_var = copy_MAT(ayb->cycle_var);
    if(NULL!=ayb->cycle_var && NULL==ayb_copy->cycle_var){ goto cleanup;}

    ayb_copy->omega = copy_MAT(ayb->omega);
    if(NULL!=ayb->omega && NULL==ayb_copy->omega){ goto cleanup;}
    
    ayb->spiked = calloc(ayb->ncluster, sizeof(bool));
    if(NULL==ayb_copy->spiked){ goto cleanup;}
    memcpy(ayb_copy->spiked, ayb->spiked, ayb->ncluster * sizeof(bool));

    return ayb_copy;

cleanup:
    free_AYB(ayb_copy);
    return NULL;
}

void show_AYB(XFILE * fp, const AYB ayb, bool showall){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    xfprintf(fp,"%u cycles from %u clusters\n",ayb->ncycle,ayb->ncluster);
//    xfputs("M:\n",fp); show_MAT(fp,ayb->M,NBASE,NBASE);
//    xfputs("P:\n",fp); show_MAT(fp,ayb->P,ayb->ncycle,ayb->ncycle);
    xfputs("N:\n",fp); show_MAT(fp,ayb->N,NBASE,ayb->ncycle);
    xfputs("At:\n",fp); show_MAT(fp,ayb->At,NBASE*ayb->ncycle,NBASE*ayb->ncycle);
    xfputs("we:\n",fp); show_MAT(fp,ayb->we,ayb->ncluster,1);
    xfputs("cycle_var:\n",fp); show_MAT(fp,ayb->cycle_var,ayb->ncycle,1);
    xfputs("omega:\n",fp); show_MAT(fp,ayb->omega,NBASE*ayb->ncycle,NBASE*ayb->ncycle);
    xfputs("lambda:\n",fp); show_MAT(fp,ayb->lambda,ayb->ncluster,1);
    if (showall) {
        xfputs("lss:\n",fp); show_MAT(fp,ayb->lss,ayb->ncluster,1);
        xfputs("spike-in:\n",fp);
        for (uint_fast32_t cl=0; cl<ayb->ncluster; cl++){
            xfputc(ayb->spiked[cl] ?'1':'0',fp);
        }
        xfputc('\n',fp); xfputc('\n',fp);
        xfputs("Bases:\n",fp); show_ARRAY(NUC)(fp,ayb->bases,"",ayb->ncycle*10);
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
real_t * get_AYB_lambdas(const AYB ayb, uint_fast32_t *num) {

    real_t * lambdas = calloc(ayb->lambda->nrow, sizeof(real_t));
    *num = 0;

    for (uint_fast32_t i = 0; i < ayb->lambda->nrow; i++) {
        if (ayb->lambda->x[i] != 0.0) {
            lambdas[(*num)++] = ayb->lambda->x[i];
        }
    }
    return lambdas;
}

/** Return the number of clusters. */
uint_fast32_t get_AYB_ncluster(const AYB ayb) {
    return ayb->ncluster;
}

/** Return the number of cycles. */
uint_fast32_t get_AYB_ncycle(const AYB ayb) {
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
void show_AYB_bases(XFILE * fp, const AYB ayb, const uint_fast32_t cl) {

    const uint_fast32_t ncycle = ayb->ncycle;

    for (uint_fast32_t cy = 0; cy < ncycle; cy++){
        show_NUC(fp, ayb->bases.elt[cl * ncycle + cy]);
    }
}

/** Show qualities in a single line. */
void show_AYB_quals(XFILE * fp, const AYB ayb, const uint_fast32_t cl) {

    const uint_fast32_t ncycle = ayb->ncycle;

    for (uint_fast32_t cy = 0; cy < ncycle; cy++){
        show_PHREDCHAR(fp, ayb->quals.elt[cl * ncycle + cy]);
    }
}


/**
 * Calculate covariance of (processed) residuals.
 * Uses the "fast approach", working with processed intensities.
 * Returns full covariance matrix.
 */
MAT calculate_covariance(AYB ayb, const bool do_full){

    validate(NULL != ayb, NULL);
    unsigned int ncluster = ayb->ncluster;  // need unsigned int for array_from_LIST
    const uint_fast32_t ncycle = ayb->ncycle;

    MAT Vsum = NULL;                    // memory allocated in multi-thread accumulate
    real_t wesum = 0.0;
    bool ok = true;

    struct structLU AtLU = LUdecomposition(ayb->At);

    /* declare variables for multi-threading */
    LIST(CLUSTER) * nodearry = NULL;
    int_fast32_t cl;
    NUC * cl_bases = NULL;
    int th_id;                              // thread number
    const int ncpu = omp_get_max_threads();
    MAT pcl_int[ncpu];                      // Shell for processed intensities
    MAT V[ncpu];                            // memory allocated in accumulate
    real_t wei[ncpu];                       // accumulate by thread determines sum order
    for (int i = 0; i < ncpu; i++) {
        pcl_int[i] = NULL;
        V[i] = NULL;
        wei[i] = 0.0;
    }
    
    /* make an array of list pointers for multi-threading */
    nodearry = array_from_LIST(CLUSTER)(ayb->tile->clusterlist, &ncluster);
    if (NULL==nodearry) { 
        ok = false;
        goto cleanup; 
    }

    /* multi-threaded loop */
    #pragma omp parallel for \
        default(shared) private(th_id, cl, cl_bases)

    for (cl = 0; cl < ncluster; cl++){
        th_id = omp_get_thread_num();

        cl_bases = ayb->bases.elt + cl * ncycle;
        pcl_int[th_id] = processNew(AtLU, ayb->N, nodearry[cl]->elt->signals, pcl_int[th_id]);
        if (NULL == pcl_int[th_id]) {
            ok = false;
        }
        else {

            /* add this cluster values */
            V[th_id] = accumulate_covariance(ayb->we->x[cl], pcl_int[th_id], ayb->lambda->x[cl], cl_bases, do_full, V[th_id]);
            if (NULL == V[th_id]) { 
                ok = false; 
            }

            /* sum denominator */
            wei[th_id] += ayb->we->x[cl];
        }
    }
    
    if (ok) {
        /* Accumulate from multi-thread */ 
        const int lda = V[0]->nrow;
        Vsum = new_MAT(lda, lda);
        if (NULL == Vsum) {
            ok = false;
        }
        else {
            for (int i = 0; i < ncpu; i++) {
                for (int j = 0; j < lda * lda; j++) {
                    Vsum->x[j] += V[i]->x[j];
                }
                wesum += wei[i];
            }
	    // V's accumulated are lower triangular. Make symmetric.
	    symmeteriseL2U(Vsum);

            /* scale sum of squares to make covariance */
            scale_MAT(Vsum, 1.0/wesum);
        }
    }

cleanup:
    if (!ok) {  
        Vsum = free_MAT(Vsum);
    }
    
    free_array_LIST(CLUSTER)(nodearry);
    for (int i = 0; i < ncpu; i++) {
        free_MAT(pcl_int[i]);
        free_MAT(V[i]);
    }
    free_MAT(AtLU.mat);
    xfree(AtLU.piv);
    return Vsum;
}

/**
 * Call bases. Includes calculate covariance and estimate lambda.
 * Last iteration parameter for final working and qualities.
 * Returns number of zero lambdas or data error indication.
 */
int estimate_bases(AYB ayb, const int blk, const bool lastiter, const bool showdebug) {

    validate(NULL != ayb, 0);
    unsigned int ncluster = ayb->ncluster;  // need unsigned int for array_from_LIST
    const uint_fast32_t ncycle = ayb->ncycle;

    MAT V_part = NULL;                  // Partial covariance matrix
    QSPIKEPTR qspike = NULL;
    WORKPTR work = NULL;
    real_t effDF = NBASE * ncycle;
    int ret_count = 0;
    bool show_processed = (ShowWorking && lastiter);

    struct structLU AtLU = LUdecomposition(ayb->At);

    /* declare multi-threading variables required before any goto */
    LIST(CLUSTER) * nodearry = NULL;
    const int ncpu = omp_get_max_threads();
    MAT pcl_int[ncpu];                      // Shell for processed intensities
    for (int i = 0; i < ncpu; i++) {
        pcl_int[i] = NULL;
    }
    int zero_lam[ncpu];
    memset(zero_lam, 0, ncpu * sizeof(int));

#ifndef NDEBUG
    XFILE *fpi2 = NULL;
    XFILE *fplam = NULL;
    XFILE *fpout = NULL;
    if (showdebug) {
        fpout = open_output("ayb2");
        if (!xfisnull(fpout)) {
            show_AYB(fpout, ayb, false);
        }
        fpout = xfclose(fpout);
    }
#endif

    /* calculate partial covariance */
    V_part = calculate_covariance(ayb,false);

    if (V_part == NULL) {
        /* set calls to null and terminate processing */
        ret_count = DATA_ERR;
        goto cleanup;
    }

#ifndef NDEBUG
    if (showdebug) {
        fpout = open_output("covpart");
        if (!xfisnull(fpout)) {
            xfputs("covariance parital:\n", fpout);
            show_MAT(fpout, V_part, 0, 0);
        }
        fpout = xfclose(fpout);
    }
#endif

    /* scale is variance of residuals; get from V full matrix */
    for (uint_fast32_t cy = 0; cy < ncycle; cy++){
        ayb->cycle_var->x[cy] = 0.;
        for (uint_fast32_t b = 0; b < NBASE; b++){
            uint_fast32_t offset = cy * NBASE + b;
            ayb->cycle_var->x[cy] += V_part->x[offset * ncycle * NBASE + offset];
        }
    }

    /* calculate restricted fitted V inverse */
	ayb->omega = fit_omega(V_part, ayb->omega);
    if (ayb->omega == NULL) {
        /* set calls to null and terminate processing */
        ret_count = DATA_ERR;
        goto cleanup;
    }
	
#ifndef NDEBUG
    if (showdebug) {
        fpout = open_output("omfit");
        if (!xfisnull(fpout)) {
            xfputs("omega fitted:\n", fpout);
            show_MAT(fpout, ayb->omega, 0, 0);
        }
        fpout = xfclose(fpout);
    }
#endif    	

#ifndef NDEBUG
    if (showdebug) {
        fpi2 = open_output("pi2");
        fplam = open_output("lam2");
        if (!xfisnull(fplam)) {
            xfputs("lambda:\n", fplam);
        }
    }
#endif

    if (lastiter) {
        /* Calculate the median value of the lss array before updating any of the entries */
        /* (only used on last iteration ) */
        effDF = median(ayb->lss->x, ncluster);

        if (SpikeIn) {
            /* storage for spike-in quality counts */
            qspike = calloc(MAX_QUALITY + 1, sizeof(*qspike));
        }
    }

    /* declare variables for multi-threading */
    int_fast32_t cl;
    uint_fast32_t cy;
    NUC * cl_bases = NULL;
    PHREDCHAR * cl_quals = NULL;
    int th_id;                              // thread number
    
    /* make an array of list pointers for multi-threading */
    nodearry = array_from_LIST(CLUSTER)(ayb->tile->clusterlist, &ncluster);
    if (NULL == nodearry) { 
        /* set calls to null and terminate processing */
        ret_count = DATA_ERR;
        goto cleanup; 
    }

    /* multi-threaded loop */
    #pragma omp parallel for \
        default(shared) private(th_id, cl, cy, cl_bases, cl_quals)

    /* process intensities then estimate lambda and call bases for each cluster */
    for (cl = 0; cl < ncluster; cl++){
        th_id = omp_get_thread_num();

        cl_bases = ayb->bases.elt + cl * ncycle;
        cl_quals = ayb->quals.elt + cl * ncycle;
        pcl_int[th_id] = processNew(AtLU, ayb->N, nodearry[cl]->elt->signals, pcl_int[th_id]);
        if (NULL == pcl_int[th_id]) {
            ret_count = DATA_ERR;
        }
        else {

            /* estimate lambda using Weighted Least Squares */
//            ayb->lambda->x[cl] = estimate_lambdaWLS(pcl_int, cl_bases, ayb->lambda->x[cl], ayb->cycle_var->x);
            ayb->lambda->x[cl] = estimate_lambda_A (nodearry[cl]->elt->signals, ayb->N, ayb->At, cl_bases);
            if (ayb->lambda->x[cl] == 0.0) {
                zero_lam[th_id]++;
            }
        
#ifndef NDEBUG
    if (showdebug) {
        if (!xfisnull(fpi2)) {
            xfprintf(fpi2, "cluster: %u\n", cl + 1);
            show_MAT(fpi2, pcl_int[th_id], pcl_int[th_id]->nrow, pcl_int[th_id]->ncol);
        }
        if (!xfisnull(fplam)) {
            xfprintf(fplam, "%u: %#12.6f\n", cl + 1, ayb->lambda->x[cl]);
        }
    }
#endif

            /* call bases for all cycles */
            NUC sp_bases[ncycle];                   // thread private storage for copy of spike-in sequence

            /* only calculate lss for spike-in data clusters unless last iteration */
            if (!lastiter && ayb->spiked[cl]) {
                ayb->lss->x[cl] = calculate_lss(pcl_int[th_id], ayb->lambda->x[cl], ayb->omega, cl_bases);
            }
            else {
                if (ayb->spiked[cl]) {
                    /* save spiked-in sequence for diff counts */ 
                    memcpy(sp_bases, cl_bases, ncycle * sizeof(NUC));
                }
                ayb->lss->x[cl] = call_bases(pcl_int[th_id], ayb->lambda->x[cl], ayb->omega, cl_bases);
            }
            
            /* call qualities, only needed on last iteration */
            if (lastiter) {
                real_t qual[ncycle];                // thread private storage for (real_t) quality values
                call_qualities_post(pcl_int[th_id], ayb->lambda->x[cl], ayb->omega, effDF, cl_bases, qual);

                if (SpikeIn) {
                    if (ayb->spiked[cl]) {
                        /* add obs/diffs to counts */
                        for (cy = 0; cy < ncycle; cy++) {
                            int q = qualint_from_quality(qual[cy]);
                            qspike[q].obs++;
                            if (cl_bases[cy] != sp_bases[cy]) {
                                qspike[q].diff++;
                            }
                        }
                    }
                }
                
                else {
                    /* calibrate using calibration tables */
                    calibrate_by_table(ncycle, ayb->lambda->x[cl], cl_bases, qual); 
                }

                /* convert qualities to phred char */
                for ( cy=0 ; cy<ncycle ; cy++){
                    cl_quals[cy] = phredchar_from_quality(qual[cy]);
                }
            }

            /* repeat estimate lambda with the new bases */
            /* don't do if last iteration for working values */
            if (!lastiter) {
//                ayb->lambda->x[cl] = estimate_lambdaWLS(pcl_int, cl_bases, ayb->lambda->x[cl], ayb->cycle_var->x);
                ayb->lambda->x[cl] = estimate_lambda_A (nodearry[cl]->elt->signals, ayb->N, ayb->At, cl_bases);

                /* store the least squares error */
                store_cluster_error(ayb, pcl_int[th_id], cl);
            }
        }
    }

    /* accumulate from multi-thread */ 
    if (ret_count != DATA_ERR) {
        for (int i = 0; i < ncpu; i++) {
            ret_count += zero_lam[i];
        }
    }

    /* calibrate using spike-in data */
    if (lastiter && SpikeFound) {
        calibrate_by_spikein(ayb, blk, qspike);
    }

    /* output processed intensities if requested and final iteration */
    if (show_processed) {
        work = open_processed(ayb, blk);

        for (cl = 0; cl < ncluster; cl++){
            /* processed output must be done outide of multi-thread so need to process intensities again */
            pcl_int[0] = processNew(AtLU, ayb->N, nodearry[cl]->elt->signals, pcl_int[0]);
            if (NULL != pcl_int[0]) {
                write_processed(work, ayb, nodearry[cl]->elt, cl, pcl_int[0]);
            }
        }

        work = close_processed(work, ayb, effDF, blk, ret_count);
        xfree(work);
    }

#ifndef NDEBUG
    if (showdebug) {
        xfclose(fpi2);
        xfclose(fplam);
    }
#endif

/* cleanup for success and error states */
cleanup:
    if (ret_count == DATA_ERR) {
        set_null_calls(ayb);
    }

    free_array_LIST(CLUSTER)(nodearry);
    for (int i = 0; i < ncpu; i++) {
        free_MAT(pcl_int[i]);
    }
    free_MAT(AtLU.mat);
    xfree(AtLU.piv);
    free_MAT(V_part);
    xfree(qspike);
    return ret_count;
}

/**
 * Code for parameter estimation loop.
 * Equivalent to parameter_estimation_loop in old AYB.
 * - Updates: At, N
 * - Recalculates weights
 * - Scales all lambda by factor
 */
real_t estimate_MPN(AYB ayb){
    real_t ret = NAN;
    validate(NULL!=ayb, ret);
    const uint_fast32_t ncycle = ayb->ncycle;

    /*  Calculate new weights */
    //timestamp("Updating weights\n",stderr);
    real_t sumLSS = update_cluster_weights(ayb);
    if (isnan(sumLSS)) {
        set_null_calls(ayb);
        return ret;
    }

    MAT J = NULL, K = NULL;
    MAT Sbar = NULL, Ibar = NULL;
    MAT lhs = NULL, rhs = NULL;
    real_t * tmp = NULL;
    real_t lambdaf = 1.0;

    if (!FixedParam) {
        /*  Precalculate terms for iteration */
        //timestamp("Calculating matrices\n",stderr);
        //timestamp("J\t",stderr);
        J = calculateNewJ(ayb->lambda,ayb->bases,ayb->we,ncycle,NULL);
        //timestamp("K\t",stderr);
        K = calculateNewK(ayb->lambda,ayb->bases,ayb->tile,ayb->we,ncycle,NULL);
        //timestamp("Others\n",stderr);
        Sbar = calculateSbar(ayb->lambda,ayb->we,ayb->bases,ncycle,NULL);
        Ibar = calculateIbar(ayb->tile,ayb->we,NULL);
        real_t Wbar = calculateWbar(ayb->we);
        tmp = calloc(ncycle*ncycle*NBASE*NBASE,sizeof(real_t));
    
        lhs = calculateLhs(Wbar, J, Sbar, NULL);
        rhs = calculateRhs(K, Ibar, NULL);

        if ((NULL==J) || (NULL==K) || (NULL==Sbar) || (NULL==Ibar) || (NULL==lhs) || (NULL==rhs)) { goto cleanup; }
        /* assume ayb->At and Initial_At same size so only need to check one */
        if ((rhs->nrow != Initial_At->nrow + 1) || (rhs->ncol != Initial_At->ncol)) { goto cleanup; }
        if (rhs->ncol != (ayb->N->nrow * ayb->N->ncol)) { goto cleanup; }

        // add the solver constant
        // lhs -> lhs + r Id
        // rhs -> rhs + r initialA^t
        uint_fast32_t nrow = lhs->nrow;
        uint_fast32_t rnrow = rhs->nrow;
        for (uint_fast32_t i = 0; i < nrow; i++) {
            lhs->x[i * nrow + i] += RIDGE_VAL;
        }
        nrow = Initial_At->nrow;
        for (uint_fast32_t i = 0; i < Initial_At->ncol; i++) {
            for (uint_fast32_t j = 0; j < nrow; j++) {
                rhs->x[i * rnrow + j] += RIDGE_VAL * Initial_At->x[i * nrow + j];
            }
        }

        solverSVD(lhs, rhs, tmp, DELTA_DIAG);

        /* extract new At and N */
        nrow = ayb->At->nrow;
        for (uint_fast32_t i = 0; i < rhs->ncol; i++) {
            for (uint_fast32_t j = 0; j < nrow; j++) {
                ayb->At->x[i * nrow + j] = rhs->x[i * rnrow + j];
            }
            ayb->N->x[i] = rhs->x[i * rnrow + rnrow - 1];
        }

        lambdaf = normalise_MAT(ayb->At,3e-8);
    }

    else {
        /* keep parameters fixed; use copy of A matrix to get lambda scaling factor */
        MAT At_copy = copy_MAT(ayb->At);
        lambdaf = normalise_MAT(At_copy,3e-8);
        free_MAT(At_copy);
    }
    
    // Scale lambdas by factor
    scale_MAT(ayb->lambda,lambdaf);

    ret = sumLSS;

/* cleanup for success and error states */
cleanup:
    free_MAT(lhs);
    free_MAT(rhs);
    xfree(tmp);
    free_MAT(Ibar);
    free_MAT(Sbar);
    free_MAT(K);
    free_MAT(J);

    //xfprintf(xstderr,"Initial %e\tImprovement %e\t = %e\n",sumLSS,delta,sumLSS-delta);
    //xfprintf(xstderr,"Updated weights %e\n", update_cluster_weights(ayb));
    if (isnan(ret)) {
        set_null_calls(ayb);
    }
    return ret;
}

/**
 * Set initial values for the model.
 * Returns false if one of the initial matrices is wrong dimension or process intensities fails.
 */
bool initialise_model(AYB ayb, const int blk, const bool showdebug) {

    validate(NULL != ayb, false);
    unsigned int ncluster = ayb->ncluster;  // need unsigned int for array_from_LIST
    const int lda = NBASE * ayb->ncycle;

    /* initial M, N, A */
    MAT M = new_MAT(NBASE, NBASE);
    M = init_matrix(M, E_CROSSTALK, NULL);
    if (M == NULL) {
        message (E_MATRIX_FAIL_S, MSG_ERR, MATRIX_TEXT[E_CROSSTALK]);
        return false;
    }

    Initial_At = new_MAT(lda, lda);
    Initial_At = init_matrix(Initial_At, E_PARAMA, M);
    free_MAT(M);
    if (Initial_At == NULL) {
        message (E_MATRIX_FAIL_S, MSG_ERR, MATRIX_TEXT[E_PARAMA]);
        return false;
    }

    ayb->N = init_matrix(ayb->N, E_NOISE, NULL);
    if (ayb->N == NULL) {
        message (E_MATRIX_FAIL_S, MSG_ERR, MATRIX_TEXT[E_NOISE]);
        return false;
    }

    /* store A as transpose and for first iteration */
    transpose_inplace(Initial_At);
    copyinto_MAT(ayb->At, Initial_At);

    /* Initialise call bases return values to nbase * ncycle */
    set_MAT(ayb->lss, NBASE * ayb->ncycle);

    /* Initial weights and cycles are all equal. Arbitrarily one */
    set_MAT(ayb->we, 1.0);
    set_MAT(ayb->cycle_var, 1.0);

    /* read in and store any spike-in data */
    read_spikein_data(ayb, blk);

    struct structLU AtLU = LUdecomposition(ayb->At);
    bool ret = true;

#ifndef NDEBUG
    XFILE *fpout = NULL;
    if (showdebug) {
        fpout = open_output("pi");
    }
#endif

    /* declare variables for multi-threading */
    LIST(CLUSTER) * nodearry = NULL;
    int_fast32_t cl;
    uint_fast32_t cy;
    NUC * cl_bases = NULL;
    PHREDCHAR * cl_quals = NULL;
    int th_id;                              // thread number
    const int ncpu = omp_get_max_threads();
    MAT pcl_int[ncpu];                      // Shell for processed intensities
    for (int i = 0; i < ncpu; i++) {
        pcl_int[i] = NULL;
    }
    
    /* make an array of list pointers for multi-threading */
    nodearry = array_from_LIST(CLUSTER)(ayb->tile->clusterlist, &ncluster);
    if (NULL==nodearry) { 
        ret = false;
        goto cleanup; 
    }

    /* multi-threaded loop */
    #pragma omp parallel for \
        default(shared) private(th_id, cl, cy, cl_bases, cl_quals)

    /* process intensities then call initial bases and lambda for each cluster */
    for (cl = 0; cl < ncluster; cl++){
        th_id = omp_get_thread_num();

        cl_bases = ayb->bases.elt + cl * ayb->ncycle;
        cl_quals = ayb->quals.elt + cl * ayb->ncycle;
        pcl_int[th_id] = processNew(AtLU, ayb->N, nodearry[cl]->elt->signals, pcl_int[th_id]);
        if (NULL == pcl_int[th_id]) {
            ret = false;
        }
        else {

#ifndef NDEBUG
    if (showdebug) {
        if (!xfisnull(fpout)) {
            xfprintf(fpout, "cluster: %u\n", cl + 1);
            show_MAT(fpout, pcl_int[th_id], pcl_int[th_id]->nrow, pcl_int[th_id]->ncol);
        }
    }
#endif

            /* call initial bases for each cycle */
            for ( cy = 0; cy < ayb->ncycle; cy++){

                /* skip any spike-in data clusters */
                if (!ayb->spiked[cl]) {
                
                    /* deal differently with missing data */
                    if (nodata(nodearry[cl]->elt->signals->xint + cy * NBASE, NBASE)) {
                        cl_bases[cy] = call_base_nodata();
                    }
                    else {
                        cl_bases[cy] = call_base_simple(pcl_int[th_id]->x + cy * NBASE);
                    }
                }
                cl_quals[cy] = MIN_PHRED;
            }
            /* initial lambda */
//            ayb->lambda->x[cl] = estimate_lambdaOLS(pcl_int, cl_bases);
            ayb->lambda->x[cl] = estimate_lambda_A (nodearry[cl]->elt->signals, ayb->N, ayb->At, cl_bases);

            /* store the least squares error */
            store_cluster_error(ayb, pcl_int[th_id], cl);
        }
    }

#ifndef NDEBUG
    if (showdebug) {
        xfclose(fpout);
    }
#endif

/* cleanup for success and error states */
cleanup:
    free_array_LIST(CLUSTER)(nodearry);
    for (int i = 0; i < ncpu; i++) {
        free_MAT(pcl_int[i]);
    }
    free_MAT(AtLU.mat);
    xfree(AtLU.piv);
    return ret;
}

/** Set show working flag. */
void set_show_working(void) {

    ShowWorking = true;
}

/** Set spike-in data calibration flag. */
void set_spike_calib(void) {

    SpikeCalib = true;
}

/**
 * Start up; call at program start after options.
 * Issues option info messages.
 * Returns true if any M, N, P matrix and quality table reads are successful.
 */
bool startup_ayb(void) {

    /* check if spike-in data configured */
    SpikeIn = spike_in();
    if (SpikeIn) {
        if (SpikeCalib) {    
            message(E_QUALCALIB_S, MSG_INFO, "spike-in data");
        }
        else {
            message(E_QUALCALIB_S, MSG_INFO, "post-processing with spike-in data");
        }
        /* no quality calibration table output */
        set_noqualout();
    }
    else {
        message(E_QUALCALIB_S, MSG_INFO, "calibration table");
    }    

    message(E_OPT_SELECT_SG, MSG_INFO, "Generalised error value", get_generr());

    /* read any M, N, A */
    if (!read_matrices()) {
        return false;
    }
    /* check for fixed parameter matrices */
    if (Matrix[E_NOISE] == NULL) {
        if (Matrix[E_PARAMA] != NULL) {
            message(E_BAD_MATRIXIN, MSG_FATAL);
            return false;
        }
    }
    else {
        if (Matrix[E_PARAMA] == NULL) {
            message(E_BAD_MATRIXIN, MSG_FATAL);
            return false;
        }
        else {
            FixedParam = true;
            message(E_OPT_SELECT_SS, MSG_INFO, "Fixed parameters", "A and N");
        }
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
    free_MAT(Initial_At);
}
