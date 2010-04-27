/** 
 * \file ayb_model.c
 * Top Level Modelling.
 * Used as a singleton class with local functions accessing global member data.
 *//* 
 *  Created : 14 Apr 2010
 *  Author  : Tim Massingham/Hazel Marsden
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

#include "ayb_model.h"
#include "call_bases.h"
#include "dirio.h"
#include "intensities.h"
#include "lambda.h"
#include "matrix.h"
#include "message.h"
#include "nuc.h"
#include "tile.h"


/* constants */
/* None      */

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
/* constants */
/* none */

/* members */

static unsigned long NCycle = 0;                ///< Number of cycles to analyse - move to Tile?.
static TILE Tile;                               ///< Intensities data read in as single Tile.

static MAT Matrix[E_NMATRIX];                   ///< Predetermined matrices.
static AYB Ayb = NULL;                          ///< The model data.


/* private functions */

/**
 * Initialise crosstalk (M), Phasing (P) or Noise (N) matrix using an internal method.
 * To be done.
 */
static bool standard_matrix(const IOTYPE idx) {
    bool found = false;

    switch (idx) {
        case E_CROSSTALK:
            break;

        case E_NOISE:
            break;

        case E_PHASING:
            break;

        default: ;
    }

    return found;
}
/**
 * Initialise crosstalk (M), Phasing (P) and Noise (N) matrices.
 * May be read in or initialised using an internal method.
 */
static bool init_matrix() {
    XFILE *fpmat = NULL;
    bool found = true;

    for (IOTYPE idx = 0; idx < E_NMATRIX; idx++) {
        fpmat = open_matrix(idx);
        if (fpmat == NULL) {
            /* no input file specified, initialise internally */
            found = standard_matrix(idx);
        }
        else {
            /* read in matrix */
            Matrix[idx] = read_MAT_from_column_file(fpmat);
            if (Matrix[idx] == NULL) {
                found = false;
            }
        }
        xfclose(fpmat);
        fpmat = NULL;

        /* exit if this iteration failed */
        if (!found) {
            break;
        }
    }

    return found;
}

/** Read and store an intensities input file. */
static void read_intensities(XFILE *fp) {

    unsigned int nc = NCycle;
    Tile = read_TILE(fp, &nc);
    if (nc < NCycle) {
        message(E_CYCLESIZE_DD, MSG_WARN, NCycle, nc);
        NCycle = nc;
    }
}

/** Set initial values for the model. */
static void initialise_model() {

    Ayb = new_AYB( NCycle, Tile->ncluster);
    if (Ayb == NULL){return;}

    /* initial intensities, just copy from read in tile;
     * will become more elaborate once selection of data implemented */
    /* tile has already been allocated memory in new_AYB */
    free_TILE(Ayb->tile);
    Ayb->tile = copy_TILE(Tile);

    /* initial M, P, N */
    copyinto_MAT(Ayb->M, Matrix[E_CROSSTALK]);
    copyinto_MAT(Ayb->N, Matrix[E_NOISE]);
    copyinto_MAT(Ayb->P, Matrix[E_PHASING]);

    /* Initial cross-talk from default array */
//    copyinto_MAT(ayb->M,initial_M);
    /* Initial noise is zero */
//    set_MAT(ayb->N,0.);
    /* Initial phasing */
//    #warning "Phasing not properly initialised yet"
//    ayb->P = identity_MAT(ayb->ncycle);

    /* Initial weights and cycles are all equal. Arbitrarily one */
    set_MAT(Ayb->we, 1.0);
    set_MAT(Ayb->cycle_var, 1.0);

#ifndef NDEBUG
    XFILE *fpout = NULL;
    fpout = open_output("inv");
#endif

    /* process intensities then call initial bases and lambda for each cluster */
    MAT pcl_int = NULL;                     // Shell for processed intensities
    MAT Minv_t = transpose_inplace(invert(Ayb->M));
    MAT Pinv_t = transpose_inplace(invert(Ayb->P));

#ifndef NDEBUG
    if (!xisnull_file(fpout)) {
        show_MAT(fpout, Minv_t, Minv_t->nrow, Minv_t->ncol);
        show_MAT(fpout, Pinv_t, Pinv_t->nrow, Pinv_t->ncol);
    }
    xfclose(fpout);
    fpout = NULL;
    fpout = open_output("pi");
#endif

    unsigned int cl = 0;
    LIST(CLUSTER) node = Ayb->tile->clusterlist;
    while (NULL != node && cl < Ayb->ncluster){
        /* process intensities */
        pcl_int = process_intensities(node->elt->signals, Minv_t, Pinv_t, Ayb->N, pcl_int);

#ifndef NDEBUG
        if (!xisnull_file(fpout)) {
            show_MAT(fpout, pcl_int, pcl_int->nrow, pcl_int->ncol);
        }
#endif

        /* call initial bases */
        NUC * cl_bases = Ayb->bases.elt + cl * Ayb->ncycle;
        PHREDCHAR * cl_quals = Ayb->quals.elt + cl * Ayb->ncycle;
        for ( uint32_t cy = 0; cy < Ayb->ncycle; cy++){
            cl_bases[cy] = call_base_simple(pcl_int->x + cy * NBASE);
            cl_quals[cy] = MIN_PHRED;
        }
        /* initial lambda */
        Ayb->lambda->x[cl++] = estimate_lambdaOLS(pcl_int, cl_bases);

        node = node->nxt;
    }

#ifndef NDEBUG
    xfclose(fpout);
#endif

    free_MAT(pcl_int);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
}

/**
 * Output the results of the base calling.
 */
static void output_results () {

    if (Ayb == NULL){return;}

    XFILE *fpout;
    fpout = open_output("seq");

    if (!xisnull_file(fpout)) {

        const uint32_t ncluster = Ayb->ncluster;
        const uint32_t ncycle = Ayb->ncycle;

        for (uint32_t cl = 0; cl < ncluster; cl++){
            xfprintf(fpout, ">cluster_%03u\n", cl + 1);
            for (uint32_t cy = 0; cy < ncycle; cy++){
                show_NUC(fpout, Ayb->bases.elt[cl * ncycle + cy]);
            }
//            fputs("\n+\n",fp);
            xfputs("\n", fpout);
//            for (uint32_t cy = 0; cy < ncycle; cy++){
//                show_PHREDCHAR(fp,ayb->quals.elt[cl*ncycle+cy]);
//            }
//            fputc('\n',fp);
        }
    }
    xfclose(fpout);
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

void free_AYB(AYB ayb){
    if(NULL==ayb){ return;}
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

void show_AYB(XFILE * fp, const AYB ayb){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    xfprintf(fp,"%u cycles from %u clusters\n",ayb->ncycle,ayb->ncluster);
    xfputs("M:\n",fp); show_MAT(fp,ayb->M,NBASE,NBASE);
    xfputs("P:\n",fp); show_MAT(fp,ayb->P,ayb->ncycle,ayb->ncycle);
    xfputs("N:\n",fp); show_MAT(fp,ayb->N,NBASE,8);
    xfputs("we:\n",fp); show_MAT(fp,ayb->we,8,1);
    xfputs("cycle_var:\n",fp); show_MAT(fp,ayb->cycle_var,8,1);
    xfputs("lambda:\n",fp); show_MAT(fp,ayb->lambda,8,1);
    xfputs("Bases:\n",fp); show_ARRAY(NUC)(fp,ayb->bases,"",ayb->ncycle);
    xfputc('\n',fp);
    xfputs("Quality:\n",fp); show_ARRAY(PHREDCHAR)(fp,ayb->quals,"",ayb->ncycle);
    xfputc('\n',fp);
#ifdef NDEBUG
    xfputs("Tile:\n",fp); show_TILE(fp,ayb->tile,10);
#else
    xfputs("Tile:\n",fp); show_TILE(fp,ayb->tile,ayb->ncluster);
#endif
    xfputc('\n',fp);
}

/** Analyse a single input file. File is already opened. */
void analyse_tile (XFILE *fp) {

    if (xisnull_file(fp)) {return;}

    /* read intensity data from supplied file */
    read_intensities(fp);
    if (Tile == NULL) {
        message(E_BAD_INPUT_S, MSG_ERR, get_current_file());
        return;
    }

    /* set initial model values */
    initialise_model();
    /* no longer need the raw data as read in */
    free_TILE(Tile);
    Tile = NULL;

    if (Ayb == NULL) {
        message(E_INIT_FAIL_S, MSG_ERR, get_current_file());
        return;
    }
/*
for 1 to niter
    estimate_MPC (ayb)
    estimate_Bases (ayb)
end for
dump_fastq (AYB)
*/
#ifndef NDEBUG
    XFILE *fpout;
    fpout = open_output("ayb");
    if (!xisnull_file(fpout)) {
        show_AYB(fpout, Ayb);
    }
    xfclose(fpout);
#endif

    /* output the results */
    output_results ();

    /* free the structure ready for next */
    free_AYB (Ayb);
    Ayb = NULL;
}

/** Set the number of cycles to analyse. */
void set_ncycle(const char *ncycle_str) {

    char *endptr;
    NCycle = strtoul(ncycle_str, &endptr, 0);
}

/** Start up; call at program start after options. */
bool startup_model(){

    /* check number of cycles supplied - hmhm may be determined in a different way */
    message(E_GENERIC_SD, MSG_DEBUG, "ncycle:", NCycle);
    if (NCycle == 0) {
        message(E_NOCYCLES, MSG_FATAL);
        return false;
    }
    else {
        message(E_CYCLE_SELECT_D, MSG_INFO, NCycle);

        /* initialise M, N, P */
        return init_matrix();
    }
}

/** Tidy up; call at program shutdown. */
void tidyup_model(){

    /* free memory */
    for (IOTYPE idx = 0; idx < E_NMATRIX; idx++) {
        free_MAT(Matrix[idx]);
    }
}
