/** 
 * \file cluster.c
 * Cluster Class.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
 *  Copyright (C) 2012 Disinformatics.org
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
#include "cluster.h"
#include "nuc.h"
#include "utility.h"


/* constants */
/* None      */

/* members */
/* None    */


/* private functions */

/** 
 * Read the coordinate and intensities part of a cluster line.
 * Number of cycles wanted is supplied, or zero with first means read all.
 * If first line then count available columns.
 * Returned value from matrix line is passed back as reference parameter.
 * Returns null if any problem.
 */
static CLUSTER read_CLUSTER_line(unsigned int *ncycle, char *ptr, const bool first) {
    CLUSTER cluster = NULL;
    MAT signals = NULL;
    int nc = 0, ncfound = 0;

    /* Read coordinate information */
    unsigned long int x, y;
    x = strtoul(ptr, &ptr, 0);
    if ('\t' != ptr[0]) {goto cleanup;}
    y = strtoul(ptr, &ptr, 0);

    if (first) {
        /* count number of cycles available, preserve original line pointer */
        char *ptr1 = ptr;
        ncfound = count_line_columns(NBASE, ptr1);
        /* set to read all available if not specified */
        if (*ncycle == 0) {
            *ncycle = ncfound;
        }
    }

    /* read cycle data */
    nc = *ncycle;
    signals = new_MAT_from_line(NBASE, &nc, ptr);
    if ((NULL == signals) || (nc == 0)) {goto cleanup;}

    /* store all the data in the cluster */
    cluster = new_CLUSTER();
    if (NULL == cluster) {goto cleanup;}
    cluster->x = x;
    cluster->y = y;
    cluster->signals = signals;

    if (first) {
        *ncycle = ncfound;
    }
    else {
        *ncycle = nc;
    }
    return cluster;

cleanup:
    free_CLUSTER(cluster);
    free_MAT(signals);
    return NULL;
}


/* public functions */

/*
 * Standard functions for cluster structure
 */
CLUSTER new_CLUSTER ( void){
    CLUSTER cl = calloc(1,sizeof(*cl));
    return cl;
}

CLUSTER free_CLUSTER(CLUSTER cluster){
    if(NULL==cluster){return NULL;}
    free_MAT(cluster->signals);
    xfree(cluster);
    return NULL;
}

CLUSTER copy_CLUSTER(const CLUSTER cluster){
    if(NULL==cluster){ return NULL;}
    CLUSTER cl = new_CLUSTER();
    if(NULL==cl){ return NULL;}
    cl->signals = copy_MAT(cluster->signals);
    if(NULL==cl->signals && NULL!=cluster->signals){ goto cleanup; }
    cl->x = cluster->x;
    cl->y = cluster->y;

    return cl;

cleanup:
    free_MAT(cl->signals);
    xfree(cl);
    return NULL;
}

void show_CLUSTER(XFILE * fp, const CLUSTER cluster){
    if(NULL==fp){ return;}
    if(NULL==cluster){ return;}
    xfprintf(fp,"Cluster coordinates: (%u,%u)\n",cluster->x,cluster->y);
    show_MAT(fp,cluster->signals,4,10);
}

/**
 * Create a new cluster from the supplied array of values.
 * It is the responsibility of the caller to ensure the supplied array is large enough
 * as no size checking can be done. Ignores x and y.
 * Returns a pointer to the following element in the array
 * (which may be beyond its range).
 */
CLUSTER coerce_CLUSTER_from_array(const unsigned int ncycle, int_t * x, int_t ** next){
    if(NULL==x){ return NULL;}
    CLUSTER cl = new_CLUSTER();
    if(NULL==cl){ return NULL;}
    cl->signals = coerce_MAT_from_intarray(NBASE, ncycle, x);
    if(NULL==cl->signals){ goto cleanup;}
    /* return the next element in the array */
    *next = x + NBASE * ncycle;
    return cl;

cleanup:
    free_MAT(cl->signals);
    xfree(cl);
    return NULL;
}

/** 
 * Append clustin onto clustout, selecting data columns. 
 * clustout may be null, in which case it is created using details from clustin.
 * See matrix() append_columns for error handling.
 */
CLUSTER copy_append_CLUSTER(CLUSTER clustout, const CLUSTER clustin, int colstart, int colend){

    /* validate parameters */
    if(NULL==clustin) {return clustout;}
    
    if (NULL==clustout){
        clustout = new_CLUSTER();
        if(NULL==clustout){ return NULL;}
        clustout->x = clustin->x;
        clustout->y = clustin->y;
        clustout->signals = NULL;
    }
    /* append the selected matrix columns */
    clustout->signals = append_columns(clustout->signals, clustin->signals, colstart, colend);
    
    /* hmhm any error checking required?? */
    return clustout;
}

/*
 * Input functions from data or file
 */

/**
 * Create a new cluster from the supplied cif data.
 * Index of cluster is supplied along with number of cycles wanted.
 */
CLUSTER read_cif_CLUSTER(CIFDATA cif, const unsigned int cl, unsigned int ncycle) {

    if (NULL==cif) {return NULL;}
    if ((ncycle > cif_get_ncycle(cif)) || (cl >= cif_get_ncluster(cif))) {return NULL;}

    MAT signals = new_MAT_int(NBASE, ncycle, true);
    if(signals==NULL) {return NULL;}
    CLUSTER cluster = new_CLUSTER();
    if(NULL==cluster) {goto cleanup;}

    for (int cy = 0; cy < ncycle; cy++) {
        for (int base = 0; base < NBASE; base++) {
            signals->xint[cy * NBASE + base] = clipint(cif_get_int(cif, cl, base, cy));
        }
    }
    /* x and y not available. Fake using cluster number */
    cluster->x = (int)floor( 0.5*(1+sqrt(1+8*cl)));
    cluster->y = 1 + cl - (cluster->x * (cluster->x-1))/2;
    cluster->signals = signals;
    return cluster;

cleanup:
    free_CLUSTER(cluster);
    free_MAT(signals);
    return NULL;
}

/** 
 * Read a cluster line from file pointer, including lane and tile. 
 * Format should be that of Illumina's _int.txt.
 * Number of cycles wanted is supplied, or zero means read all,
 * and actual number available returned.
 * A little validation of the file format is performed.
 * Returns null if any problem.
 */
CLUSTER read_first_CLUSTER( XFILE * fp, unsigned int *ncycle, 
                            unsigned int *lane, unsigned int *tile){
    if (NULL == fp) {return NULL;}

    /* Get line from file */
    char *line = NULL;
    size_t len = 0;
    line = xfgetln(fp, &len);
    if (NULL == line) {return NULL;}
    char *ptr = line;
    CLUSTER cluster = NULL;

    /* Read lane and tile information to return */
    *lane = (unsigned int)strtoul(ptr, &ptr, 0);
    if ('\t' != ptr[0]) {goto cleanup;}
    *tile = (unsigned int)strtoul(ptr, &ptr, 0);
    if ('\t' != ptr[0]) {goto cleanup;}

    /* read cluster information including available cycle count */
    cluster = read_CLUSTER_line(ncycle, ptr, true);

cleanup:
    xfree(line);
    return cluster;
}

/** 
 * Read a cluster line from file pointer. 
 * Format should be that of Illumina's _int.txt.
 * Number of cycles wanted is supplied and a different value may be returned if lower.
 * A little validation of the file format is performed.
 * Returns null if any problem.
 */
CLUSTER read_known_CLUSTER( XFILE * fp, unsigned int *ncycle){
    if (NULL == fp) {return NULL;}

    /* Get line from file */
    char *line = NULL;
    size_t len = 0;
    line = xfgetln(fp, &len);
    if (NULL == line) {return NULL;}
    char *ptr = line;
    CLUSTER cluster = NULL;

    /* Read and discard lane and tile information */
    strtoul(ptr, &ptr, 0); // skip value
    if ('\t' != ptr[0]) {goto cleanup;}
    strtoul(ptr, &ptr, 0); // skip value
    if ('\t' != ptr[0]) {goto cleanup;}

    /* read cluster information */
    cluster = read_CLUSTER_line(ncycle, ptr, false);

cleanup:
    xfree(line);
    return cluster;
}

// Stub function
CLUSTER read_unknown_CLUSTER( XFILE * fp){
    return NULL;
}

/**
 * Write the coordinates of a cluster to file.
 * Preceded with tabs as per Illumina int.txt format.
 */
void write_coordinates (XFILE * fp, const CLUSTER cluster) {
    xfprintf(fp, "\t%u\t%u", cluster->x, cluster->y);
}


#ifdef TEST
#include <err.h>

int main(int argc, char * argv[]){
    if(argc<3){
        /* arguments are number of cycles, input file (standard format) and optional cif input file */
        errx(EXIT_FAILURE, "Usage: test-cluster ncycle _int.txt_filename [cif_filename]");
    }

    unsigned int ncycle = 0;
    sscanf(argv[1], "%u", &ncycle);

    XFILE * fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    if (xfisnull(fp)) {
        errx(EXIT_FAILURE, "Failed to open supplied _int.txt file");
    }

    /* optional cif file testing */
    CIFDATA cif = NULL;
    bool docif = false;
    if (argc > 3) {
        XFILE * fpcif = xfopen(argv[3], XFILE_UNKNOWN, "r");
        if (xfisnull(fpcif)) {
            errx(EXIT_FAILURE, "Failed to open supplied cif file");
        }

        /* read in all the cif data */
        cif = readCIFfromStream(fpcif);
        if(NULL==cif){
            errx(EXIT_FAILURE, "Failed to read supplied cif file");
        }
        else {
            docif = true;
        }
        xfclose(fpcif);
    }

    xfputs("Read null file as first cluster\n", xstdout);
    unsigned int nc = ncycle;
    unsigned int lane = 0, tile = 0;
    CLUSTER cl = read_first_CLUSTER(NULL, &nc, &lane, &tile);
    if (cl==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        cl = free_CLUSTER(cl);
    }

    xfputs("Read all cycles first cluster from _int.txt file\n", xstdout);
    nc = 0;
    cl = read_first_CLUSTER(fp, &nc, &lane, &tile);
    if (cl==NULL) {
        errx(EXIT_FAILURE, "Failed to read supplied _int.txt file");
    }

    xfprintf(xstdout, "Cycles found: %u\n", nc);
    xfprintf(xstdout, "Lane: %u, Tile: %u\n", lane, tile);
    xfputs("Write coordinates: ", xstdout);
    write_coordinates(xstdout, cl);
    xfprintf(xstdout, "\n");
    show_CLUSTER(xstdout, cl);
    
    xfputs("Read first cluster from _int.txt file\n", xstdout);
    nc = ncycle;
    cl = free_CLUSTER(cl);
    cl = read_first_CLUSTER(fp, &nc, &lane, &tile);
    show_CLUSTER(xstdout, cl);
    
    xfputs("Read null file as next cluster\n", xstdout);
    CLUSTER cl2 = read_known_CLUSTER(NULL, &nc);
    if (cl2==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        cl2 = free_CLUSTER(cl2);
    }

    xfputs("Read next cluster from _int.txt file\n", xstdout);
    nc = ncycle;
    cl2 = read_known_CLUSTER(fp, &nc);
    show_CLUSTER(xstdout, cl2);
    xfclose(fp);
  
    xfputs("Copy null cluster\n", xstdout);
    cl2 = free_CLUSTER(cl2);
    cl2 = copy_CLUSTER(NULL);
    if (cl2==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        cl2 = free_CLUSTER(cl2);
    }
    
    xfputs("Copy cluster\n", xstdout);
    cl2 = copy_CLUSTER(cl);
    show_CLUSTER(xstdout, cl2);

    if (docif) {
        xfputs("Too many cycles cluster from cif file\n", xstdout);
        cl2 = free_CLUSTER(cl2);
        cl2 = read_cif_CLUSTER(cif, 0, cif_get_ncycle(cif) + 1);
        if (cl2==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            cl2 = free_CLUSTER(cl2);
        }

        xfputs("Invalid cluster from cif file \n", xstdout);
        cl2 = read_cif_CLUSTER(cif, cif_get_ncluster(cif), ncycle);
        if (cl2==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            cl2 = free_CLUSTER(cl2);
        }
    
        xfputs("Valid cluster from cif file\n", xstdout);
        cl2 = read_cif_CLUSTER(cif, 0, ncycle);
        show_CLUSTER(xstdout, cl2);
    }
    else {
        xfputs("Skip cif file tests\n", xstdout);
    }
    
    xfputs("Copy append cluster, from second column (1) to ncol/2\n", xstdout);
    cl = copy_append_CLUSTER(cl, cl2, 1, cl2->signals->ncol/2);
    show_CLUSTER(xstdout, cl);

    xfputs("Copy append cluster to null, from second column (1) to ncol/2\n", xstdout);
    CLUSTER cl3 = NULL;
    cl3 = copy_append_CLUSTER(cl3, cl2, 1, cl2->signals->ncol/2);
    show_CLUSTER(xstdout, cl3);

    xfputs("Create an array\n", xstdout);
    int_t arry[NBASE * ncycle];
    for (unsigned int idx = 0; idx < NBASE * ncycle; idx++) {
        arry[idx] = cl->signals->xint[idx];
    }
    xfputs("array values:", xstdout);
    for (unsigned int idx = 0; idx < NBASE * ncycle; idx++) {
        xfprintf(xstdout, INT_FORMAT, arry[idx]);
    }
    xfputs("\n", xstdout);

    xfputs("Coerce null array\n", xstdout);
    int_t *x;
    CLUSTER cl4 = coerce_CLUSTER_from_array(ncycle, NULL, &x);
    if (cl4==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
    }

    xfputs("Coerce cluster from array, ignores x,y\n", xstdout);
    x = arry;
    cl4 = coerce_CLUSTER_from_array(ncycle, x, &x);
    show_CLUSTER(xstdout, cl4);

    free_CLUSTER(cl);
    free_CLUSTER(cl2);
    free_CLUSTER(cl3);
    /* do not free cl4 as it points to arry */
    return EXIT_SUCCESS;
}
#endif

