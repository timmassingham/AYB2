/** 
 * \file tile.c
 * Tile Class.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
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

#include <stdlib.h>
#include "cif.h"
#include "tile.h"


/* constants */
/* None      */

/* members */
/* None    */


/* private functions */

/**
 * Create a tile from cif data.
 * Number of cycles required is specified, or zero means read all.
 */
static TILE create_TILE_from_cif(CIFDATA cif, unsigned int ncycle) {
    TILE tile = NULL;
    CLUSTER cl = NULL;
    LIST(CLUSTER) listtail = NULL;

    if(NULL==cif) {return NULL;}

    unsigned int cifcycle = cif_get_ncycle(cif);
    unsigned int cifcluster = cif_get_ncluster(cif);
    xfprintf(xstderr, "Read cif tile: %u cycles from %u clusters.\n", cifcycle, cifcluster);

    tile = new_TILE();
    if (NULL==tile) {return NULL;}

    if (ncycle == 0) {
        /* get all available */
        ncycle = cifcycle;
    }
    if (cifcycle < ncycle) {
        /* not enough cycles, just return the number without bothering to store the data */
        tile->ncycle = cifcycle;
    }
    else {
        tile->ncycle = ncycle;
        if (cifcycle > ncycle) {
            /* extra data */
            warnx("Intensity file contains more data than requested: additional %d cycles.", cifcycle - ncycle);
        }

        /* treat first cluster differently */
        cl = read_cif_CLUSTER(cif, tile->ncluster, ncycle);
        if (NULL==cl) {goto cleanup;}

        tile->clusterlist = cons_LIST(CLUSTER)(cl, tile->clusterlist);
        tile->ncluster++;

        /* store remaining clusters, appending to tail of cluster list */
        listtail = tile->clusterlist;
        while (tile->ncluster < cifcluster) {
            cl = read_cif_CLUSTER(cif, tile->ncluster, ncycle);
            if (NULL==cl) {break;}
            listtail = rcons_LIST(CLUSTER)(cl, listtail);
            tile->ncluster++;
        }
     }

    return tile;

cleanup:
    free_LIST(CLUSTER)(tile->clusterlist);
    xfree(tile);
    return NULL;
}


/* public functions */

/*
 * Standard functions for tile structure
 */
TILE new_TILE(void){
    TILE tile = calloc(1,sizeof(*tile));
    return tile;
}

TILE free_TILE(TILE tile){
    if(NULL==tile){ return NULL;}
    free_LIST(CLUSTER)(tile->clusterlist);
    xfree(tile);
    return NULL;
}

TILE copy_TILE(const TILE tile){
    if(NULL==tile){ return NULL;}
    TILE newtile = new_TILE();
    if(NULL==newtile){ return NULL;}
    newtile->lane = tile->lane;
    newtile->tile = tile->tile;
    newtile->ncluster = tile->ncluster;
    newtile->ncycle = tile->ncycle;
    newtile->clusterlist = copy_LIST(CLUSTER)(tile->clusterlist);
    if(NULL==newtile->clusterlist && NULL!=tile->clusterlist){ goto cleanup;}
    return newtile;

cleanup:
    free_LIST(CLUSTER)(newtile->clusterlist);
    xfree(newtile);
    return NULL;
}

void show_TILE(XFILE * fp, const TILE tile, const unsigned int n){
    if(NULL==fp){ return;}
    if(NULL==tile){ return;}
    xfprintf(fp,"Tile data structure: lane %u tile %u.\n",tile->lane, tile->tile);
    xfprintf(fp,"Number of clusters: %u.\n",tile->ncluster);
    xfprintf(fp,"Number of cycles: %u.\n",tile->ncycle);
    unsigned int ncl=0, maxcl=(n!=0 && n<tile->ncluster)?n:tile->ncluster;
    LIST(CLUSTER) node = tile->clusterlist;
    while (NULL!=node && ncl<maxcl){
        xfprintf(fp,"%d: ",ncl+1);
        show_CLUSTER(fp,node->elt);
        node = node->nxt;
        ncl++;
    }
    if( maxcl<tile->ncluster){ xfprintf(fp,"... (%u others)\n",tile->ncluster-maxcl); }
}

/**
 * Create a new tile from the supplied array of values.
 * It is the responsibility of the caller to ensure the supplied array is large enough
 * as no size checking can be done. Ignores lane and tile.
 */
TILE coerce_TILE_from_array(unsigned int ncluster, unsigned int ncycle, int_t * x){
    TILE tile = NULL;
    CLUSTER cl = NULL;
    LIST(CLUSTER) listtail = NULL;

    if(NULL==x){ return NULL;}
    tile = new_TILE();
    if(NULL==tile){ return NULL;}
    tile->ncluster = ncluster;
    tile->ncycle = ncycle;
    int_t * next_x;

    /* first cluster, create new list */
    cl = coerce_CLUSTER_from_array(ncycle, x, &next_x);
    if (NULL==cl){ goto cleanup;}
    tile->clusterlist = cons_LIST(CLUSTER)(cl, tile->clusterlist);

    /* read in remaining clusters, appending to tail of cluster list */
    listtail = tile->clusterlist;
    for (int idx = 1; idx < ncluster; idx++) {
        cl = coerce_CLUSTER_from_array(ncycle, next_x, &next_x);
        if (NULL==cl){ goto cleanup;}
        listtail = rcons_LIST(CLUSTER)(cl, listtail);
    }
    return tile;

cleanup:
    free_LIST(CLUSTER)(tile->clusterlist);
    xfree(tile);
    return NULL;
}

/** 
 * Append the clusters of tilein onto tileout, selecting data columns. 
 * tileout may be null or have an empty clusterlist, 
 * in which case new clusters are created using details from tilein.
 * See matrix() append_columns for error handling.
 */
TILE copy_append_TILE(TILE tileout, const TILE tilein, int colstart, int colend){

    /* validate parameters */
    if(NULL==tilein) {return tileout;}
    
    if (NULL==tileout){
        tileout = new_TILE();
        if(NULL==tileout) {return NULL;}
        tileout->lane = tilein->lane;
        tileout->tile = tilein->tile;
    }
     
    if (tileout->clusterlist == NULL) {   
        /* append each cluster to null and create a new list */
        LIST(CLUSTER) nodein = tilein->clusterlist;
        CLUSTER clustout = NULL;
        LIST(CLUSTER) listtail = NULL;;
        while (nodein != NULL) {
            clustout = copy_append_CLUSTER(clustout, nodein->elt, colstart, colend);
            if (tileout->clusterlist == NULL) {
                tileout->clusterlist = cons_LIST(CLUSTER)(clustout, tileout->clusterlist);
                listtail = tileout->clusterlist;
            }
            else {
                listtail = rcons_LIST(CLUSTER)(clustout, listtail);
            }
            tileout->ncluster++;
            nodein = nodein->nxt;
            clustout = NULL;
        }    
    }
    else {
        /* append each cluster in list */
        LIST(CLUSTER) nodein = tilein->clusterlist;
        LIST(CLUSTER) nodeout = tileout->clusterlist;
        while ((nodein != NULL) && (nodeout != NULL)) {
            nodeout->elt = copy_append_CLUSTER(nodeout->elt, nodein->elt, colstart, colend);
            nodein = nodein->nxt;
            nodeout = nodeout->nxt;
        }
    }
    tileout->ncycle += colend - colstart + 1;
    
    /* hmhm any error checking required?? */
    return tileout;
}

/**
 * Read a tile from a cif file.
 * Returns a new TILE containing a list of clusters, in the same order as file.
 */
TILE read_cif_TILE(XFILE * fp, unsigned int ncycle) {
    CIFDATA cif = NULL;
    TILE tile = NULL;

    if(NULL==fp) {return NULL;}

    /* read in all the cif data */
    cif = readCIFfromStream(fp);

    if(NULL==cif){
        warnx("Failed to read cif tile.");
    }
    else {
        tile = create_TILE_from_cif(cif, ncycle);
    }

    free_cif(cif);
    return tile;
}

/**
 * Read a tile from a cif run-folder.
 * Returns a new TILE containing a list of clusters, in the same order as files.
 */
TILE read_folder_TILE(const char *root, const unsigned int nlane, const unsigned int ntile, unsigned int ncycle) {
    CIFDATA cif = NULL;
    TILE tile = NULL;

    if (NULL==root) {return NULL;}
    
    cif = readCIFfromDir(root, nlane, ntile, XFILE_RAW);

    if(NULL==cif){
        warnx("Failed to read tile from run-folder; lane number %u tile number %u.", nlane, ntile);
    }
    else {
        xfprintf(xstderr, "Information: Run-folder data found; lane number %u tile number %u\n", nlane, ntile);
        tile = create_TILE_from_cif(cif, ncycle);
        /* store lane and tile since we have them */
        tile->lane = nlane;
        tile->tile = ntile;
    }

    free_cif(cif);
    return tile;
}

/**
 * Read a tile from an Illumina int.txt file.
 * Returns a list of clusters, in reverse order compared to file.
 * \n read_known_TILE is deprecated in favour of read_TILE, which
 * preserves the order of the clusters.
 */
TILE read_known_TILE(XFILE * fp, unsigned int ncycle){
    if(NULL==fp){return NULL;}
    warnx("%s is for demonstration and is not fully functional.",__func__);
    TILE tile = new_TILE();
    CLUSTER cl = NULL;
    while(  cl = read_known_CLUSTER(fp,&ncycle), NULL!=cl ){
        tile->clusterlist = cons_LIST(CLUSTER)(cl,tile->clusterlist);
        tile->ncluster++;
    }
    tile->ncycle = ncycle;
    return tile;
}

/**
 * Read a tile from an Illumina int.txt file.
 * Returns a new TILE containing a list of clusters, in the same order as file.
 * Number of cycles required is specified, or zero means read all.
 * Number read stored in tile structure.
 */
TILE read_TILE(XFILE * fp, unsigned int ncycle){
    TILE tile = NULL;
    CLUSTER cl = NULL;
    LIST(CLUSTER) listtail = NULL;

    if(NULL==fp){return NULL;}
    tile = new_TILE();
    if (NULL==tile) { return NULL;}

    /* Treat first cluster differently to get lane, tile and number of available cycles */
    unsigned int nc = ncycle;
    cl = read_first_CLUSTER(fp, &nc, &(tile->lane), &(tile->tile));
    if ( NULL==cl){ goto cleanup;}

    if (ncycle == 0) {
        /* get all available */
        ncycle = nc;
    }
    if (nc < ncycle) {
        /* not enough cycles, just return the number without bothering to store the data */
        tile->ncycle = nc;
        free_CLUSTER(cl);
    }
    else {
        tile->ncycle = ncycle;
        if (nc > ncycle) {
            /* extra data */
            warnx("Intensity file contains more data than requested: additional %d cycles.", nc - ncycle);
            nc = ncycle;
        }

        tile->clusterlist = cons_LIST(CLUSTER)(cl,tile->clusterlist);
        tile->ncluster++;

        /* Read in remaining clusters, appending to tail of cluster list */
        listtail = tile->clusterlist;
        while(  cl = read_known_CLUSTER(fp, &nc), NULL!=cl ){
            /* check later data */
            if (nc < ncycle) { goto cleanup;}
            listtail = rcons_LIST(CLUSTER)(cl,listtail);
            tile->ncluster++;
        }
    }
    return tile;

cleanup:
    free_LIST(CLUSTER)(tile->clusterlist);
    xfree(tile);
    return NULL;
}

/**
 * Write the lane and tile numbers to file.
 * Separated with a tab as per Illumina int.txt format.
 */
void write_lane_tile(XFILE * fp, const TILE tile) {
    xfprintf(fp, "%u\t%u", tile->lane, tile->tile);
}


#ifdef TEST
#include <err.h>
#include "dirio.h"          // for lanetile

static const unsigned int NXCAT = 4;
static const unsigned int NYCAT = 4;
static int bds[4] = {0,1000,0,1000};


static bool pick_spot( const CLUSTER cl, const void * bds){
    int * ibds = (int *) bds;
    if( cl->x>=ibds[0] && cl->x<ibds[1] && cl->y>=ibds[2] && cl->y<ibds[3] )
        return true;
    return false;
}

static unsigned int whichquad( const CLUSTER cl, const void * dat){
    unsigned int i= (unsigned int)( NXCAT * ((float)cl->x / 1794.0));
    unsigned int j= (unsigned int)( NYCAT * ((float)cl->y / 2048.0));
    return i*NYCAT+j;
}


int main ( int argc, char * argv[]){
    if(argc<3){
        /* arguments are number of cycles, input file (standard format), 
           optional cif input file, and optional cif run-folder and lane tile range */
        errx(EXIT_FAILURE, "Usage: test-tile ncycle _int.txt_filename [cif_filename run-folder lane_tile_range]");
    }

    unsigned int ncycle = 0;
    sscanf(argv[1],"%u", &ncycle);

    XFILE * fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    if (xfisnull(fp)) {
        errx(EXIT_FAILURE, "Failed to open supplied _int.txt file");
    }

    xfputs("Read tile in reverse order\n", xstdout);
    TILE tile1 = read_known_TILE(fp, ncycle);
    xfclose(fp);
    show_TILE(xstdout, tile1, 10);

    xfputs("Reverse list inplace\n", xstdout);
    tile1->clusterlist = reverse_inplace_list_CLUSTER(tile1->clusterlist);
    show_TILE(xstdout, tile1, 10);

    xfputs("Reverse list\n", xstdout);
    LIST(CLUSTER) newrcl = reverse_list_CLUSTER(tile1->clusterlist);
    show_LIST(CLUSTER)(xstdout, newrcl, 10);
    free_LIST(CLUSTER)(newrcl);

    xfputs("Read null file\n", xstdout);
    TILE tile_fwd = read_TILE(NULL, ncycle);
    if (tile_fwd==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        tile_fwd = free_TILE(tile_fwd);
    }

    xfputs("Read all cycles tile in normal order\n", xstdout);
    fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    tile_fwd = read_TILE(fp, 0);
    xfclose(fp);
    if (tile_fwd==NULL) {
        errx(EXIT_FAILURE, "Failed to read supplied _int.txt file");
    }

    xfputs("Write lane and tile: ", xstdout);
    write_lane_tile(xstdout, tile_fwd);
    xfprintf(xstdout, "\n");
    show_TILE(xstdout, tile_fwd, 10);
    tile_fwd = free_TILE(tile_fwd);
    
    xfputs("Read tile in normal order\n", xstdout);
    fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    tile_fwd = read_TILE(fp, ncycle);
    xfclose(fp);
    show_TILE(xstdout, tile_fwd, 10);

    xfputs("Copy null tile\n", xstdout);
    TILE tile2 = copy_TILE(NULL);
    if (tile2==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        tile2 = free_TILE(tile2);
    }
    
    xfputs("Copy tile\n", xstdout);
    tile2 = copy_TILE(tile_fwd);
    show_TILE(xstdout, tile2, 10);
    free_TILE(tile2);

    /* optional cif file testing */
    if (argc > 3) {
        XFILE * fpcif = xfopen(argv[3], XFILE_UNKNOWN, "r");
        if (xfisnull(fpcif)) {
            errx(EXIT_FAILURE, "Failed to open supplied cif file");
        }

        xfputs("Read null cif file\n", xstdout);
        TILE tile_cif = read_cif_TILE(NULL, ncycle);
        if (tile_cif==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            tile_cif = free_TILE(tile_cif);
        }
    
        xfputs("Read not a cif file\n", xstdout);
        fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
        tile_cif = read_cif_TILE(fp, ncycle);
        xfclose(fp);
        if (tile_cif==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            tile_cif = free_TILE(tile_cif);
        }
    
        xfputs("Too many cycles tile from cif file\n", xstdout);
        tile_cif = read_cif_TILE(fpcif, 9999);
        xfclose(fpcif);
        if (tile_cif==NULL) {
            errx(EXIT_FAILURE, "Failed to read supplied cif file");
        }

        xfprintf(xstdout, "Available cycles: %u\n", tile_cif->ncycle);
        tile_cif = free_TILE(tile_cif);
    
        xfputs("Read all cycles tile from cif file\n", xstdout);
        fpcif = xfopen(argv[3], XFILE_UNKNOWN, "r");
        tile_cif = read_cif_TILE(fpcif, 0);
        xfclose(fpcif);
        xfprintf(xstdout, "Available cycles: %u\n", tile_cif->ncycle);
        tile_cif = free_TILE(tile_cif);

        xfputs("Read tile from cif file\n", xstdout);
        fpcif = xfopen(argv[3], XFILE_UNKNOWN, "r");
        tile_cif = read_cif_TILE(fpcif, ncycle);
        xfclose(fpcif);
        show_TILE(xstdout, tile_cif, 10);
        free_TILE(tile_cif);
    }
    else {
        xfputs("Skip cif file tests\n", xstdout);
    }

    /* optional cif run-folder testing */
    if (argc > 4) {
        xfputs("Read null cif run-folder\n", xstdout);
        TILE tile_fol = read_folder_TILE(NULL, 1, 1, ncycle);
        if (tile_fol==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            tile_fol = free_TILE(tile_fol);
        }
    
        xfputs("Read not a cif run-folder\n", xstdout);
        tile_fol = read_folder_TILE("xxxx", 1, 1, ncycle);
        if (tile_fol==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            tile_fol = free_TILE(tile_fol);
        }
    
        xfputs("Invalid lane from cif run-folder (> 9)\n", xstdout);
        tile_fol = read_folder_TILE("xxxx", 10, 1, ncycle);
        if (tile_fol==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            tile_fol = free_TILE(tile_fol);
        }
    
        xfputs("Invalid tile from cif run-folder (> 9999)\n", xstdout);
        tile_fol = read_folder_TILE("xxxx", 1, 10000, ncycle);
        if (tile_fol==NULL) {
            xfputs("Return value null, ok\n", xstdout);
        }
        else {
            xfputs("Return value not null, not ok\n", xstdout);
            tile_fol = free_TILE(tile_fol);
        }

        if (argc > 5) {
            /* need dirio to assume run-folder */
            set_run_folder();

            if (set_pattern(argv[5])) {
                xfputs("Read tiles from cif run-folder\n", xstdout);
                LANETILE lanetile = {0, 0};
                bool more = true;
                while (more) {
                    lanetile = get_next_lanetile();
                    if (lanetile_isnull(lanetile)) {
                        more = false;
                    }
                    else {
                        tile_fol = read_folder_TILE(argv[4], lanetile.lane, lanetile.tile, ncycle);
                        if (tile_fol==NULL) {
                            xfprintf(xstdout, "Lane %u tile %u not available\n", lanetile.lane, lanetile.tile);
                        }
                        else {
                            show_TILE(xstdout, tile_fol, 10);
                            tile_fol = free_TILE(tile_fol);
                        }
                    }
                }
            }
            else {
                xfputs("Invalid run-folder lane tile range supplied\n", xstdout);
            }
        }
        else {
            xfputs("No run-folder lane tile range supplied\n", xstdout);
        }
    }
    else {
        xfputs("Skip cif run-folder tests\n", xstdout);
    }

    fputs("Copy append tile, from second column (1) to ncol/2\n",stdout);
    tile_fwd = copy_append_TILE(tile_fwd, tile1, 1, ncycle/2);    
    show_TILE(xstdout,tile_fwd,10);
    
    fputs("Copy append to tile with no cluster list, from second column (1) to ncol/2\n",stdout);
    free_LIST(CLUSTER)(tile_fwd->clusterlist);
    tile_fwd->clusterlist = NULL;
    tile_fwd->ncluster = 0;
    tile_fwd->ncycle = 0;
    tile_fwd = copy_append_TILE(tile_fwd, tile1, 1, ncycle/2);    
    show_TILE(xstdout,tile_fwd,10);

    fputs("Copy append to null, from second column (1) to ncol/2\n",stdout);
    tile_fwd = free_TILE(tile_fwd);
    tile_fwd = copy_append_TILE(tile_fwd, tile1, 1, ncycle/2);    
    show_TILE(xstdout,tile_fwd,10);

    xfputs("Create an array\n", xstdout);
    unsigned int ncl = 0;
    unsigned int matsize = tile1->clusterlist->elt->signals->nrow * ncycle;
    unsigned int ncluster = tile1->ncluster;
    if (ncluster > 10) {ncluster = 10;}
    int_t arry[matsize * ncluster];

    LIST(CLUSTER) node = tile1->clusterlist;
    while (NULL!=node && ncl < ncluster){
        for (unsigned int idx = 0; idx < matsize; idx++) {
            arry[ncl * matsize + idx] = node->elt->signals->xint[idx];
        }
        node = node->nxt;
        ncl++;
    }

    xfputs("array values:", xstdout);
    for (unsigned int idx = 0; idx < matsize * ncluster; idx++) {
        xfprintf(xstdout, INT_FORMAT, arry[idx]);
    }
    xfputs("\n", xstdout);

    xfputs("Coerce null array\n", xstdout);
    int_t *x;
    TILE tile_ary = coerce_TILE_from_array(ncluster, ncycle, NULL);
    if (tile_ary==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
    }

    xfputs("Coerce tile from array, ignores x,y\n", xstdout);
    x = arry;
    tile_ary = coerce_TILE_from_array(ncluster, ncycle, x);
    show_TILE(xstdout, tile_ary, 10);
    /* do not free tile_ary as it points to arry */
    
    fputs("Filter list\n",stdout);
    LIST(CLUSTER) filteredlist = filter_list_CLUSTER(pick_spot,tile1->clusterlist,(void *)bds);
    show_LIST(CLUSTER)(xstdout,filteredlist,10);
    shallow_free_list_CLUSTER(filteredlist);

    LIST(CLUSTER) * spots = split_list_CLUSTER(whichquad,tile1->clusterlist,NXCAT*NYCAT,NULL);
    for ( int i=0 ; i<NXCAT ; i++){
        for ( int j=0 ; j<NYCAT ; j++){
            fprintf(stdout,"(%d,%d) has %u elements\n",i+1,j+1,length_LIST(CLUSTER)(spots[i*NYCAT+j]));
            //show_LIST(CLUSTER)(xstdout,spots[i*NYCAT+j],10);
            shallow_free_list_CLUSTER(spots[i*NYCAT+j]);
        }
    }

    free_TILE(tile1);
    free_TILE(tile_fwd);
    return EXIT_SUCCESS;
}

#endif

