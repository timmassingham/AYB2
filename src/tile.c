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
/* None */

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
TILE coerce_TILE_from_array(unsigned int ncluster, unsigned int ncycle, real_t * x){
    TILE tile = NULL;
    CLUSTER cl = NULL;
    LIST(CLUSTER) listtail = NULL;

    if(NULL==x){ return NULL;}
    tile = new_TILE();
    if(NULL==tile){ return NULL;}
    tile->ncluster = ncluster;
    tile->ncycle = ncycle;
    real_t * next_x;

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
TILE read_cif_TILE (XFILE * fp, unsigned int ncycle) {
    CIFDATA cif = NULL;
    CLUSTER cl = NULL;
    LIST(CLUSTER) listtail = NULL;
    unsigned int cifcycle=0, cifcluster=0;

    if(NULL==fp) {return NULL;}
    TILE tile = new_TILE();
    if (NULL==tile) {return NULL;}

    /* read in all the cif data */
    cif = NULL;
    cif = readCIFfromStream(fp);
    if(NULL==cif){
        warnx("Failed to read tile.");
        goto cleanup;
    }

    cifcycle = cif_get_ncycle(cif);
    cifcluster = cif_get_ncluster(cif);
    xfprintf(xstderr, "Read cif tile: %u cycles from %u clusters.\n", cifcycle, cifcluster);

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

    free_cif(cif);
    return tile;

cleanup:
    free_cif(cif);
    free_LIST(CLUSTER)(tile->clusterlist);
    xfree(tile);
    return NULL;
}

/**
 * Read a tile from an Illumina int.txt file.
 * Returns a list of clusters, in reverse order compared to file.
 * \n read_known_TILE is deprecated in favour of read_TILE, which
 * preserves the order of the clusters.
 */
TILE read_known_TILE( XFILE * fp, unsigned int ncycle){
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
 * Number of cycles required is specified.
 * Number read stored in tile structure.
 */
TILE read_TILE( XFILE * fp, unsigned int ncycle){
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
    if (nc > ncycle) {
        /* extra data */
        warnx("Intensity file contains more data than requested: additional %d cycles.", nc - ncycle);
        /* ignore extra for rest of read */
        nc = ncycle;
    }
    else {
        /* set to read what is available */
        ncycle = nc;
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
    tile->ncycle = ncycle;
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
void write_lane_tile (XFILE * fp, const TILE tile) {
    xfprintf(fp, "%u\t%u", tile->lane, tile->tile);
}


#ifdef TEST
//#include <stdio.h>
//#include <err.h>

const unsigned int nxcat = 4;
const unsigned int nycat = 4;


bool pick_spot( const CLUSTER cl, const void * bds){
    int * ibds = (int *) bds;
    if( cl->x>=ibds[0] && cl->x<ibds[1] && cl->y>=ibds[2] && cl->y<ibds[3] )
        return true;
    return false;
}

unsigned int whichquad( const CLUSTER cl, const void * dat){
    unsigned int i= (unsigned int)( nxcat * ((float)cl->x / 1794.0));
    unsigned int j= (unsigned int)( nycat * ((float)cl->y / 2048.0));
    return i*nycat+j;
}


int main ( int argc, char * argv[]){
    int bds[4] = {0,1000,0,1000};
    if(argc!=3){
        xfputs("Usage: test ncycle filename\n", xstdout);
        return EXIT_FAILURE;
    }

    unsigned int ncycle = 0;
    sscanf(argv[1],"%u", &ncycle);

    XFILE * fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    TILE tile = read_known_TILE(fp, ncycle);
    xfclose(fp);
    show_TILE(xstdout, tile, 10);

    xfputs("Reversing list inplace\n", xstdout);
    tile->clusterlist = reverse_inplace_list_CLUSTER(tile->clusterlist);
    show_TILE(xstdout,tile,10);

    xfputs("Reversing list\n",xstdout);
    LIST(CLUSTER) newrcl = reverse_list_CLUSTER(tile->clusterlist);
    show_LIST(CLUSTER)(xstdout,newrcl,10);
    free_LIST(CLUSTER)(newrcl);

    xfputs("Reading tile in normal order\n",xstdout);
    fp = xfopen(argv[2],XFILE_UNKNOWN,"r");
    TILE tile_fwd = read_TILE(fp,ncycle);
    xfclose(fp);
    show_TILE(xstdout,tile_fwd,10);
    
    fputs("Appending tile, from second column to half way\n",stdout);
    tile_fwd = copy_append_TILE(tile_fwd, tile, 1, ncycle/2);    
    show_TILE(xstdout,tile_fwd,10);
    
    fputs("Appending to tile with no cluster list, from second column to half way\n",stdout);
    free_LIST(CLUSTER)(tile_fwd->clusterlist);
    tile_fwd->clusterlist = NULL;
    tile_fwd->ncluster = 0;
    tile_fwd->ncycle = 0;
    tile_fwd = copy_append_TILE(tile_fwd, tile, 1, ncycle/2);    
    show_TILE(xstdout,tile_fwd,10);

    fputs("Appending to null, from second column to half way\n",stdout);
    tile_fwd = free_TILE(tile_fwd);
    tile_fwd = copy_append_TILE(tile_fwd, tile, 1, ncycle/2);    
    show_TILE(xstdout,tile_fwd,10);
    free_TILE(tile_fwd);

    xfputs("Create an array\n", xstdout);
    unsigned int ncl=0;
    unsigned int matsize = tile->clusterlist->elt->signals->nrow * ncycle;
    unsigned int ncluster = tile->ncluster;
    real_t arry[matsize * ncluster];

    LIST(CLUSTER) node = tile->clusterlist;
    while (NULL!=node && ncl < ncluster){
        for (unsigned int idx = 0; idx < matsize; idx++) {
            arry[ncl * matsize + idx] = node->elt->signals->x[idx];
        }
        node = node->nxt;
        ncl++;
    }

    xfputs("array values:", xstdout);
    for (unsigned int idx = 0; idx < matsize * ncluster; idx++) {
        xfprintf(xstdout, " %#8.2f", arry[idx]);
    }
    xfputs("\n", xstdout);

    xfputs("Coerce from array, ignores x,y\n", xstdout);
    real_t *x = arry;
    TILE tile_ary = coerce_TILE_from_array(ncluster, ncycle, x);
    show_TILE(xstdout, tile_ary, 10);
    /* do not free tile_ary as it points to arry */

    fputs("Filtering list\n",stdout);
    LIST(CLUSTER) filteredlist = filter_list_CLUSTER(pick_spot,tile->clusterlist,(void *)bds);
    show_LIST(CLUSTER)(xstdout,filteredlist,10);
    shallow_free_list_CLUSTER(filteredlist);

    LIST(CLUSTER) * spots = split_list_CLUSTER(whichquad,tile->clusterlist,nxcat*nycat,NULL);
    for ( int i=0 ; i<nxcat ; i++){
        for ( int j=0 ; j<nycat ; j++){
            fprintf(stdout,"(%d,%d) has %u elements\n",i+1,j+1,length_LIST(CLUSTER)(spots[i*nycat+j]));
            //show_LIST(CLUSTER)(xstdout,spots[i*nycat+j],10);
            shallow_free_list_CLUSTER(spots[i*nycat+j]);
        }
    }

    free_TILE(tile);
    return EXIT_SUCCESS;
}

#endif

