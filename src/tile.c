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
#include "tile.h"


/* constants */
/* None      */

/* members */
/* below */


/* private functions */
/* undetermined */

/* public functions */
/* undetermined */


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
    while(  cl = read_known_CLUSTER(fp,&ncycle, false), NULL!=cl ){
        tile->clusterlist = cons_LIST(CLUSTER)(cl,tile->clusterlist);
        tile->ncluster++;
    }
    tile->ncycle = ncycle;
    return tile;
}

/**
 * Read a tile from an Illumina int.txt file.
 * Returns a new TILE containing a list of clusters, in the same order as file.
 */
TILE read_TILE( XFILE * fp, unsigned int ncycle){
    if(NULL==fp){return NULL;}
    warnx("%s is for demonstration and is not fully functional.",__func__);
    TILE tile = new_TILE();
    if (NULL==tile) { return NULL;}

    // Treat first cluster differently
    unsigned int nc = ncycle;
    CLUSTER cl = read_known_CLUSTER(fp, &nc, true);
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

    // Read in remaining clusters, appending to tail of cluster list
    LIST(CLUSTER) listtail = tile->clusterlist;
    while(  cl = read_known_CLUSTER(fp, &nc, false), NULL!=cl ){
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


#ifdef TEST
#include <stdio.h>
#include <err.h>

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
        fputs("Usage: test ncycle filename\n",stdout);
        return EXIT_FAILURE;
    }

    unsigned int ncycle = 0;
    sscanf(argv[1],"%u",&ncycle);

    XFILE * fp = xfopen(argv[2],XFILE_UNKNOWN,"r");
    TILE tile = read_known_TILE(fp,ncycle);
    xfclose(fp);
    show_TILE(xstdout,tile,10);

    fputs("Reversing list inplace\n",stdout);
        tile->clusterlist = reverse_inplace_list_CLUSTER(tile->clusterlist);
    show_TILE(xstdout,tile,10);

    fputs("Reversing list\n",stdout);
    LIST(CLUSTER) newrcl = reverse_list_CLUSTER(tile->clusterlist);
    show_LIST(CLUSTER)(xstdout,newrcl,10);
    free_LIST(CLUSTER)(newrcl);

    fputs("Reading tile in normal order\n",stdout);
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

