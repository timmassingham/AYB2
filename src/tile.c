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
 * Read a tile from an Illumina int.txt file.
 * Returns a list of clusters, in reverse order compared to file.
 * \n read_known_TILE is deprecated in favour of read_TILE, which
 * preserves the order of the clusters.
 */
/*hmhm*/
//TILE read_known_TILE( XFILE * fp, unsigned int ncycle){
TILE read_known_TILE( XFILE * fp, unsigned int *ncycle){
    if(NULL==fp){return NULL;}
    warnx("%s is for demonstration and is not fully functional.",__func__);
    TILE tile = new_TILE();
    CLUSTER cl = NULL;
    while(  cl = read_known_CLUSTER(fp,ncycle), NULL!=cl ){
        tile->clusterlist = cons_LIST(CLUSTER)(cl,tile->clusterlist);
        tile->ncluster++;
    }
    return tile;
}

/**
 * Read a tile from an Illumina int.txt file.
 * Returns a new TILE containing a list of clusters, in the same order as file.
 */
/*hmhm*/
//TILE read_TILE( XFILE * fp, unsigned int ncycle){
TILE read_TILE( XFILE * fp, unsigned int *ncycle){
    if(NULL==fp){return NULL;}
    warnx("%s is for demonstration and is not fully functional.",__func__);
    TILE tile = new_TILE();

    // Treat first cluster differently
    CLUSTER cl = read_known_CLUSTER(fp,ncycle);
    if ( NULL!=cl){
        tile->clusterlist = cons_LIST(CLUSTER)(cl,tile->clusterlist);
        tile->ncluster++;
    }
    // Read in remaining clusters, appending to tail of cluster list
    LIST(CLUSTER) listtail = tile->clusterlist;
    while(  cl = read_known_CLUSTER(fp,ncycle), NULL!=cl ){
        listtail = rcons_LIST(CLUSTER)(cl,listtail);
        tile->ncluster++;
    }
    return tile;
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
    TILE tile = read_known_TILE(fp,&ncycle);
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
    TILE tile_fwd = read_TILE(fp,&ncycle);
    xfclose(fp);
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

