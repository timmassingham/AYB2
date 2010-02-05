/*
 *  Copyright (C) 2010 by Tim Massingham
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

#include <stdlib.h>
#include "tile.h"


/*
 * Standard functions for tile structure
 */


TILE new_TILE(void){
	TILE tile = calloc(1,sizeof(*tile));
	return tile;
}

void free_TILE(TILE tile){
	if(NULL==tile){ return;}
	free_LIST(CLUSTER)(tile->clusterlist);
	xfree(tile);
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
		show_CLUSTER(fp,node->elt);
		node = node->nxt;
		ncl++;
	}
}


/*
 * Read a tile from a Illumina int.txt file
 */
TILE read_known_TILE( XFILE * fp, unsigned int ncycle){
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


#ifdef TEST
#include <stdio.h>
#include <err.h>
int main ( int argc, char * argv[]){
	if(argc!=3){
		fputs("Usage: test ncycle filename\n",stdout);
		return EXIT_FAILURE;
	}
	
	unsigned int ncycle = 0;
	sscanf(argv[1],"%u",&ncycle);

	XFILE * fp = xfopen(argv[2],XFILE_UNKNOWN,"r");
	TILE tile = read_known_TILE(fp,ncycle);
	show_TILE(xstdout,tile,10);
	free_TILE(tile);
	xfclose(fp);
	return EXIT_SUCCESS;
}

#endif


