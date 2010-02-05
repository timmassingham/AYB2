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

#include "cluster.h"
#include "utility.h"


/*
 * Standard functions for cluster structure
 */
CLUSTER new_CLUSTER ( void){
	CLUSTER cl = calloc(1,sizeof(*cl));
	return cl;
}

void free_CLUSTER(CLUSTER cluster){
	if(NULL==cluster){return;}
	free_MAT(cluster->signals);
	xfree(cluster);
}
	
CLUSTER copy_CLUSTER(const CLUSTER cluster){
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
	show_MAT(fp,cluster->signals,4,5);
}

/* 
 * Basic input functions from file
 */
/*  Reads a cluster from file pointer, assuming that the number of cycles is known.
 * Format should be that of Illumina's _int.txt
 * A little validation of the file format is performed.
 */
CLUSTER read_known_CLUSTER( XFILE * fp, const unsigned int ncycle){
    char * line = NULL;
    MAT signals = NULL;
    CLUSTER cluster = NULL;
    
	if(NULL==fp){return NULL;}
	// Get line from file
	size_t len = 0;
	line = xfgetln(fp,&len);
	char * ptr = line;
	
	/* Read lane, tile and coordinate information */
	unsigned long int lane,tile,x,y;
	lane = strtoul(ptr,&ptr,0);
	if('\t'!=ptr[0]){ goto cleanup;}
	tile = strtoul(ptr,&ptr,0);
	if('\t'!=ptr[0]){ goto cleanup;}
	x = strtoul(ptr,&ptr,0);
	if('\t'!=ptr[0]){ goto cleanup;}
	y = strtoul(ptr,&ptr,0);
	/* Read signals, tab separated quads of numbers. */
	signals = new_MAT(NBASE,ncycle);
	if(NULL==signals){ goto cleanup;}
	for ( unsigned int i=0 ; i<ncycle ; i++){
		if('\t'!=ptr[0]){ goto cleanup;}
		signals->x[i*NBASE+0] = strtor(ptr,&ptr);
		signals->x[i*NBASE+1] = strtor(ptr,&ptr);
		signals->x[i*NBASE+2] = strtor(ptr,&ptr);
		signals->x[i*NBASE+3] = strtor(ptr,&ptr);
	}	

	cluster = new_CLUSTER();
	if(NULL==cluster){goto cleanup;}
	cluster->x = x;
	cluster->y = y;
	cluster->signals = signals;

	xfree(line);	
	return cluster;

cleanup:
	free_CLUSTER(cluster);
    free_MAT(signals);
	xfree(line);
	return NULL;
}

// Stub function
CLUSTER read_unknown_CLUSTER( XFILE * fp){
	return NULL;
}

#ifdef TEST
#include <stdio.h>
int main(int argc , char * argv[]){
	if(argc!=3){
		fputs("Usage: test ncycle filename\n",stdout);
		return EXIT_FAILURE;
	}
	
	unsigned int ncycle = 0;
	sscanf(argv[1],"%u",&ncycle);

	XFILE * fp = xfopen(argv[2],XFILE_UNKNOWN,"r");
	CLUSTER cl = read_known_CLUSTER(fp,ncycle);
	show_CLUSTER(xstdout,cl);
	free_CLUSTER(cl);
	xfclose(fp);
	return EXIT_SUCCESS;
}
#endif

