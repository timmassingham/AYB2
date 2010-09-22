/** 
 * \file cif.c
 * Interface to CIF format intensity files.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden

 *  Copyright (C) 2008,2009 by Tim Massingham
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the ciftool base-calling software.
 *
 *  ciftool is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ciftool is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ciftool.  If not, see <http://www.gnu.org/licenses/>.
 */

/*  Code to read and write CIF version 1
 *  Integers in file are little endian, this is assumed in code (eg. intel architecture)
 *
 *  Format:
 *  "CIF" version(u1) datasize(u1) firstcycle(u2) #cycles(2) #clusters(4)
 *  Repeat: cycle, channel, cluster
 *
 * version: currently 1
 * datasize: number of bytes used for floats
 * firstcycle: offset for cycles
 * #cycles: number of cycles
 * #clusters: number of clusters
 * floats are signed and truncated to nearest integer, then rounded into range
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <glob.h>
#include <math.h>
#include "cif.h"

#define NCHANNEL    4

struct cifData {
    uint8_t version, datasize;
    uint16_t firstcycle, ncycle;
    uint32_t ncluster;
    encInt intensity;
};


/* constants */
/* none      */

/* members */
/* none */


/* private functions */
/* undetermined */

/* public functions */
/* undetermined */

// Accessor functions
uint8_t cif_get_version ( const CIFDATA cif ){ return cif->version;}
uint8_t cif_get_datasize ( const CIFDATA cif ){ return cif->datasize;}
uint16_t cif_get_firstcycle ( const CIFDATA cif ){ return cif->firstcycle;}
uint16_t cif_get_ncycle ( const CIFDATA cif ){ return cif->ncycle;}
uint32_t cif_get_ncluster ( const CIFDATA cif ){ return cif->ncluster;}
encInt cif_get_const_intensities ( const CIFDATA cif ){ return cif->intensity; }

/** Return a single channel value as a real. */
real_t cif_get_real (const CIFDATA cif, const uint32_t cl, const uint32_t base, const uint32_t cy) {

    if (NULL==cif) {return NAN;}
    if (NULL==cif->intensity.i8) {return NAN;}
    if ((cl >= cif->ncluster) || (base >= NCHANNEL)|| (cy >= cif->ncycle) ) {return NAN;}

    switch(cif->datasize) {
        case 1: return (real_t)cif->intensity.i8 [(cy * NCHANNEL + base) * cif->ncluster + cl]; break;
        case 2: return (real_t)cif->intensity.i16[(cy * NCHANNEL + base) * cif->ncluster + cl]; break;
        case 4: return (real_t)cif->intensity.i32[(cy * NCHANNEL + base) * cif->ncluster + cl]; break;
        default: return NAN;
    }
}

/** Set a single channel value from a real. */
void cif_set_from_real (CIFDATA cif, const uint32_t cl, const uint32_t base, const uint32_t cy, real_t x) {

    if (NULL==cif) {return;}
    if ((cl >= cif->ncluster) || (base >= NCHANNEL)|| (cy >= cif->ncycle) ) {return;}
    if (NULL==cif->intensity.i8) {
        cif->intensity.i8 = calloc(cif->ncycle * cif->ncluster * NCHANNEL, cif->datasize);
        if (NULL==cif->intensity.i8) {return;}
    }

    /* round properly first */
    x = roundr(x);

    switch(cif->datasize) {
        case 1:
            if ((x >= INT8_MIN) && (x <= INT8_MAX)) {
                cif->intensity.i8 [(cy * NCHANNEL + base) * cif->ncluster + cl] = (int8_t)x; break;
            }
        case 2:
            if ((x >= INT16_MIN) && (x <= INT16_MAX)) {
                cif->intensity.i16[(cy * NCHANNEL + base) * cif->ncluster + cl] = (int16_t)x; break;
            }
        case 4:
            if ((x >= INT32_MIN) && (x <= INT32_MAX)) {
                cif->intensity.i32[(cy * NCHANNEL + base) * cif->ncluster + cl] = (int32_t)x; break;
            }
        default: ;
    }
}

CIFDATA new_cif ( void ){
   CIFDATA cif = malloc(sizeof(*cif));
   if(NULL==cif){ return NULL; }
   cif->version = 1;
   cif->datasize = 2;
   cif->firstcycle = 1;
   cif->ncycle = 0;
   cif->ncluster = 0;
   cif->intensity.i8 = NULL;
   return cif;
}

/** Create a new empty cif structure ready for data of a known size. */
CIFDATA create_cif (const uint32_t ncycle, const uint32_t ncluster) {
    CIFDATA cif = new_cif();
    cif->ncycle = ncycle;
    cif->ncluster = ncluster;
    return cif;
}

// Delete
void free_cif ( CIFDATA cif ){
    if ( NULL==cif ) return;
    if ( NULL!=cif->intensity.i8){ free(cif->intensity.i8); }
    free(cif);
}


bool __attribute__((const)) isCifAllowedDatasize ( const uint8_t datasize );
CIFDATA readCifHeader (XFILE * ayb_fp);
encInt readCifIntensities ( XFILE * ayb_fp , const CIFDATA const header, encInt intensties );
encInt readEncodedFloats ( XFILE  * ayb_fp, const uint32_t nfloat, const uint8_t nbyte, encInt  tmp_mem );
bool writeCifHeader ( XFILE * ayb_fp, const CIFDATA const header);
bool writeCifIntensities ( XFILE * ayb_fp , const CIFDATA const header,
                           const encInt intensities );
bool writeEncodedFloats ( XFILE * ayb_fp , const uint32_t nfloat , const uint8_t nbyte,
                          const encInt floats );


bool __attribute__((const)) isCifAllowedDatasize ( const uint8_t datasize ){
    if ( 1==datasize || 2==datasize || 4==datasize ){ return true;}
    return false;
}

/*  Routines to convert between CIF format values (cycle,channel,cluster)
 * and AYB format (cluster, cycle, channel)
 *  If we reshape as a ncluster x (nchannel*ncycle) matrix, CIF is the transpose
 * of AYB
 */
/*INTENSITY_FLOAT_TYPE * transpose ( const INTENSITY_FLOAT_TYPE * restrict m, const uint32_t nr, const uint32_t nc, INTENSITY_FLOAT_TYPE * restrict intensities){
    if (NULL==intensities){ intensities = malloc( nr*nc*sizeof(INTENSITY_FLOAT_TYPE) ); }
    if (NULL==intensities){ return NULL;}
    for ( uint32_t c=0 ; c<nc ; c++){
        for ( uint32_t r=0 ; r<nr ; r++){
            intensities[r*nc+c] = m[c*nr+r];
        }
    }
    return intensities;
}
INTENSITY_FLOAT_TYPE * cif2ayb ( const INTENSITY_FLOAT_TYPE * restrict cif_ints, const uint32_t ncluster, const uint32_t ncycle, INTENSITY_FLOAT_TYPE * restrict intensities){
    assert(NULL!=cif_ints);
    return transpose ( cif_ints, ncluster, ncycle*NCHANNEL, intensities );
}
INTENSITY_FLOAT_TYPE * ayb2cif ( const INTENSITY_FLOAT_TYPE * restrict ayb_ints, const uint32_t ncluster, const uint32_t ncycle, INTENSITY_FLOAT_TYPE * restrict intensities){
    assert(NULL!=ayb_ints);
    return transpose ( ayb_ints, ncycle*NCHANNEL, ncluster, intensities );
}*/


/* Read and check header */
CIFDATA readCifHeader (XFILE * ayb_fp){
    CIFDATA header = malloc(sizeof(struct cifData));
    if(NULL==header){return NULL;}
    char str[4] = "\0\0\0\0";
    const char * hd_str = "CIF";

    // Check if this is a CIF format file
    // "Magic number" is three bytes long
    xfread(&str,3,1,ayb_fp);
    // and should be the string "CIF"
    if ( strncmp(str,hd_str,3) ){
        free(header);
        return NULL;
    }
    // Rest of header.
    xfread(&header->version,   1,1,ayb_fp);
    xfread(&header->datasize,  1,1,ayb_fp);
    xfread(&header->firstcycle,2,1,ayb_fp);
    xfread(&header->ncycle,    2,1,ayb_fp);
    xfread(&header->ncluster,  4,1,ayb_fp);

    header->intensity.i32 = NULL;

    assert( 1==header->version);
    assert( isCifAllowedDatasize(header->datasize) );

    return header;
}

/* Read intensities */
encInt readCifIntensities ( XFILE * ayb_fp , const CIFDATA header, encInt  intensities ){
    const uint32_t size = NCHANNEL * header->ncluster * header->ncycle;
    intensities = readEncodedFloats(ayb_fp,size,header->datasize,intensities);
    return intensities;
}

/* Read an array of encoded floats and return array */
encInt readEncodedFloats ( XFILE  * ayb_fp, const uint32_t nfloat, const uint8_t nbyte, encInt intensities ){
    assert(1==nbyte || 2==nbyte || 4==nbyte);
    assert(nfloat>0);

    if ( NULL==intensities.i32 ){
        intensities.i8 = calloc(nfloat,nbyte);
    }
    xfread(intensities.i8,nbyte,nfloat,ayb_fp);

    return intensities;
}


/* Write header for CIF file */
bool writeCifHeader ( XFILE * ayb_fp, const CIFDATA const header){
    assert( isCifAllowedDatasize(header->datasize) );
    assert( 1==header->version );

    xfputs("CIF",ayb_fp);
    xfwrite(&header->version,   1,1,ayb_fp);
    xfwrite(&header->datasize,  1,1,ayb_fp);
    xfwrite(&header->firstcycle,2,1,ayb_fp);
    xfwrite(&header->ncycle,    2,1,ayb_fp);
    xfwrite(&header->ncluster,  4,1,ayb_fp);

    return true;
}

/*  Write intensities in encoded format
 *  Balance between speed and memory consumption could be changed by writing
 * out in more chunks.
 */
bool writeCifIntensities ( XFILE * ayb_fp , const CIFDATA  header,
                           const encInt  intensities ){
    const uint32_t nfloat = NCHANNEL * header->ncycle * header->ncluster;
    return writeEncodedFloats ( ayb_fp, nfloat, header->datasize, intensities );
}


int32_t __attribute__((const)) getMax (const uint8_t nbyte){
    switch(nbyte){
        case 1: return INT8_MAX;
        case 2: return INT16_MAX;
        case 4: return INT32_MAX;
        default:
            err(EINVAL,"Unimplemented bytesize %u in %s:%d\n",nbyte,__FILE__,__LINE__);
    }
    return 0.; // Never reach here
}

int32_t __attribute__((const)) getMin (const uint8_t nbyte){
    switch(nbyte){
        case 1: return INT8_MIN;
        case 2: return INT16_MIN;
        case 4: return INT32_MIN;
        default:
            err(EINVAL,"Unimplemented bytesize %u in %s:%d\n",nbyte,__FILE__,__LINE__);
    }
    return 0.; // Never reach here
}

/*void * encodeFloats ( const void * restrict floats, const uint32_t nfloat , const uint8_t nbyte){
    assert(NULL!=floats);
    assert(isCifAllowedDatasize(nbyte));

    int32_t enc_max = getMax(nbyte);
    int32_t enc_min = getMin(nbyte);
    int8_t * mem8=NULL;
    int16_t * mem16=NULL;
    int32_t * mem32=NULL;
    void * ret_mem;
    INTENSITY_FLOAT_TYPE f;
    switch(nbyte){
        case 1:
            mem8 = malloc(nfloat*sizeof(int8_t));
            for ( uint32_t i=0; i<nfloat ; i++){
                f = floats[i];
                f = (f<enc_max)?f:enc_max;
                f = (f>enc_min)?f:enc_min;
                mem8[i] = (int8_t) f;
            }
            ret_mem = (void *) mem8;
            break;
        case 2:
            mem16 = malloc(nfloat*sizeof(int16_t));
            for ( uint32_t i=0; i<nfloat ; i++){
                f = floats[i];
                f = (f<enc_max)?f:enc_max;
                f = (f>enc_min)?f:enc_min;
                mem16[i] = (int16_t) f;
            }
            ret_mem = (void *) mem16;
            break;
        case 4:
            mem32 = malloc(nfloat*sizeof(int32_t));
            for ( uint32_t i=0; i<nfloat ; i++){
                f = floats[i];
                f = (f<enc_max)?f:enc_max;
                f = (f>enc_min)?f:enc_min;
                mem32[i] = (int32_t) f;
            }
            ret_mem = (void *) mem32;
            break;
        default:
            err(EINVAL,"Unimplemented bytesize %u in %s:%d\n",nbyte,__FILE__,__LINE__);
    }

    return ret_mem;
}*/

/*  Write floats in encoded format
 */
bool writeEncodedFloats ( XFILE * ayb_fp , const uint32_t nfloat , const uint8_t nbyte,
                          const encInt  intensities ){
    xfwrite( intensities.i8, nbyte, nfloat, ayb_fp);
    return true;
}



/*  Routines to read and write CIF files  */

bool write2CIFstream ( XFILE * ayb_fp, const encInt  intensities, const uint16_t firstcycle, const uint32_t ncycle, const uint32_t ncluster, const uint8_t nbyte){
    assert(NULL!=ayb_fp);
    assert(NULL!=intensities.i8);
    assert( isCifAllowedDatasize(nbyte) );

    struct cifData header = { 1, nbyte, firstcycle, ncycle, ncluster, {.i8=NULL} };
    writeCifHeader(ayb_fp,&header);
    writeEncodedFloats(ayb_fp,ncluster*ncycle*NCHANNEL,nbyte,intensities);

    return true;
}


bool write2CIFfile ( const char * fn, const XFILE_MODE mode, const encInt  intensities, const uint16_t firstcycle, const uint32_t ncycle, const uint32_t ncluster, const uint8_t nbyte){
    if(NULL==fn){ return false;}
    XFILE * ayb_fp = xfopen(fn,mode,"wb");
    if ( NULL==ayb_fp){ return false;}
    bool ret = write2CIFstream(ayb_fp,intensities,firstcycle,ncycle,ncluster,nbyte);
    xfclose(ayb_fp);
    return ret;
}

bool writeCIFtoFile ( const CIFDATA const cif, const char * fn, const XFILE_MODE mode){
    if(NULL==cif){ return false;}
    if(NULL==fn){ return false;}
    return write2CIFfile(fn,mode,cif->intensity,cif->firstcycle,cif->ncycle,cif->ncluster,cif->datasize);
}

bool writeCIFtoStream ( const CIFDATA  cif, XFILE * ayb_fp){
    if(NULL==cif){ return false;}
    if(NULL==ayb_fp){ return false;}
    bool ret = write2CIFstream(ayb_fp,cif->intensity,cif->firstcycle,cif->ncycle,cif->ncluster,cif->datasize);
    return ret;
}


CIFDATA readCIFfromStream ( XFILE * ayb_fp ){
    CIFDATA cif = readCifHeader(ayb_fp);
    encInt e = {.i8=NULL};
    cif->intensity = readCifIntensities ( ayb_fp , cif, e );
    return cif;
}

CIFDATA readCIFfromFile ( const char * fn, const XFILE_MODE mode){
    XFILE * ayb_fp = xfopen(fn,mode,"rb");
    if ( NULL==ayb_fp){ return NULL;}
    CIFDATA cif = readCIFfromStream(ayb_fp);
    xfclose(ayb_fp);
    return cif;
}

bool consistent_cif_headers( const CIFDATA cif1, const CIFDATA cif2 ){
   if ( cif1->version  != cif2->version  ){ return false; } // Check might be relaxed later
   if ( cif1->datasize != cif2->datasize ){ return false; } // Check might be relaxed later
   if ( cif1->ncluster != cif2->ncluster ){ return false; }
   return true;
}

CIFDATA cif_add_file( const char * fn, const XFILE_MODE mode, CIFDATA cif ){
   if ( NULL==fn){ goto cif_add_error;}
   XFILE * ayb_fp = xfopen(fn,mode,"rb");
   if ( NULL==ayb_fp){ goto cif_add_error;}

   CIFDATA newheader = readCifHeader(ayb_fp);
   if ( NULL==newheader ){ goto cif_add_error;}
   if ( NULL==cif->intensity.i8 ){
      cif->ncluster = newheader->ncluster;
      // First file read. Allocated memory needed
      cif->intensity.i8 = calloc(cif->ncycle*cif->ncluster*NCHANNEL,cif->datasize);
      if ( NULL==cif->intensity.i8 ){ goto cif_add_error;}
   }
   if ( ! consistent_cif_headers(cif,newheader) ){ goto cif_add_error;}
   const uint32_t offset = (newheader->firstcycle - 1) * cif->ncluster * NCHANNEL;
   encInt mem = {.i8=NULL};
   switch(cif->datasize){
       case 1: mem.i8 = cif->intensity.i8 + offset; break;
       case 2: mem.i16 = cif->intensity.i16 + offset; break;
       case 3: mem.i32 = cif->intensity.i32 + offset; break;
       default: errx(EXIT_FAILURE,"Incorrect datasize in %s (%s:%d)\n",__func__,__FILE__,__LINE__);
   }
   readCifIntensities(ayb_fp,newheader,mem);

   free_cif(newheader);
   return cif;

cif_add_error:
   free_cif(newheader);
   free_cif(cif);
   return NULL;
}


char * cif_create_cifglob ( const char * root, const uint32_t lane, const uint32_t tile ){
   if(NULL==root){ return NULL;}
   if(lane>9){
      warn("Assumption that lane numbering is less than 10 violated (asked for %u).\n",lane);
      return NULL;
   }
   if(tile>9999){
      warn("Assumption that tile numbering is less than 9999 violated (asked for %u).\n",tile);
      return NULL;
   }

   char * cif_glob = calloc(strlen(root)+41,sizeof(char));
   if(NULL==cif_glob){
      return NULL;
   }
   strcpy(cif_glob,root);
   uint32_t offset = strlen(cif_glob);
   strcpy(cif_glob+offset,"/Data/Intensities/L00");
   offset += 21;
   cif_glob[offset] = lane + 48;
   offset++;
   strcpy(cif_glob+offset,"/C*.1/s_X_");
   cif_glob[offset+8] = lane + 48;
   offset += 10;
   offset += sprintf(cif_glob+offset,"%u",tile);
   strcpy(cif_glob+offset,".cif");

   return cif_glob;
}

CIFDATA spliceCIF(const CIFDATA cif, uint32_t ncycle, uint32_t offset){
    if(NULL==cif){return NULL;}
    if(offset+ncycle>cif->ncycle){ return NULL;}

    CIFDATA newcif = new_cif();
    memcpy(newcif,cif,sizeof(*newcif));
    newcif->firstcycle = 1;
    newcif->ncycle = ncycle;

    const uint32_t nobs = NCHANNEL*newcif->ncluster*newcif->ncycle;
    newcif->intensity.i8 = calloc(nobs,newcif->datasize);
    if(NULL==newcif->intensity.i8){ goto clean;}

    uint32_t offset8 = offset*NCHANNEL*newcif->ncluster*cif->datasize;
    memcpy(newcif->intensity.i8,cif->intensity.i8+offset8,nobs*newcif->datasize);

    return newcif;

clean:
    free(newcif);
    return NULL;
}

/* Read an entire run from a run directory */
CIFDATA readCIFfromDir ( const char * root, const uint32_t lane, const uint32_t tile, const XFILE_MODE mode){
   CIFDATA cif = NULL;

   if(NULL==root){ return NULL;}
   if(lane>9){
      warn("Assumption that lane numbering is less than 10 violated (asked for %u).\n",lane);
      return NULL;
   }
   if(tile>9999){
      warn("Assumption that tile numbering is less than 9999 violated (asked for %u).\n",tile);
      return NULL;
   }

   // Find matching files
   glob_t g;
   char * cif_glob = cif_create_cifglob(root,lane,tile);
   int ret = glob(cif_glob,0,NULL,&g);
   if(0!=ret){ goto readCIF_error; }
   free(cif_glob);

   const uint32_t ncycle = g.gl_pathc;

   cif = new_cif();
   cif->ncycle = ncycle;
   for ( uint32_t i=0 ; i<ncycle ; i++){
      cif = cif_add_file(g.gl_pathv[i],XFILE_RAW,cif);
      if(NULL==cif){ fprintf(stderr,"Problem reading CIF \"%s\"\n",g.gl_pathv[i]); }
   }

   globfree(&g);
   return cif;

readCIF_error:
   globfree(&g);
   return cif;
}



/*  Print routine for CIF structure */
void showCIF ( XFILE * ayb_fp, const CIFDATA const cif, bool showall){
    static const char * basechar = "ACGT";
    if ( NULL==ayb_fp) return;
    if ( NULL==cif) return;

    xfprintf( ayb_fp, "cifData version = %u\n", cif->version );
    xfprintf( ayb_fp, "datasize = %u bytes\n", cif->datasize );
    xfprintf( ayb_fp, "ncycles = %u, of which the first is cycle number %u\n", cif->ncycle,cif->firstcycle);
    xfprintf( ayb_fp, "nclusters = %u\n", cif->ncluster);

    uint32_t mcluster, mcycle;
    if (showall) {
        mcluster = cif->ncluster;
        mcycle   = cif->ncycle;
    }
    else {
        mcluster = (cif->ncluster>5)?5:cif->ncluster;
        mcycle   = (cif->ncycle>5)?5:cif->ncycle;
    }

    for ( int cluster=0 ; cluster<mcluster ; cluster++){
        xfprintf( ayb_fp, "Cluster %d\n", cluster+1 );
        for ( int base=0 ; base<4 ; base++){
            xfputc (basechar[base],ayb_fp);
            for ( int cycle=0 ; cycle<mcycle ; cycle++){
                float f;
                switch(cif->datasize){
                    case 1: f = (float)cif->intensity.i8[(cycle*NCHANNEL+base)*cif->ncluster+cluster]; break;
                    case 2: f = (float)cif->intensity.i16[(cycle*NCHANNEL+base)*cif->ncluster+cluster]; break;
                    case 4: f = (float)cif->intensity.i32[(cycle*NCHANNEL+base)*cif->ncluster+cluster]; break;
                    default: f = NAN;
                }
                xfprintf( ayb_fp, " %5.0f", f);
            }
            xfputc('\n',ayb_fp);
        }
    }
    if( mcluster!=cif->ncluster){
        xfprintf(ayb_fp,"%u clusters omitted. ", cif->ncluster - 5 );
    }
    if( mcycle!=cif->ncycle ){
        xfprintf(ayb_fp,"%u cycles omitted. ", cif->ncycle - 5 );
    }
    xfputc('\n',ayb_fp);
}


#ifdef TEST
#include <time.h>
static void timestamp(const char * str, FILE * fp){
    time_t t = time(&t);
    char * c = ctime(&t);
    c[24]='\t';
    fprintf(fp,"%s%s",c,str);
}

int main ( int argc, char * argv[] ){
    if ( argc!=5 ){
        fputs("./a.out lane tile in_cif_filename out_cif_filename\n",stderr);
        return EXIT_FAILURE;
    }

   int lane,tile;
   sscanf(argv[1],"%d",&lane);
   sscanf(argv[2],"%d",&tile);
timestamp("Starting\n",stderr);
   CIFDATA cif = readCIFfromDir(argv[3],lane,tile,XFILE_RAW);
timestamp("Read\n",stderr);
    showCIF(xstdout,cif,false);
timestamp("Splitting\n",stderr);
    CIFDATA newcif = spliceCIF(cif, cif->ncycle/2, 0);
timestamp("Writing\n",stderr);
    writeCIFtoFile(newcif,argv[4],XFILE_RAW);
    free_cif(newcif);
    free_cif(cif);
timestamp("Reading new\n",stderr);
    cif = readCIFfromFile (argv[4],XFILE_RAW);
timestamp("Done\n",stderr);
    showCIF(xstdout,cif,false);
    free_cif(cif);

    return EXIT_SUCCESS;
}
#endif

