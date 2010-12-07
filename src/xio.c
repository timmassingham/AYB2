/**
 * \file xio.c
 * Generic File Access including Compressed.
 * Implements an XFILE type that can represent a variety of compression modes.
 * None, gzip (zlib) and bzip2 (bzlib) are currently supported.
 *//*
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2008-2010 by Tim Massingham
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
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <bzlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include "xio.h"


#ifndef HAS_REALLOCF
/** Reallocate with free memory and set pointer to null if extend fails. */
void * reallocf(void * ptr, size_t t){
	void * ptr2;
	ptr2 = realloc(ptr,t);
	if(ptr2==NULL && ptr!=NULL){
	    free(ptr);
	}
	return ptr2;
}
#endif

/** File pointer of variable mode. */
typedef union { BZFILE * bzfh; gzFile zfh; FILE * fh;} XFILE_TYPE;

/**
 * XFILE structure contains the mode and the appropriate file pointer.
 * Definition of Mode and XFILE typedef in header.
 */
struct _xfile_struct {
    XFILE_MODE mode;
    XFILE_TYPE ptr;
};

/* constants */

const char * GZ = "gz";
const char * BZ2 = "bz2";

/* members */

XFILE _xstdout, _xstderr;
XFILE * xstdout = &_xstdout;            ///< Standard output as an XFILE.
XFILE * xstderr = &_xstderr;            ///< Standard error as an XFILE.
static bool _xinit = false;


/* private functions */

/** Create standard output and error as an XFILE. */
static void initialise_std ( void ){
   _xstdout.mode = XFILE_RAW; _xstdout.ptr.fh = stdout;
   _xstderr.mode = XFILE_RAW; _xstderr.ptr.fh = stderr;
   _xinit = true;
}

/** Return a pointer to the final suffix of the supplied file name. */
static const char * find_suffix ( const char * fn ){
	const size_t len = strlen(fn);

	for ( size_t i=len-1 ; i>0 ; i-- ){
		if ( fn[i] == '.') return fn+i+1;
	}
	if(fn[0]=='.') return fn+1;
	
	return fn+len; // Pointer to '\0' at end of string
}

/** gzip printf. */
static int gzvprintf ( gzFile zfp, const char * fmt, va_list args ){
    int ret;
    char * buf;

    vasprintf(&buf,fmt,args);
    ret = gzputs(zfp,buf);
    free(buf);
    return ret;
}

/** bzip2 printf. */
static int BZ2_bzvprintf ( BZFILE * bzfp, const char * fmt, va_list args ){
    int ret,len;
    char * buf;

    vasprintf(&buf,fmt,args);
    len = strlen(buf);              
    ret = BZ2_bzwrite( bzfp,buf,len*sizeof(char));
    free(buf);                                      
    return ret;                                             
}

/** Return true if selected file pointer is not null. */
static int xnotnull_file(XFILE * fp){
    switch( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW:   if(NULL==fp->ptr.fh){ return 0;} break;
      case XFILE_GZIP:  if(NULL==fp->ptr.zfh){ return 0;} break;
      case XFILE_BZIP2: if(NULL==fp->ptr.bzfh){ return 0;} break;
    }

    return 1;
}


/* public functions */

/** Return true if XFILE is null. Checks structure and contents. */
int xfisnull(XFILE * fp) {
    if (NULL==fp) { return 1;};
    switch( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW:   if(NULL==fp->ptr.fh){ return 1;} break;
      case XFILE_GZIP:  if(NULL==fp->ptr.zfh){ return 1;} break;
      case XFILE_BZIP2: if(NULL==fp->ptr.bzfh){ return 1;} break;
    }

    return 0;
}

/** Attempt to guess the mode of file compression from the suffix of the supplied filename. */
XFILE_MODE guess_mode_from_filename ( const char * fn ){
    const char * suffix = find_suffix (fn);
    if ( strcmp(suffix, GZ) == 0 ){ return XFILE_GZIP;}
    if ( strcmp(suffix, BZ2) == 0 ){ return XFILE_BZIP2;}

    return XFILE_RAW;
}

/** Close an XFILE. Closes selected file and frees structure memory. */
XFILE * xfclose(XFILE * fp){
    if(!_xinit){initialise_std();}
    if (NULL==fp) {return NULL;}

    if(xnotnull_file(fp)) {
        switch( fp->mode ){
          case XFILE_UNKNOWN:
          case XFILE_RAW:   fclose(fp->ptr.fh); break;
          case XFILE_GZIP:  gzclose(fp->ptr.zfh); break;
          case XFILE_BZIP2: BZ2_bzclose(fp->ptr.bzfh); break;
        }
    }
    free(fp);
    return NULL;
}

/**
 * Open an XFILE. Allocates structure memory then guesses the compression mode if not supplied.
 * mode_str selects read/write etc as for normal file stream.
 * Opens the file in the appropriate mode.
 * If fails to open then frees structure memory and returns a null pointer.
 * */
XFILE * xfopen(const char * restrict fn, const XFILE_MODE mode, const char * mode_str){
    if(!_xinit){initialise_std();}
    XFILE * fp = malloc(sizeof(XFILE));
    int fail=0;

    fp->mode = mode;
    if ( XFILE_UNKNOWN==mode){ fp->mode = guess_mode_from_filename(fn);}
	
    switch ( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW:  
            fp->ptr.fh = fopen(fn,mode_str); 
            if(NULL==fp->ptr.fh){fail=1;} 
            break;
      case XFILE_GZIP:
            fp->ptr.zfh = gzopen(fn,mode_str); 
            if(NULL==fp->ptr.zfh){fail=1;} 
            break;
      case XFILE_BZIP2:
            fp->ptr.bzfh = BZ2_bzopen(fn,mode_str);
	    if(NULL==fp->ptr.zfh){fail=1;}
	    break;
      default: fail=1;
    }

    if (fail){
        perror(fn);
        free(fp);
        fp = NULL;
    }

    return fp;
}

/**
 * Read a block of data, up to nmemb elements of the given size.
 * Result placed in buffer *ptr which must be large enough to accommodate.
 * Returns number of bytes read or number of elements (mode RAW)?
 */
size_t xfread(void *ptr, size_t size, size_t nmemb, XFILE *fp){
    if(!_xinit){initialise_std();}
    size_t ret = 0;

    switch( fp->mode ){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fread(ptr,size,nmemb,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzread(fp->ptr.zfh,ptr,size*nmemb); break;
        case XFILE_BZIP2: ret = BZ2_bzread(fp->ptr.bzfh,ptr,size*nmemb); break;
    }

    return ret; 
}

/**
 * Write a block of data, up to nmemb elements of the given size.
 * Returns number of bytes written or number of elements (mode RAW)?
 */
size_t xfwrite(const void * restrict ptr, const size_t size, const size_t nmemb, XFILE * fp){
    if(!_xinit){initialise_std();}
    size_t ret = 0;

    switch( fp->mode ){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fwrite(ptr,size,nmemb,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzwrite(fp->ptr.zfh,ptr,size*nmemb); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,(void*)ptr,size*nmemb); break;
    }       

    return ret;
}

/** Write a character. */
int xfputc ( int c, XFILE * fp){
    if(!_xinit){initialise_std();}
    int ret = EOF;
    switch(fp->mode){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fputc(c,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzputc(fp->ptr.zfh,c); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,&c,sizeof(char)); break;
    }
    return ret;
}

/** Write a string. Does not append the terminating null or a new line. */
int xfputs ( const char * restrict str, XFILE * fp){
    if(!_xinit){initialise_std();}
    int ret = EOF;
    switch(fp->mode){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fputs(str,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzputs(fp->ptr.zfh,str); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,(void*)str,strlen(str)*sizeof(char)); break;
    }
    return ret;
}

/**
 * XFILE printf.
 * Returns number of characters output or negative value on error.
 */
int xfprintf( XFILE * fp, const char * fmt, ... ){
    int ret=EOF;
    va_list args;

    if(!_xinit){initialise_std();}

    va_start(args,fmt);
    switch( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW:   ret=vfprintf(fp->ptr.fh,fmt,args); break;
      case XFILE_GZIP:  ret=gzvprintf(fp->ptr.zfh,fmt,args); break;
      case XFILE_BZIP2: ret=BZ2_bzvprintf(fp->ptr.bzfh,fmt,args); break;
    }
    va_end(args);
    return ret;
}

/** Read and return a single character. Returns EOF if fails to read. */
int xfgetc(XFILE * fp){
    char c=0;
    int ret = xfread(&c,sizeof(char),1,fp);
    return (ret!=0)?c:EOF;
}

/**
 * Read a string of length n. Appends a null terminator.
 * Result placed in buffer *s which must be large enough to accommodate.
 */
char * xfgets( char * restrict s, int n, XFILE * restrict fp){
    int ret = xfread(s,sizeof(char),n-1,fp);
    s[ret] = '\0';
    return s;
}

/**
 * Get a line from a file, allocating necessary memory.
 * Calling function is responsible for freeing allocated memory.
 * Returns pointer to line, and length in len. Appends a null terminator.
 * Returns NULL if fails to allocate memory or other error.
 * Does not deal with MS-DOS style line-feeds gracefully.
 */
char * xfgetln( XFILE * fp, size_t * len ){
    char * str = NULL;
    int size = 0;
    int c = 0;
    *len = 0;
    while ( c=xfgetc(fp), (c!=EOF) && (c!='\n') && (c!='\r') ){
    	if(size<=*len){
    	    size += 80;
            str = reallocf(str, size);
            if(NULL==str){ return NULL;}
        }
        str[*len] = c;
        (*len)++;
	}
	// Make sure there is sufficient memory for terminating '\0'
	if(size<=*len){
    	size += 1;
        str = reallocf(str, size);
        if(NULL==str){ return NULL;}
    }
    str[*len] = '\0';

	
	return str;
}
