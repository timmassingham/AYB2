/*
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

#ifndef HAS_REALLOCF
void * reallocf(void * ptr, size_t t){
	void * ptr2;
	ptr2 = realloc(ptr,t);
	if(ptr2==NULL && ptr!=NULL){
	    free(ptr);
	}
	return ptr2;
}
#endif

typedef enum { XFILE_UNKNOWN, XFILE_RAW, XFILE_GZIP, XFILE_BZIP2 } XFILE_MODE;
typedef union { BZFILE * bzfh; gzFile zfh; FILE * fh;} XFILE_TYPE;

struct _xfile_struct {
    XFILE_MODE mode;
    XFILE_TYPE ptr;
};

typedef struct _xfile_struct XFILE;

XFILE _xstdout, _xstderr;
XFILE * xstdout = &_xstdout;
XFILE * xstderr = &_xstderr;
static bool _xinit = false;

void initialise_aybstd ( void ){
   _xstdout.mode = XFILE_RAW; _xstdout.ptr.fh = stdout;
   _xstderr.mode = XFILE_RAW; _xstderr.ptr.fh = stderr;
   _xinit = true;
}



const char * find_suffix ( const char * fn ){
	const size_t len = strlen(fn);

	for ( size_t i=len-1 ; i>0 ; i-- ){
		if ( fn[i] == '.') return fn+i+1;
	}
	if(fn[0]=='.') return fn+1;
	
	return fn+len; // Pointer to '\0' at end of string
}


XFILE_MODE guess_mode_from_filename ( const char * fn ){
	const char * suffix = find_suffix (fn);
	if ( strcmp(suffix,"gz") == 0 ){ return XFILE_GZIP;}
	if ( strcmp(suffix,"bz2") == 0 ){ return XFILE_BZIP2;}
	
	return XFILE_RAW;
}



int gzvprintf ( gzFile zfp, const char * fmt, va_list args ){
    int ret,len;
    char * buf;

    vasprintf(&buf,fmt,args);
    len = strlen(buf);
    ret = gzputs(zfp,buf);
    free(buf);
    return ret;
}

int BZ2_bzvprintf ( BZFILE * bzfp, const char * fmt, va_list args ){
    int ret,len;
    char * buf;

    vasprintf(&buf,fmt,args);
    len = strlen(buf);              
    ret = BZ2_bzwrite( bzfp,buf,len*sizeof(char));
    free(buf);                                      
    return ret;                                             
}



int xfprintf( XFILE * fp, const char * fmt, ... ){
    int ret=EOF;
    va_list args;
    
    if(!_xinit){initialise_aybstd();}

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

void xnull_file(XFILE * fp){
    switch( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW: fp->ptr.fh = NULL; break;
      case XFILE_GZIP: fp->ptr.zfh = NULL; break;
      case XFILE_BZIP2: fp->ptr.bzfh = NULL; break;
    }
}

int xnotnull_file(XFILE * fp){
    switch( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW:   if(NULL==fp->ptr.fh){ return 0;} break;
      case XFILE_GZIP:  if(NULL==fp->ptr.zfh){ return 0;} break;
      case XFILE_BZIP2: if(NULL==fp->ptr.bzfh){ return 0;} break;
    }

    return 1;
}

void xfclose(XFILE * fp){
    if(!_xinit){initialise_aybstd();}
    if( ! xnotnull_file(fp) ){return;}

    switch( fp->mode ){
	  case XFILE_UNKNOWN:
      case XFILE_RAW:     fclose(fp->ptr.fh); break;
      case XFILE_GZIP:  gzclose(fp->ptr.zfh); break;
      case XFILE_BZIP2: BZ2_bzclose(fp->ptr.bzfh); break;
    }
    free(fp);
}

XFILE * xfopen(const char * restrict fn, const XFILE_MODE mode, const char * mode_str){
    if(!_xinit){initialise_aybstd();}
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
      default: xnull_file(fp); fail=1;
    }

    if (fail){
        perror(fn);
        xnull_file(fp);
    }

    return fp;
}


size_t xfread(void *ptr, size_t size, size_t nmemb, XFILE *fp){
    if(!_xinit){initialise_aybstd();}
    size_t ret = 0;

    switch( fp->mode ){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fread(ptr,size,nmemb,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzread(fp->ptr.zfh,ptr,size*nmemb); break;
        case XFILE_BZIP2: ret = BZ2_bzread(fp->ptr.bzfh,ptr,size*nmemb); break;
    }

    return ret; 
}

size_t xfwrite(const void * restrict ptr, const size_t size, const size_t nmemb, XFILE * fp){
    if(!_xinit){initialise_aybstd();}
    size_t ret = 0;

    switch( fp->mode ){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fwrite(ptr,size,nmemb,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzwrite(fp->ptr.zfh,ptr,size*nmemb); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,ptr,size*nmemb); break;
    }       

    return ret;
}

int xfputc ( int c, XFILE * fp){
    if(!_xinit){initialise_aybstd();}
    int ret = EOF;
    switch(fp->mode){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fputc(c,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzputc(fp->ptr.zfh,c); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,&c,sizeof(char)); break;
    }
    return ret;
}

int xfputs ( const char * restrict str, XFILE * fp){
    if(!_xinit){initialise_aybstd();}
    int ret = EOF;
    switch(fp->mode){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fputs(str,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzputs(fp->ptr.zfh,str); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,str,strlen(str)*sizeof(char)); break;
    }
    return ret;
}


int xfgetc(XFILE * fp){
    char c=0;
    int ret = xfread(&c,sizeof(char),1,fp);
    return (ret!=0)?c:EOF;
}

char * xfgets( char * restrict s, int n, XFILE * restrict fp){
    int ret = xfread(s,sizeof(char),n-1,fp);
    s[ret] = '\0';
    return s;
}

/*  Gets a line from the file, allocating necessary memory.
 *  Returns pointer to line.
 *  Returns NULL if failed to allocate memory or other error
 *  Does not deal with MS-DOS style line-feeds gracefully.
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

