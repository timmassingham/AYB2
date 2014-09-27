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

XFILE _xstdin, _xstdout, _xstderr;
XFILE * xstdin  = &_xstdin;	            ///< Standard input as an XFILE.
XFILE * xstdout = &_xstdout;            ///< Standard output as an XFILE.
XFILE * xstderr = &_xstderr;            ///< Standard error as an XFILE.
static bool _xinit = false;


/* private functions */

/** Create standard output and error as an XFILE. */
static void initialise_std ( void ){
   _xstdin.mode = XFILE_RAW; _xstdin.ptr.fh = stdin;
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
static int __gzvprintf ( gzFile zfp, const char * fmt, va_list args ){
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
	if (NULL==fn){ return XFILE_UNKNOWN;}
    const char * suffix = find_suffix (fn);
    if ( strcmp(suffix, GZ) == 0 ){ return XFILE_GZIP;}
    if ( strcmp(suffix, BZ2) == 0 ){ return XFILE_BZIP2;}

    return XFILE_RAW;
}

/** Close an XFILE. Closes selected file and frees structure memory. */
XFILE * xfclose(XFILE * fp){
    if(!_xinit){initialise_std();}
    if (NULL==fp) { return NULL;}

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
 * Opens the file in the appropriate mode.
 * mode_str selects read/write etc as for normal file stream
 * (note: gzip does not support '+' and bzip2 is only known to support 'r' and 'w'.)
 * If fails to open then frees structure memory and returns a null pointer.
 */
XFILE * xfopen(const char * restrict fn, const XFILE_MODE mode, const char * mode_str){
    if(!_xinit){initialise_std();}
	if (NULL==fn){ return NULL;}
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
 * Returns number of bytes read (or number of elements for mode RAW).
 */
size_t xfread(void *ptr, size_t size, size_t nmemb, XFILE *fp){
    if(!_xinit){initialise_std();}
	if (NULL==fp) { return 0;}
	if (NULL==ptr) { return 0;}
    size_t ret = 0;
    int retz = 0;

    switch( fp->mode ){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fread(ptr,size,nmemb,fp->ptr.fh); break;
        case XFILE_GZIP:  retz = gzread(fp->ptr.zfh,ptr,size*nmemb);
                          if (retz>0) {ret = retz;} break;
        case XFILE_BZIP2: retz = BZ2_bzread(fp->ptr.bzfh,ptr,size*nmemb);
                          if (retz>0) {ret = retz;} break;
    }

    return ret; 
}

/**
 * Write a block of data, up to nmemb elements of the given size.
 * Returns number of bytes written (or number of elements for mode RAW).
 */
size_t xfwrite(const void * restrict ptr, const size_t size, const size_t nmemb, XFILE * fp){
    if(!_xinit){initialise_std();}
	if (NULL==fp) { return 0;}
	if (NULL==ptr) { return 0;}
    size_t ret = 0;
    int retz = 0;

    switch( fp->mode ){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fwrite(ptr,size,nmemb,fp->ptr.fh); break;
        case XFILE_GZIP:  retz = gzwrite(fp->ptr.zfh,ptr,size*nmemb);
                          if (retz>0) {ret = retz;} break;
        case XFILE_BZIP2: retz = BZ2_bzwrite(fp->ptr.bzfh,(void*)ptr,size*nmemb);
                          if (retz>0) {ret = retz;} break;
    }       

    return ret;
}

/** Write a character. Returns negative value if fails to write. */
int xfputc ( int c, XFILE * fp){
    if(!_xinit){initialise_std();}
	if (NULL==fp) { return EOF;}
    int ret = EOF;
    switch(fp->mode){
    	case XFILE_UNKNOWN:
        case XFILE_RAW:   ret = fputc(c,fp->ptr.fh); break;
        case XFILE_GZIP:  ret = gzputc(fp->ptr.zfh,c); break;
        case XFILE_BZIP2: ret = BZ2_bzwrite(fp->ptr.bzfh,&c,sizeof(char)); break;
    }
    return ret;
}

/** 
 * Write a string. Does not append the terminating null or a new line. 
 * Returns negative value if fails to write. 
 */
int xfputs ( const char * restrict str, XFILE * fp){
    if(!_xinit){initialise_std();}
	if (NULL==fp) { return EOF;}
	if (NULL==str) { return EOF;}
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
	if (NULL==fp) { return EOF;}
	if (NULL==fmt) { return EOF;}

    va_start(args,fmt);
    switch( fp->mode ){
      case XFILE_UNKNOWN:
      case XFILE_RAW:   ret=vfprintf(fp->ptr.fh,fmt,args); break;
      case XFILE_GZIP:  ret=__gzvprintf(fp->ptr.zfh,fmt,args); break;
      case XFILE_BZIP2: ret=BZ2_bzvprintf(fp->ptr.bzfh,fmt,args); break;
    }
    va_end(args);
    return ret;
}

/** Read and return a single character. Returns EOF if fails to read. */
int xfgetc(XFILE * fp){
	if (NULL==fp) { return EOF;}
    char c=0;
    int ret = xfread(&c,sizeof(char),1,fp);
    return (ret>0)?c:EOF;
}

/**
 * Read a string of length n. Appends a null terminator.
 * Result placed in buffer *s which must be large enough to accommodate.
 * Returns NULL if fails to read.
 */
char * xfgets( char * restrict s, int n, XFILE * restrict fp){
	if (NULL==fp) { return NULL;}
	if (NULL==s) { return NULL;}
    int ret = xfread(s,sizeof(char),n-1,fp);
    if (ret<=0) { return NULL;}
    s[ret] = '\0';
    return s;
}

/**
 * Get a line from a file, allocating necessary memory.
 * Calling function is responsible for freeing allocated memory.
 * Returns pointer to line, and length in len. Appends a null terminator.
 * Returns NULL and zero len if fails to read anything.
 * Returns NULL if fails to allocate memory or other error.
 * Does not deal with MS-DOS style line-feeds gracefully.
 */
char * xfgetln( XFILE * fp, size_t * len ){
	return xfgettok(fp,len,"\n");
}

/**
 * Read token from a file until sep(arator) or EOF is found, allocating 
 * necessary memory.
 * Calling function is responsible for freeing allocated memory.
 * Returns pointer to line, and length in len. Appends a null terminator.
 * Returns NULL and zero len if fails to read anything.
 * Returns NULL if fails to allocate memory or other error.
 * Does not deal with MS-DOS style line-feeds gracefully.
 */

char * xfgettok( XFILE * fp, size_t * len, const char * sep){
	*len = 0;
	if (NULL==fp) { return NULL;}
	if (NULL==sep) { return NULL;}
	size_t seplen = strlen(sep);

	char * str = NULL;
	int size = 0;
	int c = 0;
	int sepidx = 0;
	while ( c=xfgetc(fp), (c!=EOF) ){
		if(size<=*len){
			size += 80;
			str = reallocf(str, size);
			if(NULL==str){ return NULL;}
		}
		str[*len] = c;
		(*len)++;
		if(c==sep[sepidx]){
			sepidx++;
			// Has end of token been reached?
			if(sepidx==seplen){
				break;
			}
		} else {
			// Reset counter
			sepidx = 0;
		}
	}
	// If nothing read, return.
	if(EOF==c && NULL==str){ return NULL;}
	// Remove any trailing bit of separator
	*len = *len - sepidx;
	// Make sure there is sufficient memory for terminating '\0'
	if(size <= *len){
		size += 1;
		str = reallocf(str, size);
		if(NULL==str){ return NULL;}
	}
	str[*len] = '\0';

	return str;
}


#ifdef TEST
#include <err.h>
#include <stdint.h>

static const char CH1 = 'x';
static const char ALINE[] = "..a line of text..";
static const char ASWAS[] = "aswas";
static const char DATA[] = "written_as_block";
static const char FLINE[] = "A formatted line; string %s; char %c; int %d; float %f\n";
static const char * DOT = ".";
static const char * STR1 = "xxxx";
static const char * SUFF[] = {"", "txt", "txt.gz", "txt.bz2"};
static const int INT1 = 99;
static const double FL1 = 99.99;

/* for writing as block */
typedef uint16_t int_t;
size_t IntSize = sizeof(int_t);
size_t DataLen;

/* Add file type suffix to name */
static char * add_suff(const char * name, const char * suff) {
    if (NULL==name) { return NULL;}
    if (NULL==suff) { return NULL;}
    char * fullname = NULL;
    fullname = calloc(strlen(name) + strlen(DOT) + strlen(suff) + 1, sizeof(char));
    strcpy(fullname, name);
    strcat(fullname, DOT);
    strcat(fullname, suff);
    return fullname;
}

/* used for null, empty and open for write file */
static void read_bad_file(XFILE * fp, const char *desc, int_t *intdata) {
    fprintf(stdout, "Read from %s file\n", desc);
    int ch = xfgetc(fp);
    if (ch==EOF) {
        fputs("Return value xfgetc EOF, ok\n", stdout);
    }
    else {
        fputs("Return value xfgetc not EOF, not ok\n", stdout);
    }
    char * line = malloc(9);
    strcpy(line, ASWAS);
    char * res = xfgets(line, 9, fp);
    if (res==NULL) {
        fprintf(stdout, "Return value xfgets null, ok; string: %s\n", line);
    }
    else {
        fputs("Return value xfgets not null, not ok\n", stdout);
    }
    free(line);
    size_t lenz = 99;
    line = xfgetln(fp, &lenz);
    if ((lenz==0) && (line==NULL)) {
        fputs("Return value xfgetln null and zero, ok\n", stdout);
    }
    else {
        fputs("Return value xfgetln not null and zero, not ok\n", stdout);
    }
    if (line!=NULL) { free(line); }
    lenz = 99;
    lenz = xfread(intdata, IntSize, DataLen, fp);
    if (lenz==0) {
        fputs("Return value xfread zero, ok\n", stdout);
    }
    else {
        fputs("Return value xfread not zero, not ok\n", stdout);
    }
}

int main ( int argc, char * argv[]){
	if(argc<4){
	    /* arguments are separator, input file and output file */
		errx(EXIT_FAILURE, "Usage: test-xio separator in_filename out_filename");
	}

    XFILE * fp = NULL;
    int num[3] = {0, 0, 0};
    size_t lenz = 0;
    XFILE_MODE mode;
    DataLen = (strlen(DATA) + IntSize - 1) / IntSize;    // round up
    int_t intdata[DataLen];

    fputs("Guess mode from null filename\n", stdout);
    mode = guess_mode_from_filename (NULL);
    if (mode==XFILE_UNKNOWN) {
        fputs("Return value unknown, ok\n", stdout);
    }
    else {
        fputs("Return value not unknown, not ok\n", stdout);
    }

    fputs("Guess mode from empty filename\n", stdout);
    mode = guess_mode_from_filename ("");
    if (mode==XFILE_RAW) {
        fputs("Return value raw, ok\n", stdout);
    }
    else {
        fputs("Return value not raw, not ok\n", stdout);
    }

    fputs("Guess mode from arbitrary filename\n", stdout);
    mode = guess_mode_from_filename (STR1);
    if (mode==XFILE_RAW) {
        fputs("Return value raw, ok\n", stdout);
    }
    else {
        fputs("Return value not raw, not ok\n", stdout);
    }

    fputs("Open null filename\n", stdout);
    fp = xfopen(NULL, XFILE_UNKNOWN, "r");
    if (fp==NULL) {
        fputs("Return value null, ok\n", stdout);
    }
    else {
        fputs("Return value not null, not ok\n", stdout);
    	fp = xfclose(fp);
    }
    
    fputs("Open empty filename\n", stdout);
    fp = xfopen("", XFILE_UNKNOWN, "r");
    if (fp==NULL) {
        fputs("Return value null, ok\n", stdout);
    }
    else {
        fputs("Return value not null, not ok\n", stdout);
    	fp = xfclose(fp);
    }
    
    fputs("Close null file\n", stdout);
	fp = xfclose(NULL);
    if (fp==NULL) {
        fputs("Return value null, ok\n", stdout);
    }
    else {
        fputs("Return value not null, not ok\n", stdout);
    }

    fputs("Write to null file\n", stdout);
    memset(num, 9, 3 * sizeof(int));
    num[0] = xfputc(CH1, NULL);
    num[1] = xfputs(ALINE, NULL);
    num[2] = xfprintf(NULL, FLINE, STR1, CH1, INT1, FL1);
    lenz = 99;
    lenz = xfwrite(intdata, IntSize, DataLen, NULL);
    fprintf(stdout, "Values returned: %d, %d, %d, %lu\n", num[0], num[1], num[2], lenz);

    /* read from null file */
    read_bad_file(NULL, "null", intdata);
    
    fputs("Test each file type...\n", stdout);
    
    for (XFILE_MODE modi = XFILE_RAW; modi <= XFILE_BZIP2; modi++) {
        char * line = NULL;
        char * res = NULL;
        char * fullname = NULL;
        
        fputs("\n", stdout);
        fullname = add_suff(STR1, SUFF[modi]);
        if (NULL==fullname) {
		    errx(EXIT_FAILURE, "Failed to create filename");
        }
        
        fputs("Open invalid filename for read\n", stdout);
        fp = xfopen(STR1, modi, "r");
        if (xfisnull(fp)) {
            fputs("Return value null, ok\n", stdout);
        }
        else {
            fputs("Return value not null, not ok\n", stdout);
        	fp = xfclose(fp);
        }
        free(fullname);

        fullname = add_suff(argv[3], SUFF[modi]);
        if (NULL==fullname) {
		    errx(EXIT_FAILURE, "Failed to create filename");
        }

        fprintf(stdout, "Guess mode from filename: %s\n", fullname);
        mode = guess_mode_from_filename (fullname);
        if (mode==modi) {
            fputs("Return value correct, ok\n", stdout);
        }
        else {
            fputs("Return value incorrect, not ok\n", stdout);
        }

        fputs("Open invalid mode string\n", stdout);
        if (modi < XFILE_BZIP2) {
            fp = xfopen(fullname, modi, "x");
            if (xfisnull(fp)) {
                fputs("Return value null, ok\n", stdout);
            }
            else {
                fputs("Return value not null, not ok\n", stdout);
            	fp = xfclose(fp);
            }
        }
        else {
            fputs("(Skip for bzip2; fails test and appears to open read)\n", stdout);
        }

        fprintf(stdout, "Open new file for write, known mode, suffix: %s\n", SUFF[modi]);
        fp = xfopen(fullname, modi, "w");
        if (xfisnull(fp)) {
		    errx(EXIT_FAILURE, "Failed to open %s for output", fullname);
        }
        
        /* read from open for write file */
        read_bad_file(fp, "open for write", intdata);

	    fputs("Close file\n", stdout);
	    fp = xfclose(fp);
        if (fp==NULL) {
            fputs("Return value null, ok\n", stdout);
        }
        else {
            fputs("Return value not null, not ok\n", stdout);
        }

        fprintf(stdout, "Re-open empty file for read, known mode, suffix: %s\n", SUFF[modi]);
        fp = xfopen(fullname, modi, "r");
        if (xfisnull(fp)) {
		    errx(EXIT_FAILURE, "Failed to open %s for input", fullname);
        }
        
        /* read from empty file */
        read_bad_file(fp, "empty", intdata);
        
        fputs("Write to open for read file, xfputc, xfputs, xfprintf, xfwrite\n", stdout);
        memset(num, 9, 3 * sizeof(int));
        num[0] = xfputc(CH1, fp);
        num[1] = xfputs(ALINE, fp);
        num[2] = xfprintf(fp, FLINE, STR1, CH1, INT1, FL1);
        lenz = 99;
        lenz = xfwrite(intdata, IntSize, DataLen, fp);
        fprintf(stdout, "Values returned: %d, %d, %d, %lu\n", num[0], num[1], num[2], lenz);

	    fp = xfclose(fp);

        fprintf(stdout, "Open new file for write, unknown mode, suffix: %s\n", SUFF[modi]);
        fp = xfopen(fullname, XFILE_UNKNOWN, "w");
        if (xfisnull(fp)) {
		    errx(EXIT_FAILURE, "Failed to open %s for output", fullname);
        }

        fputs("Write with null data pointer, xfputs, xfprintf, xfwrite\n", stdout);
        memset(num, 9, 3 * sizeof(int));
        num[1] = xfputs(NULL, fp);
        num[2] = xfprintf(fp, NULL, STR1, CH1, INT1, FL1);
        lenz = 99;
        lenz = xfwrite(NULL, IntSize, DataLen, fp);
        fprintf(stdout, "Values returned: %d, %d, %lu\n", num[1], num[2], lenz);

        fputs("Write with xfputc, xfputs, xfprintf, xfwrite\n", stdout);
        memset(num, 0, 3 * sizeof(int));
        num[0] = xfputc(CH1, fp);
        num[1] = xfputs(ALINE, fp);
        xfputc('\n', fp);   // add a new line
        num[2] = xfprintf(fp, FLINE, STR1, CH1, INT1, FL1);

        memset(intdata, 0, DataLen * IntSize);
        memcpy(intdata, DATA, strlen(DATA));
        lenz = 0;
        lenz = xfwrite(intdata, IntSize, DataLen, fp);
        fprintf(stdout, "Values returned: %d, %d, %d, %lu\n", num[0], num[1], num[2], lenz);
      
        fprintf(stdout, "Close and re-open file for read, unknown mode\n");
	    fp = xfclose(fp);
        fp = xfopen(fullname, XFILE_UNKNOWN, "r");
        if (xfisnull(fp)) {
		    errx(EXIT_FAILURE, "Failed to open %s for input", fullname);
        }

        fputs("Read with null data pointer\n", stdout);
        line = malloc(9);
        strcpy(line, ASWAS);
        res = xfgets(NULL, 9, fp);
        if (res==NULL) {
            fprintf(stdout, "Return value xfgets null, ok; string: %s\n", line);
        }
        else {
            fputs("Return value xfgets not null, not ok\n", stdout);
        }
        free(line);
        lenz = 99;
        lenz = xfread(NULL, IntSize, DataLen, fp);
        if (lenz==0) {
            fputs("Return value xfread zero, ok\n", stdout);
        }
        else {
            fputs("Return value xfread not zero, not ok\n", stdout);
        }

        fputs("Read with xfgetc, xfgets, xfgetln, xfread\n", stdout);
        int ch = xfgetc(fp);
        lenz = strlen(ALINE);
        line = malloc(lenz + 1);
        strcpy(line, ASWAS);
        res = xfgets(line, (int)lenz + 1, fp);
        if (res==NULL) {
            fputs("Failed to read with xfgets\n", stdout);
        }
        fprintf(stdout, "%c%s\n", ch, line);
        free(line);
        xfgetc(fp);         // skip the newline
        lenz = 0;
        line = xfgetln(fp, &lenz);
        if ((lenz > 0) && (line!=NULL)) {
            fprintf(stdout, "%s (length %lu)\n", line, lenz); 
        }
        else {
            fputs("Failed to read with xfgetln\n", stdout);
        }
        if (line!=NULL) { free(line); }
        
        memset(intdata, 0, DataLen * IntSize);
        lenz = 0;
        lenz = xfread(intdata, IntSize, DataLen, fp);
        if (lenz > 0) {
            line = malloc(DataLen * IntSize + 1);
            memcpy(line, intdata, DataLen * IntSize);
            line[DataLen * IntSize] = 0;
            fprintf(stdout, "%s (returned %lu)\n", line, lenz);
            free(line);
        }
        else {
            fputs("Failed to read with xfread\n", stdout);
        }

	    fp = xfclose(fp);
        free(fullname);
    }
    
    /* test read by token */
    fputs("\n", stdout);
	char * sep = argv[1];

	XFILE * xfp = xfopen(argv[2], XFILE_UNKNOWN, "r");
	if(NULL==xfp){
		errx(EXIT_FAILURE, "Failed to open %s for input", argv[2]);
	}

	size_t len = 0;
	int ntok = 0;
	char * tok = NULL;
	
    fputs("Read null file by token\n", stdout);
	tok = xfgettok(NULL, &len, sep);
	if ((tok==NULL) && (len==0)) {
        fputs("Return values null and zero, ok\n", stdout);
    }
    else {
        fputs("Return values not null and zero, not ok\n", stdout);
     	if (tok!=NULL) {free(tok);}
    }

    fputs("Read file by token, null separator\n", stdout);
	tok = xfgettok(xfp, &len, NULL);
	if ((tok==NULL) && (len==0)) {
        fputs("Return values null and zero, ok\n", stdout);
    }
    else {
        fputs("Return values not null and zero, not ok\n", stdout);
     	if (tok!=NULL) {free(tok);}
    }
	xfp = xfclose(xfp);

    fputs("Read file by token, empty separator\n", stdout);
	xfp = xfopen(argv[2], XFILE_UNKNOWN, "r");
	while ((tok = xfgettok(xfp, &len, "")), (len!=0) ) {
		ntok++;
		fprintf(stdout, "Token %d chars %lu: [%s]\n", ntok, len, tok);
     	free(tok);
	}
	if (ntok==1) {
        fputs("Single token returned, ok\n", stdout);
	}
	else {
        fputs("More than one token returned, not ok\n", stdout);
	}
	if (tok==NULL) {
        fputs("Final Return value null, ok\n", stdout);
    }
    else {
        fputs("Final return value not null, not ok\n", stdout);
     	free(tok);
    }
	xfp = xfclose(xfp);

    fprintf(stdout, "Read file by token \'%s\'\n", sep);
	xfp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    ntok = 0;
	while ((tok = xfgettok(xfp, &len, sep)), (len!=0) ) {
		ntok++;
		fprintf(stdout, "Token %d chars %lu: [%s]\n", ntok, len, tok);
     	free(tok);
	}
	if (tok==NULL) {
        fputs("Final Return value null, ok\n", stdout);
    }
    else {
        fputs("Final return value not null, not ok\n", stdout);
     	free(tok);
    }
	xfp = xfclose(xfp);

    fputs("Read file by line\n", stdout);
	xfp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    ntok = 0;
	while ((tok = xfgetln(xfp, &len)), (len!=0) ) {
		ntok++;
		fprintf(stdout, "Line %d chars %lu: [%s]\n", ntok, len, tok);
    	free(tok);
	}
	if (tok==NULL) {
        fputs("Final Return value null, ok\n", stdout);
    }
    else {
        fputs("Final return value not null, not ok\n", stdout);
     	free(tok);
    }
	xfp = xfclose(xfp);

	return EXIT_SUCCESS;
}

#endif /* TEST */
