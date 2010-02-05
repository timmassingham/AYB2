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

#ifndef _XIO_H
#define _XIO_H

#include <stdlib.h>

typedef enum { XFILE_UNKNOWN, XFILE_RAW, XFILE_GZIP, XFILE_BZIP2 } XFILE_MODE;
typedef struct _xfile_struct XFILE;

extern XFILE * xstdout;
extern XFILE * xstderr;

// File open and close routines
XFILE * xfopen(const char * restrict fn, const XFILE_MODE mode, const char * mode_str);
void xfclose(XFILE * fp);
void xnull_file(XFILE * fp);


// Helper routine to guess type of file from suffix
XFILE_MODE guess_mode_from_filename ( const char * restrict fn );


// Functions to read / write as binary
size_t xfread(void *ptr, size_t size, size_t nmemb, XFILE *fp);
size_t xfwrite(const void *ptr, size_t size, size_t nmemb, XFILE * fp);


// Basic string functions
int xfputc ( int c, XFILE * fp);
int xfputs ( const char * restrict str, XFILE * fp);
int xfprintf( XFILE * fp, const char * fmt, ... );

int xfgetc (XFILE * fp);
char * xfgets( char * restrict s, int n, XFILE * restrict fp);
char * xfgetln( XFILE * fp, size_t * len);



#endif /* _XIO_H */

