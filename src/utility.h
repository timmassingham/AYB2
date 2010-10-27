/**
 * \file utility.h
 * Public parts of General Utilities.
 *   - global constants
 *   - xfree
 *   - validate
 *   - USEFLOAT switch
 *   - CSTRING; a simple string type
 *//*
 *  Created : 16 Mar 2010
 *  Authors : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
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
 
#ifndef _UTILITY_H
#define _UTILITY_H
 
#include <stdlib.h>
#include <stdio.h>

/** Number of sequence bases. */
//#define NBASE 4

/** Safe Free memory. Checks supplied pointer not null. */
static inline void xfree( void * ptr){ if(NULL!=ptr){ free(ptr);} }

/**
 * Generic validation.
 * Checks supplied boolean and forces a function return with any second parameter if evaluates to false.
 */
#ifdef FAILEARLY
    #define validate(A,B) { if( !(A) ){ fprintf(stderr,"Validation failure for %s in %s at %s:%d",#A,__func__,__FILE__,__LINE__); abort();} }
#elif defined(SKIPVALIDATE)
    #define validate(A,B)
#else
    #define validate(A,B) { if( !(A) ){ return B; } }
#endif

/** Define size of real data type to use. May be set to float or double.*/
#ifdef USEFLOAT
    typedef float real_t;
    /** Input format string for selected real_t. */
    #define REAL_FORMAT_IN "%f"
    /** String to real converter for selected real_t. */
    #define strtor strtof
    /** Rounding function for selected real_t sends halfway cases away from zero. */
    #define roundr roundf
    /** Infinity */
    #define HUGE_VALR HUGE_VALF
#else
    typedef double real_t;
    /** Input format string for selected real_t. */
    #define REAL_FORMAT_IN "%lf"
    /** String to real converter for selected real_t. */
    #define strtor strtod
    /** Rounding function for selected real_t sends halfway cases away from zero. */
    #define roundr round
    /** Infinity */
    #define HUGE_VALR HUGE_VAL
#endif

/** Simple string type. */
typedef char * CSTRING;
CSTRING new_CSTRING(const size_t len);
CSTRING free_CSTRING(CSTRING cstr);
CSTRING copy_CSTRING(const CSTRING cstr);
void extend_CSTRING(const CSTRING c, const size_t len);
void show_CSTRING(FILE *fp, const CSTRING cstr);
//CSTRING read_CSTRING(FILE *fp);

/* General string utilities */
int match_string(const char *string, const char *match[], int num);

#endif /* _UTILITY_H */

