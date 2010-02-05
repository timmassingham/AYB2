/*
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the simNGS software for simulating likelihoods
 *  for next-generation sequencing machines.
 *
 *  simNGS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  simNGS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with simNGS.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef _UTILITY_H
#define _UTILITY_H
 
#include <stdlib.h>

#define NBASE 4
static inline void xfree( void * ptr){ if(NULL!=ptr){ free(ptr);} }

#ifdef FAILEARLY
    #define validate(A,B) { if( !(A) ){ fprintf(stderr,"Validation failure for %s in %s at %s:%d",#A,__func__,__FILE__,__LINE__); abort();} }
#elif defined(SKIPVALIDATE)
    #define validate(A,B)
#else
    #define validate(A,B) { if( !(A) ){ return B; } }
#endif


#ifdef USEFLOAT
    typedef float real_t;
    #define strtor strtof
#else
    typedef double real_t;
    #define strtor strtod
#endif
 
#endif /* _UTILITY_H */

