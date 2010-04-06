/**
 * \file utility.c
 * General Utilities.
 *   - CSTRING; a simple string type.
 *//*
 *  Created : 16 Mar 2010
 *  Author  : Hazel Marsden
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
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
#include <stdio.h>
#include <string.h>
#include "utility.h"


/* constants */
/* none      */

/* members */
/* none    */

/* private functions */
/* none              */

/* public functions */

/* CSTRING */

/** Create a new string of specified length. Allocates memory adding space for terminator. */
CSTRING new_CSTRING(const size_t len) {
    return calloc(len+1, sizeof(char));
}

/** Free memory allocated to a string. */
void free_CSTRING(CSTRING c) {
    xfree(c);
}

/** Copy a string to a new one, allocating memory. */
CSTRING copy_CSTRING(const CSTRING c) {
    validate(NULL!=c, NULL);
    CSTRING n = new_CSTRING(strlen(c));
    validate(NULL!=n, NULL);
    strcpy(n, c);
    return n;
}

/** Increase the length of a string. No action if already at least specified length. */
void extend_CSTRING(CSTRING c, const size_t len) {
    if (NULL==c) {
        /* nothing to copy */
        c = new_CSTRING(len);
    }
    else {
        if (strlen(c) < len) {
            /* make longer */
            CSTRING n = new_CSTRING(len);
            validate(NULL!=n,);
            strcpy(n, c);
            free_CSTRING(c);
            c = n;
        }
        else {
            /* long enough already, do nothing */
        }
    }
}

/** Output contents of string to specified file. */
void show_CSTRING(FILE *fp, const CSTRING c) {
    validate(NULL!=fp,);
    validate(NULL!=c,);
    fputs(c, fp);
}
