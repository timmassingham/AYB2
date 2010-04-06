/**
 * \file handler.c
 * Signal Handler.
 *//*
 *  Created : 10 Mar 2010
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

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <ctype.h>          // for toupper
#include "ayb_options.h"    // for PROGNAME


/* constants */
/* none      */

/* members */

volatile sig_atomic_t User_Interrupt = 0;       ///< Indicates user request to terminate.


/* private functions */

/** Check if any signals have arrived. */
void checksignals() {

    /* check for user interrupt */
//    printf("checksignals\n");
    if (User_Interrupt) {
        /* ask for confirmation before proceeding */
        fprintf(stdout, "  Do you really want to quit program " PROGNAME " ? [y/n] ");
        char c = getchar();
//        printf ("return from getchar: %d\n", c);
        if (toupper(c) == 'Y') {
            exit(EXIT_FAILURE);
        }
        else {
            /* flush input buffer */
            while ((c != '\n') && (c != EOF)){
            c = getchar();
//            printf ("return from flush getchar: %d\n", c);
            }
            /* reset flag */
            User_Interrupt = 0;
        }
    }
}


/* public functions */

/** Handle floating point exception signal. */
void  FPEhandler(int sig) {

     /* disable interrupt */
     signal(sig, SIG_IGN);

     /* message and exit */
     fprintf(stdout, "Quitting because of floating point exception\n");
     exit(EXIT_FAILURE);
}

/** Handle the SIGINT (Ctrl-C) signal. */
void  INThandler(int sig) {

    /* set flag and re-install */
    User_Interrupt = 1;
    signal(sig, INThandler);
}

