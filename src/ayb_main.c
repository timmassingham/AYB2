/**
 * \file ayb_main.c
 * Main AYB module.
 *//*
 *  Created : 23 Feb 2010
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
#include <stdbool.h>
#include <signal.h>
#include "ayb_options.h"
#include "dirio.h"
#include "handler.h"
#include "message.h"
#include "tile.h"
#include "xio.h"

#include <unistd.h>         // temp for sleep


/* constants */

static const char *LOG_PREFIX = "ayb_";         ///< Log file prefix, pass to Message startup.

/* members */
/* none    */


/* private functions */

/** Tidy up before exit. Include all module tidyup routines here. */
static void tidyup() {
    /* message files file */
    tidyup_dirio();
    tidyup_message();
}

/* public functions */

/** Main AYB routine. What does it do? */
int main(int argc, char **argv) {
    AYBOPT *p_opt;

    /* install signal handler   */
    signal(SIGINT, INThandler);
    signal(SIGFPE, FPEhandler);

    /* read program options */
    if (!read_options(argc, argv)) {
        return EXIT_FAILURE;
    }
    /* get a pointer to the program options */
    p_opt = myopt();

    /* check number of cycles supplied - hmhm move to module when decide where ncycle to live */
    message(E_GENERIC_SD, MSG_DEBUG, "ncycle:", p_opt->ncycle);
    if (p_opt->ncycle <= 0) {
        message(E_NOCYCLES, MSG_FATAL);;
        return EXIT_FAILURE;
    }
    message(E_CYCLE_SELECT_D, MSG_INFO, p_opt->ncycle);

    /* create a message log */
    startup_message(LOG_PREFIX);

    /* scan the input directory */
    if (!startup_dirio()) {
        return EXIT_FAILURE;
    }

    /**** test file access ****/

    XFILE *fpin = NULL;
    XFILE *fpout = NULL;
    bool found = true;
//    int c;

    while (found) {
        fpin = open_next(fpin);

        if (xisnull_file(fpin)) {
            found = false;
        }
        else {
            /* create an output file from the input */
            fpout = open_output("out");
            if (xisnull_file(fpout)) {
                found = false;
            }
            else {
                /* read and store input */
                TILE tile;
                unsigned int nc = p_opt->ncycle;
                tile = read_known_TILE(fpin, &nc);
                if (nc < p_opt->ncycle) {
                    message(E_CYCLESIZE_DD, MSG_WARN, p_opt->ncycle, nc);
                }
                if (tile != NULL) {
                    /* output from store */
                    show_TILE(fpout, tile, 20);
                }

//                while ((c = xfgetc(fpin)) != EOF) {
//                    xfputc(c, fpout);
//                }
                free_TILE(tile);
                xfclose(fpout);
            }
        }
    }


    /********/

    /**** this part is for testing of log and signal handling ****/
    /* convert boolean value to string */
//#define BOOLSTR(b) b ? "true" : "false"
    //  fprintf(stdout, "A value argument is %d\n", p_opt->aval);
    //  fprintf(stdout, "A flag argument is %s\n", BOOLSTR(p_opt->aflag));

    /*
    int anint = 3;
    float afloat = 123.45;

    for (int i = 0; i < 5; i++) {
        printf("Sleep for %d seconds; divisor = %d\n", p_opt->aval, anint);
        sleep(p_opt->aval);
        afloat /= anint--;
        printf("float = %0.2f\n", afloat);
        checksignals();
    }
*/
    /********/

    /* tidy up before exit */
    tidyup();

    fprintf(stdout, "%s\n", "End of AYB Main");
    return EXIT_SUCCESS;
}
