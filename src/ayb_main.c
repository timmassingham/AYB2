/**
 * \file ayb_main.c
 * Main AYB module.
 * Sets up environment before passing control to ayb_model for each intensity file.
 * Tidy up to finish.
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
#include "ayb_model.h"
#include "ayb_options.h"
#include "dirio.h"
#include "handler.h"
#include "message.h"
#include "xio.h"


/* constants */

static const char *LOG_PREFIX = "ayb_";         ///< Log file prefix, pass to Message startup.

/* members */
/* none    */


/* private functions */

/** Tidy up before exit. Include all module tidyup routines here. */
static void tidyup() {
    tidyup_model();
    tidyup_dirio();
    tidyup_message();
}

/* public functions */

/**
 * Main AYB routine.
 * Read program options and set up environment.
 * Then process each input file before tidy up and exit.
 */
int main(int argc, char **argv) {

    /* install signal handler   */
    signal(SIGINT, INThandler);
    signal(SIGFPE, FPEhandler);

    /* read program options */
    if (!read_options(argc, argv)) {
        return EXIT_FAILURE;
    }

    /* create a message log */
    startup_message(LOG_PREFIX);

    /* scan the input directory */
    if (!startup_dirio()) {
        return EXIT_FAILURE;
    }

    if (!startup_model()) {
        return EXIT_FAILURE;
    }

    /* process each intensity file */
    XFILE *fp = NULL;
    bool found = true;

    while (found) {
        fp = open_next(fp);

        if (xfisnull(fp)) {
            found = false;
        }
        else {
            /* analyse this input file */
            analyse_tile(fp);
        }
    }

    /* tidy up before exit */
    tidyup();

    fprintf(stdout, "%s\n", "End of AYB Main");
    return EXIT_SUCCESS;
}
