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
#include "datablock.h"
#include "dirio.h"
#include "handler.h"
#include "message.h"
#include "xio.h"


/* constants */
/* none    */

/* members */
/* none    */


/* private functions */

/** Tidy up before exit. Include all module tidyup routines here. */
static void tidyup() {
    tidyup_model();
    tidyup_dirio();
    tidyup_datablock();
    tidyup_message();
}

/* public functions */

/**
 * Main AYB routine.
 * Read program options and set up environment.
 * Then process each input file before tidy up and exit.
 */
int main(int argc, char **argv) {

    int ret = EXIT_SUCCESS;
    XFILE *fp = NULL;
    bool more = true;


    /* install signal handler   */
    signal(SIGINT, INThandler);
    signal(SIGFPE, FPEhandler);

    /* read program options */
    OPTRET optret = read_options(argc, argv);
    switch (optret) {
        case E_FAIL:
            fprintf(stderr, "AYB failed during option read\n");
            ret = EXIT_FAILURE;
            goto cleanup;
            break;

        case E_STOP:
            ret = EXIT_SUCCESS;
            goto cleanup;
            break;

        /* nothing to do otherwise */
        case E_CONTINUE:;
        default:;
    }

    /* create a message log, name includes the prefix argument */
    if (!startup_message(get_pattern())) {
        fprintf(stderr, "AYB failed during message startup\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* scan the input directory */
    if (!startup_dirio()) {
        fprintf(stderr, "AYB failed during file search\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (!startup_model()) {
        fprintf(stderr, "AYB failed during model initialisation\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* process each intensity file until no more or a no continue error */

    while (more) {
        fp = open_next(fp);

        if (xfisnull(fp)) {
            more = false;
        }
        else {
            /* analyse this input file */
            more = analyse_tile(argc, argv, fp);
        }
    }
    fprintf(stdout, "End of AYB\n");

cleanup:
    /* tidy up before exit */
    tidyup();
    return ret;
}
