/**
 * \file ayb_main.c
 * Main AYB module.
 * Sets up environment before looping through each supplied pattern,
 * passing control to ayb_model for each intensity file found.
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
#include <string.h>
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
static void tidyup(void) {
    tidyup_model();
    tidyup_dirio();
    tidyup_datablock();
    tidyup_message();
}

/* public functions */

/**
 * Main AYB routine.
 * Read program options and set up environment.
 * Then process each input file or run-folder lane/tile before tidy up and exit.
 */
int main(int argc, char **argv) {

    int ret = EXIT_SUCCESS;
    int nextarg = 0;
    LANETILE lanetile = {0, 0};
    XFILE *fp = NULL;

    /* install signal handler - interrupt does not work as written, do not enable */
//    signal(SIGINT, INThandler);
    signal(SIGFPE, FPEhandler);

    /* read program options */
    RETOPT status = read_options(argc, argv, &nextarg);
    switch (status) {
        case E_FAIL:
            fprintf(stderr, "AYB failed during option read\n");
            ret = EXIT_FAILURE;
            goto cleanup;
            break;

        case E_STOP:
            ret = EXIT_SUCCESS;
            goto cleanup;
            break;

        case E_CONTINUE:
            /* check if any non-option pattern match arguments supplied */
            if (nextarg >= argc) {
                message(E_NOPATTERN, MSG_ERR);
                ret = EXIT_FAILURE;
                goto cleanup;
                break;
            }

        default:;
    }

    /* create a message log, stderr or name as argument */
    if (!startup_message()) {
        fprintf(stderr, "AYB failed during message startup\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* check the specified i/o locations */
    if (!startup_dirio()) {
        fprintf(stderr, "AYB failed during i/o check\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* model related startup */
    if (!startup_model()) {
        fprintf(stderr, "AYB failed during model initialisation\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    const CSTRING input_path = get_input_path();
    const unsigned int totalcycle = get_totalcycle();

    /* process each prefix or lane/tile range supplied as non-option argument */
    for (int i = nextarg; i < argc; i++) {
        if (set_pattern(argv[i])) {
            /* process each intensity file or run-folder lane/tile until no more or a no continue error */
            status = E_CONTINUE;

            while (status == E_CONTINUE) {

                if (run_folder()) {
                    /* get lane/tile numbers */
                    lanetile = get_next_lanetile();
                    if (lanetile_isnull(lanetile)) {
                        /* next prefix */
                        status = E_FAIL;
                    }
                }
                else {
                    /* open the intensities file */
                    fp = open_next(fp);
                    if (xfisnull(fp)) {
                        /* next prefix */
                        status = E_FAIL;
                    }
                }

                if (status == E_CONTINUE) {
                    if (run_folder()) {
                        /* read intensities data from run-folder */
                        read_intensities_folder(input_path, lanetile, totalcycle);
                    }
                    else {
                        /* read intensities data from supplied file */
                        read_intensities_file(fp, totalcycle);
                    }

                    /* analyse the stored tile */
                    status = analyse_tile(argc, argv);
                }
            }
            /* analysis return may indicate stop program */
            if (status == E_STOP) {break;}
        }
    }

    /* get program executable name */
    char *pname = strrchr(argv[0], PATH_DELIM);
    if (pname == NULL) {
        pname = argv[0];
    }
    else {
        pname++;
    }
    fprintf(stdout, "End of %s\n", pname);

cleanup:
    /* tidy up before exit */
    tidyup();
    return ret;
}
