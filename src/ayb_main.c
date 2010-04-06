/*
 *  File    : ayb_main.c
 *  Created : 23 Feb 2010
 *  Author  : Hazel Marsden
 *  Purpose : Main AYB module
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
#include "handler.h"
#include "message.h"

#include <unistd.h>         // temp for sleep


/* prefix required on log file */
static const char *LOG_PREFIX = "ayb_";


/* tidy up before exit */
void tidyup() {
    /* message files file */
    tidyup_message();
}

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

    /* create a message log */
    start_message(LOG_PREFIX);

    /**** this part is for testing of log and signal handling ****/
    /* convert boolean value to string */
//#define BOOLSTR(b) b ? "true" : "false"
    //  fprintf(stdout, "A value argument is %d\n", p_opt->aval);
    //  fprintf(stdout, "A flag argument is %s\n", BOOLSTR(p_opt->aflag));

    int anint = 3;
    float afloat = 123.45;
    char astring[] = "Hi There";

    MSG1(E_STRING_S, MSG_INFO, astring);
    MSG1(E_FLOAT_F, MSG_WARN, afloat);
    MSG1(E_INT_D, MSG_ERR, anint);
    MSG2(E_INT_DD, MSG_FATAL, anint, 22);

    for (int i = 0; i < 5; i++) {
        printf("Sleep for %d seconds; divisor = %d\n", p_opt->aval, anint);
        sleep(p_opt->aval);
        afloat /= anint--;
        printf("float = %0.2f\n", afloat);
        checksignals();
    }
    /********/

    /* tidy up before exit */
    tidyup();

    fprintf(stdout, "%s\n", "End of AYB Main");
    return EXIT_SUCCESS;
}
