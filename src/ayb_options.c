/*
 *  File    : ayb_options.2
 *  Created : 23 Feb 2010
 *  Author  : Hazel Marsden
 *  Purpose : AYB specific Options
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
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include "ayb_options.h"
#include "ayb_version.h"
#include "message.h"        // message file location and level


/* print lines of help */
void print_help(FILE *fp) {
    //    validate(NULL!=fp,);
    fputs(""
#include "ayb_help.h"
          , fp);
}

void print_licence(FILE *fp) {
    //    validate(NULL!=fp,);
    fputs("\n"
PROGNAME " Advanced Base Calling for Next-Generation Sequencing Machines\n"
#include "copyright.h"
          , fp);
}

/* print lines of help */
void print_usage(FILE *fp) {
    //    validate(NULL!=fp,);
    fputs(""
#include "ayb_usage.h"
          , fp);
}


/* create options structure */
static AYBOPT Options;

/* set default values for ayb options defined here */
void init_options() {
    Options.aflag = false;
    Options.aval = 2;
}


/* options with no short form */
enum {OPT_HELP, OPT_LICENCE, OPT_VERSION};

/* option structure for getopt_long */
static struct option Longopts[] = {
    {"aval",        required_argument,  NULL, 'a'},
    {"aflag",       no_argument,        NULL, 'f'},
    {"logfile",     required_argument,  NULL, 'e'},
    {"loglevel",    required_argument,  NULL, 'l'},
    {"help",        no_argument,        NULL, OPT_HELP },
    {"licence",     no_argument,        NULL, OPT_LICENCE },
    {"version",     no_argument,        NULL, OPT_VERSION },
    {0,0,0,0}
};

/* public functions */

/* read options from command line arguments */
bool read_options(const int argc, char ** const argv) {
    bool carryon = true;

    /* set default values */
    init_options();

    /* act on each option in turn */
    int ch;

    while ((ch = getopt_long(argc, argv, "a:fe:l:", Longopts, NULL)) != -1){

        switch(ch){
            case 'a':
//                printf("option -a with value `%s'\n", optarg);
                Options.aval = atoi(optarg);
                break;

            case 'f':
//                printf("option -f\n");
                Options.aflag = true;
                break;

            case 'e':
                /* message file location */
                set_message_path(optarg);
                break;

            case 'l':
                /* message output level */
                if (!set_message_level(optarg)) {
                    fprintf(stderr, "Unrecognised error level option: \'%s\'\n\n", optarg);
                     carryon = false;
                }
                break;

            case OPT_HELP:
                print_usage(stderr);
                print_help(stderr);
            	carryon = false;
                break;

            case OPT_LICENCE:
                print_licence(stderr);
                carryon = false;
                break;

            case OPT_VERSION:
                fprintf( stderr, "\n" PROGNAME " Version %0.2f  %s\n\n", Version, VersionDate);

                carryon = false;
                break;

            default:
            	// getopt_long outputs an error message
                print_usage(stderr);
            	carryon = false;
        }
    }

    return carryon;
}

/* return a pointer to the options structure */
AYBOPT *myopt() {
    return &Options;
}

