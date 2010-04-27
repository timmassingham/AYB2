/**
 * \file ayb_options.c
 * AYB specific Options.
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
#include <string.h>
#include <getopt.h>
#include "ayb_model.h"
#include "ayb_options.h"
#include "ayb_version.h"
#include "dirio.h"          // I/O for this run hm??
#include "message.h"        // message file location and level


/* constants */
/* none      */

/* private functions that output bulk text */

/** Print help information. Includes text from ayb_help.h. */
void print_help(FILE *fp) {
    validate(NULL!=fp,);
    fputs(""
#include "ayb_help.h"
          , fp);
}

/** Print licence information. Includes text from copyright.h. */
void print_licence(FILE *fp) {
    validate(NULL!=fp,);
    fputs("\n"
PROGNAME " Advanced Base Calling for Next-Generation Sequencing Machines\n"
#include "copyright.h"
          , fp);
}

/** Print usage information. Includes text from ayb_usage.h. */
void print_usage(FILE *fp) {
    validate(NULL!=fp,);
    fputs(""
#include "ayb_usage.h"
          , fp);
}


/* members */

/** AYB options structure; any here hmhm?. */
static AYBOPT Options;

/** Options with no short form. */
enum {OPT_HELP, OPT_LICENCE, OPT_VERSION};

/** Long option structure used by getopt_long. */
static struct option Longopts[] = {
    {"aval",        required_argument,  NULL, 'a'},
    {"aflag",       no_argument,        NULL, 'f'},
    {"ncycles",     required_argument,  NULL, 'n'},
    {"prefix",      required_argument,  NULL, 'x'},
    {"input",       required_argument,  NULL, 'i'},
    {"output",      required_argument,  NULL, 'o'},
    {"logfile",     required_argument,  NULL, 'e'},
    {"loglevel",    required_argument,  NULL, 'l'},
    {"M",           required_argument,  NULL, 'M'},
    {"N",           required_argument,  NULL, 'N'},
    {"P",           required_argument,  NULL, 'P'},
    {"help",        no_argument,        NULL, OPT_HELP },
    {"licence",     no_argument,        NULL, OPT_LICENCE },
    {"version",     no_argument,        NULL, OPT_VERSION },
    {0,0,0,0}
};


/* private functions */

/** Set default values for ayb options defined in this module. */
static void init_options() {
    Options.aflag = false;
    Options.aval = 2;
}


/* public functions */

/** Read options from command line arguments. Uses getopt_long to allow long and short forms. */
bool read_options(const int argc, char ** const argv) {
    bool carryon = true;

    /* set default values */
    init_options();

    /* act on each option in turn */
    int ch;

    while ((ch = getopt_long(argc, argv, "a:fn:x:i:o:e:l:M:N:P:", Longopts, NULL)) != -1){

        switch(ch){
            case 'a':
//                printf("option -a with value `%s'\n", optarg);
                Options.aval = atoi(optarg);
                break;

            case 'f':
//                printf("option -f\n");
                Options.aflag = true;
                break;

            case 'n':
                /* number of cycles */
                set_ncycle(optarg);
                break;

            case 'x':
                /* file pattern match */
                set_pattern(optarg);
                break;

            case 'i':
                /* input file location */
                set_location(optarg, E_INPUT);
                break;

            case 'o':
                /* output file location */
                set_location(optarg, E_OUTPUT);
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

            case 'M':
                /* crosstalk file name */
                set_location(optarg, E_CROSSTALK);
                break;

            case 'N':
                /* crosstalk file name */
                set_location(optarg, E_NOISE);
                break;

            case 'P':
                /* crosstalk file name */
                set_location(optarg, E_PHASING);
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

/** Return a pointer to the options structure. */
AYBOPT *myopt() {
    return &Options;
}

