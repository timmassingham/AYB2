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
#include "ayb.h"
#include "ayb_model.h"
#include "ayb_options.h"
#include "ayb_version.h"
#include "call_bases.h"
#include "datablock.h"
#include "dirio.h"
#include "message.h"
#include "qual_table.h"


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


/* constants */
/* none */

/* members */
static int NThread = 1;


/** Options with no short form. */
enum {OPT_HELP, OPT_LICENCE, OPT_VERSION};

/** Long option structure used by getopt_long. */
static struct option Longopts[] = {
    {"simdata",     required_argument,  NULL, 's'},   // Note!! index identified as E_SIMDATA = 0 in header file
    {"blockstring", required_argument,  NULL, 'b'},
    {"dataformat",  required_argument,  NULL, 'd'},
    {"generr",      required_argument,  NULL, 'g'},
    {"logfile",     required_argument,  NULL, 'e'},
    {"format",      required_argument,  NULL, 'f'},
    {"input",       required_argument,  NULL, 'i'},
    {"spikeuse",    no_argument,        NULL, 'k'},
    {"loglevel",    required_argument,  NULL, 'l'},
    {"mu",          required_argument,  NULL, 'm'},
    {"niter",       required_argument,  NULL, 'n'},
    {"output",      required_argument,  NULL, 'o'},
    {"parallel",    required_argument,  NULL, 'p'},
    {"noqualout",   no_argument,        NULL, 'q'},
    {"runfolder",   no_argument,        NULL, 'r'},
    {"working",     required_argument,  NULL, 'w'},
    {"A",           required_argument,  NULL, 'A'},
    {"spikein",     required_argument,  NULL, 'K'},
    {"M",           required_argument,  NULL, 'M'},
    {"N",           required_argument,  NULL, 'N'},
    {"qualtab",     required_argument,  NULL, 'Q'},
    {"help",        no_argument,        NULL, OPT_HELP },
    {"licence",     no_argument,        NULL, OPT_LICENCE },
    {"license",     no_argument,        NULL, OPT_LICENCE },
    {"version",     no_argument,        NULL, OPT_VERSION },
    {0,0,0,0}
};


/* private functions */

/** Set default values for ayb options defined in this module. */
static void init_options(void) {
}

/** 
 * Set the requested number of parallel threads.
 * Do not allow to be invalid.
 */
static void set_nthread(const char *n_str) {

    char *endptr;
    long n = strtol(n_str, &endptr, 0);
    if (n > 0) {
        NThread = n;
    }
    else {
        fprintf(stderr, "Warning: Invalid number of threads (\'%s\') supplied; defaulting to %d\n", n_str, NThread);
    }
}


/* public functions */

/**
 * Read options from command line arguments. Uses getopt_long to allow long and short forms.
 * Returns index to first non-option argument as reference parameter.
 * Returns whether to continue, stop and indicate error or just stop.
 */
RETOPT read_options(const int argc, char ** const argv, int *nextarg) {
    RETOPT status = E_CONTINUE;

    /* set default values */
    init_options();

    /* act on each option in turn */
    int ch;

    while ((ch = getopt_long(argc, argv, "s:b:d:e:f:g:i:kl:m:n:o:p:qrw:A:K:M:N:Q:", Longopts, NULL)) != -1){

        switch(ch){
            case 's':
                 /* output simulation data */
                 set_simdata(optarg);
                 break;

            case 'b':
                /* pattern of data blocks */
                if (!parse_blockopt(optarg)) {
                    status = E_FAIL;
                }
                break;

            case 'd':
                /* input format */
                if (!set_input_format(optarg)) {
                    fprintf(stderr, "Fatal: Unrecognised input format option: \'%s\'\n\n", optarg);
                    status = E_FAIL;
                }
                break;

            case 'e':
                /* message file location */
                set_message_path(optarg);
                break;

            case 'f':
                /* output format */
                if (!set_output_format(optarg)) {
                    fprintf(stderr, "Fatal: Unrecognised output format option: \'%s\'\n\n", optarg);
                    status = E_FAIL;
                }
                break;

            case 'g':
                /* new quality calculation */
                if (!set_generr(optarg)) {
                    fprintf(stderr, "Fatal: Generalised error must be a positive value; \'%s\' supplied\n\n", optarg);
                    status = E_FAIL;
                }
                break;

            case 'i':
                /* input file location */
                set_location(optarg, E_INPUT);
                break;

            case 'k':
                /* spike-in data calibration flag */
                set_spike_calib();
                break;

            case 'l':
                /* message output level */
                if (!set_message_level(optarg)) {
                    fprintf(stderr, "Fatal: Unrecognised error level option: \'%s\'\n\n", optarg);
                    status = E_FAIL;
                }
                break;

            case 'm':
                /* phredchar calculation */
                if (!set_mu(optarg)) {
                    fprintf(stderr, "Fatal: Mu must be a positive value; \'%s\' supplied\n\n", optarg);
                    status = E_FAIL;
                }
                break;

            case 'n':
                /* number of base call iterations */
                set_niter(optarg);
                break;

            case 'o':
                /* output file location */
                set_location(optarg, E_OUTPUT);
                break;

            case 'p':
                /* requested number of parallel threads */
                set_nthread(optarg);
                break;

            case 'q':
                /* no quality calibration table output */
                set_noqualout();
                break;

            case 'r':
                /* input from run-folder */
                set_run_folder();
                break;

            case 'w':
                /* show working flag */
                set_show_working(optarg);
                break;

            case 'A':
                /* param A file name */
                set_location(optarg, E_PARAMA);
                break;

            case 'K':
                /* location of spike-in data */
                set_location(optarg, E_SPIKEIN);
                break;

            case 'M':
                /* initial crosstalk file name */
                set_location(optarg, E_CROSSTALK);
                break;

            case 'N':
                /* noise file name */
                set_location(optarg, E_NOISE);
                break;

            case 'Q':
                /* quality calibration conversion table file location */
                set_location(optarg, E_QUALTAB);
                break;

            case OPT_HELP:
                print_usage(stderr);
                print_help(stderr);
            	status = E_STOP;
                break;

            case OPT_LICENCE:
                print_licence(stderr);
                status = E_STOP;
                break;

            case OPT_VERSION:
                fprintf(stderr, "\n" PROGNAME " Version %0.2f  %u\n\n", get_version(), get_version_date());

                status = E_STOP;
                break;

            default:
            	// getopt_long outputs an error message
                print_usage(stderr);
            	status = E_FAIL;
        }
    }

    /* return index to non-option arguments */
    *nextarg = optind;
    return status;
}

/** Return the requested number of parallel threads. */
int get_nthread(void) {

    return NThread;
}

/**
 * Return true if supplied string matches long or short form of supplied option structure index.
 * Only some indexes identified in OptIndexT enum.
 */
bool match_option(const char *string, const OPTINDEX index) {

    bool ret = false;

    /* must start with option indicator */
    if (*string == '-') {
        if (*(string + 1) == '-') {
            /* long form */
            if (strcmp(string + 2, Longopts[index].name) == 0) {
                ret = true;
            }
        }
        else {
            /* short form */
            if (*(string + 1) == Longopts[index].val) {
                ret = true;
            }
        }
    }
    return ret;
}
