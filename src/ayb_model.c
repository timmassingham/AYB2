/** 
 * \file ayb_model.c
 * Setup, read intensities, loop and output.
 * Used as a singleton class with local functions accessing global member data.
 *//* 
 *  Created : 14 Apr 2010
 *  Author  : Tim Massingham/Hazel Marsden
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

#include <math.h>
#include <string.h>
#include "ayb.h"
#include "ayb_model.h"
#include "ayb_options.h"
#include "ayb_version.h"
#include "datablock.h"
#include "matrix.h"
#include "message.h"
#include "mixnormal.h"
#include "statistics.h"
#include "tile.h"
#include "weibull.h"


/* constants */

static const unsigned int MIN_CYCLE = 2;        ///< Minimum cycles for modelling.

/** Possible output format text. Match to OUTFORM enum. Used to match program argument and also as file extension. */
static const char *OUTFORM_TEXT[] = {"fasta", "fastq"};
/** New cluster symbol in sequence file. Match to OUTFORM enum. */
static const int OUT_SYMBOL[] = {'>', '@'};

/* members */

/** Possible output formats, and number. */
typedef enum OutFormT {E_FASTA, E_FASTQ, E_OUTFORM_NUM} OUTFORM;
static OUTFORM OutputFormat  = E_FASTQ;         ///< Selected output format.

static unsigned int NIter = 5;                  ///< Number of iterations in base call loop.
static unsigned int *ZeroLambda = NULL;         ///< Count of zero lambdas before base call, per iteration.
static bool SimData = false;                    ///< Set to output simulation data.
static CSTRING SimText = NULL;                  ///< Header text for simulation data file.
static TILE MainTile = NULL;                    ///< Tile data from file or run-folder.

/* Additional data size constraint for debug output, set in analyse_tile */
static bool ShowDebug = false;


/* private functions */

/** Create the sub-tile datablocks to be analysed. */
static TILE * create_datablocks(const TILE maintile, const unsigned int numblock) {

    TILE * tileblock = NULL;
    if (get_defaultblock()) {
        /* copy all into single block */
        tileblock = calloc(1, sizeof(*tileblock));
        if(tileblock == NULL) {return NULL;}
        tileblock[0] = copy_TILE(maintile);
    }

    else {
        /* create the sub-tiles with an array of pointers */
        tileblock = calloc(numblock, sizeof(*tileblock));
        if(tileblock == NULL) {return NULL;}

        int blk = 0;
        int colstart = 0;
        int colend = 0;

        /* block specification already decoded */
        DATABLOCK datablock = get_next_block();
        while (datablock != NULL) {
            colend = colstart + datablock->num - 1;
            switch (datablock->type) {
            case E_READ :
                /* new block if not first */
                if (tileblock[blk] != NULL) {
                    blk++;
                }

            /* no break; fall through to case CONCAT */
            case E_CONCAT :
                tileblock[blk] = copy_append_TILE(tileblock[blk], maintile, colstart, colend);
                break;

            case E_IGNORE :
                /* just increment column pointers, done outside of switch */
                break;

            default :;
            }

            colstart = colend + 1;
            datablock = get_next_block();
        }
    }

    return tileblock;
}

/** Output message with counts if any zero lambdas. */
static void output_zero_lambdas(void) {

static const int MAX_ZEROS = 1e6 - 1;       // Up to 6 digits
static const int MAX_NUMLEN = 9;            // Enough for "BigNum, \0" or "999999, \0"
static const char *BIG_NUM = "BigNum";      // Use if number too big for buffer

    /* check if any recorded */
    bool any = false;
    for (int i = 0; i < NIter; i++) {
        if (ZeroLambda[i] > 0) {
            any = true;
            break;
        }
    }

    if (any) {
        /* create the message as a string to allow variable iterations */
        char numstring[MAX_NUMLEN];
        char msgstring[MAX_NUMLEN * NIter];

        /* first one has no preceding comma */
        if (ZeroLambda[0] > MAX_ZEROS) {
            sprintf(msgstring, "%s", BIG_NUM);
        }
        else {
            sprintf(msgstring, "%d", ZeroLambda[0]);
        }

        for (int i = 1; i < NIter; i++) {
            if (ZeroLambda[i] > MAX_ZEROS) {
                sprintf(numstring, ", %s", BIG_NUM);
            }
            else {
                sprintf(numstring, ", %d", ZeroLambda[i]);
            }
            strcat(msgstring, numstring);
        }
        message(E_ZERO_LAMBDA_S, MSG_INFO, msgstring);
    }
}

/**
 * Output the results of the base calling.
 * Returns true if output file opened ok.
 */
static RETOPT output_results (const AYB ayb, const int blk) {

    if (ayb == NULL) {return E_FAIL;}

    XFILE *fpout = NULL;
    /* different rules for varying input formats */
    switch (get_input_format()) {
        case E_TXT:
            fpout = open_output_blk("seq", blk);
            break;

        case E_CIF:
            fpout = open_output_blk((CSTRING)OUTFORM_TEXT[OutputFormat], blk);
            break;

        default: ;
    }

    if (xfisnull(fpout)) {return E_STOP;}

    const uint32_t ncluster = get_AYB_ncluster(ayb);

    for (uint32_t cl = 0; cl < ncluster; cl++){
        xfprintf(fpout, "%ccluster_%u\n", OUT_SYMBOL[OutputFormat], cl + 1);
        show_AYB_bases(fpout, ayb, cl);
        /* quality score */
        if (OutputFormat == E_FASTQ) {
            xfputs("\n+\n", fpout);
            show_AYB_quals(fpout, ayb, cl);
        }
        xfputc('\n', fpout);
    }
    xfclose(fpout);
    return E_CONTINUE;
}

/** Return true if the string contains any spaces. */
static bool has_whitespace(const char * str) {

    size_t len = strlen(str);
    for (int i = 0; i < len; i++) {
        if (str[i] == ' ') {
            return true;
        }
    }
    return false;
}

/**
 * Make any header formatting adjustments.
 * Check for any new lines and ensure they are followed by a hash.
 */
static CSTRING format_header (CSTRING simtext) {

    /* count newlines so can resize string */
    size_t oldlen = strlen(simtext);
    size_t count = 0;
    for (int i = 0; i < oldlen; i++) {
        if (simtext[i] == '\n') {
            count++;
        }
    }

    if (count > 0) {
        CSTRING new = new_CSTRING(oldlen + count);

        int idx = 0;
        char *token = strtok(simtext, "\n");
        while ((token != NULL) && (idx <= count)) {
            if (idx++ == 0) {
                strcpy(new, token);
            }
            else {
                strcat(new, "\n#");
                strcat(new, token);
            }
            token = strtok(NULL, "\n");
        }

        free_CSTRING(simtext);
        simtext = new;
    }
    return simtext;
}

/** Output data for use by simulator. */
static void output_simdata(AYB ayb, const int argc, char ** const argv, const int blk) {

    enum FitT {E_LOGISTIC, E_MIXED, E_NORMAL, E_WEIBULL};       // possible distributions to fit
    const enum FitT fitdist = E_LOGISTIC;                       // fixed within this version of AYB
    const unsigned int MIX_ITER = 100;
    const unsigned int MIX_NUM = 3;
    const unsigned int SIM_VERSION = 5;

    XFILE *fpsim = open_output_blk("runfile", blk);
    if (xfisnull(fpsim)) {return;}

    /* header required if a single or first block */
    if (blk != BLK_APPEND) {
        /* header text followed by AYB version */
        SimText = format_header(SimText);
        xfprintf(fpsim, "# %s\n", SimText);
        xfprintf(fpsim, "# AYB Version %0.2f  %u\n", get_version(), get_version_date());

        /* add the command line */
        xfputs("# ", fpsim);
        xfputs(argv[0], fpsim);
        bool sfound = false;
        for (int i = 1; i < argc; i++){
            xfputc(' ', fpsim);
            /* substitute for header as already printed above */
            if (sfound) {
                xfputs("\"header\"", fpsim);
                sfound = false;
            }
            else {
                /* add quotes if would have been required for whitespace */
                if (has_whitespace(argv[i])) {
                    xfputc('\"', fpsim);
                    xfputs(argv[i], fpsim);
                    xfputc('\"', fpsim);
                }
                else {
                    xfputs(argv[i], fpsim);
                }

                /* detect simdata option and flag that next is the header */
                if (match_option(argv[i], E_SIMDATA)) {
                    sfound = true;
                }
            }
        }
        xfputs("\n", fpsim);
        
        /* sim version */
        xfprintf(fpsim, "Version %u\n", SIM_VERSION);
    }

    /* get non-zero lambdas as weibull uses log */
    uint32_t num = 0;
    real_t * lambdas = get_AYB_lambdas(ayb, &num);

    /* get parameters for fitted lambda distribution according to fixed selection */
    char fitc = 'X';
    pair_real lambdafit = {NAN, NAN};
    NormMixParam mparam = NULL;
    
    switch (fitdist) {
        case E_LOGISTIC:
            fitc = 'L';
            lambdafit.e1 = mean(lambdas, num);
            /* scale stdev by root 3 over pi */
            lambdafit.e2 = sqrt(variance(lambdas, num) * 3) / M_PI;
            break;

        case E_MIXED:
            fitc = 'M';
            /* mix and iterations fixed */
            mparam = fit_mixnormal(lambdas, num, MIX_NUM, MIX_ITER);
            break;

        case E_NORMAL:
            fitc = 'N';
            lambdafit.e1 = mean(lambdas, num);
            lambdafit.e2 = sqrt(variance(lambdas, num));
            break;

        case E_WEIBULL:
            fitc = 'W';
            lambdafit = fit_weibull(lambdas, num);
            break;

        default: ;
    }

    /* number of cycles and lambda fit parameters */
    if (fitdist == E_MIXED) {
        xfprintf(fpsim, "%u %c %u", get_AYB_ncycle(ayb), fitc, MIX_NUM);
        if (mparam == NULL) {
            message(E_NOCREATE_S, MSG_ERR, "mixed normal distribution");
        }
        else {
            for (int i = 0; i < MIX_NUM; i++) {
                xfprintf(fpsim, " %f %f %f", mparam->prob[i], mparam->mean[i], mparam->sd[i]);
            }
        }
        xfputs("\n", fpsim);
        free_NormMixParam(mparam); 
    }
    else {
        xfprintf(fpsim, "%u %c %f %f\n", get_AYB_ncycle(ayb), fitc, lambdafit.e1, lambdafit.e2);
    }
    xfree(lambdas);

    /* calculate and output all covariance */
    MAT * V = calculate_covariance(ayb, true);
    
    if (V == NULL) {
        message(E_NOCREATE_S, MSG_ERR, "full covariance");
    }
    else {
        show_MAT_rownum(fpsim, V[0], 0, 0, false);
        free_MAT(V[0]);
        xfree(V);
    }

    xfclose(fpsim);
}


/* public functions */

/**
 * Analyse a single tile. Intensities data already read in and stored in MainTile.
 * Returns continue if analysis should continue to next file,
 * else fail if to continue to next prefix, else stop.
 */
RETOPT analyse_tile (const int argc, char ** const argv) {

    if (MainTile == NULL) {
        if (!run_folder()) {
            /* if not a run-folder then problem caused by bad input file */
            message(E_BAD_INPUT_S, MSG_ERR, get_current_file());
        }
        return E_CONTINUE;
    }

    if (MainTile->ncycle < get_totalcycle()) {
        /* not enough data */
        message(E_CYCLESIZE_DD, MSG_ERR, MainTile->ncycle, get_totalcycle());
        MainTile = free_TILE(MainTile);
        return E_FAIL;
    }
    else if (MainTile->ncycle < MIN_CYCLE) {
        message(E_CYCLESIZE_D, MSG_ERR, MainTile->ncycle);
        MainTile = free_TILE(MainTile);
        return E_FAIL;
    }
    else {
        message(E_TILESIZE_DD, MSG_INFO, MainTile->ncluster, MainTile->ncycle);

        ShowDebug = false;
#ifndef NDEBUG
    /* set additional debug output, can overload the process if too much data; vary as required */
    ShowDebug = ((NIter == 1) && (MainTile->ncycle <= 20) && (MainTile->ncluster <= 100));
#endif
    }

    const unsigned int ncluster = MainTile->ncluster;
    const unsigned int numblock =  get_defaultblock() ? 1 : get_numblock();
    RETOPT status = E_CONTINUE;

    /* put the data into distinct blocks */
    TILE * tileblock = NULL;
    tileblock = create_datablocks(MainTile, numblock);
    /* no longer need the raw data as read in */
    MainTile = free_TILE(MainTile);

    if (tileblock == NULL) {
        message(E_DATABLOCK_FAIL_S, MSG_FATAL, get_current_file());
        return E_FAIL;
    }

    /* analyse each tile block separately */
    AYB ayb = NULL;
    for (int blk = 0; blk < numblock; blk++) {

        ayb = new_AYB(tileblock[blk]->ncycle, ncluster);
        if (ayb == NULL) {
            message(E_NOMEM_S, MSG_FATAL, "model structure creation");
            message(E_INIT_FAIL_DD, MSG_ERR, blk + 1, tileblock[blk]->ncycle);
            status = E_FAIL;
            goto cleanup;
        }

        /* get next tile block of raw intensities */
        ayb = replace_AYB_tile(ayb, tileblock[blk]);

        /* set initial model values */
        if (initialise_model(ayb, ShowDebug)) {
            message(E_PROCESS_DD, MSG_INFO, blk + 1, get_AYB_ncycle(ayb));

#ifndef NDEBUG
    if (ShowDebug) {
        XFILE *fpout = NULL;
        fpout = open_output_blk("ayb1", (numblock > 1) ? blk : BLK_SINGLE);
        if (!xfisnull(fpout)) {
            show_AYB(fpout, ayb, true);
        }
        xfclose(fpout);
    }
#endif

            /* base calling loop */
            int res;
            for (int i = 0; i < NIter; i++){
                xfprintf(xstdout, "Iteration: %d\n", i+1);

                estimate_MPN(ayb);

                /* parameters to estimate bases are block index and flag to indicate last iteration */
                /* return is number of zero lambdas or error */
                res = estimate_bases(ayb, (numblock > 1) ? blk : BLK_SINGLE, (i == (NIter - 1)), ShowDebug);

                if (res == DATA_ERR) {
                    /* terminate processing */
                    message(E_PROCESS_FAIL_D, MSG_ERR, i + 1);
                    break;
                }
                else {
                    ZeroLambda[i] = res;
                }
            }

            /* output any zero lambdas */
            output_zero_lambdas();

            /* output the results */
            status = output_results(ayb, (numblock > 1) ? blk : BLK_SINGLE);
            
            /* output simulation data if requested */
            if (SimData) {
                /* block indicator is append if not first of multiple blocks otherwise single */
                output_simdata(ayb, argc, argv, ((numblock > 1) && (blk > 0)) ? BLK_APPEND : BLK_SINGLE);
            }
        }

        else {
            message(E_INIT_FAIL_DD, MSG_ERR, blk + 1, get_AYB_ncycle(ayb));
            status = E_FAIL;
        }

        /* free the structure ready for next */
        ayb = free_AYB(ayb);
        if (status != E_CONTINUE) {break;}
    }

cleanup:
    if (tileblock != NULL) {
        for (int blk = 0; blk < numblock; blk++)  {
            free_TILE(tileblock[blk]);
        }
        xfree(tileblock);
    }
    return status;
}

/**
 * Read and store a single intensities input file.
 */
void read_intensities_file(XFILE *fp, unsigned int ncycle) {

    /* should always be null on entry, but check anyway */
    if (MainTile != NULL) {
        MainTile = free_TILE(MainTile);
    }

    switch (get_input_format()) {
        case E_TXT:
            MainTile = read_TILE(fp, ncycle);
            break;

        case E_CIF:
            MainTile = read_cif_TILE (fp, ncycle);
            break;

        default: ;
    }
}

/**
 * Read and store a single lane/tile of intensities from a run-folder.
 */
void read_intensities_folder(const char *root, const LANETILE lanetile, unsigned int ncycle) {

    /* should always be null on entry, but check anyway */
    if (MainTile != NULL) {
        MainTile = free_TILE(MainTile);
    }

    MainTile = read_folder_TILE(root, lanetile.lane, lanetile.tile, ncycle);
}

/** Set the number of base call iterations. */
void set_niter(const char *niter_str) {

    char *endptr;
    NIter = strtoul(niter_str, &endptr, 0);
}

/**
 * Set the output format. Text must match one of the output format text list. Ignores case.
 * Returns true if match found.
 */
bool set_output_format(const char *outform_str) {

    /* match to one of the possible options */
    int matchidx = match_string(outform_str, OUTFORM_TEXT, E_OUTFORM_NUM);
    if (matchidx >= 0) {
        OutputFormat = (OUTFORM)matchidx;
        return true;
    }
    else {
        return false;
    }
}

/** Set simdata flag and text for file. */
void set_simdata(const CSTRING simdata_str) {
    SimData = true;
    SimText = copy_CSTRING(simdata_str);
}

/**
 * Start up; call at program start after options.
 * Issues option info messages.
 * Returns true if cycle blocks and iterations parameters are ok
 * and ayb startup is successful.
 */
bool startup_model(void) {

    message(E_OPT_SELECT_SS, MSG_INFO, "Output format" ,OUTFORM_TEXT[OutputFormat]);

    /* check number of cycles and data blocks supplied */
    const unsigned int totalcycle = get_totalcycle();
    const unsigned int numblock = get_numblock();

    if (get_defaultblock()) {
        message(E_DEFAULTBLOCK, MSG_INFO);
    }
    else {
        if ((totalcycle == 0) || (numblock == 0)) {
            message(E_NOBLOCKS, MSG_FATAL);
            return false;
        }
        if (totalcycle < MIN_CYCLE) {
            message(E_CYCLESIZE_D, MSG_FATAL, totalcycle);
            return false;
        }
        message(E_OPT_SELECT_SD, MSG_INFO, "cycles total", totalcycle);
        message(E_OPT_SELECT_SD, MSG_INFO, "distinct data blocks", numblock);
    }

    /* check number of iterations supplied - may be left at default */
    if (NIter == 0) {
        message(E_BAD_ITER, MSG_FATAL);
        return false;
    }

    /* storage for zero lambda count */
    ZeroLambda = calloc(NIter, sizeof(int));

    message(E_OPT_SELECT_SD, MSG_INFO, "iterations", NIter);

    return startup_ayb();
}

/** Tidy up; call at program shutdown. */
void tidyup_model(void) {

    /* free memory */
    SimText = free_CSTRING(SimText);
    xfree(ZeroLambda);
    tidyup_ayb();
}
