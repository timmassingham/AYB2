/**
 * \file message.c
 * General Messaging Utility.
 * This provides a central messaging system for the output of program information at various levels up to debugging.
 * The level of message output is selected using a program option (default Warning).
 * 
 * Program messages are output to stderr which can be redirected in a run script.
 * Alternatively a file path can be specified as a program option and a message file will be created.
 * 
 * A program message is implemented by calling the message() function with a type and severity parameter,
 * followed by a variable number of additional parameters of varying type corresponding to the selected message type.
 * The type enumerations include a suffix to indicate the parameters required.
 * The varying parameters are handled using vfprintf which takes a va_list argument. 
 * 
 * New message types can be added as needed by extending the MsgTypeT enumeration (#MSGTYPE) and adding the appropriate message text. 
 * Message text is stored in the #MSG_TEXT char array which is formatted for printf output.
 * The text is matched to the correct message type enumeration by the order. 
 *//*
 *  Created : 26 Feb 2010
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
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "dirio.h"
#include "message.h"


/* constants */

/**
 * Formatted text for each message type. Used directly in a printf.
 * Ensure list matches MsgTypeT enum.
 */
static const char *MSG_TEXT[] = {
        "All available cycles in one block selected\n",                         // E_DEFAULTBLOCK
        "Blockstring option contains no datablocks\n",                          // E_NOBLOCKS
        "No file pattern match supplied\n",                                     // E_NOPATTERN
        "Number of model iterations incorrectly supplied\n",                    // E_BAD_ITER
        "Need to supply neither or both parameter matrices (A and N)\n",        // E_BAD_MATRIXIN
        "Run folder option invalid with txt input format\n",                    // E_BAD_RUNOPT
        "",                                                                     // E_END_NONE
        "Memory allocation failed during %s\n",                                 // E_NOMEM_S
        "Log message output level: %s\n",                                       // E_MSG_LEVEL_S
        "Input from directory: %s\n",                                           // E_INPUT_DIR_S
        "Output to directory: %s\n",                                            // E_OUTPUT_DIR_S
        "Input file found: %s\n",                                               // E_INPUT_FOUND_S
        "No filename supplied in pattern: \'%s\'\n",                            // E_NOPATTERN_FILE_S
        "Failed to read input file: %s\n",                                      // E_BAD_INPUT_S
        "Failed to create data blocks for input file: %s\n",                    // E_DATABLOCK_FAIL_S
        "Failed to initialise %s matrix\n",                                     // E_MATRIX_FAIL_S
        "Failed to create %s\n",                                                // E_NOCREATE_S
        "Zero lambdas per iteration: %s\n",                                     // E_ZERO_LAMBDA_S
        "",                                                                     // E_END_S
        "Supplied %s location parameter \'%s\' is not a directory\n",           // E_BAD_DIR_SS
        "Failed to create new %s directory \'%s\'\n",                           // E_NOCREATE_DIR_SS
        "Created new %s directory: %s\n",                                       // E_CREATED_DIR_SS
        "Supplied %s has incorrect file format near item: %s\n",                // E_BAD_INPUT_SS
        "No input files in directory \'%s\' matching pattern: \'%s\'\n",        // E_NOINPUT_SS
        "%s file failed to open: %s\n",                                         // E_OPEN_FAIL_SS
        "Lane tile range selected: lanes: %s, tiles %s\n",                      // E_LANETILE_SS
        "%s selected: %s\n",                                                    // E_OPT_SELECT_SS
        "%s error; %s\n",                                                       // E_BAD_TXT_SS
        "%s contains invalid numeric: \'%s\'\n",                                // E_BAD_NUM_SS
        "",                                                                     // E_END_SS
        "%s contains invalid character: \'%c\'\n",                              // E_BAD_CHAR_SC
        "",                                                                     // E_END_SC
        "Input file pattern match: \'%s\'; %d files found\n",                   // E_PATTERN_MATCH_SD
        "Number of %s selected: %d\n",                                          // E_OPT_SELECT_SD
        "",                                                                     // E_END_SD
        "%s matrix wrong size, need dimension %d not %d\n",                     // E_MATRIXINIT_SDD
        "",                                                                     // E_END_SDD
        "%s selected: %G\n",                                                    // E_OPT_SELECT_SG
        "",                                                                     // E_END_SG
        "Unrecognised nucleotide \'%c\'; returning NUC_AMBIG\n",                // E_BAD_NUC_C
        "",                                                                     // E_END_C
        "Processing failed at iteration %d; calls set to null\n",               // E_PROCESS_FAIL_D
        "Insufficient cycles for model; %d selected or found\n",                // E_CYCLESIZE_D
        "",                                                                     // E_END_D
        "Input file contains fewer cycles than requested; %d instead of %d\n",  // E_CYCLESIZE_DD
        "Tile data size: %d clusters of %d cycles\n",                           // E_TILESIZE_DD
        "Failed to initialise model for block %d, %d cycles\n",                 // E_INIT_FAIL_DD
        "Processing block %d, %d cycles\n",                                     // E_PROCESS_DD
        "",                                                                     // E_END_DD

        "%s %20s\n",                                                            // E_GENERIC_SS
        "%s %d\n",                                                              // E_GENERIC_SD
        "%s %u\n",                                                              // E_GENERIC_SU
        "%s %lx\n",                                                             // E_GENERIC_SX
        "%s %f\n",                                                              // E_GENERIC_SF
        "%s %0.2f %0.2f %0.2f %0.2f\n",                                         // E_GENERIC_SFFFF

        "%s (%s:%d): %s\n",                                                     // E_DEBUG_SSD_S
        "%s (%s:%d): %s %s\n",                                                  // E_DEBUG_SSD_SS
        "%s (%s:%d): %s %d\n"                                                   // E_DEBUG_SSD_SD
        };

/**
 * Message severity text. Used to match program argument and as text in log file.
 * Ensure list matches MsgSeverityT enum. 
 */
static const char *MSG_SEV_TEXT[] = {
        "None",
        "Fatal",
        "Error",
        "Information",
        "Warning",
        "Debug",
        ""};

/* members */

/** Selected level of messages, default Warning. */
static int Msg_Level = MSG_WARN;
/**
 * Selected pathname for message file.
 * Output is to stderr unless redirected by a program option.
 */
static CSTRING Msg_Path = NULL;


/* private functions */

/** Check path part of file path. Try to create if does not exist. */
static bool check_path(const CSTRING filepath) {

    if (filepath == NULL) {return false;}
    bool ret = true;

    /* find the end of any path */
    char *pdlm = strrchr(filepath, PATH_DELIM);
    if (pdlm != NULL) {
        /* path included */
        size_t len = pdlm - filepath;
        CSTRING path = new_CSTRING(len);
        strncpy(path, filepath, len);

        /* check specified path exists or can be created */
        if (!check_outdir(path, "message")) {
            ret = false;
        }
        free_CSTRING(path);
    }
    return ret;
}

/**
 * Generate a randomised virtual message name of form ayb_xxxxxx_yymmdd_hhmm.
 * Use urandom instead of rand function because typical seeding with current time
 * is insufficient to ensure different values when many runs started simultaneously.
 * Returns NULL if any error.
 */
static CSTRING generate_name(const struct tm *p_tm) {

    const char *RAND_FILE = "/dev/urandom";         // random number file
    const unsigned int DIVISOR = 1E6;               // limit random value to six digits
    const char *LOG_START = "ayb_%06u";             // includes random value template
    const size_t LEN = 10;                          // for prefix and random value
    const char *TIME_SUFFIX = "_%y%m%d_%H%M";       // gives _yymmdd_hhmm
    const size_t TIME_SUFFIX_LEN = 13;              // maximum length, including null terminator
    char timestring[TIME_SUFFIX_LEN];
    FILE *frand = NULL;

    CSTRING name = new_CSTRING(LEN + TIME_SUFFIX_LEN);
    if (name == NULL) {
        message(E_NOMEM_S, MSG_FATAL, "message name creation");
        goto cleanup;
    }

    /* start with default prefix and random number string of limited digits */
    unsigned int number = 0;
    frand = fopen(RAND_FILE, "r");
    if (frand == NULL) {
        message(E_OPEN_FAIL_SS, MSG_FATAL, "Random number generator", RAND_FILE);
        goto cleanup;
    }
    else {
        fread(&number, 1, sizeof(number), frand);
        fclose(frand);
        if (number == 0) {
            message(E_BAD_INPUT_S, MSG_FATAL, RAND_FILE);
            goto cleanup;
        }
    }

    unsigned int val = number % DIVISOR;
    sprintf(name, LOG_START, val);

    /* add the current date and time */
    strftime(timestring, TIME_SUFFIX_LEN, TIME_SUFFIX, p_tm);
    strcat(name, timestring);

    return name;

cleanup:
    free_CSTRING(name);
    return NULL;
}


/* public functions */

/** Return the selected message level. */
MSGSEV get_message_level(void) {

    return (MSGSEV)Msg_Level;
}

/** Return the selected message path. */
CSTRING get_message_path(void) {

    return Msg_Path;
}

/** 
 * Output a log message. 
 * Message type and severity parameters are required, followed by variable number of 
 * additional parameters of varying type according to the selected message type.
 * The parameters required are indicated by the type enumeration suffix.
 */
int message(MSGTYPE type, MSGSEV sev, ...) {
    va_list args;
    int ret = 0;

    /* ignore if this level not selected */
    if (sev <= Msg_Level) {
        /* the severity first */
        ret = fprintf(stderr, "%s: ", MSG_SEV_TEXT[sev]);

        /* the rest with the variable args */
        va_start(args, sev);
        ret += vfprintf(stderr, MSG_TEXT[type], args);
        va_end(args);

        /* write immediately */
        fflush(stderr);
    }

    /* returns number of characters printed */
    return ret;
}

/**
 * Set the message level. Text must match one of the severity text list. Ignores case.
 * Returns true if match found.
 */
bool set_message_level(const char *levelstr) {

    /* match to one of the possible options */
    int matchidx = match_string(levelstr, MSG_SEV_TEXT, MSG_NUM);
    if (matchidx >= 0) {
        Msg_Level = matchidx;
        return true;
    }
    else {
        return false;
    }
}

/** Set the message file name and location. */
void set_message_path(const CSTRING path) {

    /* validate later */
    Msg_Path = copy_CSTRING(path);

}

/**
 * Start up; call at program start after options.
 * Redirects stderr to log file if requested and outputs a log file header.
 * Log file name also used by program associated files so if none then
 * generates a randomised virtual name of form ayb_xxxxxx_yymmdd_hhmm.
 * Returns true unless requested log file cannot be opened or virtual name cannot be created.
 */
bool startup_message(void) {

    const char *DATE_TIME = "%d %B %Y %H:%M";           // log header date; gives 'dd mmmm yyyy hh:mm'
    const size_t DATE_TIME_LEN = 24;                    // maximum length, including null terminator

    /* get the current date and time */
    time_t lt = time(NULL);
    struct tm *p_tm = localtime(&lt);
    char timestring[DATE_TIME_LEN];

    if (Msg_Path == NULL) {
        /* create a virtual name */
        Msg_Path = generate_name(p_tm);
        if (Msg_Path == NULL) {
            /* error message already issued */
            return false;
        }
        else {
            fprintf(stdout, "AYB virtual log name is %s\n", Msg_Path);
        }
    }

    else {
        /* error file specified, redirect stderr from within program */
        /* validate any path supplied */
        if (!check_path(Msg_Path)) {
            return false;
        }

        /* check file can be created */
        FILE *fp;
        if ((fp = fopen(Msg_Path, "w")) == NULL) {
            message(E_OPEN_FAIL_SS, MSG_FATAL, "Message", Msg_Path);
            return false;
        }

        fclose(fp);
        freopen(Msg_Path, "w", stderr);
        fprintf(stdout, "AYB message log is %s\n", Msg_Path);
    }

    /* create a time string for the log file header */
    strftime(timestring, DATE_TIME_LEN, DATE_TIME, p_tm);

    /* create a log file header */
    fprintf(stderr, "AYB Message Log;\tCreated by %s;\t%s\n\n", getenv("USER"), timestring);
    message(E_MSG_LEVEL_S, MSG_INFO, MSG_SEV_TEXT[Msg_Level]);
    return true;
}

/** Tidy up; call at program shutdown. */
void tidyup_message (void) {

    /* free string memory */
    free_CSTRING(Msg_Path);
}


#ifdef TEST
#include <err.h>

static const char * SETMSGLEV = "error";
static const char *STRING1 = "xxx1";
static const char *STRING2 = "xxx2";
static const char CHAR1 = 'x';
static const int INT1 = 9;
static const int INT2 = 99;
static const real_t EXP1 = 1e-5;
static const real_t FLOAT1 = 99.99;

int main ( int argc, char * argv[]){
    if(argc<3){
        /* arguments are log path and whether to test anything else in addition to path selection */
        errx(EXIT_FAILURE,"Usage: test-message log_path more_flag");
    }

    bool more = false;
    if ((*argv[2] == 'T') || (*argv[2] == 't')) {
        more = true;
    }

    if (more) {
        fputs("Startup no path:\n", stdout);
        startup_message();
        free_CSTRING(Msg_Path);
    }

    fprintf(stdout, "Startup with path: %s\n", argv[1]);
    set_message_path((const CSTRING)argv[1]);
    startup_message();
    
    if (more) {
        fputs("Generate Randomised Name:\n", stdout);
        const int NUM = 10;
        time_t lt = time(NULL);
        struct tm *p_tm = localtime(&lt);
        CSTRING names[NUM];

        for (int i = 0; i < NUM; i++) {
            names[i] = generate_name(p_tm);
        }
        for (int i = 0; i < NUM; i++) {
            fprintf(stdout, "  Name %d: %s\n", i, names[i]);
            free_CSTRING(names[i]);
        }

        fprintf(stdout, "Set invalid message level, %s:\n", STRING1);
        bool ret = set_message_level(STRING1);
        if (ret==false) {
            xfputs("Return value false, ok\n", xstdout);
        }
        else {
            xfputs("Return value not false, not ok\n", xstdout);
        }
        
        fprintf(stdout, "Set valid message level, %s:\n", SETMSGLEV);
        ret = set_message_level(SETMSGLEV);
        if (ret==true) {
            xfputs("Return value true, ok\n", xstdout);
        }
        else {
            xfputs("Return value not true, not ok\n", xstdout);
        }

        fputs("Get message path and level:\n", stdout);
        fprintf(stdout, "Message path; level: %s; %s\n", get_message_path(), MSG_SEV_TEXT[get_message_level()]);

        fputs("Suppress messages by level:\n", stdout);
        for (MSGSEV sevi = (MSGSEV)0; sevi < MSG_NUM; sevi++) {
            message((MSGTYPE)0, sevi);
        }

        fputs("For message output see message log:\n", stdout);
        /* maximum level ensures output */
        Msg_Level = MSG_NUM - 1;
        MSGSEV sev = (MSGSEV)0;
        for (MSGTYPE typi = 0; typi < E_END_NONE; typi++) {
            message(typi, sev);
        }
        /* vary the severity */
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_NONE + 1; typi < E_END_S; typi++) {
            message(typi, sev, STRING1);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_S + 1; typi < E_END_SS; typi++) {
            message(typi, sev, STRING1, STRING2);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_SS + 1; typi < E_END_SC; typi++) {
            message(typi, sev, STRING1, CHAR1);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_SC + 1; typi < E_END_SD; typi++) {
            message(typi, sev, STRING1, INT1);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_SD + 1; typi < E_END_SDD; typi++) {
            message(typi, sev, STRING1, INT1, INT2);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_SDD + 1; typi < E_END_SG; typi++) {
            /* test a variety of numbers with g specifier */
            message(typi, sev, STRING1, EXP1);
            message(typi, sev, STRING1, FLOAT1);
            message(typi, sev, STRING1, (real_t) INT1);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_SG + 1; typi < E_END_C; typi++) {
            message(typi, sev, CHAR1);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_C + 1; typi < E_END_D; typi++) {
            message(typi, sev, INT1);
        }
        sev = (sev + 1) % MSG_NUM;
        for (MSGTYPE typi = E_END_D + 1; typi < E_END_DD; typi++) {
            message(typi, sev, INT1, INT2);
        }
    }
    
    tidyup_message();
    return EXIT_SUCCESS;
}

#endif /* TEST */
