/**
 * \file message.c
 * General Messaging Utility.
 * This provides a central messaging system for the output of program information at various levels up to debugging.
 * The level of message output is selected using a program option (default Warning).
 * \n\n Program messages are output to stderr which can be redirected in a run script.
 * Alternatively a location can be specified as a program option and a message file with a generated unique name will be created.
 * \n\n A program message is implemented by calling the message() function with a type and severity parameter,
 * followed by a variable number of additional parameters of varying type corresponding to the selected message type.
 * The type enumerations include a suffix to indicate the parameters required.
 * The varying parameters are handled using vfprintf which takes a va_list argument. 
 * \n\n New message types can be added as needed by extending the MsgTypeT enumeration (#MSGTYPE) and adding the appropriate message text. 
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
        "Number of cycles to analyse must be supplied as a positive integer\n", // E_NOCYCLES
        "Blockstring option not supplied or contains no datablocks\n",          // E_NOBLOCKS
        "No file pattern match supplied\n",                                     // E_NOPATTERN
        "Number of model iterations incorrectly supplied\n",                    // E_BADITER
        "Memory allocation failed during %s\n",                                 // E_NOMEM_S
        "Log message output level: %s\n",                                       // E_MSG_LEVEL_S
        "Input from directory: %s\n",                                           // E_INPUT_DIR_S
        "Output to directory: %s\n",                                            // E_OUTPUT_DIR_S
        "Input Format selected: %s\n",                                          // E_INPUT_FORM_S
        "Output Format selected: %s\n",                                         // E_OUTPUT_FORM_S
        "Input file found: %s\n",                                               // E_INPUT_FOUND_S
        "Failed to read input file: %s\n",                                      // E_BAD_INPUT_S
        "Failed to create data blocks for input file: %s\n",                    // E_DATABLOCK_FAIL_S
        "Failed to initialise %s matrix\n",                                     // E_MATRIX_FAIL_S
        "Zero lambdas per iteration: %s\n",                                     // E_ZERO_LAMBDA_S
        "Supplied %s location parameter \'%s\' is not a directory\n",           // E_BAD_DIR_SS
        "Failed to create new %s directory \'%s\'\n",                           // E_NOCREATE_DIR_SS
        "Created new %s directory: %s\n",                                       // E_CREATED_DIR_SS
        "No input files in directory \'%s\' matching pattern: \'%s\'\n",        // E_NOINPUT_SS
        "%s file failed to open: %s\n",                                         // E_OPEN_FAIL_SS
        "Input file pattern match: \'%s...\'; %d files found\n",                // E_PATTERN_MATCH_SD
        "Number of %s selected: %d\n",                                          // E_OPT_SELECT_SD
        "%s matrix wrong size, need dimension %d not %d\n",                     // E_MATRIXINIT_SDD
        "Unrecognised nucleotide \'%c\'; returning NUC_AMBIG\n",                // E_BAD_NUC_C
        "Processing failed at iteration %d; calls set to null\n",               // E_PROCESS_FAIL_D
        "Intensity file contains less data than requested; %d instead of %d\n", // E_CYCLESIZE_DD
        "Tile data size: %d clusters of %d cycles\n",                           // E_TILESIZE_DD
        "Failed to initialise model for block %d, %d cycles\n",                 // E_INIT_FAIL_DD
        "Processing block %d, %d cycles\n",                                     // E_PROCESS_DD

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

/** Message severity text. Used to match program argument and as text in log file.
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

/* location and name of message file */
//static const char *DEFAULT_PATH = "./";            ///< Use current directory if none supplied.
static const char PATH_DELIM = '/';                 ///< Path delimiter. For linux, other OS?
static const char *NAME_EXT = ".log";               ///< Extension for log file.
static const char *TIME_SUFFIX = "_%y%m%d_%H%M";    ///< Log file suffix; gives _yymmdd_hhmm.
static const size_t TIME_SUFFIX_LEN = 13;           ///< Length of log file suffix, including null terminator.
static const char *DATE_TIME = "%d %B %Y %H:%M";    ///< Log header date; gives 'dd mmmm yyyy hh:mm'.
static const size_t DATE_TIME_LEN = 24;             ///< Maximum length of log header date, including null terminator.
#define FILENAME_LEN 80

/* members */

/** Selected level of messages, default Warning. */
static int Msg_Level = MSG_WARN;
/** Selected location for message file. 
 * Output is to stderr unless redirected by a program option, 
 * which defaults to the current directory if the supplied path does not exist.
 */
static char Msg_Path[FILENAME_LEN] = "";


/* private functions */

/** 
 * Generate a unique message file name including date and time and add to path.
 * Filename is <prefix>_<yymmdd>_<hhmm>.log.
 * Returns true unless problem with log directory.
 */
static bool create_filename(const char *prefix, char *timestring, char *name) {

    /* check specified path exists */
    if (!check_outdir(Msg_Path, "message")) {
        return false;
    }

    strcpy(name, Msg_Path);
    size_t len = strlen(name);

    /* add path delimiter if not supplied */
    if (*(name + len - 1) != PATH_DELIM) {
        *(name + len) = PATH_DELIM;
        *(name + len + 1) = '\0';
    }

    /* add the default filename prefix */
    strcat(name, prefix);

    /* add the date and time to make it unique given it will take some time to run */
    strcat(name, timestring);

    /* add the extension */
    strcat(name, NAME_EXT);
    printf("name: %s\n", name);
    return true;
}


/* public functions */

/** Output a log message. */
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

/** Set the message file location. */
void set_message_path(const char *path) {

    /* validate later */
    strcpy(Msg_Path, path);
}

/**
 * Start up; call at program start after options.
 * Redirects stderr to log file if requested and outputs a log file header.
 * Returns true unless requested log file cannot be opened.
 */
bool startup_message(const char *prefix) {

    /* get the current date and time */
    time_t lt = time(NULL);
    struct tm *p_tm = localtime(&lt);
    char timestring[DATE_TIME_LEN];

    if (strlen(Msg_Path) > 0) {
        /* error file specified, redirect stderr from within program */
        char filename[FILENAME_LEN];
        /* create a unique name including current date/time and add path */
        strftime(timestring, TIME_SUFFIX_LEN, TIME_SUFFIX, p_tm);
        if (!create_filename(prefix, timestring, (char*)&filename)) {
            return false;
        }

        /* check file can be created */
        FILE *fp;
        if ((fp = fopen(filename, "w")) == NULL) {
            message(E_OPEN_FAIL_SS, MSG_FATAL, "Message", filename);
            return false;
        }

        fclose(fp);
        freopen(filename, "w", stderr);
    }

    /* create a time string for the log file header */
    strftime(timestring, DATE_TIME_LEN, DATE_TIME, p_tm);

    /* create a log file header */
    fprintf(stderr, "AYB Message Log;\tCreated by %s;\t%s\n\n", getenv("USER"), timestring);
    message(E_MSG_LEVEL_S, MSG_INFO, MSG_SEV_TEXT[Msg_Level]);
    return true;
}

/** Tidy up; call at program shutdown. */
void tidyup_message () {

    /* close the message file */
    fclose(stderr);
}

/* delete all message files; run on supplied parameter */
void tidy_message() {
    //needed??
}
