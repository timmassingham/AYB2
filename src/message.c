/**
 * \file message.c
 * General Messaging Utility.
 * Use function message (<message type>, <message severity>, <parameters...>) to output a message.
 * Add new message types as required.
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
#include <sys/stat.h>
#include "message.h"


/* constants */

/**
 * Formatted text for each message type. Used directly in a printf.
 * Ensure list matches MsgTypeT enum.
 */
static const char *MSG_TEXT[] = {
        "Number of cycles to analyse must be supplied as a positive integer\n", // E_NOCYCLES
        "No file pattern match supplied\n",                                     // E_NOPATTERN
        "File pattern selected: %s...\n",                                       // E_PATTERN_SELECT_S
        "Memory allocation failed during %s\n",                                 // E_NOMEM_S
        "No input files located matching pattern: \'%s\'\n",                    // E_NOINPUT_S
        "Input file found: %s\n",                                               // E_INPUT_FOUND_S
        "Failed to read input file: %s\n",                                      // E_BAD_INPUT_S
        "Failed to initialise model for input file: %s\n",                      // E_INIT_FAIL_S
        "%s directory \'%s\' not found\n",                                      // E_NODIR_SS
        "%s file failed to open: %s\n",                                         // E_OPEN_FAIL_SS
        "Unrecognised nucleotide \'%c\'; returning NUC_AMBIG\n",                // E_BAD_NUC_C
        "Number of cycles selected: %d\n",                                      // E_CYCLE_SELECT_D
        "Intensity file contains less data than requested; "
            "number of cycles changed from %d to %d\n",                         // E_CYCLESIZE_DD
        "Matrix size incorrectly specified: read in as %d by %d\n",             // E_BAD_MATSIZE_DD
        "Insufficient data or incorrect file format; "
             "expected %d %s but found only %d\n",                              // E_READ_ERR_DSD

        "%s %20s\n",                                                            // E_GENERIC_SS
        "%s %d\n",                                                              // E_GENERIC_SD
        "%s %u\n",                                                              // E_GENERIC_SU
        "%s %f\n",                                                              // E_GENERIC_SF
        "%s %0.2f %0.2f %0.2f %0.2f\n",                                         // E_GENERIC_SFFFF

        "%s (%s:%d): %s\n",                                                     // E_DEBUG_SSD_S
        "%s (%s:%d): %s %s\n",                                                  // E_DEBUG_SSD_SS
        "%s (%s:%d): %s %d\n"                                                   // E_DEBUG_SSD_SD
        };

/** Message severity text. Ensure list matches to MsgSeverityT enum. */
static const char *MSG_SEV_TEXT[] = {
        "None",
        "Fatal",
        "Error",
        "Warning",
        "Information",
        "Debug",
        ""};

/* location and name of message file */
static const char *DEFAULT_PATH = "./";         ///< Use current directory if none supplied.
static const char PATH_DELIM = '/';             ///< Path delimiter. For linux, other OS?
static const char *NAME_EXT = ".log";           ///< Extension for log file.
static const char *TIME_SUFFIX = "%y%m%d_%H%M"; ///< Log file suffix; gives yymmdd_hhmm.
static const size_t TIME_SUFFIX_LEN = 12;       ///< Length of log file suffix.
#define FILENAME_LEN 80

/* members */

static int Msg_Level = MSG_WARN;                ///< Selected level of messages.
static char Msg_Path[FILENAME_LEN] = "";        ///< Selected path for message file.


/* private functions */

/** Generate a unique message file name including date and time. */
static void create_filename(const char* prefix, char *name) {

    /* check specified path exists */
    struct stat st;

    if(stat(Msg_Path, &st) == 0) {
        strcpy(name, Msg_Path);
        size_t len = strlen(name);

        /* add path delimiter if not supplied */
        if (*(name + len - 1) != PATH_DELIM) {
            *(name + len) = PATH_DELIM;
            *(name + len + 1) = '\0';
        }
    }
    else {
        fprintf(stderr, "Message file directory: \'%s\' does not exist; Using default location\n\n", Msg_Path);
        /* use default path */
        strcpy(name, DEFAULT_PATH);
    }

    /* add the default filename prefix */
    strcat(name, prefix);

    /* add the date and time to make it unique given it will take some time to run */
    char timestring[TIME_SUFFIX_LEN + 1];

    time_t lt = time(NULL);
    struct tm *p_tm = localtime(&lt);
    strftime(timestring, TIME_SUFFIX_LEN, TIME_SUFFIX, p_tm);
    strcat(name, timestring);

    /* add the extension */
    strcat(name, NAME_EXT);
    printf("name: %s\n", name);
}

/** Match a string to one of a list. Returns index of match or -1 if none. */
static int match_string(const char *string, const char *match[], int num) {

    int result = -1;

    for (int idx = 0; idx < num; idx++) {
        if (strcasecmp(string, match[idx]) == 0) {
            result = idx;
            break;
        }
    }
    return result;
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

/** Set the message level. Text must match one of the severity text list. Ignores case. */
bool set_message_level(const char *levelstr) {

    /* match to one of the possible options */
    int found = match_string(levelstr, MSG_SEV_TEXT, MSG_NUM);
    if (found >= 0) {
        Msg_Level = found;
        return true;
    }
    else {
        return false;
    }
}

/** Set the message file path. */
void set_message_path(const char *path) {

    /* validate later */
    strcpy(Msg_Path, path);
}

/** Start up; call at program start after options. */
void startup_message(const char *prefix) {

    if (strlen(Msg_Path) > 0) {
        /* error file specified, redirect stderr from within program*/
        /* create a unique name */
        char filename[FILENAME_LEN];
        create_filename(prefix, (char*)&filename);

        /* check file can be created */
        FILE *fp;
        if ((fp = fopen(filename, "w")) == NULL) {
            fprintf(stderr, "Cannot open message file: \'%s\'; Messages to standard error and level set to Fatal\n\n", filename);
            Msg_Level = MSG_FATAL;        }
        else {
            fclose(fp);
            freopen(filename, "w", stderr);
        }
    }
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
