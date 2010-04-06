/*
 *  File    : message.c
 *  Created : 26 Feb 2010
 *  Author  : Hazel Marsden
 *  Purpose : General Messaging Utility
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
#include <time.h>
#include <sys/stat.h>
#include "message.h"


/*
 * formatted text for each message type, to be used directly in a printf;
 * first parameter is severity; suffix indicates type of arguments required;
 * match list to MsgTypeT enum
 */
const char *MSG_TEXT[] = {
        "%s: An integer (%d)\n",                //E_INT_D
        "%s: A float (%0.2f)\n",                //E_FLOAT_F
        "%s: A string (%s)\n",                  //E_STRING_S
        "%s: Give 2 integers: %d %d\n" };       //E_INT_DD

/* message severity text; match list to MsgSeverityT enum */
const char *MSG_SEV_TEXT[] = {
        "None",
        "Fatal",
        "Error",
        "Warning",
        "Information",
        ""};

/* selected level of messages */
int Msg_Level = MSG_WARN;

/* location and name of message file */
static const char *DEFAULT_PATH = "./";
static const char PATH_DELIM = '/';
static const char *NAME_EXT = ".log";
static const char *TIME_SUFFIX = "%y%m%d_%H%M";
static const size_t TIME_SUFFIX_LEN = 12;
#define FILENAME_LEN 80
char Msg_Path[FILENAME_LEN] = "";


/* generate a unique message file name including date and time */
void create_filename(const char* prefix, char *name) {

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

/* match a string to one of a list */
/* return index of match or -1 if none match */
int match_string(const char *string, const char *match[], int num) {

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

/* set the message level; must match one of the options */
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

/* set the message path */
void set_message_path(const char *path) {

    /* validate later */
    strcpy(Msg_Path, path);
}

/* start a message file; call at program start */
void start_message(const char *prefix) {

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

/* tidy message file; call at program shutdown */
void tidyup_message () {

    /* close the message file */
    fclose(stderr);
}

/* delete all message files; run on supplied parameter */
void tidy_message() {
    //needed??
}
