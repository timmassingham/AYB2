/**
 * \file dirio.c
 * I/O Environment.
 *   - Directory search.
 *   - Input file pattern match.
 *   - File open/close
 *//*
 *  Created : 16 Mar 2010
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

//#define __USE_GNU
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>          // for toupper
#include <dirent.h>
#include "dirio.h"
#include "message.h"


/* constants */

static const char *DEFAULT_PATH = "./";         ///< Use current directory if none supplied.
static const char *PATH_DELIMSTR = "/";         ///< Path delimiter. For linux, other OS?
static const char DOT = '.';                    ///< Before file extension, used for tag location.
static const char DELIM = '_';                  ///< Before tag, used for tag location.
static const char *INT_TAG = "int";             ///< Intensities tag. Should this be an argument hmhm?
static const int INT_TAG_LEN = 3;               ///< Length of Intensities tag.

/* members */

/* the I/O settings for this run */
static CSTRING Input_Path = NULL;               ///< Input path, default or program argument.
static CSTRING Output_Path = NULL;              ///< Output path, default or program argument.
static CSTRING Pattern = NULL;                  ///< File pattern match, currently a prefix, program argument.
static size_t Pattern_Len = 0;                  ///< Length of pattern match, program exits if not > 0.


static struct dirent **Dir_List = NULL;         ///< The pattern matched directory list.
static int Dir_Num = 0;                         ///< Number of pattern matched files found.
static int Index = -1;                          ///< Current index to pattern matched files.

static CSTRING Current;                         ///< Current input file, also used to create output filename.


/* private functions */

/**
 * Create a full path name from a directory and filename. Places a path delimiter between them;
 * delimiter may change according to target system.
 */
static bool full_path(const CSTRING dir, const CSTRING filename, CSTRING *filepath) {
/*  Parameters: dir - the directory, may be null
                filename - the filename
                filepath - full pathname (returned)
    Returns:    whether succeeded
*/
    /* determine length required, allow for delimiter */
    size_t len = strlen(filename);
    if (dir != NULL) {
        len += strlen(dir) + strlen(PATH_DELIMSTR);
    }

    *filepath = new_CSTRING(len);
    if (*filepath == NULL) {
        message(E_NOMEM_S, MSG_FATAL, " file path creation");
        return false;
    }
    else {
        if (dir != NULL) {
            strcpy(*filepath, dir);
            strcat(*filepath, PATH_DELIMSTR);
        }
        strcat(*filepath, filename);
        message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, *filepath);
        return true;
    }
}

//#ifndef __USE_GNU
/* substring search ignoring case, needed hmhm? */
char *strcasestr(const char *s1, const char *s2) {

    CSTRING lookin = copy_CSTRING((CSTRING)s1);
    CSTRING lookfor = copy_CSTRING((CSTRING)s2);
    char *ret;
    int pos;

    for (pos = 0; pos < strlen(lookin); pos++) {
        lookin[pos] = toupper(lookin[pos]);
    }
    for (pos = 0; pos < strlen(lookfor); pos++) {
        lookfor[pos] = toupper(lookfor[pos]);
    }

    ret = strstr(lookin, lookfor);
    if (ret != NULL) {
        /* move pointer to original string */
        ret = (char*)s1 + (int)(ret - lookin);
    }
    free_CSTRING(lookin);
    free_CSTRING(lookfor);

    return ret;
}
//#endif


/** Selector function for scandir. Matches to a prefix and includes a fixed tag. */
static int match_pattern(const struct dirent *list) {

    int ret = 0;

    /* check if prefix matches */
    int diff = strncasecmp(list->d_name, Pattern, Pattern_Len);
    if (diff == 0) {

        if (INT_TAG_LEN > 0) {
            /* check if contains input tag, surrounded by delim/dot */
            CSTRING compare_tag = new_CSTRING(INT_TAG_LEN + 2);
            const char *pnext;
            int pos = 0;
            compare_tag[pos++] = DELIM;
            for (pnext = INT_TAG; pnext < INT_TAG + INT_TAG_LEN; pnext++) {
                compare_tag[pos++] = *pnext;
            }
            compare_tag[pos++] = DOT;
            compare_tag[pos] = '\0';

            char *it;
            it = strcasestr(list->d_name, compare_tag);
            if(it != NULL) {
                ret = 1;
            }
            free_CSTRING(compare_tag);
        }
        else {
            /* no tag to match */
            ret = 1;
        }
    }

    return ret;
}

/** Create a new file name from an original. Replaces the part beyond the last delimiter with a new tag. */
CSTRING output_name(const CSTRING oldname, const CSTRING tag) {

    CSTRING newname = NULL;

    /* find last delimiter and first dot */
    char *pdlm = strrchr(oldname, DELIM);
    char *pdot = strchr(oldname, DOT);

    size_t taglen = strlen(tag);
    size_t oldlen = strlen(oldname);

    /* get a string of size for new name */
    newname = new_CSTRING(oldlen + taglen - (pdot - pdlm - 1));

    if (newname == NULL) {
        message(E_NOMEM_S, MSG_FATAL, " output name creation");
    }
    else {
        /* copy the component parts to the new name */
        char *pnext;
        int pos = 0;
        for (pnext = oldname; pnext <= pdlm; pnext++) {
            newname[pos++] = *pnext;
        }
        for (pnext = tag; pnext < tag + taglen; pnext++) {
            newname[pos++] = *pnext;
        }
        for (pnext = pdot; pnext < oldname + oldlen; pnext++) {
            newname[pos++] = *pnext;
        }
        /* finish with string terminator */
        newname[pos] = '\0';
    }
    message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, newname);
    return newname;
}


/** Scan the input directory for any files that match the specified pattern. */
static int scan_inputs() {

    int num = scandir (Input_Path, &Dir_List, match_pattern, alphasort);

    if (num < 0) {
        /* hmhm new message func to return a message string */
        char msg[80];
        sprintf(msg, "Couldn't open the directory: \'%s\'", Input_Path);
        perror (msg);
    }
    else {
        for (int cnt = 0; cnt < num; ++cnt) {
//            message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, Dir_List[cnt]->d_name);
        }
    }

    return num;
}


/* public functions */

/** Open the next input file in the directory. Return the file handle. Also closes previous file if any. */
XFILE * open_next(XFILE *fplast) {

    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    if (fplast != NULL) {
        xfclose(fplast);
        free_CSTRING(Current);
    }

    if (++Index < Dir_Num) {
        Current = copy_CSTRING(Dir_List[Index]->d_name);
        message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, Current);

        /* make a full path name */
        if (full_path(Input_Path, Current, &filepath)) {
           fp =  xfopen(filepath, XFILE_UNKNOWN, "r" );

           if (xisnull_file(fp)) {
               xfclose(fp);
//               fp = NULL;
           }
           else {
               message(E_INPUTFOUND_S, MSG_INFO, Current );
           }
        }
        free_CSTRING(filepath);
    }

    return fp;
}

/** Open an output file corresponding to current input file with supplied tag. Return the file handle. */
XFILE * open_output(const CSTRING tag) {

    CSTRING filename = NULL;
    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    /* create output file name */
    filename = output_name(Current, tag);

    if (filename != NULL) {
        if (full_path(Output_Path, filename, &filepath)) {
           fp =  xfopen(filepath, XFILE_UNKNOWN, "w" );

           if (xisnull_file(fp)) {
               xfclose(fp);
           }
           else {
               message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, filename);
           }
        }
    }
    free_CSTRING(filename);
    free_CSTRING(filepath);

    return fp;
}

/** Set the input/output path. */
void set_path(const CSTRING path, IOMODE mode){

    switch(mode) {
        case E_INPUT:
            Input_Path = copy_CSTRING(path);
            break;
        case E_OUTPUT:
            Output_Path = copy_CSTRING(path);
            break;
        default: ;
    }
}

/** Set the input filename pattern to match to. */
void set_pattern(const CSTRING pattern) {

    Pattern = copy_CSTRING(pattern);
}

/** Start up; call at program start after options. */
bool startup_dirio() {

    /* check pattern match supplied */
    if (Pattern != NULL) {
        Pattern_Len = strlen(Pattern);
    }

    if (Pattern_Len == 0) {
        message(E_NOPATTERN, MSG_FATAL);
        return false;
    }
    else {
        message(E_PATTERN_SELECT_S, MSG_INFO, Pattern);

        /* set default i/o paths if not supplied */
        if (Input_Path == NULL) {
            Input_Path = copy_CSTRING((CSTRING)DEFAULT_PATH);
        }
        if (Output_Path == NULL) {
            Output_Path = copy_CSTRING((CSTRING)DEFAULT_PATH);
        }

        /* scan for matching input files */
        Dir_Num = scan_inputs();
        if (Dir_Num <= 0) {
            message (E_NOINPUT_S, MSG_FATAL, Pattern);
            return false;
        }
        else {
            /* at least one input file */
            return true;
        }
    }
}

/** Tidy up; call at program shutdown. */
void tidyup_dirio() {

    /* free string memory */
    free_CSTRING(Input_Path);
    free_CSTRING(Output_Path);
    free_CSTRING(Pattern);
}
