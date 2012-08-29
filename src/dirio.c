/**
 * \file dirio.c
 * I/O Environment.
 * Handles directory search and file open and close and file errors.
 * 
 * The input and output file locations are specified as program options (default current directory).
 * Currently allows (by program option) for input files of type cif or txt (standard illumina format).
 * Input intensity files are selected by matching with a prefix (non-option program argument) and substring (fixed).
 * Standard function scandir is used to search the input directory for files matching the required pattern.
 * On each request for the next input any previous file is closed, the next one found is opened and the file handle returned.
 *
 * If cif input from a run-folder is selected (by program option) then the prefix is replaced by
 * a lane and tile range string of the form Ln[-n]Tn[-n].
 * Each request for the next input returns the next lane and tile in the range.
 * A virtual intensities filename s_L_TTTT is created.
 *
 * Other files such as predetermined input matrix and spike-in data files are also opened here.
 * 
 * Output file names are generated from the intensities file name by replacing the 'tag' with a new one.
 * The tag is the file suffix for cif files or between the last delimiter and the first dot for txt files. 
 * Any txt compression suffix is also removed as output is always uncompressed.
 * An attempt is made to create output directories if they do not exist.
 * 
 * Used as a singleton class with local functions accessing global member data.
 *//*
 *  Created : 16 Mar 2010
 *  Author  : Hazel Marsden
 *  Author  : Tim Massingham
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
 *  Copyright (C) 2012 Disinformatics.org
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
#include <strings.h>
#include <ctype.h>          // for toupper
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <err.h>
#include "dirio.h"
#include "message.h"


/* constants */

static const char *DEFAULT_PATH = "./";         ///< Use current directory if none supplied.
static const char *PATH_DELIMSTR = "/";         ///< Path delimiter as string. For linux, other OS?
static const char PREFIXCHAR = '+';             ///< Indicates pattern to be treated as a prefix.
static const char RANGECHAR = '-';              ///< Used to separate the min/max in a range of values.
static const char DOT = '.';                    ///< Before file extension, used for tag location.
static const char DELIM = '_';                  ///< Before tag, used for tag location.
static const char BLOCKCHAR = 'a';              ///< Start for additional block suffix.

/**
 * Possible input format text. Match to INFORM enum. Used to match program argument.
 * - Illumina files match name template {prefix}[*]_int.txt*[.{zip ext}]
 * - cif files match name template {prefix}[*].cif
 */
static const char *INFORM_TEXT[] = {"TXT", "CIF"};
/** Text for input format messages. Match to INFORM enum. */
static const char *INFORM_MESS_TEXT[] = {"standard illumina txt", "cif"};
static const char *INTEN_TAG[] = {"int", ""};           ///< Fixed Intensities file tags.
static const char *INTEN_SUF[] = {"txt", "cif"};        ///< Fixed Intensities file suffixes.
static const char *SPIKEIN_TAG = "spike";               ///< Spike-in file tag.

static const char *LTMESS_TEXT = "Lane tile string";    ///< Lane Tile parameter name for messages.

/**
 * Permission flags for a directory; owner/group all plus other read/execute.
 * Used by output directory create but seems to be ignored in favour of parent adoption.
 */
static const mode_t DIR_MODE = S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH;

/** Directory check possible outcomes. */
typedef enum DirOptT {E_ISDIR, E_NODIR, E_NOEXIST} DIROPT;

/* members */

/* the I/O settings for this run */
static INFORM Input_Format = E_CIF;             ///< Selected input format.
static CSTRING Input_Path = NULL;               ///< Input path, default or program argument.
static CSTRING Output_Path = NULL;              ///< Output path, default or program argument.
static CSTRING Spikearg_Path = NULL;            ///< Spike-in path, program argument.
static CSTRING IntenSubstr = NULL;              ///< Substring an intensities file must contain.
static CSTRING Sample_Name = "Sample";		///< Sample for these tiles.

/**
 * File pattern match. Currently a prefix or the whole part of the filename before the substring.
 * Supplied as non-option program argument. Any path parts moved to pattern path.
 */
static CSTRING Pattern = NULL;
static CSTRING Pattern_Path = NULL;             ///< Input path and any pattern path parts.
static size_t Pattern_Len = 0;                  ///< Length of pattern match, program exits if not > 0.
static CSTRING Spikein_Path;                    ///< Spike-in path to use.
static CSTRING Matrix[E_NMATRIX];               ///< Predetermined matrix input file locations.

static struct dirent **Dir_List = NULL;         ///< The pattern matched directory list.
static int Dir_Num = 0;                         ///< Number of pattern matched files found.
static int Index = -1;                          ///< Current index to pattern matched files.

static bool RunFolder = false;                  ///< Set to read intensities from a run-folder.
static LANETILE LTMin = {0, 0};                 ///< Selected minimum run-folder lane and tile.
static LANETILE LTMax = {0, 0};                 ///< Selected maximum run-folder lane and tile.
static LANETILE LTCurrent = {0, 0};             ///< Current run-folder lane and tile.

static CSTRING Current = NULL;                  ///< Current intensities file, used to create other filenames.
static CSTRING Last_Input = NULL;               ///< Last non-intensities input file opened.
static CSTRING Last_Spikein = NULL;             ///< Last spike-in file opened.


/* private functions */

/** Return whether specified directory exists. */
static DIROPT check_dir(const CSTRING dirname) {

    struct stat st;
    if (stat(dirname, &st) == 0) {
        /* check it is a directory; ISDIR returns non-zero for a directory */
        if (S_ISDIR (st.st_mode) != 0) {
            return E_ISDIR;
        }
        else {
            return E_NODIR;
        }
    }
    else {
        return E_NOEXIST;
    }
}

/** Clear structures and values associated with a pattern. */
static void clear_pattern(void) {

    Pattern = free_CSTRING(Pattern);
    Pattern_Path = free_CSTRING(Pattern_Path);
    Pattern_Len = 0;
    Index = -1;

    if ((Dir_Num > 0) && (Dir_List != NULL)) {
        for (int i = 0; i < Dir_Num; i++) {
            free(Dir_List[i]);
        }
        free(Dir_List);
        Dir_List = NULL;
        Dir_Num = 0;
    }
}

/** Return a string containing a single value or a range. */
static CSTRING format_range(unsigned int min, unsigned int max) {

    /* sufficient for MAXINT */
    static const size_t LEN = 21;
    CSTRING str = new_CSTRING(LEN);

    if (min == max) {
        sprintf(str, "%u", min);
    }
    else {
        sprintf(str, "%u%c%u", min, RANGECHAR, max);
    }

    return str;
}

/**
 * Create a full path name from a directory and filename.
 * If filename contains full path then use alone,
 * otherwise combine, placing a path delimiter between them;
 * delimiter may change according to target system.
 * Return false if supplied filename is null or cannot allocate memory for name string.
 */
static bool full_path(const CSTRING dir, const CSTRING filename, CSTRING *filepath) {
/*  Parameters: dir - the directory, may be null
                filename - the filename
                filepath - full pathname (returned)
    Returns:    whether succeeded
*/
    if (filename == NULL) {return false;}

    if (filename[0] == PATH_DELIM) {
        /* full path given in filename, return unchanged ignoring any dir */
        *filepath = copy_CSTRING(filename);
        if (*filepath == NULL) {
            message(E_NOMEM_S, MSG_FATAL, "file path creation");
            return false;
        }
        else {
            return true;
        }
    }

    /* determine length required, allow for delimiter */
    size_t len = strlen(filename);
    if (dir != NULL) {
        len += strlen(dir) + strlen(PATH_DELIMSTR);
    }

    *filepath = new_CSTRING(len);
    if (*filepath == NULL) {
        message(E_NOMEM_S, MSG_FATAL, "file path creation");
        return false;
    }
    else {
        if (dir != NULL) {
            strcpy(*filepath, dir);
            strcat(*filepath, PATH_DELIMSTR);
        }
        strcat(*filepath, filename);
//        message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, *filepath);
        return true;
    }
}

/**
 * Extract a number from supplied string and store as lane or tile.
 * Returns pointer to next element in string.
 * Returns false if an error.
 */
static bool get_lanetile(const char token, const char *restrict ptr, char **restrict endptr, LANETILE *lanetile) {

    bool ok = true;

    int num = strtol(ptr, endptr, 0);
    if (num <= 0) {
        message(E_BAD_NUM_SS, MSG_ERR, LTMESS_TEXT, ptr);
        ok = false;
    }
    else {
        switch (token) {
            case 'L':
            case 'l':
                lanetile->lane = num;
                break;

            case 'T':
            case 't':
                lanetile->tile = num;
                break;

            default :
                message(E_BAD_CHAR_SC, MSG_ERR, LTMESS_TEXT, token);
                ok = false;
        }
    }

    return ok;
}

/**
 * Extract lane or tile range from lane and tile string.
 * Expected string is Xn[-n].
 * Returns pointer to next element in string.
 * Returns false if an error.
 */
static bool get_lanetile_range(const char *restrict ptr, char **restrict endptr) {

    /* note if lane or tile */
    char token = ptr[0];

    /* starts with minimum from second character onwards */
    bool ok = get_lanetile(token, ++ptr, endptr, &LTMin);

    if (ok) {
        /* check if a range */
        if (*endptr[0] == RANGECHAR) {
            /* find first non-range char */
            while (*endptr[0] == RANGECHAR) {
                (*endptr)++;
            }
            /* finish with maximum */
            ok = get_lanetile(token, *endptr, endptr, &LTMax);
        }
        else {
            /* single number */
            switch (token) {
                case 'L':
                case 'l':
                    LTMax.lane = LTMin.lane;
                    break;

                case 'T':
                case 't':
                    LTMax.tile = LTMin.tile;
                    break;

                default : ;
            }
        }
    }

    return ok;
}

/** Create the fixed tag and suffix string that an intensities file must contain. */
static void make_substring(void) {

    /* check if anything to match */
    size_t taglen = strlen(INTEN_TAG[Input_Format]);
    size_t suflen = strlen(INTEN_SUF[Input_Format]);
    size_t len = 0;

    if (taglen > 0) {
        len += taglen + 1;
    }
    if (suflen > 0) {
        len += suflen + 1;
    }
    if (len > 0) {
        /* add room for terminator */
        IntenSubstr = new_CSTRING(++len);
        const char *pnext;
        int pos = 0;

        if (taglen > 0) {
            /* add delimiter and tag */
            IntenSubstr[pos++] = DELIM;
            for (pnext = INTEN_TAG[Input_Format]; pnext < INTEN_TAG[Input_Format] + taglen; pnext++) {
                IntenSubstr[pos++] = *pnext;
            }
        }
        if (suflen > 0) {
            /* add dot and suffix */
            IntenSubstr[pos++] = DOT;
            for (pnext = INTEN_SUF[Input_Format]; pnext < INTEN_SUF[Input_Format] + suflen; pnext++) {
                IntenSubstr[pos++] = *pnext;
            }
        }
        IntenSubstr[pos] = '\0';
    }
}

//#ifndef __USE_GNU
/* substring search ignoring case, needed hmhm? */
char *strcasestr(const char *s1, const char *s2) {

    if ((s1 == NULL) || (s2 == NULL)) {return NULL;}

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


/**
 * Selector function for scandir. Matches to a prefix and fixed tag and suffix.
 * If pattern ends with prefix indicator then allows inclusion of additional characters between.
 * Assumes match substring already set up.
 * Returns non-zero value if match found as required by scandir.
 */
static int match_pattern(const struct dirent *list) {

    int ret = 0;

    if (Pattern[Pattern_Len - 1] == PREFIXCHAR) {
        /* check if prefix matches beginning */
        if (strncasecmp(list->d_name, Pattern, Pattern_Len - 1) == 0) {
            /* check if contains _input tag and suffix  */
            if (IntenSubstr != NULL) {
                if (strcasestr(list->d_name, IntenSubstr) != NULL) {
                    ret = 1;
                }
            }
            else {
                /* nothing to match */
                ret = 1;
            }
        }
    }

    else {
        /* create whole filename */
        size_t len = Pattern_Len;
        if (IntenSubstr != NULL) {
            len += strlen(IntenSubstr);
        }
        CSTRING filename = new_CSTRING(len);
        strcpy(filename, Pattern);
        if (IntenSubstr != NULL) {
            strcat(filename, IntenSubstr);
        }
        /* check for whole filename match, still allows for compression suffix */
        if (strncasecmp(list->d_name, filename, len) == 0) {
            ret = 1;
        }
        free_CSTRING(filename);
    }

    return ret;
}

/**
 * Move any path parts from filename to file path.
 * If filename contains full path then use alone,
 * otherwise add any partial path from filename to supplied path.
 * In any case remove path parts from filename.
 * Return a new file path.
 */
static CSTRING move_path(const CSTRING filepath, CSTRING *filename) {

    if (filepath == NULL) {return NULL;}
    if (*filename == NULL) {return NULL;}

    CSTRING newpath = NULL;

    /* find the beginning of the name */
    char *pdlm = strrchr(*filename, PATH_DELIM);
    if (pdlm == NULL) {
        /* no path in filename, return supplied path and filename unchanged */
        newpath = copy_CSTRING(filepath);
    }

    else {
        size_t sublen = pdlm - *filename;
        if (*filename[0] == PATH_DELIM) {
            /* full path given in filename, extract and use alone */
            newpath = new_CSTRING(sublen);
            strncpy(newpath, *filename, sublen);
        }
        else {
            /* partial path given in filename, add to supplied path */
            size_t oldlen = strlen(filepath);
            newpath = new_CSTRING(oldlen + sublen + strlen(PATH_DELIMSTR));
            strcpy(newpath, filepath);
            if (newpath[oldlen - 1] != PATH_DELIM) {
                strcat(newpath, PATH_DELIMSTR);
            }
            strncat(newpath, *filename, sublen);
        }

        /* remove the path parts from the filename */
        CSTRING old = *filename;
        *filename = copy_CSTRING(old + sublen + 1);
        free_CSTRING(old);
    }

    return newpath;
}

/** Return just the filename part of a path. */
static CSTRING name_only(const CSTRING filepath) {

    if (filepath == NULL) {return NULL;}
    CSTRING filename = NULL;

    /* find the beginning of the name */
    char *pdlm = strrchr(filepath, PATH_DELIM);
    if (pdlm == NULL) {
        /* no path in filename, return filename unchanged */
        filename = copy_CSTRING(filepath);
    }
    else {
        filename = copy_CSTRING(pdlm + 1);
    }

    return filename;
}

/**
 * Create a new file name by replacing the part between the last delimiter and the first dot with a new tag.
 * Add a block suffix to the name if non-negative blk supplied.
 * Also removes any compression suffix.
 * Used for standard illumina (txt) outputs and additional inputs.
 */
static CSTRING new_name_body(const CSTRING oldname, const CSTRING tag, int blk) {

    if ((oldname == NULL) || (tag == NULL)) {return NULL;}

    CSTRING newname = NULL;
    size_t oldlen = strlen(oldname);
    size_t taglen = strlen(tag);

    /* find last delimiter */
    char *pdlm = strrchr(oldname, DELIM);
    if (pdlm == NULL) {
        /* no delimiter, replace whole name */
        pdlm = oldname - 1;
    }

    /* find first dot */
    char *pdot = strchr(oldname, DOT);
    if (pdot == NULL) {
        /* no dot, replace to end of name */
        pdot = oldname + oldlen;
    }

    /* get a string of size for new name */
    newname = new_CSTRING(oldlen - (pdot - pdlm - 1) + taglen + ((blk >= 0)?1:0));

    if (newname == NULL) {
        message(E_NOMEM_S, MSG_FATAL, " new name creation");
    }
    else {
        /* copy the component parts to the new name */
        char *pnext;
        int pos = 0;
        for (pnext = oldname; pnext < pdlm; pnext++) {
            newname[pos++] = *pnext;
        }
        /* add block suffix before delimiter */
        if (blk >= 0) {
            newname[pos++] = BLOCKCHAR + blk;
        }
        newname[pos++] = *pdlm;
        for (pnext = tag; pnext < tag + taglen; pnext++) {
            newname[pos++] = *pnext;
        }
        for (pnext = pdot; pnext < oldname + oldlen; pnext++) {
            newname[pos++] = *pnext;
        }
        /* finish with string terminator */
        newname[pos] = '\0';

        /* Remove any compression suffix */
        XFILE_MODE mode = guess_mode_from_filename (newname);
        if (mode != XFILE_RAW) {
            char *psuf = strrchr(newname, DOT);
            *psuf = '\0';
        }
    }
//    message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, newname);
    return newname;
}

/**
 * Create a new file name by replacing the suffix with a new tag.
 * Add a block suffix to the body if non-negative blk supplied.
 * Used for cif and program outputs and additional inputs.
 */
static CSTRING new_name_suffix(const CSTRING oldname, const CSTRING tag, int blk) {

    if ((oldname == NULL) || (tag == NULL)) {return NULL;}

    CSTRING newname = NULL;
    size_t oldlen = strlen(oldname);
    size_t taglen = strlen(tag);

    /* find the suffix */
    char *pdot = strrchr(oldname, DOT);
    if (pdot == NULL) {
        /* no dot, add to name */
        pdot = oldname + oldlen;
    }

    /* get a string of size for new name */
    newname = new_CSTRING(pdot - oldname + taglen + 1 + ((blk >= 0)?1:0));

    if (newname == NULL) {
        message(E_NOMEM_S, MSG_FATAL, " new name creation");
    }
    else {
        /* copy the component parts to the new name */
        char *pnext;
        int pos = 0;
        for (pnext = oldname; pnext < pdot; pnext++) {
            newname[pos++] = *pnext;
        }
        /* add block suffix before dot */
        if (blk >= 0) {
            newname[pos++] = BLOCKCHAR + blk;
        }
        newname[pos++] = DOT;
        for (pnext = tag; pnext < tag + taglen; pnext++) {
            newname[pos++] = *pnext;
        }
        /* finish with string terminator */
        newname[pos] = '\0';
    }

//    message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, newname);
    return newname;
}

/**
 * Open an input file corresponding to current intensities file with supplied location and tag.
 * A non-negative blk indicates a block suffix should be added to the name.
 * BLK_SINGLE indicates no block suffix.
 * Return the file handle or NULL if failed to open.
 */
static XFILE * open_input_blk(const CSTRING location, const CSTRING tag, int blk) {
    CSTRING filename = NULL;
    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    if (Current == NULL) {
        /* use the tag on its own */
        filename = copy_CSTRING(tag);
    }
    else {
        /* create input file name from current intensities name */
        switch (Input_Format) {
            case E_TXT:
                filename = new_name_body(Current, tag, blk);
                break;

            case E_CIF:
                filename = new_name_suffix(Current, tag, blk);
                break;

            default: ;
        }
    }

    if (filename != NULL) {
        if (full_path(location, filename, &filepath)) {
            fp = xfopen(filepath, XFILE_UNKNOWN, "r" );

            if (xfisnull(fp)) {
                message(E_OPEN_FAIL_SS, MSG_ERR, "Input", filepath);
                fp = xfclose(fp);
            }
            else {
                message(E_INPUT_FOUND_SS, MSG_INFO, tag, filename);
            }
        }
    }
    
    if (fp != NULL) {
        Last_Input = free_CSTRING(Last_Input);
        Last_Input = copy_CSTRING(filename);
    }

    free_CSTRING(filename);
    free_CSTRING(filepath);
    
    return fp;
}

/**
 * Scan the input directory for any files that match the specified pattern.
 * Result placed in Dir_List. Return the number found.
 */
static int scan_inputs(void) {

    /* hmhm new message func to return a message string */
    static const char *ERRMESS = "Error: Couldn't open the directory: \'%s\'";

    int num = scandir (Pattern_Path, &Dir_List, match_pattern, alphasort);

    if (num < 0) {
        CSTRING msg = new_CSTRING(strlen(ERRMESS) + strlen(Pattern_Path));
        sprintf(msg, ERRMESS, Pattern_Path);
        perror (msg);
        free_CSTRING(msg);
    }
    else {
        for (int cnt = 0; cnt < num; ++cnt) {
//            message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, Dir_List[cnt]->d_name);
        }
    }

    return num;
}

/**
 * Set up next lane and tile range from supplied lane and tile string.
 * Parse and store the lane and tile string, expected format is Ln[-n]Tn[-n].
 * Result stored in LTMin/Max.
 * Returns false if a problem with the string.
 */
static bool set_lanetile(const CSTRING lanetilestr) {

    char *ptr = lanetilestr;
    bool ok = true;

    /* initialise lane and tile values */
    LTCurrent.lane = 0;
    LTCurrent.tile = 0;
    LTMin = LTMax = LTCurrent;

    /* decode and store lane and tile */
    while ((ptr[0] != 0) && ok) {
        ok = get_lanetile_range(ptr, &ptr);
    }

    if (ok) {
        if (lanetile_isnull(LTMin)) {
            message(E_BAD_TXT_SS, MSG_ERR, LTMESS_TEXT, "need lane and tile");
            ok = false;
        }
        else {
            CSTRING lanestr = format_range(LTMin.lane, LTMax.lane);
            CSTRING tilestr = format_range(LTMin.tile, LTMax.tile);
            message(E_LANETILE_SS, MSG_INFO, lanestr, tilestr);
            free_CSTRING(lanestr);
            free_CSTRING(tilestr);
        }
    }

    return ok;
}

/**
 * Set the spikein file location. May be full path or relative to input path.
 * Returns true if path supplied and exists.
 */
static bool set_spikein(void) {

    if (Spikearg_Path == NULL) { return false; }
    
    if (Spikearg_Path[0] == PATH_DELIM) {
        /* full path given, use as is */
        Spikein_Path = copy_CSTRING(Spikearg_Path);
    }
    else {
        /* path relative to general input */
        size_t inputlen = strlen(Input_Path);
        size_t spikelen = strlen(Spikearg_Path);
        Spikein_Path = new_CSTRING(inputlen + spikelen + strlen(PATH_DELIMSTR));
        strcpy(Spikein_Path, Input_Path);
        if (Spikein_Path[inputlen - 1] != PATH_DELIM) {
            strcat(Spikein_Path, PATH_DELIMSTR);
        }
        strncat(Spikein_Path, Spikearg_Path, spikelen);
    }
    message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, Spikein_Path);

    /* check specified path exists */
    if (check_dir(Spikein_Path) != E_ISDIR) {
        message(E_BAD_DIR_SS, MSG_FATAL, "spike-in data", Spikein_Path);
        return false;
    }
    else {
        return true;
    }
}


/* public functions */

/** Return whether specified output directory exists or can be created. */
bool check_outdir(const CSTRING dirname, const char *type_str) {

    DIROPT ret = check_dir(dirname);

    if (ret == E_NODIR) {
        message(E_BAD_DIR_SS, MSG_FATAL, type_str, dirname);
        return false;
    }
    else if (ret == E_NOEXIST) {
        /* try to create it; mkdir returns zero on success */
        if (mkdir (dirname, DIR_MODE) != 0){
            message(E_NOCREATE_DIR_SS, MSG_FATAL, type_str, dirname);
            return false;
        }
        else{
            message(E_CREATED_DIR_SS, MSG_INFO, type_str, dirname);
        }
    }
    return true;
}

/** Return the name of the current intensities file. */
CSTRING get_current_file(void) {

    if (Current == NULL) {
        return "";
    }
    else {
        return Current;
    }
}

/** Return the selected input format. */
INFORM get_input_format(void) {

    return Input_Format;
}

/** Return the selected input path. */
CSTRING get_input_path(void) {

    return Input_Path;
}

/** Return the name of the last opened spike-in file. */
CSTRING get_last_spikein(void) {

    if (Last_Spikein == NULL) {
        return "";
    }
    else {
        return Last_Spikein;
    }
}

/**
 * Return the next run-folder lane and tile.
 * Create a current filename for use in output: s_L_TTTT.
 */
LANETILE get_next_lanetile(void) {

    static const char *NAME = "s_%u_%04u";
    static const size_t LEN = 8;

    if (lanetile_isnull(LTCurrent)) {
        /* first lane/tile */
        LTCurrent = LTMin;
    }
    else {
        if (LTCurrent.tile < LTMax.tile) {
            /* next tile */
            LTCurrent.tile++;
        }
        else {
            if (LTCurrent.lane < LTMax.lane) {
                /* next lane */
                LTCurrent.lane++;
                LTCurrent.tile = LTMin.tile;
            }
            else {
                /* run out, set to null */
                LTCurrent.lane = 0;
                LTCurrent.tile = 0;
            }
        }
    }

    Current = free_CSTRING(Current);
    if (!lanetile_isnull(LTCurrent)) {
        /* create a virtual lane/tile filename */
        Current = new_CSTRING(LEN);
        sprintf(Current, NAME, LTCurrent.lane, LTCurrent.tile);
    }

    return LTCurrent;
}

/** Return true if LANETILE is null. */
bool lanetile_isnull(const LANETILE lanetile) {

    return (lanetile.lane == 0) || (lanetile.tile == 0);
}

/** Return if a predetermined matrix input file is specified. */
bool matrix_from_file(IOTYPE idx) {

    if (idx < E_NMATRIX) {
        return (Matrix[idx] != NULL);
    }
    else {
        return false;
    }
}

/**
 * Open a predetermined input matrix file.
 * Return the file handle or NULL if failed to open.
 */
XFILE * open_matrix(IOTYPE idx) {

    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    /* make a full path name */
    if (Matrix[idx] == NULL) {
        message(E_DEBUG_SSD_SD, MSG_DEBUG, __func__, __FILE__, __LINE__, "No Matrix file at position:", idx);
    }
    else {
        if (full_path(Input_Path, Matrix[idx], &filepath)) {
           fp = xfopen(filepath, XFILE_UNKNOWN, "r" );

           if (xfisnull(fp)) {
               message(E_OPEN_FAIL_SS, MSG_ERR, "Input matrix", filepath);
               fp = xfclose(fp);
           }
           else {
               message(E_INPUT_FOUND_SS, MSG_INFO, "matrix", Matrix[idx] );
           }
        }
        free_CSTRING(filepath);
    }

    return fp;
}

/**
 * Open the next intensities file in the directory.
 * Output error and go to next if fails to open.
 * Return the file handle or NULL if no more files.
 * Also closes previous intensities file if any.
 */
XFILE * open_next(XFILE *fplast) {

    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    if (fplast != NULL) {
        xfclose(fplast);
        Current = free_CSTRING(Current);
    }

    while ((fp == NULL) && (++Index < Dir_Num)) {
        Current = copy_CSTRING(Dir_List[Index]->d_name);
	LTCurrent = parse_lanetile_from_filename(Current);
        message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, Current);

        /* make a full path name */
        if (full_path(Pattern_Path, Current, &filepath)) {
           fp = xfopen(filepath, XFILE_UNKNOWN, "r" );

           if (xfisnull(fp)) {
               message(E_OPEN_FAIL_SS, MSG_ERR, "Input", filepath);
               fp =xfclose(fp);
           }
           else {
               message(E_INPUT_FOUND_SS, MSG_INFO, "intensities", Current );
           }
        }
        free_CSTRING(filepath);
    }

    if ((fp == NULL) && (Current != NULL)) {
        Current = free_CSTRING(Current);
    }
    return fp;
}

/** Open an output file with no block suffix. */
XFILE * open_output(const CSTRING tag) {
    return open_output_blk(tag, BLK_SINGLE);
}

/**
 * Open an output file corresponding to current intensities file with supplied tag.
 * A non-negative blk indicates a block suffix should be added to the name.
 * BLK_SINGLE indicates no block suffix.
 * BLK_APPEND indicates open in append mode.
 * Return the file handle or NULL if failed to open.
 */
XFILE * open_output_blk(const CSTRING tag, int blk) {

    CSTRING filename = NULL;
    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    if (Current == NULL) {
        /* use the tag on its own */
        filename = copy_CSTRING(tag);
    }
    else {
        /* create output file name from current intensities */
        switch (Input_Format) {
            case E_TXT:
                filename = new_name_body(Current, tag, blk);
                break;

            case E_CIF:
                filename = new_name_suffix(Current, tag, blk);
                break;

            default: ;
        }
    }

    if (filename != NULL) {
        if (full_path(Output_Path, filename, &filepath)) {
            const char *mode_str = ((blk == BLK_APPEND) ? "a" : "w");
            fp =  xfopen(filepath, XFILE_RAW, mode_str );

            if (xfisnull(fp)) {
                message(E_OPEN_FAIL_SS, MSG_ERR, "Output", filepath);
                fp = xfclose(fp);
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

/**
 * Open an output file associated with a program run with supplied tag.
 * Return the file handle or NULL if failed to open.
 */
XFILE * open_run_output(const CSTRING tag) {

    CSTRING filename = NULL;
    CSTRING filepath = NULL;
    XFILE *fp = NULL;

    /* get the log file path and extract the file name */
    CSTRING logname = name_only(get_message_path());

    if (logname != NULL) {
        /* create new name from log name */
        filename = new_name_suffix(logname, tag, BLK_SINGLE);

        if (filename != NULL) {
            if (full_path(Output_Path, filename, &filepath)) {
                fp =  xfopen(filepath, XFILE_RAW, "w" );

                if (xfisnull(fp)) {
                    message(E_OPEN_FAIL_SS, MSG_ERR, "Output", filepath);
                    fp = xfclose(fp);
                }
                else {
                    message(E_DEBUG_SSD_S, MSG_DEBUG, __func__, __FILE__, __LINE__, filename);
                }
            }
        }
    }

    free_CSTRING(filename);
    free_CSTRING(filepath);
    free_CSTRING(logname);

    return fp;
}

/** Open a spike-in data file if any. */
XFILE * open_spikein(int blk) {
    
    if (Spikein_Path == NULL) {
        return NULL;
    }
    else {
        free_CSTRING(Last_Spikein);
        XFILE *fp = open_input_blk(Spikein_Path, (CSTRING)SPIKEIN_TAG, blk);
        if (fp != NULL) {
            Last_Spikein = copy_CSTRING(Last_Input);
        }
        return fp;
    }
}

/** Return if run-folder selected. */
bool run_folder(void) {

    return RunFolder;
}

/**
 * Set the input format. Text must match one of the input format text list. Ignores case.
 * Returns true if match found.
 */
bool set_input_format(const char *inform_str) {

    /* match to one of the possible options */
    int matchidx = match_string(inform_str, INFORM_TEXT, E_INFORM_NUM);
    if (matchidx >= 0) {
        Input_Format = (INFORM)matchidx;
        return true;
    }
    else {
        return false;
    }
}

/** Set file location information. */
void set_location(const CSTRING path, IOTYPE idx){

    switch(idx) {
        case E_INPUT:
            Input_Path = copy_CSTRING(path);
            break;
        case E_OUTPUT:
            Output_Path = copy_CSTRING(path);
            break;
        case E_SPIKEIN:
            Spikearg_Path = copy_CSTRING(path);
            break;
        case E_CROSSTALK:
        case E_NOISE:
        case E_PARAMA:
        case E_QUALTAB:
            Matrix[idx] = copy_CSTRING(path);
            break;
        default: ;
    }
}

/**
 * Set the sample name for tiles to be run.
 */
void set_sample_name(const CSTRING sample_name){
	if(NULL==sample_name){
	}
	Sample_Name = copy_CSTRING(sample_name);
}

/**
 * Set the input filename pattern to match. Moves any path parts to pattern path.
 * Checks pattern argument supplied, and at least one input file found.
 * If a run-folder then set lane and tile range instead.
 */
bool set_pattern(const CSTRING pattern) {

    if (RunFolder) {
        /* treat pattern as a supplied lane and tile string */
        return set_lanetile(pattern);
    }

    clear_pattern();
    Pattern = copy_CSTRING(pattern);

    /* ensure any path parts are in the pattern path */
    Pattern_Path = move_path(Input_Path, &Pattern);

    /* check pattern match supplied */
    if (Pattern != NULL) {
        Pattern_Len = strlen(Pattern);
    }

    if ((Pattern_Len == 0) || (Pattern_Path == NULL)) {
        message(E_NOPATTERN_FILE_S, MSG_ERR, pattern);
        return false;
    }

    /* scan for matching input files */
    Dir_Num = scan_inputs();
    if (Dir_Num < 0) {
        /* failed to access directory */
        return false;
    }

    if (Dir_Num == 0) {
        message (E_NOINPUT_SS, MSG_ERR, Pattern_Path, Pattern);
        return false;
    }
    else {
        /* at least one input file */
        message(E_PATTERN_MATCH_SD, MSG_INFO, pattern, Dir_Num);
        return true;
    }
}

/** Set run-folder flag. */
void set_run_folder(void) {

    RunFolder = true;
}

/** Return if spike-in data selected. */
bool spike_in(void) {
    return (Spikein_Path != NULL);
}

/**
 * Start up; call at program start after options.
 * Checks input and output directories exists and creates the match substring.
 * Return true if no errors.
 */
bool startup_dirio(void) {

    /* set default i/o paths if not supplied */
    if (Input_Path == NULL) {
        Input_Path = copy_CSTRING((CSTRING)DEFAULT_PATH);
    }
    if (Output_Path == NULL) {
        Output_Path = copy_CSTRING((CSTRING)DEFAULT_PATH);
    }

    /* check specified input path exists */
    if (check_dir(Input_Path) != E_ISDIR) {
        message(E_BAD_DIR_SS, MSG_FATAL, "input", Input_Path);
        return false;
    }

    /* check specified output path exists */
    if (!check_outdir(Output_Path, "output")) {
         return false;
    }
    /* create the input filename match substring */
    make_substring();

    message(E_INPUT_DIR_SS, MSG_INFO, "Input", Input_Path);
    message(E_OPT_SELECT_SS, MSG_INFO, "Input format" ,INFORM_MESS_TEXT[Input_Format]);

    /* check for run-folder */
    if (RunFolder) {
        if (Input_Format == E_CIF) {
            message(E_OPT_SELECT_SS, MSG_INFO, "Run-folder", "");
        }
        else {
            /* invalid with txt */
            message(E_BAD_RUNOPT, MSG_FATAL);
            return false;
        }
    }

    /* check for spike-in data */
    if (Spikearg_Path != NULL) {
        if (set_spikein()) {
            message(E_INPUT_DIR_SS, MSG_INFO, "Spike-in data", Spikein_Path);
        }
        else {
            /* error message already given */
            return false;
        }
    }

    message(E_OUTPUT_DIR_S, MSG_INFO, Output_Path);
    return true;
}

/** Tidy up; call at program shutdown. */
void tidyup_dirio(void) {

    /* clear any pattern info */
    clear_pattern();

    /* free string memory */
    Input_Path = free_CSTRING(Input_Path);
    Output_Path = free_CSTRING(Output_Path);
    IntenSubstr = free_CSTRING(IntenSubstr);
    Spikearg_Path = free_CSTRING(Spikearg_Path);
    Spikein_Path = free_CSTRING(Spikein_Path);
    Last_Input = free_CSTRING(Last_Input);
    Last_Spikein = free_CSTRING(Last_Spikein);
    for (IOTYPE idx = (IOTYPE)0; idx < E_NMATRIX; idx++) {
        Matrix[idx] = free_CSTRING(Matrix[idx]);
    }
}

LANETILE parse_lanetile_from_filename ( const char * fn){
	LANETILE lt;
	int ret = sscanf(fn,"%*[^0123456789]%u_%u",&lt.lane,&lt.tile);
	if(2!=ret){
		warnx("Failed to parse lane and tile numbers from %s.\nWill use default if no others are found",fn);
	}
	return lt;
}

LANETILE get_current_lanetile ( void ){
	return LTCurrent;
}

CSTRING get_sample_name ( void ){
	return Sample_Name;
}

