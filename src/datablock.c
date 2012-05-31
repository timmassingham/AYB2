/** 
 * \file datablock.c
 * Data Block Class.
 * Decodes and stores how the cycle data in an intensity file should be grouped for analysis.
 
 * The blockstring (b) program option has the form RnInCn, decoded as:
 * - R => read
 * - I => ignore
 * - C => concatenate to previous block (first R must precede first C)
 * 
 * The default is all cycles in one block.
 *
 * The parsing function checks the blockstring for syntax and decodes it, 
 * storing the results in a list of datablocks. 
 * Get next block then returns each datablock in turn until no more.
 *//* 
 *  Created : 28 May 2010
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

#include <stdlib.h>
#include "datablock.h"
#include "message.h"
#include "utility.h"

#define X(A) A ## DATABLOCK
    #include "list.def"
#undef X

/* constants */

static const char *MESS_TEXT = "Blockstring option";    ///< Blockstring option name for messages.

/* members */

static bool DefaultBlock = true;                        ///< Default to all available cycles in one block if no argument.
static unsigned int TotalCycle = 0;                     ///< Total number of cycles to analyse.
static unsigned int NumBlock = 0;                       ///< Number of distinct blocks to analyse.

/* doxygen confused by LIST */
/** BlockList = NULL; List of data blocks.
 * Also static LIST(DATABLOCK) Current = NULL; 
 *
 * Current data block.
 */
static LIST(DATABLOCK) BlockList = NULL;
static LIST(DATABLOCK) Current = NULL;


/* private functions */

/* helps doxygen confused by LIST */
static void __attribute__((__used__)) null1(void){}

/** Return the type of block string token. */
static BLOCKTYPE decode_token(const char ch) {

    switch (ch) {
        case 'R':
        case 'r': return E_READ;
        case 'I':
        case 'i': return E_IGNORE;
        case 'C':
        case 'c': return E_CONCAT;

        default : return E_ERR;
    }
}

/** Free the list of blocks. */
static void free_blocklist(void) {
    free_LIST(DATABLOCK)(BlockList);
    BlockList = NULL;
    TotalCycle = 0;
    NumBlock = 0;
}


/* public functions */

/* standard functions */
DATABLOCK new_DATABLOCK(void) {
    DATABLOCK datablock = calloc(1, sizeof(*datablock));
    return datablock;
}

DATABLOCK free_DATABLOCK(DATABLOCK datablock) {
    if(NULL == datablock) {return NULL;}
    xfree(datablock);
    return NULL;
}

DATABLOCK copy_DATABLOCK(const DATABLOCK datablock) {
    if (NULL == datablock) {return NULL;}
    DATABLOCK newblock = new_DATABLOCK();
    if(NULL == newblock) {return NULL;}
    newblock->type = datablock->type;
    newblock->num = datablock->num;
    return newblock;
}

void show_DATABLOCK(XFILE * fp, const DATABLOCK datablock) {
    if(NULL == fp) {return;}
    if(NULL == datablock) {return;}
    xfprintf(fp, "Block type: %d; Cycles: %u\n", datablock->type, datablock->num);
}

/** Return if default block selected. */
bool get_defaultblock(void) {
    return DefaultBlock;
}

/** Return the next data block. Returns NULL if are none or no more. */
DATABLOCK get_next_block(void) {
    if (BlockList == NULL) {return NULL;}
    if (Current == NULL) {
        Current = BlockList;
    }
    else {
        Current = Current->nxt;
    }
    if (Current == NULL) {return NULL;}
    return Current->elt;
}

/** Return the number of distinct blocks to analyse. */
unsigned int get_numblock(void) {
    return NumBlock;
}

/** Return the total number of cycles to be analysed. */
unsigned int get_totalcycle(void) {
    return TotalCycle;
}

/**
 * Parse the blockstring option to find the pattern of data blocks.
 * Stores the result in the BlockList.
 * Returns false if a problem with the option string.
 */
bool parse_blockopt(const char *blockstr) {

    DATABLOCK newblock = NULL;
    LIST(DATABLOCK) tail = NULL;
    BLOCKTYPE type;
    const char *ch = NULL;
    char *endptr = NULL;
    int cycles;
    bool ok = true;

    if (blockstr == NULL) {
        message(E_BAD_TXT_SS, MSG_FATAL, MESS_TEXT, "no blockstring supplied");
        return false;
    }
    /* argument supplied, turn default off even if an error detected */
    DefaultBlock = false;

    /* parse the string from the beginning */
    ch = blockstr;
    while (ok && (ch != NULL) && (ch[0] != 0)) {
        /* create a new data block */
        newblock = new_DATABLOCK();
        
        /* decode the blockstring token */
        type = decode_token(ch[0]);

        if (type == E_ERR) {
            message(E_BAD_CHAR_SC, MSG_FATAL, MESS_TEXT, ch[0]);
            ok = false;
        }
        else {
            /* check for concat before first read */
            if ((type == E_CONCAT) && (NumBlock == 0)) {
                message(E_BAD_TXT_SS, MSG_FATAL, MESS_TEXT, "concatenate before first read");
                ok = false;
            }
            else {
                /* extract cycles from second character onwards */
                cycles = strtol(++ch, &endptr, 0);
                if (cycles <= 0) {
                    message(E_BAD_NUM_SS, MSG_FATAL, MESS_TEXT, ch);
                    ok = false;
                }
                else {
                    newblock->type = type;
                    newblock->num = cycles;
                    /* increment cycle count and block count if new block */
                    TotalCycle += cycles;
                    if (type == E_READ) {
                        NumBlock++;
                    }
                }
            }
        }
        
        if (ok) {
            /* add new data block to list */
            if (BlockList == NULL) {
                BlockList = cons_LIST(DATABLOCK)(newblock, BlockList);
                tail = BlockList;
            }
            else {
                tail = rcons_LIST(DATABLOCK)(newblock, tail);
            }
        }
        else {
            /* free the new data block */
            xfree(newblock);
        }

        /* go to next non-numeric character */
        ch = endptr;
    }

    if (!ok) {
        /* clear the structure */
        free_blocklist();
    }

    return ok;
}

/** Tidy up; call at program shutdown. */
void tidyup_datablock(void) {
    /* clear the structure */
    free_blocklist();
}
