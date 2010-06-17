/** 
 * \file datablock.h
 * Public parts of Data Block Class.
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

#ifndef DATABLOCK_H_
#define DATABLOCK_H_

#include <stdbool.h>
#include "xio.h"


/** Possible types in a blockstring option. */
typedef enum BlockTypeT {E_READ, E_IGNORE, E_CONCAT, E_ERR} BLOCKTYPE;

/** Data for a single defined block. */
struct DataBlockT {
    BLOCKTYPE type;
    unsigned int num;
};

typedef struct DataBlockT * DATABLOCK;


/* function prototypes */

// Standard functions
DATABLOCK new_DATABLOCK();
DATABLOCK free_DATABLOCK(DATABLOCK datablock);
DATABLOCK copy_DATABLOCK(const DATABLOCK datablock);
void show_DATABLOCK(XFILE * fp, const DATABLOCK datablock);

DATABLOCK get_next_block();
unsigned int get_numblock() ;
unsigned int get_totalcycle();
bool parse_blockopt(const char *blockstr);
void tidyup_datablock();


#endif /* DATABLOCK_H_ */
