/** 
 * \file spikein.h
 * Public parts of spike-in data.
 *//* 
 *  Created : 19 Jan 2012
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

#ifndef SPIKEIN_H_
#define SPIKEIN_H_

#include <stdbool.h>
#include <stdint.h>
#include "nuc.h"
#include "xio.h"


/** Spike-in data structure. */
struct SpikeInT {
    uint_fast32_t knum;
    ARRAY(NUC) kbases;    
};

typedef struct SpikeInT * SPIKEIN;


/* function prototypes */

// Standard functions
SPIKEIN new_SPIKEIN(const uint_fast32_t ncycle);
SPIKEIN free_SPIKEIN(SPIKEIN spikein);
SPIKEIN copy_SPIKEIN(const SPIKEIN spikein);
void show_SPIKEIN(XFILE * fp, const SPIKEIN spikein);
SPIKEIN read_SPIKEIN(XFILE * fp, const uint_fast32_t ncycle, bool *err);

SPIKEIN get_next_spikein(void);
bool read_spikein_file(const uint_fast32_t ncycle, const int blk);
void tidyup_spikein(void);

#endif /* SPIKEIN_H_ */
