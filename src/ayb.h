/** 
 * \file ayb.h
 * Public parts of Ayb Class.
 *//* 
 *  Created : 16 Mar 2011
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

#ifndef AYB_H_
#define AYB_H_

#include <stdbool.h>
#include <stdint.h>
#include "matrix.h"
#include "tile.h"
#include "utility.h"
#include "xio.h"


/** Indicates an error in the data being processed, usually an overflow or memory allocation. */
static const int DATA_ERR = -1;

/** AYB defined as a hidden data structure. Access via structure pointer. */
typedef struct AybT * AYB;


/* function prototypes */

/* standard functions */
AYB new_AYB(const uint_fast32_t ncycle, const uint_fast32_t ncluster);
AYB free_AYB(AYB ayb);
AYB copy_AYB(const AYB ayb);
void show_AYB(XFILE * fp, const AYB ayb, bool showall);

/* access functions */
real_t * get_AYB_lambdas(AYB ayb, uint_fast32_t *num);
uint_fast32_t get_AYB_ncluster(AYB ayb);
uint_fast32_t get_AYB_ncycle(AYB ayb);
AYB replace_AYB_tile(AYB ayb, const TILE tile);
void show_AYB_bases(XFILE * fp, const AYB ayb, const uint_fast32_t cl);
void show_AYB_quals(XFILE * fp, const AYB ayb, const uint_fast32_t cl);

MAT calculate_covariance(AYB ayb, const bool do_full);
int estimate_bases(AYB ayb, const int blk, const bool lastiter, const bool showdebug);
real_t estimate_MPN(AYB ayb);
bool initialise_model(AYB ayb, const int blk, const bool showdebug);

void set_show_working(const CSTRING optarg);
void set_spike_calib(void);
bool startup_ayb(void);
void tidyup_ayb(void);

#endif /* AYB_H_ */
