/** 
 * \file qual_table.h
 * Public parts of Quality Calibration Table.
 * Table structure matches ayb_recal version of this file.
 *//* 
 *  Created : 01 Feb 2011
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

#ifndef QUAL_TABLE_H_
#define QUAL_TABLE_H_

#include "nuc.h"
#include "utility.h"

/** Calibration table structure. */
struct qual_table_str {
    char *filename, *comment;
    real_t slope, intcept;
    real_t *prior, *next, *triplet;
};

/** Pointer to calibration table structure. */
typedef struct qual_table_str * QUALTAB;

/* function prototypes */

QUALTAB free_qualtab(QUALTAB qtab);
QUALTAB new_qualtab(void);
void show_qualtab(XFILE * fp, const QUALTAB qtab);

real_t adjust_quality(const real_t qual, const NUC prior, const NUC base, const NUC next);

/**  Adjust quality score for first base by setting prior to NUC_AMBIG. */
static inline real_t adjust_first_quality(const real_t qual, const NUC base, const NUC next){
    return adjust_quality(qual, NUC_AMBIG, base, next);
}
/**  Adjust quality score for last base by setting next to NUC_AMBIG. */
static inline real_t adjust_last_quality(const real_t qual, const NUC prior, const NUC base){
    return adjust_quality(qual, prior, base, NUC_AMBIG);
}

void output_quality_table(void);
bool read_quality_table(void);
void set_noqualout(void);
void tidyup_qual_table(void);

#endif /* QUAL_TABLE_H_ */
