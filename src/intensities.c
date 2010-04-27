/** 
 * \file intensities.c
 * Intensities Processing.
 *//* 
 *  Created : 19 Apr 2010
 *  Author  : Tim Massingham/Hazel Marsden
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

#include <string.h>
#include "intensities.h"
#include "nuc.h"

/* constants */
/* None      */


/* members */
/* None    */


/* private functions */
/* None              */


/* public functions */

/**
 * Process intensities.
 * ip = Minv %*% (Intensities-N) %*% Pinv
 *   - Uses identity: Vec(ip) = ( Pinv^t kronecker Minv) Vec(Intensities-N)
 *   - Storing Intensities-N as an intermediate saved < 3%
 *   - Calculating ip^t rather than ip (pcol loop is over minor index) made no difference
 *   - Using Pinv rather than Pinv^t makes little appreciable difference
 */
MAT process_intensities(const MAT intensities,
                        const MAT Minv_t, const MAT Pinv_t, const MAT N, MAT ip) {

    validate(NULL != intensities, NULL);
    validate(NULL != Minv_t, NULL);
    validate(NULL != Pinv_t, NULL);
    validate(NULL != N, NULL);

    const uint32_t ncycle = Pinv_t->nrow;
    if (NULL==ip){
        ip = new_MAT(NBASE, ncycle);
        validate(NULL != ip, NULL);
    }
//    bzero(p->x,p->nrow * p->ncol * sizeof(real_t));
    memset(ip->x, 0, ip->nrow * ip->ncol * sizeof(real_t));

    for (uint32_t icol = 0; icol < ncycle; icol++) {    // Columns of Intensity
        for (uint32_t base = 0; base < NBASE; base++) { // Bases (rows of Minv, cols of Minv_t)
            real_t dp = 0;
            for (uint32_t chan = 0; chan < NBASE; chan++) {  // Channels
                dp += Minv_t->x[base * NBASE + chan] *
                        (intensities->x[icol * NBASE + chan] - N->x[icol * NBASE + chan]);
            }
            for (uint32_t pcol = 0; pcol < ncycle; pcol++) { // Columns of ip
                ip->x[pcol * NBASE + base] += Pinv_t->x[icol * ncycle + pcol] * dp;
            }
        }
    }

    return ip;
}
