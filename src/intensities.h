/** 
 * \file intensities.h
 * Public parts of Intensities Calculation.
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


#ifndef INTENSITIES_H_
#define INTENSITIES_H_

#include "matrix.h"
#include "nuc.h"

/* function prototypes */
MAT process_intensities(const MAT intensities,
                        const MAT Minv_t, const MAT Pinv_t, const MAT N, MAT ip);
MAT expected_intensities(const real_t lambda, const NUC * bases,
                         const MAT M, const MAT P, const MAT N, MAT e);
MAT processNew(const struct structLU AtLU, const MAT N, const MAT intensities, MAT p);


#endif /* INTENSITIES_H_ */
