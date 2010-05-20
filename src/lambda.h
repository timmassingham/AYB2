/** 
 * \file lambda.h
 * Public parts of Lambda Calculation.
 *//* 
 *  Created : 23 Apr 2010
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

#ifndef LAMBDA_H_
#define LAMBDA_H_

#include "matrix.h"
#include "nuc.h"
#include "utility.h"

/* function prototypes */

real_t estimate_lambdaOLS( const MAT p, const NUC * base);
real_t estimate_lambdaWLS( const MAT p, const NUC * base, const real_t oldlambda, const real_t * v);

#endif /* LAMBDA_H_ */
