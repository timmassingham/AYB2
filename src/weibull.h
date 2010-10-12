/** 
 * \file weibull.h
 * Public parts of routines for Weibull distribution
 *//* 
 *  Created : 2010
 *  Author : Tim Massingham
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
 *
 *  This file is part of the AYB base-calling software.
 *
 *  AYB is free software: you can redistribute it and/or modify
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

#ifndef _WEIBULL_H
#define _WEIBULL_H

#include <stdbool.h>
#include <stdint.h>
#include "utility.h"

typedef struct {
   real_t e1,e2;
} pair_real;

real_t pweibull(  real_t x,  real_t shape,  real_t scale,  bool tail,  bool logp);
real_t qweibull(  real_t p,  real_t shape,  real_t scale,  bool tail,  bool logp);
real_t dweibull(  real_t x,  real_t shape,  real_t scale,  bool logd );
pair_real fit_weibull( const real_t * x_orig, const uint32_t n );

#endif
