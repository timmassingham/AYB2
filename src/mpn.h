/** 
 * \file mpn.h
 * Public parts of Calculations of terms for Parameter Estimation.
 *//* 
 *  Created : 9 Jun 2010
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

#ifndef MPN_H_
#define MPN_H_

#include "matrix.h"
#include "nuc.h"
#include "tile.h"


/* function prototypes */

//MAT calculateIbar( const ARRAY(int16_t) intmat, const MAT we, MAT Ibar){
MAT calculateIbar( const TILE tile, const MAT we, MAT Ibar);
MAT calculateSbar( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const uint_fast32_t ncycle, MAT Sbar);
MAT calculateWe( const MAT lssi, MAT we);
real_t calculateWbar( const MAT we);
MAT calculateJ( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const uint_fast32_t ncycle, MAT J);
//MAT calculateK( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const ARRAY(int16_t) ints, const uint_fast32_t ncycle, MAT K);
MAT calculateK( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const TILE tile, const uint_fast32_t ncycle, MAT K);
MAT calculateMlhs( const MAT var, const real_t wbar, const MAT SbarT, const MAT P, const MAT Jt, real_t * tmp, MAT lhs);
MAT calculateMrhs( const MAT var, const MAT IbarT, const MAT P, const MAT Kt, real_t * tmp, MAT rhs);
MAT calculatePlhs( const real_t wbar, const MAT Sbar, const MAT Mt, const MAT J, real_t * tmp, MAT lhs);
MAT calculatePrhs( const MAT Ibar, const MAT Mt, const MAT Sbar, const MAT N, const MAT K, real_t * tmp, MAT rhs);
real_t calculateDeltaLSE(const MAT Mt, const MAT P, const MAT N, const MAT J, const MAT K, real_t * tmp);

MAT calculateNewJ(const MAT lambda, const ARRAY(NUC) bases, const MAT we, const int ncycle, MAT newJ);
//MAT calculateNewK(const MAT lambda, const ARRAY(NUC) bases,const ARRAY(int16_t) intmat, const MAT we, const int ncycle, MAT newK);
MAT calculateNewK(const MAT lambda, const ARRAY(NUC) bases, const TILE tile, const MAT we, const int ncycle, MAT newK);
MAT calculateLhs( const real_t wbar,const MAT J, const MAT Ibar, MAT lhs);
MAT calculateRhs( const MAT K, const MAT Sbar, MAT rhs);

int solverChol( MAT lhs, MAT rhs, real_t * tmp);
int solverSVD(MAT lhs, MAT rhs, real_t * tmp, const real_t delta_diag);
int solverZeroSVD(MAT lhs, MAT rhs, real_t * tmp, const real_t delta_diag);
int solverNNLS(MAT lhs, MAT rhs, real_t * tmp, const real_t delta_diag);

#endif /* MPN_H_ */
