/**
 * \file statistics.h
 * Public parts of Statistical Functions.
 *//*
 *  Created : 2010
 *  Author : Tim Massingham
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
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
 
#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <stdint.h>
#include "utility.h"

// Compensated means and variances
real_t mean( const real_t * x, const uint32_t n);
real_t variance( const real_t * x, const uint32_t n);
real_t wmean( const real_t * w, const real_t * x, const uint32_t n);
real_t wvariance( const real_t * w, const real_t * x, const uint32_t n);

// Linear regression
real_t * linearRegression( const real_t * x, const real_t * y, const uint32_t n, real_t *res);
real_t * linearResiduals( const real_t * x, const real_t * y, const real_t * param, const uint32_t nobs, real_t * resid);
real_t * wLinearRegression( const real_t * w, const real_t * x, const real_t * y, const uint32_t n, real_t * res);

// Iteratively Rewweighted Least Squares
real_t __attribute__((const)) tukey_biweight( const real_t xsqr, const real_t v);
real_t __attribute__((const)) cauchy( const real_t xsqr, const real_t v);
real_t * iwlsLinearRegression(real_t (*f)(const real_t,const real_t), const real_t * x, const real_t * y, const uint32_t niter, const uint32_t n, real_t * res);

#endif /* STATISTICS_H_ */
