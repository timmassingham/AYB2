/** 
 * \file calibration.h
 * Default calibration table for AYB.
 * Table values in calibration source file generated by ayb_recal.
 *//*
 *  Created : 03 Feb 2011
 *  Author  : Tim Massingham
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
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

#ifndef CALIBRATION_H_
#define CALIBRATION_H_

#include "qual_table.h"

/**
 * Pointer to quality calibration table.
 * Initially points to default. Redirect if new table read in.
 */
extern QUALTAB QTable;

#endif /* CALIBRATION_H_ */
