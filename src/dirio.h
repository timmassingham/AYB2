/**
 * \file dirio.h
 * Public parts of I/O Environment.
 *   - Directory search.
 *   - Input file pattern match.
 *   - File open/close
  *//*
 *  Created : 16 Mar 2010
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

#ifndef DIRIO_H_
#define DIRIO_H_

#include <stdbool.h>
#include "utility.h"            // for CSTRING
#include "xio.h"                // for XFILE

/**
 * Types of file location information. Also used as index into predetermined input matrices.
 * E_NMATRIX indicates number of such matrices.
 */
typedef enum IOTypeT {E_CROSSTALK, E_NOISE, E_PHASING, E_INPUT, E_OUTPUT, E_NMATRIX = 3} IOTYPE;


/* function prototypes */

bool check_outdir(const CSTRING dirname, const char * typestr);
CSTRING get_current_file();
CSTRING get_pattern();

XFILE * open_matrix(IOTYPE mode);
XFILE * open_next(XFILE *fplast);
XFILE * open_output(CSTRING tag);

void set_location(const CSTRING path, IOTYPE mode);
void set_pattern(const CSTRING pattern);

bool startup_dirio();
void tidyup_dirio();

#endif /* DIRIO_H_ */
