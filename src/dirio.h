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

/** Select input or output. Used in set_path. */
typedef enum IOModeT {E_INPUT, E_OUTPUT} IOMODE;


/* function prototypes */

/* open the next input file in the directory */
XFILE * open_next(XFILE *fplast);

/* open an output file corresponding to current input file with supplied suffix */
XFILE * open_output(CSTRING tag);

/* set the input/output paths */
void set_path(const CSTRING path, IOMODE mode);

/* set the filename pattern to match to */
void set_pattern(const CSTRING pattern);

/* start up; call at program start after options */
bool startup_dirio();

/* tidy up; call at program shutdown */
void tidyup_dirio();

#endif /* DIRIO_H_ */
