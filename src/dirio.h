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


/** Path delimiter character. For linux, other OS? */
static const char PATH_DELIM = '/';

/** Possible input formats. */
typedef enum InFormT {E_TXT, E_CIF, E_INFORM_NUM} INFORM;

/**
 * Types of file location information. Also used as index into predetermined input matrices.
 * E_NMATRIX indicates total number of such matrices, and E_MNP the modelling ones.
 */
typedef enum IOTypeT {E_CROSSTALK, E_NOISE, E_PHASING, E_QUALTAB, E_INPUT, E_OUTPUT, E_MNP = 3, E_NMATRIX = 4} IOTYPE;

/**
 * Open output file special block options.
 * Indicate append to existing file, if any, or a single data block (no block suffix).
 */
enum NoBlkT {BLK_APPEND = -2, BLK_SINGLE = -1};

/** Pair type allows storage and return of lane and tile. */
typedef struct LaneTileT {unsigned int lane; unsigned int tile;} LANETILE;


/* function prototypes */

bool check_outdir(const CSTRING dirname, const char * type_str);
CSTRING get_current_file(void);
INFORM get_input_format(void);
CSTRING get_input_path(void);
LANETILE get_next_lanetile(void);
bool lanetile_isnull(const LANETILE lanetile);
bool matrix_from_file(IOTYPE idx);

XFILE * open_matrix(IOTYPE mode);
XFILE * open_next(XFILE *fplast);
XFILE * open_output(CSTRING tag);
XFILE * open_output_blk(CSTRING tag, int blk);

bool run_folder(void);
bool set_input_format(const char *inform_str);
void set_location(const CSTRING path, IOTYPE mode);
bool set_pattern(const CSTRING pattern);
void set_run_folder(void);

bool startup_dirio(void);
void tidyup_dirio(void);

#endif /* DIRIO_H_ */
