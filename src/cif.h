/** 
 * \file cif.h
 * Public parts of Interface to CIF format intensity files.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2008,2009 by Tim Massingham
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the ciftool base-calling software.
 *
 *  ciftool is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ciftool is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ciftool.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CIF_H_
#define CIF_H_

#include <stdbool.h>
#include <stdint.h>
#include "utility.h"
#include "xio.h"

typedef union { int8_t * i8; int16_t * i16; int32_t * i32;} encInt;
typedef struct cifData * CIFDATA;

// Access
uint8_t cif_get_version( const CIFDATA cif );
uint8_t cif_get_datasize( const CIFDATA cif );
uint16_t cif_get_firstcycle( const CIFDATA cif );
uint16_t cif_get_ncycle( const CIFDATA cif );
uint32_t cif_get_ncluster( const CIFDATA cif );
encInt cif_get_const_intensities( const CIFDATA cif);

real_t cif_get_real (const CIFDATA cif, const uint32_t cl, const uint32_t base, const uint32_t cy);
void cif_set_from_real (CIFDATA cif, const uint32_t cl, const uint32_t base, const uint32_t cy, real_t x );

// Creation/Deletion
CIFDATA create_cif (const uint32_t ncycle, const uint32_t ncluster);
void free_cif ( CIFDATA cif);

// Other
CIFDATA readCIFfromFile ( const char * fn, const XFILE_MODE mode);
CIFDATA readCIFfromStream ( XFILE * ayb_fp );
CIFDATA readCIFfromDir ( const char * fn, const uint32_t lane, const uint32_t tile, const XFILE_MODE mode);
bool writeCIFtoFile ( const CIFDATA  cif, const char * fn, const XFILE_MODE mode);
bool writeCIFtoStream ( const CIFDATA  cif, XFILE * ayb_fp);
bool write2CIFfile ( const char * fn, const XFILE_MODE mode, const encInt  intensities, const uint16_t firstcycle, const uint32_t ncycle, const uint32_t ncluster, const uint8_t nbyte);
void showCIF ( XFILE * ayb_fp, const CIFDATA cif, uint32_t mcluster, uint32_t mcycle);
CIFDATA spliceCIF(const CIFDATA cif, uint32_t ncycle, uint32_t offset);

#endif /* CIF_H_ */
