/** 
 * \file ayb_model.h
 * Public parts of Top Level Modelling.
 *//* 
 *  Created : 14 Apr 2010
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

#ifndef AYB_MODEL_H_
#define AYB_MODEL_H_

#include <stdbool.h>
#include "dirio.h"
#include "utility.h"
#include "xio.h"


/* function prototypes */

RETOPT analyse_tile (const int argc, char ** const argv);
void read_intensities_file(XFILE *fp, unsigned int ncycle);
void read_intensities_folder(const char *root, LANETILE lanetile, unsigned int ncycle);
void set_niter(const char *niter_str);
bool set_output_format(const char *outform_str);
void set_simdata(const CSTRING simdata_str);
bool startup_model(void);
void tidyup_model(void);

#endif /* AYB_MODEL_H_ */
