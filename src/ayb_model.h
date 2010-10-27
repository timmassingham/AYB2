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
#include <stdint.h>
#include "utility.h"
#include "xio.h"

/** AYB defined as a hidden data structure. Access via structure pointer. */
typedef struct AybT * AYB;

/* function prototypes */

/* standard functions */
AYB new_AYB(const uint32_t ncycle, const uint32_t ncluster);
AYB free_AYB(AYB ayb);
AYB copy_AYB(const AYB ayb);
void show_AYB(XFILE * fp, const AYB ayb, bool showall);

bool analyse_tile (const int argc, char ** const argv, XFILE *fp);
bool set_composition(const char *comp_str);
void set_niter(const char *niter_str);
bool set_output_format(const char *outform_str);
void set_show_working(void);
void set_simdata(const CSTRING simdata_str);
bool set_solver(const char *solver_str);
bool startup_model(void);
void tidyup_model(void);

#endif /* AYB_MODEL_H_ */
