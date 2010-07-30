/** 
 * \file cluster.h
 * Public parts of Cluster Class.
 *//* 
 *  Created : 2010
 *  Authors : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk

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

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "cif.h"
#include "matrix.h"
#include "xio.h"

struct _struct_cluster {
    unsigned long int x,y;
    MAT signals;
};
typedef struct _struct_cluster * CLUSTER;

// Standard functions
CLUSTER new_CLUSTER();
CLUSTER free_CLUSTER(CLUSTER cluster);
CLUSTER copy_CLUSTER(const CLUSTER cluster);
void show_CLUSTER(XFILE * fp, const CLUSTER cluster);

// standard variations
CLUSTER coerce_CLUSTER_from_array(const unsigned int ncycle, real_t * x, real_t ** next);
CLUSTER copy_append_CLUSTER(CLUSTER clustout, const CLUSTER clustin, int colstart, int colend);

// Input
CLUSTER read_cif_CLUSTER(CIFDATA cif, const unsigned int cl, unsigned int ncycle);
CLUSTER read_known_CLUSTER( XFILE * fp, unsigned int *ncycle, bool moredata);
CLUSTER read_unknown_CLUSTER( XFILE * fp);

#endif /* CLUSTER_H_ */
