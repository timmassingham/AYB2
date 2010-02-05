/*
 *  Copyright (C) 2010 by Tim Massingham
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

#ifndef _CLUSTER_H
#define _CLUSTER_H

#include "xio.h"
#include "matrix.h"

struct _struct_cluster {
	unsigned long int x,y;
	MAT signals;
};
typedef struct _struct_cluster * CLUSTER;

// Standard functions
CLUSTER new_CLUSTER();
void free_CLUSTER(CLUSTER cluster);
CLUSTER copy_CLUSTER(const CLUSTER cluster);
void show_CLUSTER(XFILE * fp, const CLUSTER cluster);

// Input
CLUSTER read_known_CLUSTER( XFILE * fp, const unsigned int ncycle);
CLUSTER read_unknown_CLUSTER( XFILE * fp);

#endif /* _CLUSTER_H */