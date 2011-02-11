/** 
 * \file qual_table.c
 * Quality Calibration Table.
 * Includes read in, output and quality score adjustment.
 * Table structure matches ayb_recal version of this file.
 *//* 
 *  Created : 01 Feb 2011
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

#include <stdio.h>
#include <ctype.h>          // for isblank
#include "calibration.h"    // for QTable external
#include "dirio.h"
#include "matrix.h"
#include "message.h"
#include "qual_table.h"
#include "xio.h"


/* constants */

/** Quality calibration items. */
typedef enum QCCIdT {E_INTERCEPT, E_SLOPE, E_PRIOR, E_NEXT, E_TRIPLET} QCCID;
/** Name text for quality calibration table messages and output. */
static const char *MESS_TEXT[] = {"Intercept", 
                                  "Slope", 
                                  "Prior base adjustment", 
                                  "Next base adjustment", 
                                  "Triplet base adjustment"};

/* members */

static bool NewQTable = false;                  ///< Indicates a quality calibration table has been read in.
static bool QualOut = true;                     ///< Unset to turn off quality calibration table output.


/* private functions */

/**
 * Read a single value from a column matrix file.
 * Returns false if fails to read.
 */
static bool get_matrix_single(XFILE * fp, real_t *value) {

    MAT mat = read_MAT_from_column_file(fp);
    if (mat == NULL) {return false;}

    /* matrix read returns at least one element */
    *value = mat->x[0];
    free_MAT(mat);
    return true;
}

/**
 * Read an array of values from a column matrix file.
 * User needs to ensure supplied return array is large enough.
 * Returns false if fails to read or matrix wrong size. 
 */
static bool get_matrix_array(XFILE * fp, const int nrow, const int ncol, real_t *values) {

    if(values == NULL) {return false;}

    MAT mat = read_MAT_from_column_file(fp);
    if (mat == NULL) {return false;}
    if ((mat->nrow != nrow) || (mat->ncol != ncol)) {
        free_MAT(mat);
        return false;
    }

    /* Store values into array */
    memcpy(values, mat->x, nrow * ncol * sizeof(real_t));

    free_MAT(mat);
    return true;
}

/** Output a set of array values as a list of columns. */
static void show_array_bycol (XFILE * fp, real_t *values, const int nrow, const int ncol, const char *name) {

    xfprintf(fp, "# %s\n", name);
    MAT mat = coerce_MAT_from_array(nrow, ncol, values);
    write_MAT_to_column_file (fp, mat, true);

    /* only free top structure as values refer to array */
    xfree(mat);
}

/** Output a set of array values as a list of rows. */
static void show_array_byrow (XFILE * fp, real_t *values, const int nrow, const int ncol, const char *name) {

    xfprintf(fp, "%s:\n", name);
    MAT mat = coerce_MAT_from_array(nrow, ncol, values);
    show_MAT_rownum(fp, mat, nrow, ncol, false);

    /* only free top structure as values refer to array */
    xfree(mat);
}


/* public functions */

QUALTAB free_qualtab(QUALTAB qtab){
    if(NULL==qtab){ return NULL; }
    xfree(qtab->filename);
    xfree(qtab->comment);
    xfree(qtab->prior);
    xfree(qtab->next);
    xfree(qtab->triplet);
    xfree(qtab);
    return NULL;
}

QUALTAB new_qualtab(void){
	// Note: only triplets of bases implemented
    static const char *NOCALIB_TEXT = "uncalibrated table";

	QUALTAB qtab = calloc(1,sizeof(*qtab));
	if(NULL==qtab){ return NULL;}

	qtab->filename = calloc(strlen(NOCALIB_TEXT)+1, sizeof(char));
    if(NULL==qtab->filename) {goto cleanup;}
	strcpy(qtab->filename, NOCALIB_TEXT);
	qtab->comment = NULL;
	qtab->slope = 1.0;
	qtab->intcept = 0.0;
	qtab->prior = calloc(NBASE*NBASE,sizeof(real_t));
	qtab->next = calloc(NBASE*NBASE,sizeof(real_t));
	qtab->triplet  = calloc(NBASE*NBASE*NBASE,sizeof(real_t));
	if(NULL==qtab->prior || NULL==qtab->next || NULL==qtab->triplet){
		goto cleanup;
	}
	return qtab;

cleanup:
	free_qualtab(qtab);
	return NULL;
}

void show_qualtab(XFILE * fp, const QUALTAB qtab) {
    if(NULL==fp){ return;}
    if(NULL==qtab){ return;}

    xfprintf(fp, "Filename: %s\n", qtab->filename);
    xfprintf(fp, "Comment: %s\n", qtab->comment);
    xfprintf(fp, "%s: %f\n", MESS_TEXT[E_SLOPE], qtab->slope);
    xfprintf(fp, "%s: %f\n", MESS_TEXT[E_INTERCEPT], qtab->intcept);

    show_array_byrow (fp, qtab->prior, NBASE, NBASE, MESS_TEXT[E_PRIOR]);
    show_array_byrow (fp, qtab->next, NBASE, NBASE, MESS_TEXT[E_NEXT]);
    show_array_byrow (fp, qtab->triplet, NBASE, NBASE * NBASE, MESS_TEXT[E_TRIPLET]);
}

/**
 * Adjust quality score for base using a linear calibration and neighbours.
 * Vectors used are defined in calibration table, from default or optionally read in.
 * First and last bases of a read are special cases, dealt with by setting the
 * prior or next base (respectively) to be NUC_AMBIG.
 */
real_t adjust_quality(const real_t qual, const NUC prior, const NUC base, const NUC next){

    if(isambig(base)){ return MIN_QUALITY; }
    real_t new_qual = QTable->intcept + QTable->slope * qual;

    if(!isambig(prior)){
        new_qual += QTable->prior[prior*NBASE+base];
    }
    if(!isambig(next)){
        new_qual += QTable->next[next*NBASE+base];
    }
    if(!isambig(next) && !isambig(prior)){
        new_qual += QTable->triplet[(next*NBASE+prior)*NBASE+base];
    }
    return new_qual;
}

/**
 * Output quality calibration table.
 */
void output_quality_table(void) {

    if (QualOut) {
        XFILE *fpout = open_run_output("tab");
        if (!xfisnull(fpout)) {

            /* start with any header, followed by a blank line */
            if (QTable->comment == NULL) {
                xfprintf(fpout, "\n");
            }
            else {
                xfprintf(fpout, "%s\n\n", QTable->comment);
            }

            /* show_array needs an array */
            real_t value[1];
            value[0] = QTable->intcept;
            show_array_bycol (fpout, value, 1, 1, MESS_TEXT[E_INTERCEPT]);
            value[0] = QTable->slope;
            show_array_bycol (fpout, value, 1, 1, MESS_TEXT[E_SLOPE]);
            show_array_bycol (fpout, QTable->prior, NBASE, NBASE, MESS_TEXT[E_PRIOR]);
            show_array_bycol (fpout, QTable->next, NBASE, NBASE, MESS_TEXT[E_NEXT]);
            show_array_bycol (fpout, QTable->triplet, NBASE, NBASE * NBASE, MESS_TEXT[E_TRIPLET]);
        }
        xfclose(fpout);
    }
}

/**
 * Read in quality calibration table if supplied and use to replace the default values.
 * Returns false if file supplied but failed to read.
 */
bool read_quality_table(void) {

    XFILE *fp = NULL;
    QUALTAB qtable = NULL;
    int index = 0;

    if (matrix_from_file(E_QUALTAB)) {
        fp = open_matrix(E_QUALTAB);
        if (fp == NULL) {return false;}

        qtable = new_qualtab();
        if (qtable == NULL) {goto cleanup;}

        /* read the header */
        size_t comlen = 0;
        /* default comment is null so nothing to free */
        qtable->comment = xfgettok(fp, &comlen, "\n\n");

        /* read intercept and slope */
        index = E_INTERCEPT;
        if (!get_matrix_single(fp, &qtable->intcept)) {goto cleanup;}
        index = E_SLOPE;
        if (!get_matrix_single(fp, &qtable->slope)) {goto cleanup;}

        /* read adjustments */
        index = E_PRIOR;
        if(!get_matrix_array(fp, NBASE, NBASE, qtable->prior)) {goto cleanup;}
        index = E_NEXT;
        if(!get_matrix_array(fp, NBASE, NBASE, qtable->next)) {goto cleanup;}
        index = E_TRIPLET;
        if(!get_matrix_array(fp, NBASE, NBASE * NBASE, qtable->triplet)) {goto cleanup;}

        /* use new table */
        QTable = qtable;
        NewQTable = true;
        xfclose(fp);
    }

    return true;

 cleanup:
    message(E_BAD_INPUT_SS, MSG_FATAL, "quality calibration table", MESS_TEXT[index]);
    free_qualtab(qtable);
    xfclose(fp);
    return false;
}

/** Unset show quality calibration table flag. */
void set_noqualout(void) {

    QualOut = false;
}

/** Tidy up; call at program shutdown. */
void tidyup_qual_table(void) {

    if (NewQTable) {
        free_qualtab(QTable);
    }
}
