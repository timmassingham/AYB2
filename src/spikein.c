/** 
 * \file spikein.c
 * Spike-in data Class.
 * Reads in and stores spike-in data from a tab separated file containing cluster numbers and sequence.
 * Validated results are stored in a list of spikeins. 
 * Get next spikein then returns each spikein in turn until no more.
 *//* 
 *  Created : 19 Jan 2012
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
#include "dirio.h"
#include "message.h"
#include "spikein.h"

#define X(A) A ## SPIKEIN
    #include "list.def"
#undef X

/* constants */

/* members */

/* doxygen confused by LIST */
/** SpikeList = NULL; List of spike-in data.
 * Also static LIST(SPIKEIN) Current = NULL; 
 *
 * Current spike-in data block.
 */
static LIST(SPIKEIN) SpikeList = NULL;
static LIST(SPIKEIN) Current = NULL;

/* helps doxygen confused by LIST */
static void __attribute__((__used__)) null1(void){}


/* private functions */

/** Free spike-in data list. */
static LIST(SPIKEIN) free_SPIKEIN_list(LIST(SPIKEIN) spikelist) {
    free_LIST(SPIKEIN)(spikelist);
    return NULL;
}

/** Show first n blocks of spike-in data list.  */
static void show_SPIKEIN_list(XFILE * fp, LIST(SPIKEIN) spikelist, const uint32_t nblock) __attribute__((used));
static void show_SPIKEIN_list(XFILE * fp, LIST(SPIKEIN) spikelist, const uint32_t nblock) {

    if (NULL == fp) { return; }
    if (NULL == spikelist) { return; }
    
    LIST(SPIKEIN) node = spikelist;
    unsigned int n = 0;
    while ((NULL != node) && (nblock == 0 || n < nblock)) {
        show_SPIKEIN(fp, node->elt);
        node = node->nxt;
        n++;
    }
}

/** Read in spike-in data list. */
static LIST(SPIKEIN) read_SPIKEIN_list(XFILE * fp, const uint32_t ncycle) {

    if (NULL == fp) { return NULL; }

    LIST(SPIKEIN) spikelist = NULL;
    SPIKEIN newblock = NULL;
    LIST(SPIKEIN) tail = NULL;
    bool err = false;

    while(newblock = read_SPIKEIN(fp, ncycle, &err), NULL != newblock ) {

        /* add new data block to list */
        if (spikelist == NULL) {
            spikelist = cons_LIST(SPIKEIN)(newblock, spikelist);
            tail = spikelist;
        }
        else {
            tail = rcons_LIST(SPIKEIN)(newblock, tail);
        }
   	}
    
    if (err) {
        spikelist = free_SPIKEIN_list(spikelist);
    }
    
    return spikelist;
}


/* public functions */

/* standard functions */
SPIKEIN new_SPIKEIN(const uint32_t ncycle) {
    SPIKEIN spikein = calloc(1, sizeof(*spikein));
    if (NULL == spikein) { return NULL;}
    spikein->kbases = new_ARRAY(NUC)(ncycle);
    return spikein;
}

SPIKEIN free_SPIKEIN(SPIKEIN spikein) {
    if(NULL == spikein) { return NULL; }
    free_ARRAY(NUC)(spikein->kbases);
    xfree(spikein);
    return NULL;
}

SPIKEIN copy_SPIKEIN(const SPIKEIN spikein) {
    if (NULL == spikein) { return NULL; }
    SPIKEIN spike_copy = calloc(1, sizeof(*spikein));
    if(NULL == spike_copy) { return NULL; }

    spike_copy->knum = spikein->knum;
    spike_copy->kbases = copy_ARRAY(NUC)(spikein->kbases);
    if(NULL != spikein->kbases.elt && NULL == spike_copy->kbases.elt){ 
        free_SPIKEIN(spike_copy);
        return NULL;
    }

    return spike_copy;
}

void show_SPIKEIN(XFILE * fp, const SPIKEIN spikein) {
    if(NULL == fp) { return; }
    if(NULL == spikein) { return; }
    
    /* convert from 0-based for cluster loop to 1-based for file */
    xfprintf(fp, "%d\t", spikein->knum + 1);
    show_ARRAY(NUC)(fp, spikein->kbases, "", 0);
}

/** 
 * Read a single line of spike-in data.
 * Expects tab separated cluster number and sequence.
 * Sequence must be at least ncycle long.
 */
SPIKEIN read_SPIKEIN(XFILE * fp, const uint32_t ncycle, bool *err) {
    if (NULL == fp) { return NULL; }
    
    /* get a line from file */
    char *line = NULL;
    size_t len = 0;
    line = xfgetln(fp, &len);
    if (NULL == line) { return NULL; }

    char msg[80];
    char * ptr = line;
    SPIKEIN spikein = new_SPIKEIN(ncycle);
    if (NULL == spikein) {
        strcpy(msg, "memory allocation");
        goto cleanup; 
    }

    /* line starts with the cluster number */
    long n = strtol(ptr, &ptr, 0);
    if (n <= 0) {
        strcpy(msg, "invalid cluster number");
        goto cleanup; 
    }
    /* convert from 1-based in file to 0-based to match cluster loop */
    spikein->knum = n - 1;
    
    /* check for tab and skip if found */
    if ('\t' != ptr[0]) {
        strcpy(msg, "file not tab separated");
        goto cleanup; 
    }
    ptr++;
    
    /* check length of base string and store */
    if (strlen(ptr) < ncycle) {
        strcpy(msg, "insufficient cycles");
        goto cleanup; 
    }

    for( uint32_t i = 0; i < ncycle; i++){
        spikein->kbases.elt[i] = nuc_from_char(ptr[i]);
    }
    
    /* check for any erroneous bases, will read as an ambiguous base */
    if (has_ambiguous_base(spikein->kbases.elt, ncycle)) {
        strcpy(msg, "invalid base");
        goto cleanup; 
    }

    xfree(line);
    return spikein;
    
cleanup:
    message(E_BAD_TXT_SS, MSG_ERR, "Read spike-in", msg);
    *err = true;
    xfree(line);
    free_SPIKEIN(spikein);
    return NULL;
}


/** Return the next spike-in cluster. Returns NULL if are none or no more. */
SPIKEIN get_next_spikein(void) {

    if (SpikeList == NULL) { return NULL; }
    if (Current == NULL) {
        Current = SpikeList;
    }
    else {
        Current = Current->nxt;
    }
    if (Current == NULL) { return NULL; }
    return Current->elt;
}

/** 
 * Read the spike-in data file.
 * Returns true if successfully read.
 */
bool read_spikein_file(const uint32_t ncycle, const int blk) {

    bool res = false;

    XFILE *fp = open_spikein(blk);
    if (!xfisnull(fp)) {
        SpikeList = read_SPIKEIN_list(fp, ncycle);
        if (NULL == SpikeList) {
            message(E_BAD_INPUT_S, MSG_ERR, get_last_spikein());
        }
        else {
            res = true;
        }
    }
    xfclose(fp);
    
    return res;
}

/** Tidy up; call at program shutdown (or before). */
void tidyup_spikein(void) {

    /* clear the structure */
    SpikeList = free_SPIKEIN_list(SpikeList);
    Current = NULL;
}


#ifdef TEST
#include <err.h>

static const char BADTAB[] = "9 ACGTACGTAC";
static const char BADNUM[] = "-1\tACGTACGTAC";
static const char BADBASE[] = "9\tACGXACGTAC";
static const int NBAD = 10;

/* Read invalid spikein data. */
void read_spikein_error(XFILE * fp , char message[]) {
    bool err = false;

    xfprintf(xstdout, "Read spike-in with %s\n", message);
    SPIKEIN spikein = read_SPIKEIN(fp, NBAD, &err);
    if ((spikein==NULL) && (err==true)) {
        xfputs("Return value null and error, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null or no error, not ok\n", xstdout);
        spikein = free_SPIKEIN(spikein);
    }
}

int main ( int argc, char * argv[]){
    if(argc<3){
        /* arguments are input file and output file */
        errx(EXIT_FAILURE, "Usage: test-spikein ncycle in_filename out_filename");
    }
    
    unsigned int ncycle = 0;
    sscanf(argv[1],"%u", &ncycle);

    /* null file */
    xfputs("Read null file\n", xstdout);
    LIST(SPIKEIN) spikelist = read_SPIKEIN_list(NULL, ncycle);
    if (spikelist==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        spikelist = free_SPIKEIN_list(spikelist);
    }
    
    xfputs("Write null file\n", xstdout);
    show_SPIKEIN_list(NULL, spikelist, 0);
    xfputs("No result, ok\n", xstdout);

    /* empty list */
    xfputs("Write empty spike-in list\n", xstdout);
    XFILE * fp = xfopen(argv[3], XFILE_RAW, "w");
    if (xfisnull(fp)) {
	    errx(EXIT_FAILURE, "Failed to open %s for output", argv[3]);
    }

    show_SPIKEIN_list(fp, spikelist, 0);
    fp = xfclose(fp);

    xfputs("Read empty spike-in file\n", xstdout);
    fp = xfopen(argv[3], XFILE_UNKNOWN, "r");
    if (xfisnull(fp)) {
        errx(EXIT_FAILURE, "Failed to re-open output file");
    }

    spikelist = read_SPIKEIN_list(fp, ncycle);
    if (spikelist==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        spikelist = free_SPIKEIN_list(spikelist);
    }
    fp = xfclose(fp);

    xfputs("Get next from empty spike-in list\n", xstdout);
    SpikeList = spikelist;
    SPIKEIN spikein = get_next_spikein();
    if (spikein==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        spikein = free_SPIKEIN(spikein);
    }

    xfputs("Free empty spike-in list\n", xstdout);
    spikelist = free_SPIKEIN_list(spikelist);
    if (spikelist==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
    }

    /* invalid file contents */
    xfputs("(Writing invalid spike-in data to file)\n", xstdout);
    fp = xfopen(argv[3], XFILE_RAW, "w");
    if (xfisnull(fp)) {
	    errx(EXIT_FAILURE, "Failed to open %s for output", argv[3]);
    }

    xfprintf(fp, "%s\n%s\n%s\n", BADTAB, BADNUM, BADBASE);
    fp = xfclose(fp);

    fp = xfopen(argv[3], XFILE_UNKNOWN, "r");
    if (xfisnull(fp)) {
        errx(EXIT_FAILURE, "Failed to re-open output file");
    }
    
    read_spikein_error(fp, "missing tab");
    read_spikein_error(fp, "invalid cluster number");
    read_spikein_error(fp, "invalid base");
    fp = xfclose(fp);

    /* wrong ncycle */
    xfputs("Read spike-in file with not enough cycles\n", xstdout);
    fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    if (xfisnull(fp)) {
        errx(EXIT_FAILURE, "Failed to open supplied spike data file");
    }

    spikelist = read_SPIKEIN_list(fp, ncycle + 50);
    if (spikelist==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
        spikelist = free_SPIKEIN_list(spikelist);
    }
    fp = xfclose(fp);
    
    /* normal operation */
    xfputs("Read spike-in file\n", xstdout);
    fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
    if (xfisnull(fp)) {
        errx(EXIT_FAILURE, "Failed to open supplied spike data file");
    }

    spikelist = read_SPIKEIN_list(fp, ncycle);
    fp = xfclose(fp);

    xfputs("Write spike-in file\n", xstdout);
    fp = xfopen(argv[3], XFILE_RAW, "w");
    if (xfisnull(fp)) {
	    errx(EXIT_FAILURE, "Failed to open %s for output", argv[3]);
    }

    show_SPIKEIN_list(fp, spikelist, 0);
    xfputs("**Diff input/output files to check result\n", xstdout);
    fp = xfclose(fp);

    xfputs("Get next, first in list\n", xstdout);
    SpikeList = spikelist;
    spikein = get_next_spikein();
    show_SPIKEIN(xstdout, spikein);

    xfputs("Get next, next in list\n", xstdout);
    spikein = get_next_spikein();
    show_SPIKEIN(xstdout, spikein);

    xfputs("Copy spikein\n", xstdout);
    SPIKEIN spike2 = copy_SPIKEIN(spikein);
    show_SPIKEIN(xstdout, spike2);
    free_SPIKEIN(spike2);

    /* count number of blocks */
    uint32_t cnt = 1;
    while (spikein != NULL) {
        cnt++;
        spikein = get_next_spikein();
    }
    xfputs("Write fewer spike-in blocks\n", xstdout);
    show_SPIKEIN_list(xstdout, spikelist, 2);

    xfputs("Write more spike-in blocks\n", xstdout);
    show_SPIKEIN_list(xstdout, spikelist, cnt + 2);

    xfputs("Free spike-in list\n", xstdout);
    spikelist = free_SPIKEIN_list(spikelist);
    if (spikelist==NULL) {
        xfputs("Return value null, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null, not ok\n", xstdout);
    }

    return EXIT_SUCCESS;
}
#endif
