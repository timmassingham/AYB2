/** 
 * \file nuc.c
 * Nucleotide Base and Phred Quality Score Class.
 *//* 
 *  Created : 2010
 *  Author : Tim Massingham/Hazel Marsden
 *
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
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

#include <tgmath.h>
#include "message.h"
#include "nuc.h"
#include "utility.h"


/* constants */
/* None      */

/* members */
/* None    */


/* private functions */

/*
static bool isprob( const real_t p){
    if(p>1.0){ return false; }
    if(p<0.0){ return false; }
    return true;
}
*/

/* public functions */

void show_NUC(XFILE * fp, const NUC nuc){
    validate(NULL!=fp,);
    xfputc(char_from_nuc(nuc),fp);
}

void show_PHREDCHAR(XFILE *fp, const PHREDCHAR pc){
    validate(NULL!=fp,);
    /* replace with minimum if not a printable character */
    if(pc<MIN_PHRED || pc>MAX_PHRED) {
        xfputc(MIN_PHRED,fp);
    }
    else {
        xfputc(pc,fp);
    }
}

NUC read_NUC(XFILE *fp){
    return nuc_from_char(xfgetc(fp));
}

PHREDCHAR read_PHREDCHAR(XFILE * fp){
    return phredchar_from_char(xfgetc(fp));
}

NUC nuc_from_char( const char c){
    switch(c){
        case 'A':
        case 'a': return NUC_A;
        case 'C':
        case 'c': return NUC_C;
        case 'G':
        case 'g': return NUC_G;
        case 'T':
        case 't': return NUC_T;
        case 'N':
        case 'n': return NUC_AMBIG;
        default:
            message (E_BAD_NUC_C, MSG_WARN, c);
    }
    return NUC_AMBIG;
}

char char_from_nuc(const NUC nuc){
    switch(nuc){
    case NUC_A: return 'A';
    case NUC_C: return 'C';
    case NUC_G: return 'G';
    case NUC_T: return 'T';
    }
    return 'N';
}

ARRAY(NUC) nucs_from_string( const char * nucstr ){
    validate(NULL!=nucstr,null_ARRAY(NUC));
    const uint32_t len = strlen(nucstr);
    ARRAY(NUC) nucs = new_ARRAY(NUC)(len);
    validate(0!=nucs.nelt,nucs);
    for( uint32_t i=0 ; i<len ; i++){
        nucs.elt[i] = nuc_from_char(nucstr[i]);
    }
    return nucs;
}

NUC complement(const NUC nuc){
    switch(nuc){
    case NUC_A: return NUC_T;
    case NUC_C: return NUC_G;
    case NUC_G: return NUC_C;
    case NUC_T: return NUC_A;
    }
    return NUC_AMBIG;
}

ARRAY(NUC) reverse_complement(const ARRAY(NUC) nucs){
    validate(NULL!=nucs.elt,null_ARRAY(NUC));
    ARRAY(NUC) new_nuc = new_ARRAY(NUC)(nucs.nelt);
    validate(NULL!=new_nuc.elt,new_nuc);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        new_nuc.elt[i] = complement(nucs.elt[nucs.nelt-i-1]);
    }
    return new_nuc;
}

/** Convert a char to PHRED-style character, ensuring within range. */
PHREDCHAR phredchar_from_char( const char c){
    validate(c>=MIN_PHRED,MIN_PHRED);
    validate(c<=MAX_PHRED,MAX_PHRED);
    return c;
}

/** Convert probability of error to PHRED-style character representation. */
PHREDCHAR phredchar_from_prob( real_t p){
//    validate(isprob(p),ERR_PHRED);
    real_t c = MIN_PHRED - 10.*log1p(-p)/log(10.);
    if(isnan(c)){ c = MIN_PHRED; }              // probability > 1
    else if(!isfinite(c)){ c = MAX_PHRED; }     // probability = 1 
    if(c<MIN_PHRED){c=MIN_PHRED;}
    if(c>MAX_PHRED){c=MAX_PHRED;}
    return (PHREDCHAR)(c+0.5);
}

/** Convert probability of error to PHRED-style quality value. */
real_t quality_from_prob( real_t p){
//    validate(isprob(p),ERR_PHRED);
	return -10.*log1p(-p)/log(10.);
}

/** Convert PHRED-style quality value to character representation. */
PHREDCHAR phredchar_from_quality( real_t qual){
    real_t c  = MIN_PHRED + qual;
    if(isnan(c)){ c = MIN_PHRED; }              // probability > 1
    else if(!isfinite(c)){ c = MAX_PHRED; }     // probability = 1 
    if(c<MIN_PHRED){c=MIN_PHRED;}
    if(c>MAX_PHRED){c=MAX_PHRED;}
    return (PHREDCHAR)(c+0.5);
}

/** Convert PHRED-style quality value to nearest valid integer. */
int qualint_from_quality( real_t qual){
    if(isnan(qual)){ qual = MIN_QUALITY; }              // probability > 1
    else if(!isfinite(qual)){ qual = MAX_QUALITY; }     // probability = 1 
    if(qual<MIN_QUALITY){ qual = MIN_QUALITY; }
    if(qual>MAX_QUALITY){ qual = MAX_QUALITY; }
    return (int)(qual+0.5);
}

/** Convert character representation to PHRED-style quality integer value. */
int qualint_from_phredchar( const PHREDCHAR pc){
    int qual = (int)pc - MIN_PHRED;
    if(qual<MIN_QUALITY){ qual = MIN_QUALITY; }
    if(qual>MAX_QUALITY){ qual = MAX_QUALITY; }    
    return qual;
}


#ifdef TEST
#include <err.h>
#include <string.h>
#include <limits.h>

static real_t probs[] = {-0.1, 0.0, 0.5, 0.8, 0.9, 0.955333, 0.99, 0.993, 0.9994, 0.99995, 0.999996, 0.9999997, 0.99999998, 1.0, 1.1};

int main ( int argc, char * argv[]){
    if(argc<2){
        /* argument is test sequence and optional sequence/phred file with a line of sequence */
        /* and up to 2 lines of potential quality characters, for valid and invalid               */
        errx(EXIT_FAILURE, "Usage: test-nuc sequence [sequence/quality filename]");
    }

    /* optional read sequence/quality testing */
    if (argc > 2) {
        XFILE * fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
        if (xfisnull(fp)) {
            errx(EXIT_FAILURE, "Failed to open supplied sequence/quality file");
        }
        /* find length of lines and output as string */
        char *line = NULL;
        size_t len[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {
            line = xfgetln(fp, &len[i]);
            xfprintf(xstdout, "Test line %d              : %s\n", i + 1, line);
            xfree(line);
        }
        xfclose(fp);

        /* reopen to read as nuc/phredchar */
        fp = xfopen(argv[2], XFILE_UNKNOWN, "r");
        xfputs("Read in sequence         : ", xstdout);
        NUC nuc;
        for (int i = 0; i < len[0]; i++) {
            nuc = read_NUC(fp);
            show_NUC(xstdout, nuc);
        }
        xfputc('\n', xstdout);
        /* skip end of line */
        xfgetc(fp);
        
        xfputs("Read in valid qualchar   : ", xstdout);
        PHREDCHAR pc;
        for (int i = 0; i < len[1]; i++) {
            pc = read_PHREDCHAR(fp);
            show_PHREDCHAR(xstdout, pc);
        }
        xfputc('\n', xstdout);
        /* skip end of line */
        xfgetc(fp);
              
        xfputs("Read in invalid qualchar : ", xstdout);
        for (int i = 0; i < len[2]; i++) {
            pc = read_PHREDCHAR(fp);
            show_PHREDCHAR(xstdout, pc);
        }
        xfprintf(xstdout, "\n(should all be %c or %c)\n\n", MIN_PHRED, MAX_PHRED);
        xfclose(fp);
    }

    xfprintf(xstdout, "Supplied sequence string : %s\n", argv[1]);
    xfprintf(xstdout, "Sequence length          : %u\n", strlen(argv[1]));
    
    xfputs("Null sequence            : ", xstdout);
    ARRAY(NUC) nucs = nucs_from_string(NULL);
    if ((nucs.elt==NULL) && (nucs.nelt==0)) {
        xfputs("Return value null array, ok\n", xstdout);
    }
    else {
        xfputs("Return value not null array, not ok\n", xstdout);
    }
    free_ARRAY(NUC)(nucs);

    nucs = nucs_from_string(argv[1]);
    xfprintf(xstdout, "Stored length            : %u\n", nucs.nelt);
    xfputs("Store/write sequence     : ", xstdout);
    for (uint32_t i = 0; i < nucs.nelt; i++){
        show_NUC(xstdout, nucs.elt[i]);
    }
    xfputc('\n', xstdout);

    ARRAY(NUC) rc = reverse_complement(nucs);
    xfputs("Reversed sequence        : ", xstdout);
    for (uint32_t i = 0; i < rc.nelt; i++){
        show_NUC(xstdout, rc.elt[i]);
    }
    xfputs("\n\n", xstdout);
    free_ARRAY(NUC)(nucs);
    free_ARRAY(NUC)(rc);
    
    xfputs("Outside range qualchars  : ", xstdout);
    for (char ch = 0; ch < MIN_PHRED; ch++) {
        show_PHREDCHAR(xstdout, ch);
    }
    xfputs("...", xstdout);
    /* use int because char will alway be <= CHAR_MAX */
    for (int ch = MAX_PHRED + 1; ch <= CHAR_MAX; ch++) {
        show_PHREDCHAR(xstdout, (char)ch);
    }
    xfprintf(xstdout, "\n(should all be %c)\n", MIN_PHRED);

    xfputs("Coerced to quality range : ", xstdout);
    for (char ch = 0; ch < MIN_PHRED; ch++) {
        show_PHREDCHAR(xstdout, phredchar_from_char(ch));
    }
    xfputs("...", xstdout);
    for (int ch = MAX_PHRED + 1; ch <= CHAR_MAX; ch++) {
        show_PHREDCHAR(xstdout, phredchar_from_char((char)ch));
    }
    xfprintf(xstdout, "\n(should all be %c or %c)\n", MIN_PHRED, MAX_PHRED);
    
    xfputs("Range of qualchars       : ", xstdout);
    for (char ch = MIN_PHRED; ch <= MAX_PHRED; ch++) {
        show_PHREDCHAR(xstdout, phredchar_from_char(ch));
    }
    xfputs("\n\n", xstdout);
    
    xfputs("Quality character from probability and qualint from qualchar:\n", xstdout);
    int nprob = sizeof(probs) / sizeof(real_t);
    PHREDCHAR pc;
    int q;
    for (int i = 0; i < nprob; i++) {
        pc = phredchar_from_prob(probs[i]);
        q = qualint_from_phredchar(pc);
        xfprintf(xstdout, "prob: % 10.8f, value: %3d, qualchar: ", probs[i], pc);
        show_PHREDCHAR(xstdout, pc);
        xfprintf(xstdout, ", qualint: %2d\n", q);
    }

    xfputs("Quality character from probability via value and qualint from value:\n", xstdout);
    real_t val;
    for (int i = 0; i < nprob; i++) {
        val = quality_from_prob(probs[i]);
        pc = phredchar_from_quality(val);
        q = qualint_from_quality(val);
        xfprintf(xstdout, "prob: % 10.8f, qual value: %7.4f, qualchar: ", probs[i], val);
        show_PHREDCHAR(xstdout, pc);
        xfprintf(xstdout, ", qualint: %2d\n", q);
    }

    return EXIT_SUCCESS;
}
#endif
