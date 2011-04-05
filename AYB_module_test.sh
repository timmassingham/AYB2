#!/bin/bash
# Run AYB module tests.
# Compares output against reference and sends a report to stdout.
# File locations assume run from AYB top directory.
# No arguments.

# Location and name of test files and other inputs
INDIR=test
OUTDIR=log
REFDIR=logref
LOGEXT=log
REFEXT=ref
ERRFILE=msgerr
NC=6
ININT=test100_int.txt
INCIF=s_2_0001.cif
INFOLDER=/nfs/research2/goldman/NextGen/Data/sample-runfolder
LANETILE=L2T8-9
INSEQ=ACGTaatgXc
INSEQPHRED=nuc_phred.txt
MATFROM=mat_from.txt
MATTO=mat_to.txt
INTOK=xio_in.txt
SEP=ow

echo "AYB module test results  " $(date +"%d %B %Y %H:%M")
echo ""

MODULE=cluster
echo "Testing $MODULE"
# arguments ncycle _int.txt_filename [cif_filename]
#bin/test-$MODULE $NC $INDIR/$ININT >$OUTDIR/$MODULE.$LOGEXT  2>$OUTDIR/$ERRFILE.$LOGEXT
bin/test-$MODULE $NC $INDIR/$ININT $INDIR/$INCIF >$OUTDIR/$MODULE.$LOGEXT  2>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=matrix
echo "Testing $MODULE"
# arguments appendto appendfrom (filenames)
bin/test-$MODULE $INDIR/$MATTO $INDIR/$MATFROM >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT 

MODULE=message
echo "Testing $MODULE"
# arguments none
bin/test-$MODULE >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=mpn
echo "Testing $MODULE"
# arguments none
bin/test-$MODULE >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=nuc
echo "Testing $MODULE"
# arguments sequence [sequence/quality filename]
#bin/test-$MODULE $INSEQ >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
bin/test-$MODULE $INSEQ $INDIR/$INSEQPHRED >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT 

MODULE=tile
echo "Testing $MODULE"
# arguments ncycle _int.txt_filename [cif_filename run-folder lane_tile_range]
#bin/test-$MODULE $NC $INDIR/$ININT >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
#bin/test-$MODULE $NC $INDIR/$ININT $INDIR/$INCIF >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
bin/test-$MODULE $NC $INDIR/$ININT $INDIR/$INCIF $INFOLDER $LANETILE >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=xio
echo "Testing $MODULE"
# arguments separator filename
bin/test-$MODULE $SEP $INDIR/$INTOK >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

# compare the error output file
echo "Checking program messages"
diff -q -s $OUTDIR/$ERRFILE.$LOGEXT $REFDIR/$ERRFILE.$REFEXT

