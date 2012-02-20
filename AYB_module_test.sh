#!/bin/bash
# Run AYB module tests.
# Compares output against reference and sends a report to stdout.
# File locations assume run from AYB top directory.
# No arguments.

# Location and name of test files and other inputs
BIN=bin
INDIR=test
OUTDIR=log
REFDIR=logref
LOGEXT=log
REFEXT=ref
ERRFILE=msgerr
NC=6
NC2=20
ININT=test100_int.txt
INCIF=s_2_0001.cif
INFOLDER=/nfs/research2/goldman/NextGen/Data/sample-runfolder
LANETILE=L2T8-9
MATFROM=mat_from.txt
MATTO=mat_to.txt
NODIR=nodir.txt
NEWDIR=newdir
MESSLOG=messagelog
NMIX=3
NITER=20
INMIXVAL=test100_lambdas.txt
INSEQ=ACGTaatgXc
INSEQPHRED=nuc_phred.txt
INSPIKE=spike_in.txt
OUTSPIKE=spike_out.txt
INTOK=xio_in.txt
OUTXIO=xio_out
SEP=ow

echo "AYB module test results  " $(date +"%d %B %Y %H:%M")
echo ""

MODULE=cluster
echo "Testing $MODULE"
# arguments ncycle _int.txt_filename [cif_filename]
#$BIN/test-$MODULE $NC $INDIR/$ININT >$OUTDIR/$MODULE.$LOGEXT  2>$OUTDIR/$ERRFILE.$LOGEXT
$BIN/test-$MODULE $NC $INDIR/$ININT $INDIR/$INCIF >$OUTDIR/$MODULE.$LOGEXT  2>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=matrix
echo "Testing $MODULE"
# arguments appendto appendfrom (filenames)
$BIN/test-$MODULE $INDIR/$MATTO $INDIR/$MATFROM >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT 

MODULE=message
echo "Testing $MODULE I/O"
# arguments log_path more_flag
# test log dir not a dir
cat > $NODIR << ENDIT
empty
ENDIT
$BIN/test-$MODULE $NODIR/notadir F >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
rm $NODIR

#test log dir not exist, creatable
if [ `find ./ -name $NEWDIR` ]; then
    chmod u+w $NEWDIR
    rm -r $NEWDIR
fi
$BIN/test-$MODULE $NEWDIR/canmake.log F >>$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
echo -e "Log file in new directory \c"
if [ `find $NEWDIR -name canmake.log` ]; then
    echo "created ok"
else
    echo "failed to create"
fi

#test log dir not writable
chmod u-w $NEWDIR
$BIN/test-$MODULE $NEWDIR/nomake.log F >>$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
#test log dir not exist, not creatable
$BIN/test-$MODULE $NEWDIR/nomakedir/nomake.log F >>$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
chmod u+w $NEWDIR
rm -r $NEWDIR

echo "Testing $MODULE"
# arguments log_path more_flag
$BIN/test-$MODULE $OUTDIR/$MESSLOG.$LOGEXT T >>$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT
diff -s $OUTDIR/$MESSLOG.$LOGEXT $REFDIR/$MESSLOG.$REFEXT

MODULE=mixnormal
echo "Testing $MODULE"
# arguments nmix niter filename
$BIN/test-$MODULE $NMIX $NITER $INDIR/$INMIXVAL >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=mpn
echo "Testing $MODULE"
# arguments none
$BIN/test-$MODULE >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=nuc
echo "Testing $MODULE"
# arguments sequence [sequence/quality filename]
#$BIN/test-$MODULE $INSEQ >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
$BIN/test-$MODULE $INSEQ $INDIR/$INSEQPHRED >$OUTDIR/$MODULE.$LOGEXT  2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT 

MODULE=spikein
echo "Testing $MODULE"
# arguments ncycle infilename outfilename
$BIN/test-$MODULE $NC2 $INDIR/$INSPIKE $OUTDIR/$OUTSPIKE >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT
# cpmpare outfile to infile - will be identical if ncycle matches infile
diff -q -s $INDIR/$INSPIKE $OUTDIR/$OUTSPIKE

MODULE=tile
echo "Testing $MODULE"
# arguments ncycle _int.txt_filename [cif_filename run-folder lane_tile_range]
#$BIN/test-$MODULE $NC $INDIR/$ININT >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
#$BIN/test-$MODULE $NC $INDIR/$ININT $INDIR/$INCIF >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
$BIN/test-$MODULE $NC $INDIR/$ININT $INDIR/$INCIF $INFOLDER $LANETILE >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

MODULE=xio
echo "Testing $MODULE"
# arguments separator infilename outfilename
$BIN/test-$MODULE $SEP $INDIR/$INTOK $OUTDIR/$OUTXIO >$OUTDIR/$MODULE.$LOGEXT 2>>$OUTDIR/$ERRFILE.$LOGEXT
diff -q -s $OUTDIR/$MODULE.$LOGEXT $REFDIR/$MODULE.$REFEXT

# compare the error output file
echo "Checking program messages"
diff -s $OUTDIR/$ERRFILE.$LOGEXT $REFDIR/$ERRFILE.$REFEXT

