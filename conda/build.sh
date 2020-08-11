export CFLAGS="${CFLAGS} -funroll-loops -DNDEBUG -std=gnu99 -fopenmp -I${BUILD_PREFIX}/include"

sed -e 's,^\(CFLAGS .*\)$,#\1,' \
    -e 's,^\(LDFLAGS .*\)$,#\1,' \
    -e 's,^\(CC .*\)$,#\1,' \
  $SRC_DIR/src/Makefile > $SRC_DIR/src/Makefile.fixed

(cd $SRC_DIR/src && make -f Makefile.fixed -j $CPU_COUNT)

install -m 755 -s $SRC_DIR/bin/AYB $PREFIX/bin

