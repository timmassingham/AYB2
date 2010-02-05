CC = gcc
CFLAGS = -Wall -O3 -funroll-loops -DNDEBUG -fomit-frame-pointer -std=gnu99
#CFLAGS = -std=c99 -Wall -g -std=gnu99
#CFLAGS = -g -pg -std=c99
LD = ld
LDFLAGS =  -lz -lbz2 -lm -lc -lblas -llapack
INCFLAGS = 
DEFINES = -DUSEFLOAT
objects =  xio.o matrix.o cluster.o tile.o

ayb: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) -o ../bin/$@ $(LDFLAGS) $(objects)

test: test-cluster test-tile

test-cluster: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o ../bin/$@ -DTEST cluster.c $(filter-out cluster.o,$(objects))

test-tile: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o ../bin/$@ -DTEST tile.c $(filter-out tile.o,$(objects))


.c.o:
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

.f.o:
	g77 -O3 -o $@ -c $<

clean:
	rm -f *~ *.o  *.obj  *.exe *.lib
	mkdir -p ../bin

archive: clean
	cd ../.. && tar zcf ayb.tgz --exclude=.git --exclude=bin --exclude=archive --exclude=test ayb
