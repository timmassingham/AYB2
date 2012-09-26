# makefile for AYB

BINDIR = ../bin
MANDIR = ../man
CC = gcc
FC = gfortran
CFLAGS = -Wall -O3 -funroll-loops -DNDEBUG -std=gnu99 -fopenmp
LDFLAGS =  -lm -lz -lbz2 -lblas -llapack
INCFLAGS = 
DEFINES =
objects =  ayb_main.o ayb_model.o ayb_options.o ayb_version.o ayb.o calibration.o call_bases.o \
           cif.o cluster.o conjugate.o datablock.o dirio.o handler.o intensities.o lambda.o matrix.o \
           message.o mixnormal.o mpn.o nuc.o qual_table.o spikein.o statistics.o tile.o utility.o weibull.o xio.o

AYB: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/$@ $(LDFLAGS) $(objects)
	
test: test-cluster test-matrix test-message test-mixnormal test-mpn test-nuc test-spikein test-tile test-xio

test-cluster: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST cluster.c $(filter-out cluster.o ayb_main.o,$(objects))

test-conjugate: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST conjugate.c $(filter-out conjugate.o ayb_main.o,$(objects))

test-matrix: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST matrix.c $(filter-out matrix.o ayb_main.o,$(objects))

test-message: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST message.c $(filter-out message.o ayb_main.o,$(objects))

test-mixnormal: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST mixnormal.c $(filter-out mixnormal.o ayb_main.o,$(objects))

test-mpn: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST mpn.c $(filter-out mpn.o ayb_main.o,$(objects))

test-nuc: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST nuc.c $(filter-out nuc.o ayb_main.o,$(objects))

test-spikein: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST spikein.c $(filter-out spikein.o ayb_main.o,$(objects))

test-statistics: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST statistics.c $(filter-out statistics.o ayb_main.o,$(objects))

test-tile: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST tile.c $(filter-out tile.o ayb_main.o,$(objects))

test-weibull: statistics.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST weibull.c $< && \
	   ln -fs $(BINDIR)/$@ $(BINDIR)/qweibull && \
	   ln -fs $(BINDIR)/$@ $(BINDIR)/pweibull && \
	   ln -fs $(BINDIR)/$@ $(BINDIR)/dweibull && \
	   ln -fs $(BINDIR)/$@ $(BINDIR)/fweibull

test-xio: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ -DTEST xio.c $(filter-out xio.o ayb_main.o,$(objects))

ayb_options.o: ayb_options.c ayb_usage.h ayb_help.h
ayb.o: ayb.c message.h
ayb_main.o: ayb_main.c message.h
ayb_model.o: ayb_model.c message.h
datablock.o: datablock.c message.h
dirio.o: dirio.c message.h
nuc.o: nuc.c message.h
qual_table.o: qual_table.c message.h
spikein.o: spikein.c message.h

.c.o:
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

.f.o:
	$(FC) -m64 -O3 -o $@ -c $<

clean:
	rm -f $(BINDIR)/* *~ *.o  *.obj  *.exe *.lib ../*.lis $(MANDIR)/*.1 $(MANDIR)/*.1.html
	find .. -name "*.log" -exec rm {} \;
	
%.1: %.1.txt
	cd $(MANDIR); a2x -d manpage -f manpage $<

%.1.html: %.1.txt
	cd $(MANDIR); asciidoc -d manpage $<

man: man_troff man_html

man_troff: $(MANDIR)/AYB.1

man_html: $(MANDIR)/AYB.1.html

archive: clean man
	cd ../.. && tar zcf AYB.tgz --exclude=archive --exclude=test --exclude=doc --exclude=.git --exclude=.gitignore AYB2

archive_test:
	cd ../.. && tar zcf AYB_test.tgz AYB2/test
