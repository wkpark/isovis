# This makefile is for NCSA Isovis.
# Uncomment the appropriate three lines (CFLAGS, LIBS, INCLUDES) for
# your system.
#
CC = cc
LFLAGS  =  

# non SGI systems
#CFLAGS   = -O
#LIBS     = -lvg -ldf -lm
#INCLUDES = -I. -I/ncsa/hdf/include

# SGI systems
CFLAGS   = -O -DSGI
LIBS     = -lvg -ldf -lgl_s -lm
INCLUDES = -I. -I/rels/shared/soft/hdf/include

# SGI - DGL systems
#CFLAGS   = -O -DSGI
#LIBS     = -lvg -ldf -ldgl -lm -lsun -lbsd
#INCLUDES = -I. -I/ncsa/hdf/include

#
all: \
	isovis
#
isovis:	file_io.o graphics.o iso.o main.o misc.o parse_opts.o
	$(CC) $(LFLAGS) -o isovis file_io.o graphics.o iso.o \
	main.o misc.o parse_opts.o $(LIBS)
#
file_io.o:	file_io.c isovis.h
	$(CC) $(CFLAGS) -c file_io.c $(INCLUDES)
#
graphics.o:	graphics.c isovis.h
	$(CC) $(CFLAGS) -c graphics.c $(INCLUDES)
#
iso.o:	iso.c isovis.h
	$(CC) $(CFLAGS) -c iso.c $(INCLUDES)
#
main.o:	main.c isovis.h
	$(CC) $(CFLAGS) -c main.c $(INCLUDES)
#
misc.o:	misc.c isovis.h
	$(CC) $(CFLAGS) -c misc.c $(INCLUDES)
#
parse_opts.o:	parse_opts.c isovis.h
	$(CC) $(CFLAGS) -c parse_opts.c $(INCLUDES)
#
clean:
	-rm -f ,*
	-rm *.o *.u
