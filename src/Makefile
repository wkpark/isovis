# This makefile is for OCLSC/NCSA Isovis.
# Uncomment the appropriate lines for your system ...
# your system.
#
#
CC = cc -cckr

# on your local cray
#CC = scc -hintrinsics,inline,aggress,report=isv

LFLAGS  =  

# set these if you have HDF and VSET libraries
#N_CFLAGS = -DVSET -DHDF
#N_LIBS = -lvg -ldf
#N_INCLUDES = -I/usr/local/include/ncsa

# non SGI systems
#CFLAGS   = -O $(N_CFLAGS)
#LIBS     = $(N_LIBS) -lm
#INCLUDES = -I. $(N_INCLUDES)

# SGI systems

# uncomment if you have/like my popup menu
#PUP = -lpup
#PUP = ../pup/libpup.a
#PFLAG = -DMYPUP

# uncomment if you like to save things in SGI image format
IFLAG = -DSGI_IMAGE
IMAGE = -limage

# uncomment if you want to use polymode() on the VGX to do hidden line
POLYMODE = -DPOLYMODE
# uncomment of you want to use linesmooth() to antialias wire frame drawings
ANTI = -DANTIALIAS

CFLAGS   = -g -DSGI $(N_CFLAGS) -float ${PFLAG} ${IFLAG} $(POLYMODE) $(ANTI)
LIBS     = $(N_LIBS) ${PUP} -lgl_s -lfastm ${IMAGE} -lm -lc_s
INCLUDES = -I. $(N_INCLUDES)

# SGI - DGL systems
#CFLAGS   = -O -DSGI $(N_CFLAGS)
#LIBS     = $(N_LIBS) -ldgl -lm -lsun -lbsd
#INCLUDES = -I. $(N_INCLUDES)

# if you want the more complicated, but twice as fast marching cubes
# implementation
ISO = niso
#if you don't
#ISO = iso

SRCS = file_io.c graphics.c ${ISO}.c main.c misc.c parse_opts.c tmesh.c \
	smooth.c lines.c compact.c
OBJS = file_io.o graphics.o ${ISO}.o main.o misc.o parse_opts.o tmesh.o \
	smooth.o lines.o compact.o

#
all: \
	isovis-O
#
# for -O3
isovis-O3 : $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -O3 -o isovis $(SRCS) $(LIBS)

# for regular compilation
isovis-O:	$(OBJS)
	$(CC) $(CFLAGS) -o isovis $(OBJS) $(LIBS)
#
file_io.o:	isovis.h
graphics.o:	isovis.h
iso.o:		isovis.h
main.o:		isovis.h
misc.o:		isovis.h
parse_opts.o:	isovis.h

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c
#
clean:
	-rm -f ,*
	-rm *.o *.u isovis.1

#
TAR_FILES = file_io.c graphics.c iso.c niso.c main.c misc.c parse_opts.c \
	tmesh.c smooth.c lines.c compact.c cell_table.h isovis.h Makefile \
	isovis.1 isovis.man README

tar : ${TAR_FILES}
	rm -f isovis.tar.Z
	tar cvLf isovis.tar ${TAR_FILES}
	compress isovis.tar

isovis.1 : isovis.man
	nroff -man isovis.man > isovis.1
