# Protomat Makefile for UNIX
#   Edit CC for your site
#   Type  "make protomat" to compile all protomat programs
#
CC = /usr/bin/gcc 
CFLAGS = -c
CL = /usr/bin/gcc 
LFLAGS = -lm -o
#
protomat:	protomot uextract motifj motomat matodat 

#	protomot - Extract groups from PROSITE
#
protomot:	protomot.o motmisc.o
		$(CL) protomot.o motmisc.o $(LFLAGS) protomot
#	uextract -  Extract sequences from SWISS-PROT
uextract:	uextract.o motmisc.o
		$(CL) uextract.o motmisc.o $(LFLAGS) uextract
#	motifj -  Find Hamilton Smith's motifs
motifj:		motifj.o motmisc.o
		$(CL) motifj.o motmisc.o $(LFLAGS) motifj
#	motomat - Assemble blocks
motomat:	motomat.o motmisc.o
		$(CL) motomat.o motmisc.o $(LFLAGS) motomat
#	matodat - Calibrate blocks
matodat:	matodat.o motmisc.o
		$(CL) matodat.o motmisc.o $(LFLAGS) matodat
#
#
protomot.o:	protomot.c motifj.h
		$(CC) $(CFLAGS) protomot.c
uextract.o:	uextract.c motifj.h
		$(CC) $(CFLAGS) uextract.c
motifj.o:	motifj.c motifj.h
		$(CC) $(CFLAGS) motifj.c
motomat.o:	motomat.c motomat2.c motifj.h
		$(CC) $(CFLAGS) motomat.c
matodat.o:	matodat.c motifj.h
		$(CC) $(CFLAGS) matodat.c
#
motmisc.o:	motmisc.c motifj.h
		$(CC) $(CFLAGS) motmisc.c
#
rmos:
	rm protomot.o
	rm uextract.o
	rm motifj.o
	rm motomat.o
	rm matodat.o
