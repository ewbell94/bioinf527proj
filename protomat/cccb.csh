#!/bin/csh
#       Compile one program with blimps routines
#	Usage to compile program.c:   cccb.csh <program>

#  NOTE:  Set BLIMPS_DIR and CC for your site

set BLIMPS_DIR = /howard/blimps/blimps-3.2
set CC = "/opt/SUNWspro/bin/cc -fast"
#set CC = "/opt/gnu/bin/gcc -g"

$CC -I$BLIMPS_DIR/include -L$BLIMPS_DIR/lib -o $1 $1.c -lblocks -lmatrix -lsequences -lfrequency -lpssm -lfiles -lmemory -loptions -lgcode -lstrutil -lsl -lm
