# Grandi model make
# 
# files: source files: main.c
#        is taken care of
# executable file: gfb
#
# first define target file: gfb
#              dependencies are the object files that build the program

# SHELL = /bin/sh

CC       = gcc
CFLAGS   = -Wall -g -fno-inline
INCLUDE  = /home/sanjay/sundials/include
MY_APP	 = ord
LIB	 = -L/home/sanjay/sundials/lib

ord:	ord.c
	${CC} ${CFLAGS} -I${INCLUDE} -c ord.c -fpic -o ord.o
	${CC} ${CFLAGS} ord.o -I${INCLUDE} -g -lm ${LIB} -lsundials_cvodes -lsundials_nvecserial -mcmodel=large -o ${MY_APP}

# file numbers, if_sa for both gsa and fsa
sa:	ord
	./${MY_APP} 21 22 0

good:	gooddata.c
	gcc -o good gooddata.c -lm
	./good

clean:
	rm -f *.o ex*.txt

veryclean:
	rm -f *.o ex*.txt *.dat *.csv *.vtk
	rm -f ${MY_APP} *~
run:
	./${MY_APP} 1 10
	./${MY_APP} 21 25
	./${MY_APP} 31 39

edit:	ord.c
	gedit ord.c Makefile &

