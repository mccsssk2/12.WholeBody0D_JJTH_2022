# makefile.
# customized for jjj/pm3
# change all DRIVER instances to whatever it needs to be.
# whoisthis=pm3user
whoisthis=kharches
#
CC       		= mpicc
CFLAGS   	= -g -O2
INCLUDE  	= /home/${whoisthis}/software/sundials/instdir/include
MY_APP	 	= hdt
LIB	 		= -L/home/${whoisthis}/software/sundials/instdir/lib

hdt:	ursino.c
	${CC} ${CFLAGS} -I${INCLUDE} -c ursino.c -o ursino.o
	${CC} ${CFLAGS} ursino.o -I${INCLUDE} -lm ${LIB} -lsundials_cvodes -lsundials_nvecserial -o ${MY_APP}

run:
	./${MY_APP}

clean:
	rm  ${MY_APP}

veryclean:
	rm -r *.dat ${MY_APP} *.o *~ *.txt *.out
