# Malekar model make
# whoisthis=pm3user
whoisthis=kharches
CC       = mpicc
CFLAGS   =
INCLUDE  = /home/${whoisthis}/software/sundials/instdir/include
MY_APP = ra
LIB	 = -L/home/${whoisthis}/software/sundials/instdir/lib

ra:	generateRandomPars.c
	${CC} ${CFLAGS} -I${INCLUDE} -c generateRandomPars.c -o generateRandomPars.o
	${CC} ${CFLAGS} generateRandomPars.o -I${INCLUDE} -lm ${LIB} -lsundials_cvodes -lsundials_nvecserial -o ${MY_APP}

run:
	./${MY_APP}

veryclean:
	rm -r *.dat ${MY_APP} ${MY_APPI} *.o *~

