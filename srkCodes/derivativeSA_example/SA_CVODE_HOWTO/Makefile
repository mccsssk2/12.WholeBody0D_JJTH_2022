# -----------------------------------------------------------------
# test the forward SA file and customise.
#
SHELL = /bin/sh

prefix       = /home/sanjay/sundials
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

CPP      = cc -E
CPPFLAGS = 
CC       = cc
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = -lm 

INCLUDES = -I${includedir}
LIBRARIES = -lsundials_cvodes -lsundials_nvecserial ${LIBS}
LIBRARIES_BL = 

EXAMPLES = cvsRoberts_FSA_dns 

OBJECTS = ${EXAMPLES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${CFLAGS} ${LDFLAGS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL}" ; \
	  ${CC} -o $${i} $${i}.o ${CFLAGS} ${LDFLAGS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL}; \
	done

run:
	./cvsRoberts_FSA_dns -sensi stg t

clean:
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}

# -----------------------------------------------------------------------------------------

