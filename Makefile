GCC = gcc
GPP = g++
OBJ_LIB = libsgp4ansi.o epoch.o coord.o vector.o
OUTPUTDIR=Plot
LIBFLAGS = -fPIC -DUSE_WGS72
CCFLAGS = -std=c11
LIBS = -lm
CXXFLAGS = -lm
RELEASEFLAGS = -O3 -flto
DEBUGFLAGS = -g3
DEBUG = 0

ifeq ($(OS),Windows_NT)
	LIB_NAME = libsgp4ansi.dll
	TEST_NAME = ansi.exe
	REF_TEST_NAME = aiaa.exe
else
	LIB_NAME = libsgp4ansi.so
	TEST_NAME = ansi
	REF_TEST_NAME = aiaa
endif

ifeq ($(DEBUG), 1)
	CCFLAGS += ${DEBUGFLAGS}
	CXXFLAGS += ${DEBUGFLAGS}
else
	CCFLAGS += ${RELEASEFLAGS}
	CXXFLAGS += ${RELEASEFLAGS}
endif

all: test ref_test

test: library ${TEST_NAME}
	cp ${TEST_NAME} ${OUTPUTDIR}/${TEST_NAME}

${TEST_NAME}:
	${GCC} ${CCFLAGS} -L${CURDIR} -o ${TEST_NAME} test.c ${LIBS} -lsgp4ansi

ref_test:
	${MAKE} -C AIAA DEBUG=${DEBUG} all

library: ${LIB_NAME}
	cp ${LIB_NAME} ${OUTPUTDIR}/${LIB_NAME}

${LIB_NAME}: libsgp4ansi.o epoch.o coord.o vector.o
	${GCC} ${CCFLAGS} ${LIBFLAGS} -shared libsgp4ansi.o epoch.o coord.o vector.o -o ${LIB_NAME} ${LIBS}

libsgp4ansi.o:
	${GCC} ${CCFLAGS} ${LIBFLAGS} -c libsgp4ansi.c ${LIBS}

epoch.o:
	${GCC} ${CCFLAGS} ${LIBFLAGS} -c epoch.c ${LIBS}

coord.o:
	${GCC} ${CCFLAGS} ${LIBFLAGS} -c coord.c ${LIBS}

vector.o:
	${GCC} ${CCFLAGS} ${LIBFLAGS} -c vector.c ${LIBS}

lib_clean:
	rm -f libsgp4ansi.o
	rm -f epoch.o
	rm -f coord.o
	rm -f vector.o
	rm -f ${LIB_NAME}
	rm -f ${OUTPUTDIR}/${LIB_NAME}

ref_test_clean:
	rm -f ${OUTPUTDIR}/${REF_TEST_NAME}
	${MAKE} -C AIAA clean

test_clean:
	rm -f ${TEST_NAME}
	rm -f ${OUTPUTDIR}/${TEST_NAME}

test_data_clean:
	rm -f ${OUTPUTDIR}/*.out

clean: lib_clean ref_test_clean test_clean test_data_clean
