GCC = gcc
GPP = g++
OBJ_LIB = libsgp4ansi.o epoch.o coord.o vector.o
OBJ_REF_TEST = $(addprefix AIAA/, sgp4ext.o sgp4io.o sgp4unit.o testcpp.o )
OUTPUTDIR=Plot
LIBFLAGS = -fPIC -DUSE_WGS72
CCFLAGS = -std=c11 -lm
CXXFLAGS = -lm
RELEASEFLAGS = -O3 -flto
DEBUGFLAGS = -g3
DEBUG = 0

ifeq ($(OS),Windows_NT)
	LIB_NAME = ($addprefix ${OUTPUTDIR}/, libsgp4ansi.dll)
	TEST_NAME = ($addprefix ${OUTPUTDIR}/, ansi.exe)
	REF_TEST_NAME = ($addprefix ${OUTPUTDIR}/, aiaa.exe)
else
	LIB_NAME = $(addprefix ${OUTPUTDIR}/, libsgp4ansi.so)
	TEST_NAME = $(addprefix ${OUTPUTDIR}/, ansi)
	REF_TEST_NAME = $(addprefix ${OUTPUTDIR}/, aiaa)
endif

ifeq ($(DEBUG), 1)
	CCFLAGS += ${DEBUGFLAGS}
	CXXFLAGS += ${DEBUGFLAGS}
else
	CCFLAGS += ${RELEASEFLAGS}
	CXXFLAGS += ${RELEASEFLAGS}
endif

all: ${LIB_NAME} ${TEST_NAME} ${REF_TEST_NAME}

${LIB_NAME}: ${OBJ_LIB}
	${GCC} ${CCFLAGS} ${LIBFLAGS} -shared ${OBJ_LIB} -o ${LIB_NAME}

${REF_TEST_NAME}: ${OBJ_REF_TEST}
	${GPP} ${CXXFLAGS} ${OBJ_REF_TEST} -o ${REF_TEST_NAME}

${TEST_NAME}:
	${GCC} ${CCFLAGS} -L${OUTPUTDIR}/ -lsgp4ansi test.c -o $@

${OBJ_REF_TEST}: %.o: %.cpp
	${GPP} ${CXXFLAGS} -c $< -o $@

${OBJ_LIB}: %.o: %.c
	${GCC} ${CCFLAGS} ${LIBFLAGS} -c $< -o $@

lib_clean:
	rm -f libsgp4ansi.o
	rm -f epoch.o
	rm -f coord.o
	rm -f vector.o
	rm -f ${LIB_NAME}

ref_test_clean:
	rm -f ${REF_TEST_NAME}
	rm -f AIAA/*.o

test_clean:
	rm -f ${TEST_NAME}

test_data_clean:
	rm -f ${OUTPUTDIR}/*.out

clean: lib_clean ref_test_clean test_clean test_data_clean
