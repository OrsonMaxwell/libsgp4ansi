GCC = gcc
GPP = g++
OBJ_LIB = libsgp4ansi.o epoch.o coord.o vector.o
OBJ_REF_TEST = sgp4ext.o sgp4io.o sgp4unit.o testcpp.o
OUTPUTDIR=Plot
LIBFLAGS = -fPIC -DUSE_WGS72
CCFLAGS = -std=c11 -lm
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

lib: ${OBJ_LIB}
	${GCC} ${CCFLAGS} ${LIBFLAGS} -shared ${OBJ_LIB} -o ./${OUTPUTDIR}/${LIB_NAME}

ref_test: ${OBJ_REF_TEST}
	cd AIAA && ${GPP} ${CXXFLAGS} ${OBJ_REF_TEST} -o ../${OUTPUTDIR}/${REF_TEST_NAME}

test: lib
	${GCC} ${CCFLAGS} test.c -o ./${OUTPUTDIR}/${TEST_NAME} -L./${OUTPUTDIR}/ -lsgp4ansi

%.o: AIAA/%.cpp
	cd AIAA && ${GPP} ${CXXFLAGS} -c $*.cpp

${OBJ_LIB}: %.o: %.c
	${GCC} ${CCFLAGS} ${LIBFLAGS} -c $< -o $@

lib_clean:
	rm -f libsgp4ansi.o
	rm -f epoch.o
	rm -f coord.o
	rm -f vector.o
	rm -f Plot/${LIB_NAME}

ref_test_clean:
	rm -f AIAA/*.o
	rm -f Plot/${REF_TEST_NAME}

test_clean:
	rm -f Plot/${TEST_NAME}

clean: lib_clean ref_test_clean test_clean
	rm -f Plot/*.out
