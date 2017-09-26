GPP = g++
OUTPUTDIR=Plot
CXXFLAGS = -lm
RELEASEFLAGS = -O3 -flto
DEBUGFLAGS = -g3
DEBUG = 0

ifeq ($(OS),Windows_NT)
	REF_TEST_NAME = aiaa.exe
else
	REF_TEST_NAME = aiaa
endif

ifeq ($(DEBUG), 1)
	CXXFLAGS += ${DEBUGFLAGS}
else
	CXXFLAGS += ${RELEASEFLAGS}
endif

all: ${REF_TEST_NAME}
	cp ${REF_TEST_NAME} ../${OUTPUTDIR}/${REF_TEST_NAME}

${REF_TEST_NAME}: sgp4ext.o sgp4io.o sgp4unit.o testcpp.o
	${GPP} ${CXXFLAGS} sgp4ext.o sgp4io.o sgp4unit.o testcpp.o -o ${REF_TEST_NAME}

sgp4ext.o:
	${GPP} ${CXXFLAGS} -c sgp4ext.cpp

sgp4io.o:
	${GPP} ${CXXFLAGS} -c sgp4io.cpp

sgp4unit.o:
	${GPP} ${CXXFLAGS} -c sgp4unit.cpp

testcpp.o:
	${GPP} ${CXXFLAGS} -c testcpp.cpp

clean:
	rm -f ${REF_TEST_NAME}
	rm -f *.o
