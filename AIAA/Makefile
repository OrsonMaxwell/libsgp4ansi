GPP = g++
OUTPUTDIR=Plot
CXXFLAGS = -lm -std=c++11
RELEASEFLAGS = -O3 -flto
DEBUGFLAGS = -g3
DEBUG = 0
MATHTRACE = 0

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

ifeq ($(MATHTRACE), 1)
	CXXFLAGS += -DMATH_TRACE
endif

all: ../${OUTPUTDIR}/${REF_TEST_NAME}
	@echo "${CURDIR}/../${OUTPUTDIR}/${REF_TEST_NAME} done"

../${OUTPUTDIR}/${REF_TEST_NAME}: sgp4ext.o sgp4io.o sgp4unit.o testcpp.o
	${GPP} ${CXXFLAGS} sgp4ext.o sgp4io.o sgp4unit.o testcpp.o -o ../${OUTPUTDIR}/${REF_TEST_NAME}

sgp4ext.o: sgp4ext.cpp sgp4ext.h
	${GPP} ${CXXFLAGS} -c sgp4ext.cpp

sgp4io.o: sgp4io.cpp sgp4io.h
	${GPP} ${CXXFLAGS} -c sgp4io.cpp

sgp4unit.o: sgp4unit.cpp sgp4unit.h
	${GPP} ${CXXFLAGS} -c sgp4unit.cpp

testcpp.o: testcpp.cpp sgp4ext.h sgp4io.h sgp4unit.h
	${GPP} ${CXXFLAGS} -c testcpp.cpp

clean:
	rm -f *.o
