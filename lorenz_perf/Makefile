# CXX = icpc
# CC = icc
# OPTIMIZER = -fast -ipo
CC = gcc
CXX = g++
OPTIMIZER = -Ofast
ANALYZER = -fno-omit-frame-pointer

CFLAGS = -I${CVODE_INCLUDE} -DNDEBUG -static ${OPTIMIZER} ${ANALYZER}
CXXFLAGS = -I${BOOST_ROOT} -DNDEBUG -static ${OPTIMIZER} ${ANALYZER}
LDFLAGS = -L${CVODE_LIB} -lsundials_cvode -lsundials_nvecserial -lm -Wl,-rpath,/home/mario/sundials/lib

lorenz_cvode: lorenz_cvode.c
	${CC} ${CFLAGS} lorenz_cvode.c ${LDFLAGS} -o lorenz_cvode

lorenz_odeint: lorenz_odeint.cpp
	${CXX} -std=c++0x ${CXXFLAGS} lorenz_odeint.cpp -o lorenz_odeint
