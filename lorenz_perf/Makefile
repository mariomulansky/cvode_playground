# CXX = icpc
# CC = icc
OPTIMIZER = -fast -ipo
CC = gcc
CXX = g++
OPTIMIZER = -O3

CFLAGS = -I${CVODE_INCLUDE} -DNDEBUG -static ${OPTIMIZER}
CXXFLAGS = -I${BOOST_ROOT} -I${ODEINT_ROOT} -DNDEBUG -static ${OPTIMIZER}
LDFLAGS = -L${CVODE_LIB} -lsundials_cvode -lsundials_nvecserial -lm -Wl,-rpath,/home/mario/sundials/lib

lorenz_cvode: lorenz_cvode.c
	${CC} ${CFLAGS} lorenz_cvode.c ${LDFLAGS} -o lorenz_cvode

lorenz_odeint: lorenz_odeint.cpp
	${CXX} -std=c++0x ${CXXFLAGS} lorenz_odeint.cpp -o lorenz_odeint
