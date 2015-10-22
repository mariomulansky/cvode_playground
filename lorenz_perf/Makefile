CXX = icpc

CXXFLAGS = -I${CVODE_INCLUDE} -I${BOOST_ROOT} -I${ODEINT_ROOT} -fast -ipo -DNDEBUG -static
LDFLAGS = -L${CVODE_LIB} -lsundials_cvode -lsundials_nvecserial -Wl,-rpath,/home/mario/sundials/lib

lorenz_cvode: lorenz_cvode.cpp
	${CXX} ${CXXFLAGS} lorenz_cvode.cpp ${LDFLAGS} -o lorenz_cvode

lorenz_odeint: lorenz_odeint.cpp
	${CXX} -std=c++0x ${CXXFLAGS} lorenz_odeint.cpp -o lorenz_odeint
