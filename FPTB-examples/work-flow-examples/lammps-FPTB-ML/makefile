
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

INC_DIR = /path/to/include/dir
LAMMPS_DIR = /path/to/lammps/with/DEEPMD
EXAMPLE5_DIR = ./ 
EIGEN_DIR =/path/to/eigen
ARPACK_DIR = /path/to/arpack

CXXFLAGS= -I${INC_DIR} -I${LAMMPS_DIR} -I${EXAMPLE5_DIR} -I${EIGEN_DIR} -I${ARPACK_DIR} -L${LAMMPS_DIR} -L${ARPACK_DIR}/lib64 -std=c++2a -llammps
DEPS = ${INC_DIR}/blot/blotio.hpp ${INC_DIR}/blot/blotmath.hpp  ${INC_DIR}/blot/tight_binding.hpp ${INC_DIR}/blot/kspace.hpp ${INC_DIR}/blot/basics.hpp ${INC_DIR}/blot/units.hpp ${INC_DIR}/blot/dynamics.hpp ${INC_DIR}/blot/lookuptable.hpp ${EXAMPLE5_DIR}/init_tb.cpp ${EXAMPLE5_DIR}/init_hoppings.cpp ${EXAMPLE5_DIR}/perovskite_tb.cpp


ehrenfest: ehrenfest.o $(DEPS)
	-${CLINKER} $(CXXFLAGS) ehrenfest.o -o ehrenfest  ${SLEPC_EPS_LIB}
	${RM} *.o 
