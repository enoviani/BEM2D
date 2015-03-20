default: interface

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

interface: petscInterface.o  chkopts
test: interface
	-${CLINKER} -o tests test_petscInterface.c petscInterface.o ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./tests
	@make clean
	@rm tests
