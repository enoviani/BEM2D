include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


default: display.o green.o emptyMatrix.o matrixOperation.o petscInterface.o  chkopts
	-${CLINKER} main.c display.o green.o emptyMatrix.o matrixOperation.o petscInterface.o ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
