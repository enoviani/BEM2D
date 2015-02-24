include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


default: petscInterface.o  chkopts
	-${CLINKER} testPetscInterface.c petscInterface.o  ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}

cleanup: 
	rm -fr *.o a.out *.dSYM *~
