default: potflow

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

potflow: display.o green.o matrixOperation.o chkopts
	-${CLINKER} main.c display.o green.o matrixOperation.o ../petscInterface/petscInterface.o -I../petscInterface/ ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}

clean::
	rm -fr *.o a.out* *~ *.txt
