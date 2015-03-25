default: potflow

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

potflow: display.o matrixOperation.o chkopts
	-${CLINKER} main.c display.o matrixOperation.o ../common/common.a -I../common/ ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}

clean::
	rm -fr *.o a.out* *~ *.txt
