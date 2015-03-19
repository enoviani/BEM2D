default: neum2diri

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

neum2diri:  green.o chkopts
	-${CLINKER} main.c green.o ../petscInterface/petscInterface.o -I../petscInterface/ ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}

convergencetest: neum2diri
	./a.out 50
	./a.out 100
	./a.out 200
	./a.out 400
	make clean

clean::
	@rm -fr *.o a.out* *~ *.txt
