default: neum2diri

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

neum2diri:  green.o chkopts
	-${CLINKER} main.c green.o ../petscInterface/petscInterface.o -I../petscInterface/ ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}

convergencetest: neum2diri
	./a.out 50 > res.txt
	./a.out 100 >> res.txt
	./a.out 200 >> res.txt
	./a.out 400 >> res.txt
	make clean

clean::
	@rm -fr *.o a.out* *~
