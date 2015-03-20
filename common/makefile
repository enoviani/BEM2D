default: common 

common: interface greenstuff
	ar -cvq common.a *.o

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

interface: petscInterface.o  chkopts
greenstuff: green.o  chkopts

test: common
	-${CLINKER} -o tests test_petscInterface.c petscInterface.o ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./tests
	@make clean
	@rm tests

clean ::
	@rm *.a
