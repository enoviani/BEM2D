default: common 

common: interface greenstuff geom disp
	ar -cvq common.a *.o

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

interface: petscInterface.o  chkopts
greenstuff: green.o  chkopts
geom: geometry.o chkopts
disp: display.o chkopts

test: clean
	make -s common
	@-${CLINKER} -o petsctests test_petscInterface.c common.a ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./petsctests
	@-${CLINKER} -o greentests test_green.c common.a ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./greentests
	@rm petsctests greentests

clean ::
	@rm -f *.a
