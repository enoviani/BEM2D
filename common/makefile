default: common 

common: interface greenstuff geom expinteg
	ar -cvq common.a *.o

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

interface: petscInterface.o  chkopts
greenstuff: green.o  chkopts
geom: geometry.o chkopts
expinteg: expE1.o chkopts

test: clean
	make -s common
	@-${CLINKER} -o petsctests test_petscInterface.c common.a ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./petsctests
	@-${CLINKER} -o greentests test_green.c common.a ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./greentests
	@-${CLINKER} -o expinttests test_expE1.c common.a ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB}
	./expinttests
	@rm petsctests greentests expinttests

clean ::
	@rm -f *.a
