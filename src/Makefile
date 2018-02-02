COMP = gfortran
#COMP = ifort 
FLAGS = -fopenmp -O3 
#FLAGS = 
#BLASLIB = -lblas
MKLLIB = -lblas -llapack
#MKLLIB = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

../bin/samsa.x : constants.o io.o data.o ints.o wavefun.o energy.o cc.o props.o trans.o occupy.o scf.o dscf.o rdmft.o samsa.o
	${COMP} ${FLAGS} constants.o io.o data.o ints.o wavefun.o energy.o cc.o props.o trans.o occupy.o scf.o dscf.o rdmft.o samsa.o ${MKLLIB} -o ../bin/samsa.x

constants.o : constants.f90
	${COMP} ${FLAGS} -c constants.f90 

io.o : io.f90
	${COMP} ${FLAGS} -c io.f90

data.o : data.f90
	${COMP} ${FLAGS} -c data.f90 

ints.o : ints.f90
	${COMP} ${FLAGS} -c ints.f90 

wavefun.o : wavefun.f90
	${COMP} ${FLAGS} -c wavefun.f90

energy.o : energy.f90
	${COMP} ${FLAGS} -c energy.f90 

cc.o : cc.f90
	${COMP} ${FLAGS} -c cc.f90 

props.o : props.f90
	${COMP} ${FLAGS} -c props.f90

trans.o : trans.f90
	${COMP} ${FLAGS} -c trans.f90

occupy.o : occupy.f90
	${COMP} ${FLAGS} -c occupy.f90

scf.o : scf.f90
	${COMP} ${FLAGS} -c scf.f90

dscf.o : dscf.f90
	${COMP} ${FLAGS} -c dscf.f90

rdmft.o : rdmft.f90
	${COMP} ${FLAGS} -c rdmft.f90

samsa.o : samsa.f90
	${COMP} ${FLAGS} -c samsa.f90 

clean:
	rm *.o *.mod

