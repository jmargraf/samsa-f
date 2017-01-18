../bin/samsa.x : constants.o data.o samsa.o
	gfortran constants.o data.o samsa.o -o ../bin/samsa.x

constants.o : constants.f90
	gfortran -c constants.f90

data.o : data.f90
	gfortran -c data.f90

samsa.o : samsa.f90
	gfortran -c samsa.f90
