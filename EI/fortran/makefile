# Start of the makefile
# Defining variables
# Makefile
programa.x: molarMass.o selectGroup.o main.o
	gfortran -o programa.x  molarMass.o selectGroup.o main.o
molarMass.o: molarMass.f90
	gfortran -c  molarMass.f90
main.o: main.f90
	gfortran -c  main.f90
selectGroup.o: selectGroup.f90
	gfortran -c  selectGroup.f90

clean:
	rm -f programa.x
	rm -f molarMass.o selectGroup.o main.o
run:
	./programa.x
