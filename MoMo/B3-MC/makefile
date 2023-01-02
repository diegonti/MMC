FC=gfortran
FLAGS=-g

# Linking
main.x: ran2.o r1279.o practice3.o
	$(FC) $(FLAGS) $^ -o $@

# Compiling
ran2.o: ran2.f
	$(FC) $(FLAGS) -c ran2.f

r1279.o: r1279.f90
	$(FC) $(FLAGS) -c r1279.f90

practice3.o: practice3.f90
	$(FC) $(FLAGS) -c practice3.f90

clean:
	rm *.o
