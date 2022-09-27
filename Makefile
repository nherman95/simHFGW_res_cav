all : harmosc_solver.mod
	 python -m numpy.f2py -c -m harmosc  harmosc.f90 --f90flags="-fopenmp" --opt="-Ofast" -lgomp --fcompiler=gfortran

harmosc_solver.mod:
	gfortran -c harmosc.f90 -fopenmp -Ofast

clean:
	rm *.o *.mod *.so