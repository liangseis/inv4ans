OBIN= ~/bin
LIB= .

surftom: surfinv.for surfpath.for Sphere.for output.for lsqr.for inversion.for Datamod.for CORRECT.F90 azindex.for
	gfortran -g -O surfinv.for surfpath.for Sphere.for output.for lsqr.for inversion.for Datamod.for CORRECT.F90 azindex.for -o surftom2_etp0.50

