FC = gfortran

Fsw = -O2  -ffree-line-length-none 
Lsw = 

F1=test_SOLVERS
F2=../src/SOLVERS
F3=../src/utilites

$(F1): $(F1).o $(F2).o  $(F3).o
	$(FC) $(F3).o $(F2).o  $(F1).o -o $(F1) $(Lsw)

$(F1).o: $(F1).f90
	$(FC) -c $(Fsw) $(F1).f90

$(F2).o: $(F2).f90
	$(FC) -c $(Fsw) $(F2).f90 -o $(F2).o

$(F3).o: $(F3).f90
	$(FC) -c $(Fsw) $(F3).f90 -o $(F3).o






