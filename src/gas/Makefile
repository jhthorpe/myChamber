bindir:= ../../bin

FC= gfortran 
FFLAGS= -fcheck=bounds

util:= ../util/util.f90

all:
	$(FC) $(FFLAGS) test.f90 gas.f90 $(util) 
