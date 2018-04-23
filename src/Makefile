FC=gfortran
FFLAGS=-fcheck=bounds -O2

bindir:= ../bin

progdirs:= parser
modules:= gas/gas.f90 util/util.f90 integral/integral.f90 pop/pop.f90 photo/photo.f90

.PHONY: all prog

all: prog
	$(FC) -o $(bindir)/myChamber $(FFLAGS) myChamber.f90 $(modules)

prog: 
		@set -e; \
        for i in $(progdirs); do \
                if [ -d $$i ]; then \
                  if [ -f $$i/Makefile ]; then \
                    $(MAKE) -C $$i all;\
                    fi; \
                  fi; \
        done;

