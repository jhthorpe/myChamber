FC=gfortran
FFLAGS=-fcheck=bounds

bindir:= ../bin

progdirs:= parser

.PHONY: all prog

all: prog
	$(FC) -o $(bindir)/myChamber $(FFLAGS) myChamber.f90

prog: 
		@set -e; \
        for i in $(progdirs); do \
                if [ -d $$i ]; then \
                  if [ -f $$i/Makefile ]; then \
                    $(MAKE) -C $$i all;\
                    fi; \
                  fi; \
        done;

