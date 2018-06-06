#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
MAKE = make
AR = ar
FC = gfortran-mp-6

all:  mstm-scan clean

mstm-scan: mpidefs-serial-scan.o intrinsics-scan.o modules-scan.o main-scan.o
	$(FC) mpidefs-serial-scan.o intrinsics-scan.o modules-scan.o main-scan.o -o mstm-scan.x 

mpidefs-serial-scan.o:
	$(FC) -c mpidefs-serial-scan.f90

intrinsics-scan.o:
	$(FC) -c intrinsics-scan.f90


modules-scan.o:
	$(FC) -c modules-scan.f90


main-scan.o:
	$(FC) -c main-scan.f90


clean:
	rm -f *.o *.mod 

install:
	mkdir -p $(HOME)/bin
	mv mstm-scan.x $(HOME)/bin
