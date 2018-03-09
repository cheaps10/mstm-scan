#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
MAKE = make
AR = ar


all:  mstm-scan

mstm-scan: mpidefs-serial-scan.o intrinsics-scan.o modules-scan.o main-scan.o
	gfortran-mp-6 mpidefs-serial-scan.o intrinsics-scan.o modules-scan.o main-scan.o -o ../run_mstm_jobs/mstm-scan.x 

mpidefs-serial-scan.o:
	gfortran-mp-6 -c mpidefs-serial-scan.f90

intrinsics-scan.o:
	gfortran-mp-6 -c intrinsics-scan.f90


modules-scan.o:
	gfortran-mp-6 -c modules-scan.f90


main-scan.o:
	gfortran-mp-6 -c main-scan.f90


clean:
	rm -f *.o *.mod 
